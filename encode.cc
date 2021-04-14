/*
OFDM TV encoder

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#include <iostream>
#include <cassert>
#include <cmath>
#include "netpbm.hh"
#include "complex.hh"
#include "utils.hh"
#include "decibel.hh"
#include "fft.hh"
#include "wav.hh"
#include "pcm.hh"
#include "mls.hh"
#include "crc.hh"
#include "galois_field.hh"
#include "bose_chaudhuri_hocquenghem_encoder.hh"

template <typename value, typename cmplx, int rate>
struct Encoder
{
	static const int symbol_len = (1280 * rate) / 8000;
	static const int guard_len = symbol_len / 8;
	static const int img_width = 320;
	static const int img_height = 240;
	static const int frame_width = 32;
	static const int mls0_len = 127;
	static const int mls0_poly = 0b10001001;
	static const int mls1_len = img_width;
	static const int mls1_poly = 0b1100110001;
	static const int mls2_poly = 0b10001000000001011;
	static const int mls3_poly = 0b10111010010000001;
	static const int mls4_len = 255;
	static const int mls4_poly = 0b100101011;
	DSP::WritePCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	CODE::CRC<uint16_t> crc;
	CODE::BoseChaudhuriHocquenghemEncoder<255, 63> bchenc;
	cmplx fdom[2 * symbol_len];
	cmplx tdom[symbol_len];
	cmplx kern0[symbol_len], kern1[symbol_len];
	cmplx guard[guard_len];
	value rgb_line[2 * 3 * img_width];
	cmplx papr_min, papr_max;
	int img_off;
	int mls0_off;
	int mls1_off;
	int mls4_off;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
	void rgb_to_yuv(value *yuv, const value *rgb)
	{
		value WR(0.299), WB(0.114), WG(1-WR-WB), UMAX(0.493), VMAX(0.877);
		yuv[0] = WR * rgb[0] + WG * rgb[1] + WB * rgb[2];
		yuv[1] = (UMAX/(1-WB)) * (rgb[2]-yuv[0]);
		yuv[2] = (VMAX/(1-WR)) * (rgb[0]-yuv[0]);
	}
	void rgb_to_cmplx(cmplx *out0, cmplx *out1, const value *rgb0, const value *rgb1)
	{
		value yuv0[3], yuv1[3];
		rgb_to_yuv(yuv0, rgb0);
		rgb_to_yuv(yuv1, rgb1);
		out0[0] = cmplx(yuv0[0]+(yuv0[2]+yuv1[2])/2, yuv0[0]-(yuv0[1]+yuv1[1])/2);
		out1[0] = cmplx((yuv0[1]+yuv1[1])/2-yuv1[0], (yuv0[2]+yuv1[2])/2-yuv1[0]);
		rgb_to_yuv(yuv0, rgb0+3);
		rgb_to_yuv(yuv1, rgb1+3);
		out0[1] = cmplx((yuv0[1]+yuv1[1])/2-yuv0[0], (yuv0[2]+yuv1[2])/2-yuv0[0]);
		out1[1] = cmplx(yuv1[0]+(yuv0[2]+yuv1[2])/2, yuv1[0]-(yuv0[1]+yuv1[1])/2);
	}
	void improve_papr(const cmplx *kern)
	{
		for (int n = 0; n < 1000; ++n) {
			int peak = 0;
			for (int i = 1; i < symbol_len; ++i)
				if (norm(tdom[peak]) < norm(tdom[i]))
					peak = i;
			cmplx orig = tdom[peak];
			for (int i = 0; i < peak; ++i)
				tdom[i] -= orig * kern[symbol_len-peak+i];
			for (int i = peak; i < symbol_len; ++i)
				tdom[i] -= orig * kern[i-peak];
		}
	}
	void symbol(const cmplx *kern = 0)
	{
		bwd(tdom, fdom);
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] /= sqrt(value(8 * symbol_len));
		if (kern)
			improve_papr(kern);
		for (int i = 0; i < guard_len; ++i) {
			value x = value(i) / value(guard_len - 1);
			x = value(0.5) * (value(1) - std::cos(DSP::Const<value>::Pi() * x));
			guard[i] = DSP::lerp(x, guard[i], tdom[i+symbol_len-guard_len]);
		}
		cmplx peak, mean;
		for (int i = 0; i < symbol_len; ++i) {
			cmplx power(tdom[i].real() * tdom[i].real(), tdom[i].imag() * tdom[i].imag());
			peak = cmplx(std::max(peak.real(), power.real()), std::max(peak.imag(), power.imag()));
			mean += power;
		}
		if (mean.real() > 0 && mean.imag() > 0) {
			cmplx papr(peak.real() / mean.real(), peak.imag() / mean.imag());
			papr *= value(symbol_len);
			papr_min = cmplx(std::min(papr_min.real(), papr.real()), std::min(papr_min.imag(), papr.imag()));
			papr_max = cmplx(std::max(papr_max.real(), papr.real()), std::max(papr_max.imag(), papr.imag()));
		}
		pcm->write(reinterpret_cast<value *>(guard), guard_len, 2);
		pcm->write(reinterpret_cast<value *>(tdom), symbol_len, 2);
		for (int i = 0; i < guard_len; ++i)
			guard[i] = tdom[i];
	}
	void pilot_block()
	{
		CODE::MLS seq1(mls1_poly);
		value mls1_fac = sqrt(value(symbol_len) / value(mls1_len));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		for (int i = mls1_off; i < mls1_off + mls1_len; ++i)
			fdom[bin(i)] = mls1_fac * (1 - 2 * seq1());
		symbol(kern1);
	}
	void schmidl_cox(bool inverted = false)
	{
		CODE::MLS seq0(mls0_poly);
		value mls0_fac = sqrt(value(symbol_len) / value(mls0_len));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		fdom[bin(mls0_off-2)] = mls0_fac;
		for (int i = 0; i < mls0_len; ++i)
			fdom[bin(2*i+mls0_off)] = (1 - 2 * (inverted^seq0()));
		for (int i = 0; i < mls0_len; ++i)
			fdom[bin(2*i+mls0_off)] *= fdom[bin(2*(i-1)+mls0_off)];
		symbol();
	}
	void meta_data(uint64_t md)
	{
		uint8_t data[8] = { 0 }, parity[24] = { 0 };
		for (int i = 0; i < 47; ++i)
			CODE::set_be_bit(data, i, (md>>i)&1);
		crc.reset();
		uint16_t cs = crc(md);
		for (int i = 0; i < 16; ++i)
			CODE::set_be_bit(data, i+47, (cs>>i)&1);
		bchenc(data, parity);
		CODE::MLS seq4(mls4_poly);
		value mls4_fac = sqrt(value(symbol_len) / value(mls4_len));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		fdom[bin(mls4_off-1)] = mls4_fac;
		for (int i = 0; i < 63; ++i)
			fdom[bin(i+mls4_off)] = (1 - 2 * (CODE::get_be_bit(data, i)^seq4()));
		for (int i = 63; i < mls4_len; ++i)
			fdom[bin(i+mls4_off)] = (1 - 2 * (CODE::get_be_bit(parity, i-63)^seq4()));
		for (int i = 0; i < mls4_len; ++i)
			fdom[bin(i+mls4_off)] *= fdom[bin(i-1+mls4_off)];
		symbol(kern0);
	}
	Encoder(DSP::WritePCM<value> *pcm, DSP::ReadPEL<value> *pel, int freq_off, uint64_t call_sign) :
		pcm(pcm), crc(0xA8F4), bchenc({
			0b100011101, 0b101110111, 0b111110011, 0b101101001, 0b110111101,
			0b111100111, 0b100101011, 0b111010111, 0b000010011, 0b101100101,
			0b110001011, 0b101100011, 0b100011011, 0b100111111, 0b110001101,
			0b100101101, 0b101011111, 0b111111001, 0b111000011, 0b100111001,
			0b110101001, 0b000011111, 0b110000111, 0b110110001, 0b101001101})
	{
		mls1_off = img_off = (freq_off * symbol_len) / rate - img_width / 2;
		mls0_off = mls1_off + 34;
		mls4_off = mls1_off + 33;
		int car_min = mls1_off - frame_width;
		int car_max = mls1_off+mls1_len + frame_width;
		papr_min = cmplx(1000, 1000), papr_max = cmplx(-1000, -1000);
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		int count = 0;
		for (int i = car_min; i <= car_max; ++i) {
			if (i < mls4_off-1 || i >= mls4_off+mls4_len) {
				fdom[bin(i)] = 1;
				++count;
			}
		}
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] /= value(10 * count);
		bwd(kern0, fdom);
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		count = 0;
		for (int i = car_min; i <= car_max; ++i) {
			if (i < mls1_off || i >= mls1_off+mls1_len) {
				fdom[bin(i)] = 1;
				++count;
			}
		}
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] /= value(10 * count);
		bwd(kern1, fdom);
		pilot_block();
		schmidl_cox(true);
		meta_data(call_sign);
		CODE::MLS seq2(mls2_poly), seq3(mls3_poly);
		value img_fac = sqrt(value(symbol_len) / value(img_width));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		for (int j = 0; j < img_height; j += 2) {
			if (j%8 == 0)
				pilot_block();
			pel->read(rgb_line, 2 * img_width);
			for (int i = 0; i < img_width; i += 2)
				rgb_to_cmplx(fdom+bin(i+img_off), fdom+bin(i+img_off)+symbol_len, rgb_line+3*i, rgb_line+3*(img_width+i));
			for (int k = 0; k < 2; ++k) {
				for (int i = 0; i < img_width; ++i)
					fdom[bin(i+img_off)] = img_fac * cmplx(
						fdom[bin(i+img_off)+symbol_len*k].real() * (1 - 2 * seq2()),
						fdom[bin(i+img_off)+symbol_len*k].imag() * (1 - 2 * seq3()));
				symbol(kern1);
			}
		}
		pilot_block();
		schmidl_cox();
		meta_data(call_sign);
		pilot_block();
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		symbol();
		std::cerr << "real PAPR: " << DSP::decibel(papr_min.real()) << " .. " << DSP::decibel(papr_max.real()) << " dB" << std::endl;
		if (pcm->channels() == 2)
			std::cerr << "imag PAPR: " << DSP::decibel(papr_min.imag()) << " .. " << DSP::decibel(papr_max.imag()) << " dB" << std::endl;
	}
};

long long int base37_encoder(const char *str)
{
	long long int acc = 0;
	for (char c = *str++; c; c = *str++) {
		acc *= 37;
		if (c >= '0' && c <= '9')
			acc += c - '0' + 1;
		else if (c >= 'a' && c <= 'z')
			acc += c - 'a' + 11;
		else if (c >= 'A' && c <= 'Z')
			acc += c - 'A' + 11;
		else if (c != ' ')
			return -1;
	}
	return acc;
}

int main(int argc, char **argv)
{
	if (argc < 6 || argc > 8) {
		std::cerr << "usage: " << argv[0] << " OUTPUT RATE BITS CHANNELS PICTURE [OFFSET] [CALLSIGN]" << std::endl;
		return 1;
	}

	const char *output_name = argv[1];
	int output_rate = std::atoi(argv[2]);
	int output_bits = std::atoi(argv[3]);
	int output_chan = std::atoi(argv[4]);
	const char *picture_name = argv[5];

	int freq_off = output_chan == 1 ? 2000 : 0;
	if (argc >= 7)
		freq_off = std::atoi(argv[6]);

	if ((output_chan == 1 && freq_off < 1200) || freq_off < 1200 - output_rate / 2 || freq_off > output_rate / 2 - 1200) {
		std::cerr << "Unsupported frequency offset." << std::endl;
		return 1;
	}

	long long int call_sign = base37_encoder("ANONYMOUS");
	if (argc >= 8)
		call_sign = base37_encoder(argv[7]);

	if (call_sign < 0 || call_sign >= 129961739795077L) {
		std::cerr << "Unsupported call sign." << std::endl;
		return 1;
	}

	typedef float value;
	typedef DSP::Complex<value> cmplx;

	DSP::ReadPNM<value> picture_file(picture_name);
	if (picture_file.width() != 320 || picture_file.height() != 240 || picture_file.mono()) {
		std::cerr << "Unsupported picture format." << std::endl;
		return 1;
	}

	DSP::WriteWAV<value> output_file(output_name, output_rate, output_bits, output_chan);
	output_file.silence(output_rate);
	switch (output_rate) {
	case 8000:
		delete new Encoder<value, cmplx, 8000>(&output_file, &picture_file, freq_off, call_sign);
		break;
	case 16000:
		delete new Encoder<value, cmplx, 16000>(&output_file, &picture_file, freq_off, call_sign);
		break;
	case 44100:
		delete new Encoder<value, cmplx, 44100>(&output_file, &picture_file, freq_off, call_sign);
		break;
	case 48000:
		delete new Encoder<value, cmplx, 48000>(&output_file, &picture_file, freq_off, call_sign);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}
	output_file.silence(output_rate);

	return 0;
}

