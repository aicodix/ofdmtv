/*
OFDM TV encoder

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#include <iostream>
#include <cmath>
#include "netpbm.hh"
#include "complex.hh"
#include "utils.hh"
#include "decibel.hh"
#include "fft.hh"
#include "wav.hh"
#include "pcm.hh"
#include "mls.hh"

template <typename value, typename cmplx, int rate>
struct Encoder
{
	static const int symbol_len = (1280 * rate) / 8000;
	static const int guard_len = symbol_len / 8;
	static const int img_off = 160;
	static const int img_width = 320;
	static const int img_height = 240;
	static const int mls0_off = 192;
	static const int mls0_len = 127;
	static const int mls0_poly = 0b10001001;
	static const int mls1_off = img_off;
	static const int mls1_len = img_width;
	static const int mls1_poly = 0b1100110001;
	DSP::WritePCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	cmplx fdom[symbol_len];
	cmplx tdom[symbol_len];
	cmplx guard[guard_len];
	value rgb_line[3 * img_width];
	cmplx papr_min, papr_max;

	void rgb_to_yuv(value *yuv, const value *rgb)
	{
		value WR(0.299), WB(0.114), WG(1-WR-WB), UMAX(0.493), VMAX(0.877);
		yuv[0] = WR * rgb[0] + WG * rgb[1] + WB * rgb[2];
		yuv[1] = (UMAX/(1-WB)) * (rgb[2]-yuv[0]);
		yuv[2] = (VMAX/(1-WR)) * (rgb[0]-yuv[0]);
	}
	void rgb_to_cmplx(cmplx *out, const value *rgb)
	{
		value yuv[6];
		rgb_to_yuv(yuv, rgb);
		rgb_to_yuv(yuv+3, rgb+3);
		out[0] = cmplx(yuv[0]+(yuv[2]+yuv[5])/2, yuv[0]-(yuv[1]+yuv[4])/2);
		out[1] = cmplx((yuv[1]+yuv[4])/2-yuv[3], (yuv[2]+yuv[5])/2-yuv[3]);
	}
	void symbol()
	{
		bwd(tdom, fdom);
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] /= sqrt(value(symbol_len));
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
		for (int i = 0; i < guard_len; ++i) {
			cmplx power(guard[i].real() * guard[i].real(), guard[i].imag() * guard[i].imag());
			peak = cmplx(std::max(peak.real(), power.real()), std::max(peak.imag(), power.imag()));
			mean += power;
		}
		if (mean.real() > 0 && mean.imag() > 0) {
			cmplx papr(peak.real() / mean.real(), peak.imag() / mean.imag());
			papr *= value(symbol_len + guard_len);
			papr_min = cmplx(std::min(papr_min.real(), papr.real()), std::min(papr_min.imag(), papr.imag()));
			papr_max = cmplx(std::max(papr_max.real(), papr.real()), std::max(papr_max.imag(), papr.imag()));
		}
		pcm->write(reinterpret_cast<value *>(guard), guard_len, 2);
		pcm->write(reinterpret_cast<value *>(tdom), symbol_len, 2);
		for (int i = 0; i < guard_len; ++i)
			guard[i] = tdom[i];
	}
	Encoder(DSP::WritePCM<value> *pcm, DSP::ReadPEL<value> *pel) : pcm(pcm)
	{
		papr_min = cmplx(1000, 1000), papr_max = cmplx(-1000, -1000);
		CODE::MLS seq0(mls0_poly), seq1(mls1_poly);
		value mls1_fac = sqrt(value(symbol_len) / value(4 * mls1_len));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		for (int i = mls1_off; i < mls1_off + mls1_len; ++i)
			fdom[i] = mls1_fac * (1 - 2 * seq1());
		symbol();
		value mls0_fac = sqrt(value(symbol_len) / value(4 * mls0_len));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		for (int i = 0; i < mls0_len; ++i)
			fdom[2*i+mls0_off] = - mls0_fac * (1 - 2 * seq0());
		seq1.reset();
		for (int i = mls1_off; i < mls1_off + mls1_len; ++i)
			fdom[i] *= 1 - 2 * seq1();
		symbol();
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		seq1.reset();
		for (int i = mls1_off; i < mls1_off + mls1_len; ++i)
			fdom[i] = mls1_fac * (1 - 2 * seq1());
		symbol();
		value img_fac = sqrt(value(symbol_len) / value(4 * img_width));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		for (int j = 0; j < img_height; ++j) {
			pel->read(rgb_line, img_width);
			for (int i = 0; i < img_width; i += 2)
				rgb_to_cmplx(fdom+i+img_off, rgb_line+3*i);
			seq1.reset();
			for (int i = 0; i < img_width; ++i)
				fdom[i+img_off] *= img_fac * (1 - 2 * seq1());
			symbol();
		}
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		seq1.reset();
		for (int i = mls1_off; i < mls1_off + mls1_len; ++i)
			fdom[i] = mls1_fac * (1 - 2 * seq1());
		symbol();
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		seq0.reset();
		for (int i = 0; i < mls0_len; ++i)
			fdom[2*i+mls0_off] = mls0_fac * (1 - 2 * seq0());
		seq1.reset();
		for (int i = mls1_off; i < mls1_off + mls1_len; ++i)
			fdom[i] *= 1 - 2 * seq1();
		symbol();
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		seq1.reset();
		for (int i = mls1_off; i < mls1_off + mls1_len; ++i)
			fdom[i] = mls1_fac * (1 - 2 * seq1());
		symbol();
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		symbol();
		std::cerr << "real PAPR: " << DSP::decibel(papr_min.real()) << " .. " << DSP::decibel(papr_max.real()) << " dB" << std::endl;
		if (pcm->channels() == 2)
			std::cerr << "imag PAPR: " << DSP::decibel(papr_min.imag()) << " .. " << DSP::decibel(papr_max.imag()) << " dB" << std::endl;
	}
};

int main(int argc, char **argv)
{
	if (argc != 6) {
		std::cerr << "usage: " << argv[0] << " OUTPUT RATE BITS CHANNELS PICTURE" << std::endl;
		return 1;
	}

	const char *output_name = argv[1];
	int output_rate = std::atoi(argv[2]);
	int output_bits = std::atoi(argv[3]);
	int output_chan = std::atoi(argv[4]);
	const char *picture_name = argv[5];

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
		delete new Encoder<value, cmplx, 8000>(&output_file, &picture_file);
		break;
	case 11025:
		delete new Encoder<value, cmplx, 11025>(&output_file, &picture_file);
		break;
	case 16000:
		delete new Encoder<value, cmplx, 16000>(&output_file, &picture_file);
		break;
	case 22050:
		delete new Encoder<value, cmplx, 22050>(&output_file, &picture_file);
		break;
	case 44100:
		delete new Encoder<value, cmplx, 44100>(&output_file, &picture_file);
		break;
	case 48000:
		delete new Encoder<value, cmplx, 48000>(&output_file, &picture_file);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}
	output_file.silence(output_rate);

	return 0;
}

