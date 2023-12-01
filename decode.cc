/*
OFDM TV decoder

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#include <iostream>
#include <cassert>
#include <cmath>
namespace DSP { using std::abs; using std::min; using std::cos; using std::sin; }
#include "bip_buffer.hh"
#include "theil_sen.hh"
#include "trigger.hh"
#include "complex.hh"
#include "blockdc.hh"
#include "hilbert.hh"
#include "phasor.hh"
#include "netpbm.hh"
#include "delay.hh"
#include "sma.hh"
#include "wav.hh"
#include "pcm.hh"
#include "fft.hh"
#include "mls.hh"
#include "crc.hh"
#include "osd.hh"

template <typename value, typename cmplx, int search_pos, int symbol_len, int guard_len>
struct SchmidlCox
{
	typedef DSP::Const<value> Const;
	static const int match_len = guard_len | 1;
	static const int match_del = (match_len - 1) / 2;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	DSP::SMA4<cmplx, value, symbol_len, false> cor;
	DSP::SMA4<value, value, 2*symbol_len, false> pwr;
	DSP::SMA4<value, value, match_len, false> match;
	DSP::Delay<value, match_del> delay;
	DSP::SchmittTrigger<value> threshold;
	DSP::FallingEdgeTrigger falling;
	cmplx tmp0[symbol_len], tmp1[symbol_len], tmp2[symbol_len];
	cmplx seq[symbol_len], kern[symbol_len];
	cmplx cmplx_shift = 0;
	value timing_max = 0;
	value phase_max = 0;
	int index_max = 0;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
public:
	int symbol_pos = 0;
	value cfo_rad = 0;
	value frac_cfo = 0;

	SchmidlCox(const cmplx *sequence) : threshold(value(0.17*match_len), value(0.19*match_len))
	{
		for (int i = 0; i < symbol_len; ++i)
			seq[i] = sequence[i];
		fwd(kern, sequence);
		for (int i = 0; i < symbol_len; ++i)
			kern[i] = conj(kern[i]) / value(symbol_len);
	}
	bool operator()(const cmplx *samples)
	{
		cmplx P = cor(samples[search_pos+symbol_len] * conj(samples[search_pos+2*symbol_len]));
		value R = value(0.5) * pwr(norm(samples[search_pos+2*symbol_len]));
		value min_R = 0.0001 * symbol_len;
		R = std::max(R, min_R);
		value timing = match(norm(P) / (R * R));
		value phase = delay(arg(P));

		bool collect = threshold(timing);
		bool process = falling(collect);

		if (!collect && !process)
			return false;

		if (timing_max < timing) {
			timing_max = timing;
			phase_max = phase;
			index_max = match_del;
		} else if (index_max < symbol_len + guard_len + match_del) {
			++index_max;
		}

		if (!process)
			return false;

		frac_cfo = phase_max / value(symbol_len);

		DSP::Phasor<cmplx> osc;
		osc.omega(frac_cfo);
		symbol_pos = search_pos - index_max;
		index_max = 0;
		timing_max = 0;
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = samples[i+symbol_pos+symbol_len] * osc();
		fwd(tmp0, tmp1);
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = 0;
		for (int i = 0; i < symbol_len; ++i)
			if (norm(tmp0[bin(i-1)]) > 0 &&
				std::min(norm(tmp0[i]), norm(tmp0[bin(i-1)])) * 2 >
				std::max(norm(tmp0[i]), norm(tmp0[bin(i-1)])))
					tmp1[i] = tmp0[i] / tmp0[bin(i-1)];
		fwd(tmp0, tmp1);
		for (int i = 0; i < symbol_len; ++i)
			tmp0[i] *= kern[i];
		bwd(tmp2, tmp0);

		int shift = 0;
		value peak = 0;
		value next = 0;
		for (int i = 0; i < symbol_len; ++i) {
			value power = norm(tmp2[i]);
			if (power > peak) {
				next = peak;
				peak = power;
				shift = i;
			} else if (power > next) {
				next = power;
			}
		}
		if (peak <= next * 4)
			return false;

		int pos_err = std::nearbyint(arg(tmp2[shift]) * symbol_len / Const::TwoPi());
		if (abs(pos_err) > guard_len / 2)
			return false;
		symbol_pos -= pos_err;

		cfo_rad = shift * (Const::TwoPi() / symbol_len) - frac_cfo;
		if (cfo_rad >= Const::Pi())
			cfo_rad -= Const::TwoPi();
		return true;
	}
};

void base37_decoder(char *str, long long int val, int len)
{
	for (int i = len-1; i >= 0; --i, val /= 37)
		str[i] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[val%37];
}

template <typename value, typename cmplx, int rate>
struct Decoder
{
	typedef DSP::Const<value> Const;
	static const int symbol_len = (1280 * rate) / 8000;
	static const int filter_len = (((21 * rate) / 8000) & ~3) | 1;
	static const int guard_len = symbol_len / 8;
	static const int img_width = 320;
	static const int img_height = 240;
	static const int teeth_count = 16;
	static const int teeth_dist = img_width / teeth_count;
	static const int teeth_off = teeth_dist / 2;
	static const int mls0_len = 127;
	static const int mls0_off = -mls0_len + 1;
	static const int mls0_poly = 0b10001001;
	static const int mls1_len = img_width + teeth_count;
	static const int mls1_off = -mls1_len / 2;
	static const int mls1_poly = 0b1100110001;
	static const int mls2_poly = 0b10001000000001011;
	static const int mls3_poly = 0b10111010010000001;
	static const int mls4_len = 255;
	static const int mls4_off = -mls4_len / 2;
	static const int mls4_poly = 0b100101011;
	static const int buffer_len = 6 * (symbol_len + guard_len);
	static const int search_pos = buffer_len - 4 * (symbol_len + guard_len);
	DSP::ReadPCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	DSP::BlockDC<value, value> blockdc;
	DSP::Hilbert<cmplx, filter_len> hilbert;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	DSP::TheilSenEstimator<value, mls1_len> tse;
	SchmidlCox<value, cmplx, search_pos, symbol_len/2, guard_len> correlator;
	CODE::CRC<uint16_t> crc;
	CODE::OrderedStatisticsDecoder<255, 71, 4> osddec;
	int8_t genmat[255*71];
	cmplx chan[mls1_len];
	cmplx fdom[2 * symbol_len], tdom[symbol_len];
	value rgb_line[2 * 3 * img_width];
	value phase[mls1_len], index[mls1_len];
	value cfo_rad, sfo_rad;
	int symbol_pos;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
	static value nrz(bool bit)
	{
		return 1 - 2 * bit;
	}
	void yuv_to_rgb(value *rgb, const value *yuv)
	{
		value WR(0.299), WB(0.114), WG(1-WR-WB), UMAX(0.493), VMAX(0.877);
		rgb[0] = yuv[0] + ((1-WR)/VMAX) * yuv[2];
		rgb[1] = yuv[0] - (WB*(1-WB)/(UMAX*WG)) * yuv[1] - (WR*(1-WR)/(VMAX*WG)) * yuv[2];
		rgb[2] = yuv[0] + ((1-WB)/UMAX) * yuv[1];
	}
	void cmplx_to_rgb(value *rgb0, value *rgb1, cmplx *inp0, cmplx *inp1)
	{
		value upv[2] = {
			inp0[0].real()-inp0[0].imag(),
			inp1[1].real()-inp1[1].imag()
		};
		value umv[2] = {
			inp1[0].real()-inp1[0].imag(),
			inp0[1].real()-inp0[1].imag()
		};
		value yuv0[6] = {
			(inp0[0].real()+inp0[0].imag()+umv[0])/2, (upv[0]+umv[0])/2, (upv[0]-umv[0])/2,
			(upv[1]-inp0[1].real()-inp0[1].imag())/2, (upv[1]+umv[1])/2, (upv[1]-umv[1])/2
		};
		value yuv1[6] = {
			(upv[0]-inp1[0].real()-inp1[0].imag())/2, (upv[0]+umv[0])/2, (upv[0]-umv[0])/2,
			(inp1[1].real()+inp1[1].imag()+umv[1])/2, (upv[1]+umv[1])/2, (upv[1]-umv[1])/2
		};
		yuv_to_rgb(rgb0, yuv0);
		yuv_to_rgb(rgb1, yuv1);
		yuv_to_rgb(rgb0+3, yuv0+3);
		yuv_to_rgb(rgb1+3, yuv1+3);
	}
	const cmplx *mls0_seq()
	{
		CODE::MLS seq0(mls0_poly);
		for (int i = 0; i < symbol_len/2; ++i)
			fdom[i] = 0;
		for (int i = 0; i < mls0_len; ++i)
			fdom[(i+mls0_off/2+symbol_len/2)%(symbol_len/2)] = nrz(seq0());
		return fdom;
	}
	const cmplx *next_sample()
	{
		cmplx tmp;
		pcm->read(reinterpret_cast<value *>(&tmp), 1);
		if (pcm->channels() == 1)
			tmp = hilbert(blockdc(tmp.real()));
		return input_hist(tmp);
	}
	Decoder(DSP::WritePEL<value> *pel, DSP::ReadPCM<value> *pcm, int skip_count) :
		pcm(pcm), correlator(mls0_seq()), crc(0xA8F4)
	{
		CODE::BoseChaudhuriHocquenghemGenerator<255, 71>::matrix(genmat, true, {
			0b100011101, 0b101110111, 0b111110011, 0b101101001,
			0b110111101, 0b111100111, 0b100101011, 0b111010111,
			0b000010011, 0b101100101, 0b110001011, 0b101100011,
			0b100011011, 0b100111111, 0b110001101, 0b100101101,
			0b101011111, 0b111111001, 0b111000011, 0b100111001,
			0b110101001, 0b000011111, 0b110000111, 0b110110001});

		DSP::Phasor<cmplx> osc;
		blockdc.samples(2*(symbol_len+guard_len));
		const cmplx *buf;
		bool okay;
		do {
			okay = false;
			do {
				if (!pcm->good())
					return;
				buf = next_sample();
			} while (!correlator(buf));

			symbol_pos = correlator.symbol_pos;
			cfo_rad = correlator.cfo_rad;
			std::cerr << "symbol pos: " << symbol_pos << std::endl;
			std::cerr << "coarse cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;

			osc.omega(-cfo_rad);
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i+symbol_pos+(symbol_len+guard_len)] * osc();
			fwd(fdom, tdom);
			CODE::MLS seq4(mls4_poly);
			for (int i = 0; i < mls4_len; ++i)
				fdom[bin(i+mls4_off)] *= nrz(seq4());
			int8_t soft[mls4_len];
			uint8_t data[(mls4_len+7)/8];
			for (int i = 0; i < mls4_len; ++i)
				soft[i] = std::min<value>(std::max<value>(
					std::nearbyint(127 * (fdom[bin(i+mls4_off)] /
					fdom[bin(i-1+mls4_off)]).real()), -128), 127);
			bool unique = osddec(data, soft, genmat);
			if (!unique) {
				std::cerr << "OSD error." << std::endl;
				return;
			}
			uint64_t md = 0;
			for (int i = 0; i < 55; ++i)
				md |= (uint64_t)CODE::get_be_bit(data, i) << i;
			uint16_t cs = 0;
			for (int i = 0; i < 16; ++i)
				cs |= (uint16_t)CODE::get_be_bit(data, i+55) << i;
			crc.reset();
			if (crc(md<<9) != cs) {
				std::cerr << "CRC error." << std::endl;
				return;
			}
			if ((md&255) != 1) {
				std::cerr << "operation mode unsupported." << std::endl;
				return;
			}
			if ((md>>8) == 0 || (md>>8) >= 129961739795077L) {
				std::cerr << "call sign unsupported." << std::endl;
				return;
			}
			char call_sign[10];
			base37_decoder(call_sign, md>>8, 9);
			call_sign[9] = 0;
			std::cerr << "call sign: " << call_sign << std::endl;
			okay = true;
		} while (skip_count--);

		if (!okay)
			return;

		for (int i = 0; i < symbol_pos+(symbol_len+guard_len); ++i)
			buf = next_sample();

		CODE::MLS seq1(mls1_poly), seq2(mls2_poly), seq3(mls3_poly);
		for (int j = 0; j < img_height; j += 2) {
			if (j%8==0) {
				for (int i = 0; i < symbol_len; ++i)
					tdom[i] = buf[i+(symbol_len+guard_len)] * osc();
				fwd(fdom, tdom);
				seq1.reset();
				for (int i = 0; i < mls1_len; ++i)
					chan[i] = nrz(seq1()) * fdom[bin(i+mls1_off)];
				int count = mls1_len - 1;
				for (int i = 0; i < count; ++i)
					phase[i] = arg(chan[i+1] / chan[i]);
				std::nth_element(phase, phase+count/2, phase+count);
				value angle_err = phase[count/2];
				int pos_err = std::nearbyint((symbol_len * angle_err) / Const::TwoPi());
				for (int i = 0; i < (symbol_len+guard_len) - pos_err; ++i)
					buf = next_sample();
				for (int i = 0; i < symbol_len; ++i)
					tdom[i] = buf[i] * osc();
				for (int i = 0; i < guard_len; ++i)
					osc();
				fwd(fdom, tdom);
				seq1.reset();
				for (int i = 0; i < mls1_len; ++i)
					chan[i] = nrz(seq1()) * fdom[bin(i+mls1_off)];
			}
			for (int k = 0; k < 2; ++k) {
				for (int i = 0; i < (symbol_len+guard_len); ++i)
					buf = next_sample();
				for (int i = 0; i < symbol_len; ++i)
					tdom[i] = buf[i] * osc();
				for (int i = 0; i < guard_len; ++i)
					osc();
				fwd(fdom+symbol_len*k, tdom);
				for (int i = teeth_off, l = 0; i < mls1_len; i += teeth_dist+1, ++l) {
					index[l] = i+mls1_off;
					phase[l] = arg(fdom[bin(i+mls1_off)+symbol_len*k] / chan[i]);
				}
				tse.compute(index, phase, teeth_count);
				for (int i = 0; i < mls1_len; ++i)
					chan[i] *= DSP::polar<value>(1, tse(i+mls1_off));
				for (int i = teeth_off; i < mls1_len; i += teeth_dist+1)
					chan[i] = fdom[bin(i+mls1_off)+symbol_len*k];
				for (int i = 0, l = 0; i < img_width; ++i, ++l) {
					if ((i + teeth_off) % teeth_dist == 0)
						++l;
					fdom[bin(l+mls1_off)+symbol_len*k] /= chan[l];
					fdom[bin(l+mls1_off)+symbol_len*k] = cmplx(
						fdom[bin(l+mls1_off)+symbol_len*k].real() * nrz(seq2()),
						fdom[bin(l+mls1_off)+symbol_len*k].imag() * nrz(seq3()));
				}
			}
			for (int i = 0, l = 0; i < img_width; i += 2, l += 2) {
				if ((i + teeth_off) % teeth_dist == 0)
					++l;
				cmplx_to_rgb(rgb_line+3*i, rgb_line+3*(i+img_width), fdom+bin(l+mls1_off), fdom+bin(l+mls1_off)+symbol_len);
			}
			pel->write(rgb_line, 2 * img_width);
		}
	}
};

int main(int argc, char **argv)
{
	if (argc < 3 || argc > 4) {
		std::cerr << "usage: " << argv[0] << " OUTPUT INPUT [SKIP]" << std::endl;
		return 1;
	}

	typedef float value;
	typedef DSP::Complex<value> cmplx;

	const char *output_name = argv[1];
	const char *input_name = argv[2];

	DSP::ReadWAV<value> input_file(input_name);

	if (input_file.channels() < 1 || input_file.channels() > 2) {
		std::cerr << "Only real or analytic signal (one or two channels) supported." << std::endl;
		return 1;
	}

	int skip_count = 0;
	if (argc > 3)
		skip_count = std::atoi(argv[3]);

	DSP::WritePNM<value> output_file(output_name, 320, 240);

	switch (input_file.rate()) {
	case 8000:
		delete new Decoder<value, cmplx, 8000>(&output_file, &input_file, skip_count);
		break;
	case 16000:
		delete new Decoder<value, cmplx, 16000>(&output_file, &input_file, skip_count);
		break;
	case 44100:
		delete new Decoder<value, cmplx, 44100>(&output_file, &input_file, skip_count);
		break;
	case 48000:
		delete new Decoder<value, cmplx, 48000>(&output_file, &input_file, skip_count);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}

	return 0;
}

