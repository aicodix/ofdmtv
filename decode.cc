/*
OFDM TV decoder

Copyright 2020 Ahmet Inan <inan@aicodix.de>
*/

#include <iostream>
#include <cmath>
namespace DSP { using std::abs; using std::min; using std::cos; using std::sin; }
#include "bip_buffer.hh"
#include "regression.hh"
#include "resampler.hh"
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

template <typename value, typename cmplx, int buffer_len, int symbol_len, int guard_len>
struct SchmidlCox
{
	typedef DSP::Const<value> Const;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	DSP::SMA4<cmplx, value, symbol_len, false> cor;
	DSP::SMA4<value, value, symbol_len, false> pwr;
	DSP::SMA4<value, value, guard_len, false> match;
	DSP::Delay<value, guard_len/2> delay;
	DSP::SchmittTrigger<value> threshold;
	DSP::FallingEdgeTrigger falling;
	cmplx tmp0[symbol_len], tmp1[symbol_len], tmp2[symbol_len];
	cmplx seq[symbol_len], kern[symbol_len];
	cmplx cmplx_shift = 0;
	value timing_max = 0;
	value phase_max = 0;
	int index_max = 0;
public:
	int symbol_pos = 0;
	value cfo_rad = 0;
	value frac_cfo = 0;

	SchmidlCox(const cmplx *sequence) : threshold(value(0.17*guard_len), value(0.19*guard_len))
	{
		for (int i = 0; i < symbol_len; ++i)
			seq[i] = sequence[i];
		fwd(kern, sequence);
		for (int i = 0; i < symbol_len; ++i)
			kern[i] = conj(kern[i]) / value(symbol_len);
	}
	bool operator()(const cmplx *samples)
	{
		cmplx P = cor(samples[buffer_len-5*symbol_len] * conj(samples[buffer_len-4*symbol_len]));
		value R = pwr(norm(samples[buffer_len-4*symbol_len]));
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
			index_max = guard_len / 2;
		} else if (index_max < 3*symbol_len) {
			++index_max;
		}

		if (!process)
			return false;

		frac_cfo = phase_max / value(symbol_len);

		DSP::Phasor<cmplx> osc;
		osc.omega(frac_cfo);
		symbol_pos = buffer_len - 6*symbol_len - index_max;
		index_max = 0;
		timing_max = 0;
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = samples[i+symbol_pos+symbol_len] * osc();
		fwd(tmp0, tmp1);
		value avg_pwr(0);
		for (int i = 0; i < symbol_len; ++i)
			avg_pwr += norm(tmp0[i]);
		avg_pwr /= value(symbol_len);
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = 0;
		for (int i = 0; i < symbol_len; ++i)
			if (norm(tmp0[i]) * 1000 > avg_pwr && norm(tmp0[i]) < avg_pwr * 1000 &&
				std::min(norm(tmp0[i]), norm(tmp0[(i-1+symbol_len)%symbol_len])) * 2 >
				std::max(norm(tmp0[i]), norm(tmp0[(i-1+symbol_len)%symbol_len])))
					tmp1[i] = tmp0[i] / tmp0[(i-1+symbol_len)%symbol_len];
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
		if (peak <= next * 8)
			return false;

		if (abs(arg(tmp2[shift])) >= Const::FourthPi())
			return false;

		cfo_rad = shift * (Const::TwoPi() / symbol_len) - frac_cfo;
		if (cfo_rad >= Const::Pi())
			cfo_rad -= Const::TwoPi();
		return true;
	}
};

template <typename value, typename cmplx, int rate>
struct Decoder
{
	typedef DSP::Const<value> Const;
	static const int symbol_len = (1280 * rate) / 8000;
	static const int filter_len = (((21 * rate) / 8000) & ~3) | 1;
	static const int guard_len = symbol_len / 8;
	static const int img_off = 160;
	static const int img_width = 320;
	static const int img_height = 240;
	static const int mls0_off = 194;
	static const int mls0_len = 127;
	static const int mls0_poly = 0b10001001;
	static const int mls1_off = img_off;
	static const int mls1_len = img_width;
	static const int mls1_poly = 0b1100110001;
	static const int mls2_poly = 0b10001000000001011;
	static const int mls3_poly = 0b10111010010000001;
	static const int buffer_len = (7 + img_height) * (symbol_len + guard_len);
	DSP::ReadPCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	DSP::BlockDC<value, value> blockdc;
	DSP::Hilbert<cmplx, filter_len> hilbert;
	DSP::Resampler<value, filter_len, 3> resample;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	SchmidlCox<value, cmplx, buffer_len, symbol_len/2, guard_len> correlator;
	cmplx head[2 * symbol_len], tail[2 * symbol_len];
	cmplx fdom[2 * symbol_len], tdom[buffer_len], resam[buffer_len];
	value rgb_line[2 * 3 * img_width];
	value phase[2 * mls1_len], index[2 * mls1_len];
	value cfo_rad, sfo_rad;
	int symbol_pos;

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
			fdom[i+mls0_off/2] = 1 - 2 * seq0();
		return fdom;
	}
	int pos_error(const cmplx *symbol)
	{
		value avg = 0;
		for (int i = 1; i < mls1_len; ++i)
			if ((symbol[i+mls1_off] / symbol[i-1+mls1_off]).real() >= 0)
				avg += phase[i] = arg(symbol[i+mls1_off] / symbol[i-1+mls1_off]);
			else
				avg += phase[i] = arg(- symbol[i+mls1_off] / symbol[i-1+mls1_off]);
		avg /= value(mls1_len-1);
		value var = 0;
		for (int i = 1; i < mls1_len; ++i)
			var += (phase[i] - avg) * (phase[i] - avg);
		value std_dev = std::sqrt(var/(mls1_len-2));
		int count = 0;
		value sum = 0;
		for (int i = 1; i < mls1_len; ++i) {
			if (2 * std::abs(phase[i] - avg) <= std_dev) {
				sum += phase[i];
				++count;
			}
		}
		return std::nearbyint(sum * symbol_len / (count * Const::TwoPi()));
	}
	int displacement(const cmplx *sym0, const cmplx *sym1)
	{
		fwd(head, sym0);
		fwd(tail, sym1);
		for (int i = 0; i < symbol_len; ++i)
			head[i] *= conj(tail[i]);
		bwd(tail, head);
		int idx = 0;
		for (int i = 0; i < symbol_len; ++i)
			if (norm(tail[i]) > norm(tail[idx]))
				idx = i;
		if (idx > symbol_len / 2)
			idx -= symbol_len;
		return -idx;
	}
	Decoder(DSP::WritePEL<value> *pel, DSP::ReadPCM<value> *pcm) : pcm(pcm), resample(rate, (rate * 19) / 40, 2), correlator(mls0_seq())
	{
		bool real = pcm->channels() == 1;
		blockdc.samples(2*(symbol_len+guard_len));
		const cmplx *buf;
		do  {
			if (!pcm->good())
				return;
			cmplx tmp;
			pcm->read(reinterpret_cast<value *>(&tmp), 1);
			if (real)
				tmp = hilbert(blockdc(tmp.real()));
			buf = input_hist(tmp);
		} while (!correlator(buf));

		symbol_pos = correlator.symbol_pos;
		cfo_rad = correlator.cfo_rad;
		int dis = displacement(buf+symbol_pos-(img_height+2)*(symbol_len+guard_len), buf+symbol_pos+(symbol_len+guard_len));
		sfo_rad = (dis * Const::TwoPi()) / ((img_height+3)*(symbol_len+guard_len));
		std::cerr << "symbol pos: " << symbol_pos << std::endl;
		std::cerr << "coarse sfo: " << 1000000 * sfo_rad / Const::TwoPi() << " ppm" << std::endl;
		std::cerr << "coarse cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;
		if (dis) {
			value diff = sfo_rad * (rate / Const::TwoPi());
			resample(resam, buf, -diff, buffer_len);
			symbol_pos = std::nearbyint(correlator.symbol_pos * (1 - sfo_rad / Const::TwoPi()));
			std::cerr << "resam pos: " << symbol_pos << std::endl;
			cmplx P;
			for (int i = 0; i < symbol_len/2; ++i)
				P += resam[i+symbol_pos] * conj(resam[i+symbol_pos+symbol_len/2]);
			cfo_rad = correlator.cfo_rad + correlator.frac_cfo - arg(P) / value(symbol_len/2);
			std::cerr << "resam cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;
		} else {
			for (int i = 0; i < buffer_len; ++i)
				resam[i] = buf[i];
		}

		for (int n = 0; n < 4; ++n) {
			DSP::Phasor<cmplx> osc;
			osc.omega(-cfo_rad);
			for (int i = 0; i < buffer_len; ++i)
				tdom[i] = resam[i] * osc();
			fwd(tail, tdom+symbol_pos+(symbol_len+guard_len));
			if (1) {
				int finer_pos = symbol_pos - pos_error(tail);
				if (finer_pos != symbol_pos && correlator.symbol_pos - guard_len / 2 < finer_pos && finer_pos < correlator.symbol_pos + guard_len / 2) {
					symbol_pos = finer_pos;
					fwd(tail, tdom+symbol_pos+(symbol_len+guard_len));
					std::cerr << "finer pos: " << symbol_pos << std::endl;
				}
			}
			int distance = 2*(symbol_len+guard_len);
			if (n > 1)
				distance = (img_height+3)*(symbol_len+guard_len);
			fwd(head, tdom+symbol_pos+(symbol_len+guard_len)-distance);
			int length = 0;
			for (int i = 0; i < mls1_len; ++i) {
				phase[length] = arg(head[i+mls1_off] / tail[i+mls1_off]);
				index[length] = i+mls1_off;
				++length;
			}
			if (n > 1) {
				int head_pos = symbol_pos-(img_height+4)*(symbol_len+guard_len);
				fwd(head+symbol_len, tdom+head_pos);
				fwd(tail+symbol_len, tdom+head_pos+distance);
				for (int i = 0; i < mls1_len; ++i) {
					phase[length] = arg(head[i+mls1_off+symbol_len] / tail[i+mls1_off+symbol_len]);
					index[length] = i+mls1_off;
					++length;
				}
			}

			DSP::SimpleLinearRegression<value> dirty(index, phase, length);
			value avg_diff = 0;
			for (int i = 0; i < length; ++i)
				avg_diff += dirty(index[i]) - phase[i];
			avg_diff /= value(length);
			value var_diff = 0;
			for (int i = 0; i < length; ++i)
				var_diff += (dirty(index[i]) - phase[i] - avg_diff) * (dirty(index[i]) - phase[i] - avg_diff);
			value std_dev = std::sqrt(var_diff/(length-1));
			int count = 0;
			for (int i = 0; i < length; ++i) {
				if (2 * std::abs(dirty(index[i])-phase[i]) < std_dev) {
					index[count] = index[i];
					phase[count] = phase[i];
					++count;
				}
			}
			DSP::SimpleLinearRegression<value> sfo_cfo(index, phase, count);
			sfo_rad += sfo_cfo.slope() * symbol_len / value(distance);
			value sfo_max = 0.0005 * Const::TwoPi();
			sfo_rad = std::min(std::max(sfo_rad, -sfo_max), sfo_max);
			cfo_rad -= sfo_cfo.yint() / value(distance);
			std::cerr << "finer sfo: " << 1000000 * sfo_rad / Const::TwoPi() << " ppm" << std::endl;
			std::cerr << "finer cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz" << std::endl;
		}
		CODE::MLS seq1(mls1_poly), seq2(mls2_poly), seq3(mls3_poly);
		for (int i = 0; i < mls1_len; ++i)
			head[i+mls1_off] *= (1 - 2 * seq1());
		seq1.reset();
		for (int i = 0; i < mls1_len; ++i)
			tail[i+mls1_off] *= (1 - 2 * seq1());
		for (int j = 0; j < img_height; j += 2) {
			for (int k = 0; k < 2; ++k) {
				fwd(fdom+symbol_len*k, tdom+symbol_pos+(j+k-img_height-1)*(symbol_len+guard_len));
				value x = value(j+k+1) / value(img_height+3);
				for (int i = 0; i < img_width; ++i)
					fdom[i+img_off+symbol_len*k] /= DSP::lerp(x, head[i+img_off], tail[i+img_off]);
				for (int i = 0; i < img_width; ++i)
					fdom[i+img_off+symbol_len*k] = cmplx(
						fdom[i+img_off+symbol_len*k].real() * (1 - 2 * seq2()),
						fdom[i+img_off+symbol_len*k].imag() * (1 - 2 * seq3()));
			}
			for (int i = 0; i < img_width; i += 2)
				cmplx_to_rgb(rgb_line+3*i, rgb_line+3*(i+img_width), fdom+i+img_off, fdom+i+img_off+symbol_len);
			pel->write(rgb_line, 2 * img_width);
		}
	}
};

int main(int argc, char **argv)
{
	if (argc != 3) {
		std::cerr << "usage: " << argv[0] << " OUTPUT INPUT" << std::endl;
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

	DSP::WritePNM<value> output_file(output_name, 320, 240);

	switch (input_file.rate()) {
	case 8000:
		delete new Decoder<value, cmplx, 8000>(&output_file, &input_file);
		break;
	case 16000:
		delete new Decoder<value, cmplx, 16000>(&output_file, &input_file);
		break;
	case 44100:
		delete new Decoder<value, cmplx, 44100>(&output_file, &input_file);
		break;
	case 48000:
		delete new Decoder<value, cmplx, 48000>(&output_file, &input_file);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}

	return 0;
}

