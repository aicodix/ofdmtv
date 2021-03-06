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
#include "sma.hh"
#include "wav.hh"
#include "pcm.hh"
#include "fft.hh"
#include "mls.hh"

template <typename value, typename cmplx, int buffer_len, int symbol_len, int guard_len>
struct SchmidlCox
{
	typedef DSP::Const<value> Const;
	static const int correlator_len = symbol_len / 2;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	DSP::SMA4<cmplx, value, correlator_len, false> cor;
	DSP::SMA4<value, value, correlator_len, false> pwr;
	DSP::SchmittTrigger<value> treshold;
	DSP::FallingEdgeTrigger falling;
	cmplx tmp0[symbol_len], tmp1[symbol_len], tmp2[symbol_len];
	cmplx kern[symbol_len];
	value phase_buf[symbol_len], timing_buf[symbol_len];
	cmplx cmplx_shift = 0;
	value timing_max = 0;
	int buf_pos = 0;
public:
	int symbol_pos = 0;
	value cfo_rad = 0;

	SchmidlCox(const cmplx *sequence) : treshold(value(0.26), value(0.29))
	{
		fwd(kern, sequence);
		for (int i = 0; i < symbol_len; ++i)
			kern[i] = conj(kern[i]) / value(symbol_len);
	}
	bool operator()(const cmplx *samples)
	{
		cmplx P = cor(samples[buffer_len-2*symbol_len-correlator_len] * conj(samples[buffer_len-2*symbol_len]));
		value R = pwr(norm(samples[buffer_len-2*symbol_len]));
		value min_R = 0.0001 * correlator_len;
		R = std::max(R, min_R);
		value timing = norm(P) / (R * R);

		bool collect = treshold(timing);
		bool process = falling(collect);

		if (!collect && !process)
			return false;

		value phase = arg(P);
		timing_max = std::max(timing_max, timing);
		if (buf_pos < symbol_len) {
			timing_buf[buf_pos] = timing;
			phase_buf[buf_pos] = phase;
			++buf_pos;
		}

		if (!process)
			return false;

		int left = 0;
		while (left < buf_pos-1 && timing_buf[left] < timing_max * 0.9)
			++left;
		int right = buf_pos-1;
		while (right > 0 && timing_buf[right] < timing_max * 0.9)
			--right;
		int middle = (left + right) / 2;
		timing_max = 0;

		value frac_cfo = phase_buf[middle] / value(correlator_len);

		DSP::Phasor<cmplx> osc;
		osc.omega(frac_cfo);
		int sample_pos = buffer_len - 3*symbol_len + middle - buf_pos;
		buf_pos = 0;
		for (int i = 0; i < symbol_len; ++i)
			tmp2[i] = samples[i+sample_pos] * osc();
		fwd(tmp0, tmp2);
		for (int i = 0; i < guard_len; ++i)
			osc();
		for (int i = 0; i < symbol_len; ++i)
			tmp2[i] = samples[i+sample_pos+symbol_len+guard_len] * osc();
		fwd(tmp1, tmp2);
		value avg_pwr(0);
		for (int i = 0; i < symbol_len; ++i)
			avg_pwr += norm(tmp1[i]);
		avg_pwr /= value(symbol_len);
		for (int i = 0; i < symbol_len; ++i)
			if (norm(tmp1[i]) <= avg_pwr * 0.001)
				tmp0[i] = 0;
			else
				tmp0[i] /= tmp1[i];
		fwd(tmp1, tmp0);
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] *= kern[i];
		bwd(tmp2, tmp1);

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

		int guard_frac = symbol_len / guard_len;
		int frac_shift = shift % guard_frac;
		value rad_shift = frac_shift * (Const::TwoPi() / guard_frac);
		cmplx cmplx_shift = DSP::polar<value>(1, rad_shift);

		if (abs(arg(tmp2[shift] * cmplx_shift)) >= Const::FourthPi())
			return false;

		symbol_pos = sample_pos;
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
	static const int mls2_poly = 0b10001000000001011;
	static const int mls3_poly = 0b10111010010000001;
	static const int buffer_len = (6 + img_height) * (symbol_len + guard_len);
	DSP::ReadPCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	DSP::BlockDC<value, value> blockdc;
	DSP::Hilbert<cmplx, 129> hilbert;
	DSP::Resampler<value, 129, 3> resample;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	SchmidlCox<value, cmplx, buffer_len, symbol_len, guard_len> correlator;
	cmplx head[symbol_len], tail[symbol_len];
	cmplx fdom[2 * symbol_len], tdom[buffer_len];
	value rgb_line[2 * 3 * img_width];
	value phase[mls1_len], index[mls1_len];
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
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		for (int i = 0; i < mls0_len; ++i)
			fdom[2*i+mls0_off] = 1 - 2 * seq0();
		return fdom;
	}
	int pos_error()
	{
		fwd(tail, tdom+symbol_pos+(symbol_len+guard_len));
		CODE::MLS seq1(mls1_poly);
		for (int i = 0; i < mls1_len; ++i)
			tail[i+mls1_off] *= (1 - 2 * seq1());
		value avg = 0;
		for (int i = 1; i < mls1_len; ++i)
			avg += phase[i] = arg(tail[i+mls1_off] / tail[i-1+mls1_off]);
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
		sfo_rad = 0;
		std::cerr << "symbol pos: " << symbol_pos << std::endl;
		std::cerr << "coarse sfo: " << 1000000 * sfo_rad / Const::TwoPi() << " ppm" << std::endl;
		std::cerr << "coarse cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;

		for (int n = 0; n < 4; ++n) {
			value diff = sfo_rad * (rate / Const::TwoPi());
			resample(tdom, buf, -diff, buffer_len);
			DSP::Phasor<cmplx> osc;
			osc.omega(-cfo_rad);
			for (int i = 0; i < buffer_len; ++i)
				tdom[i] *= osc();

			int finer_pos = symbol_pos - pos_error();
			if (correlator.symbol_pos - guard_len / 2 < finer_pos && finer_pos < correlator.symbol_pos + guard_len / 2)
				symbol_pos = finer_pos;
			std::cerr << "finer pos: " << symbol_pos << std::endl;

			int distance;
			if (n > 1) {
				fwd(head, tdom+symbol_pos-(img_height+2)*(symbol_len+guard_len));
				fwd(tail, tdom+symbol_pos+(symbol_len+guard_len));
				distance = (img_height+3)*(symbol_len+guard_len);
			} else {
				fwd(head, tdom+symbol_pos-(symbol_len+guard_len));
				fwd(tail, tdom+symbol_pos+(symbol_len+guard_len));
				distance = 2*(symbol_len+guard_len);
			}
			for (int i = 0; i < mls1_len; ++i) {
				phase[i] = arg(head[i+mls1_off] / tail[i+mls1_off]);
				index[i] = i+mls1_off;
			}

			DSP::SimpleLinearRegression<value> dirty(index, phase, mls1_len);
			value avg_diff = 0;
			for (int i = 0; i < mls1_len; ++i)
				avg_diff += dirty(index[i]) - phase[i];
			avg_diff /= value(mls1_len);
			value var_diff = 0;
			for (int i = 0; i < mls1_len; ++i)
				var_diff += (dirty(index[i]) - phase[i] - avg_diff) * (dirty(index[i]) - phase[i] - avg_diff);
			value std_dev = std::sqrt(var_diff/(mls1_len-1));
			int count = 0;
			for (int i = 0; i < mls1_len; ++i) {
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
	case 11025:
		delete new Decoder<value, cmplx, 11025>(&output_file, &input_file);
		break;
	case 16000:
		delete new Decoder<value, cmplx, 16000>(&output_file, &input_file);
		break;
	case 22050:
		delete new Decoder<value, cmplx, 22050>(&output_file, &input_file);
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

