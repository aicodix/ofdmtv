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
#include "hilbert.hh"
#include "phasor.hh"
#include "netpbm.hh"
#include "sma.hh"
#include "wav.hh"
#include "pcm.hh"
#include "fft.hh"
#include "mls.hh"

template <typename value, typename cmplx, int buffer_len, int symbol_len, int guard_len, int seq_len, int seq_off>
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
	cmplx seq[symbol_len], kern[symbol_len];
	value phase_buf[symbol_len], timing_buf[symbol_len];
	cmplx cmplx_shift = 0;
	value timing_max = 0;
	int buf_pos = 0;
public:
	int symbol_pos = 0;
	value cfo_rad = 0;
	value sfo_rad = 0;

	SchmidlCox(const cmplx *sequence) : treshold(value(0.26), value(0.29))
	{
		for (int i = 0; i < symbol_len; ++i)
			seq[i] = 0;
		for (int i = 0; i < seq_len; ++i)
			seq[2*i] = sequence[i];
		fwd(kern, seq);
		for (int i = 0; i < symbol_len; ++i)
			kern[i] = conj(kern[i]) / value(symbol_len);
	}
	bool operator()(const cmplx *samples)
	{
		cmplx P = cor(samples[buffer_len-2*symbol_len-correlator_len] * conj(samples[buffer_len-2*symbol_len]));
		value R = pwr(norm(samples[buffer_len-2*symbol_len]));
		value min_R = 0.01 * correlator_len;
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
		if (peak <= next * 8 || shift >= symbol_len - 2 * seq_len)
			return false;

		int guard_frac = symbol_len / guard_len;
		int frac_shift = shift % guard_frac;
		value rad_shift = frac_shift * (Const::TwoPi() / guard_frac);
		cmplx cmplx_shift = DSP::polar<value>(1, rad_shift);

		cmplx tmp = tmp2[shift] * cmplx_shift;
		if (tmp.real() <= 10 * std::abs(tmp.imag()))
			return false;

		symbol_pos = sample_pos;
		sfo_rad = 0;
		cfo_rad = (shift-seq_off) * (Const::TwoPi() / symbol_len) - frac_cfo;

		avg_pwr = 0;
		for (int i = 0; i < seq_len; ++i)
			avg_pwr += norm(tmp0[2*i+shift]);
		avg_pwr /= value(seq_len);
		int count = 0;
		for (int i = 0; i < seq_len; ++i) {
			value power = norm(tmp0[2*i+shift]);
			if (2 * power > avg_pwr && power < 2 * avg_pwr) {
				phase_buf[count] = arg(tmp0[2*i+shift] * seq[2*i] * cmplx_shift);
				timing_buf[count] = 2*i+shift;
				++count;
			}
		}

		DSP::SimpleLinearRegression<value> sfo_cfo(timing_buf, phase_buf, count);
		sfo_rad += sfo_cfo.slope() * symbol_len / value(symbol_len+guard_len);
		cfo_rad -= sfo_cfo.yint() / value(symbol_len+guard_len);

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
	static const int buffer_len = (6 + img_height) * (symbol_len + guard_len);
	DSP::ReadPCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	DSP::Hilbert<cmplx, 129> hilbert;
	DSP::Resampler<value, 129, 3> resample;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	SchmidlCox<value, cmplx, buffer_len, symbol_len, guard_len, mls0_len, mls0_off> correlator;
	cmplx mls0_sig[mls0_len], mls1_sig[mls1_len];
	cmplx head[symbol_len], tail[symbol_len];
	cmplx fdom[symbol_len], tdom[buffer_len];
	value rgb_line[3 * img_width];
	value phase[mls1_len], index[mls1_len];
	value cfo_rad, sfo_rad;
	int symbol_pos;

	void unwrap(value *phase, int len)
	{
		value sum0 = 0, sum1 = 0;
		for (int i = 1; i < len; ++i) {
			if (phase[i] >= Const::HalfPi() && phase[i-1] <= -Const::HalfPi())
				sum1 -= Const::TwoPi();
			if (phase[i] <= -Const::HalfPi() && phase[i-1] >= Const::HalfPi())
				sum1 += Const::TwoPi();
			phase[i-1] += sum0;
			sum0 = sum1;
		}
		phase[len-1] += sum1;
	}
	void yuv_to_rgb(value *rgb, const value *yuv)
	{
		value WR(0.299), WB(0.114), WG(1-WR-WB), UMAX(0.492), VMAX(0.877);
		rgb[0] = yuv[0] + ((1-WR)/VMAX) * yuv[2];
		rgb[1] = yuv[0] - (WB*(1-WB)/(UMAX*WG)) * yuv[1] - (WR*(1-WR)/(VMAX*WG)) * yuv[2];
		rgb[2] = yuv[0] + ((1-WB)/UMAX) * yuv[1];
	}
	void cmplx_to_rgb(value *rgb, const cmplx *inp)
	{
		value y0py1 = inp[0].real()+inp[1].real();
		value y0my1 = inp[0].imag()-inp[1].imag();
		value yuv[6] = {
			(y0py1+y0my1)/2, (inp[0].imag()+inp[1].imag()-y0py1)/2, (inp[0].real()-inp[1].real()-y0my1)/2,
			(y0py1-y0my1)/2, (inp[0].imag()+inp[1].imag()-y0py1)/2, (inp[0].real()-inp[1].real()-y0my1)/2
		};
		yuv_to_rgb(rgb, yuv);
		yuv_to_rgb(rgb+3, yuv+3);
	}
	cmplx *mls0_init()
	{
		CODE::MLS seq0(mls0_poly);
		for (int i = 0; i < mls0_len; ++i)
			mls0_sig[i] = 1 - 2 * seq0();
		return mls0_sig;
	}
	void mls1_init()
	{
		CODE::MLS seq1(mls1_poly);
		for (int i = 0; i < mls1_len; ++i)
			mls1_sig[i] = 1 - 2 * seq1();
	}
	Decoder(DSP::WritePEL<value> *pel, DSP::ReadPCM<value> *pcm) : pcm(pcm), resample(rate, (rate * 19) / 40, 2), correlator(mls0_init())
	{
		mls1_init();
		bool real = pcm->channels() == 1;
		const cmplx *buf;
		do  {
			if (!pcm->good())
				return;
			cmplx tmp;
			pcm->read(reinterpret_cast<value *>(&tmp), 1);
			if (real)
				tmp = hilbert(tmp.real());
			buf = input_hist(tmp);
		} while (!correlator(buf));

		symbol_pos = correlator.symbol_pos;
		cfo_rad = correlator.cfo_rad;
		sfo_rad = correlator.sfo_rad;
		std::cerr << "symbol pos: " << symbol_pos << std::endl;
		std::cerr << "coarse sfo: " << 1000000 * sfo_rad / Const::TwoPi() << " ppm" << std::endl;
		std::cerr << "coarse cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;

		for (int n = 0; n < 2; ++n) {
			if (n) {
				for (int i = 0; i < mls1_len; ++i)
					fdom[i+mls1_off] = head[i+mls1_off] / tail[i+mls1_off];
				value avg_pwr(0);
				for (int i = 0; i < mls1_len; ++i)
					avg_pwr += norm(fdom[i+mls1_off]);
				avg_pwr /= value(mls1_len);
				int count = 0;
				for (int i = 0; i < mls1_len; ++i) {
					value power = norm(fdom[i+mls1_off]);
					if (2 * power > avg_pwr && power < 2 * avg_pwr) {
						phase[count] = arg(fdom[i+mls1_off]);
						index[count] = i+mls1_off;
						++count;
					}
				}
				unwrap(phase, count);
				DSP::SimpleLinearRegression<value> sfo_cfo(index, phase, count);
				value slope = sfo_cfo.slope();
				value yint = sfo_cfo.yint();
				yint -= Const::TwoPi() * std::nearbyint(sfo_cfo(mls1_off+mls1_len/2) / Const::Pi());
				sfo_rad += slope * symbol_len / value((img_height+3)*(symbol_len+guard_len));
				cfo_rad -= yint / value((img_height+3)*(symbol_len+guard_len));
				std::cerr << "finer sfo: " << 1000000 * sfo_rad / Const::TwoPi() << " ppm" << std::endl;
				std::cerr << "finer cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz" << std::endl;
			}
			value diff = sfo_rad * (rate / Const::TwoPi());
			resample(tdom, buf, -diff, buffer_len);
			DSP::Phasor<cmplx> osc;
			osc.omega(-cfo_rad);
			for (int i = 0; i < buffer_len; ++i)
				tdom[i] *= osc();
			symbol_pos = correlator.symbol_pos * (1 - sfo_rad / Const::TwoPi());
			std::cerr << "resampled pos: " << symbol_pos << std::endl;
			fwd(head, tdom+symbol_pos-(img_height+2)*(symbol_len+guard_len));
			fwd(tail, tdom+symbol_pos+symbol_len+guard_len);
		}

		value img_fac = sqrt(value(symbol_len) / value(64 * img_width));
		value mls1_fac = sqrt(value(symbol_len) / value(4 * mls1_len));
		for (int i = 0; i < mls1_len; ++i)
			head[i+mls1_off] *= (img_fac / mls1_fac) * mls1_sig[i];
		for (int i = 0; i < mls1_len; ++i)
			tail[i+mls1_off] *= (img_fac / mls1_fac) * mls1_sig[i];
		for (int j = 0; j < img_height; ++j) {
			fwd(fdom, tdom+symbol_pos+(j-img_height-1)*(symbol_len+guard_len));
			value x = value(j+1) / value(img_height+1);
			for (int i = 0; i < img_width; ++i)
				fdom[i+img_off] /= DSP::lerp(x, head[i+img_off], tail[i+img_off]);
			for (int i = 0; i < img_width; i += 2)
				cmplx_to_rgb(rgb_line+3*i, fdom+i+img_off);
			pel->write(rgb_line, img_width);
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

