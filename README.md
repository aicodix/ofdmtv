
### OFDMTV

Quick start:

Encode [smpte.ppm](smpte.ppm) [PNM](https://en.wikipedia.org/wiki/Netpbm) picture file with 320x240 pixels to ```encoded.wav``` [WAV](https://en.wikipedia.org/wiki/WAV) audio file with 8000 Hz sample rate, 16 bits and only 1 (real) channel:

```
./encode encoded.wav 8000 16 1 smpte.ppm
```

Start recording to ```recorded.wav``` audio file and stop after 60 seconds:

```
arecord -c 1 -f S16_LE -r 8000 -d 60 recorded.wav
```

Start playing ```encoded.wav``` audio file:

```
aplay encoded.wav
```

Decode ```recorded.wav``` audio file to ```decoded.ppm``` picture file:

```
./decode decoded.ppm recorded.wav
```

Watch ```decoded.ppm``` picture file in [feh](https://feh.finalrewind.org/):

```
feh decoded.ppm
```

### Simulating

Prerequisite: [disorders](https://github.com/aicodix/disorders)

Encode [smpte.ppm](smpte.ppm) to [analytic](https://en.wikipedia.org/wiki/Analytic_signal) audio signal, add [multipath](https://en.wikipedia.org/wiki/Multipath_propagation), [AWGN](https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise), [SFO, CFO](https://en.wikipedia.org/wiki/Carrier_frequency_offset), decode and compare to original in [feh](https://feh.finalrewind.org/):

```
./encode ping.wav 8000 16 2 smpte.ppm && ../disorders/multipath pong.wav ping.wav ../disorders/multipath.txt 10 && ../disorders/awgn ping.wav pong.wav -50 && ../disorders/cfo pong.wav ping.wav 234.567 && ../disorders/sfo ping.wav pong.wav 40 && ./decode decoded.ppm ping.wav && feh decoded.ppm smpte.ppm
```

