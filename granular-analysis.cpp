// Assignment 3 due 01/29/21 by 12pm
// Written by Stejara Dinulescu
// Built off of stft=peaks-zoom-corrections.cpp from Assignment 2

#include <algorithm>  // std::sort
#include <complex>
#include <iostream>
#include <valarray>
#include <vector>

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
typedef std::pair<double, double>
    freq_amp;  // https://www.geeksforgeeks.org/pair-in-cpp-stl/

double dbtoa(double db) { return pow(10.0, db / 20.0); }

// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray &x) {
  const size_t N = x.size();
  if (N <= 1) return;

  // divide
  CArray even = x[std::slice(0, N / 2, 2)];
  CArray odd = x[std::slice(1, N / 2, 2)];

  // conquer
  fft(even);
  fft(odd);

  // combine
  for (size_t k = 0; k < N / 2; ++k) {
    Complex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
    x[k] = even[k] + t;
    x[k + N / 2] = even[k] - t;
  }
}

void load(std::vector<float> &input, const char *filePath) {
  unsigned int channels;
  unsigned int sampleRate;
  drwav_uint64 totalPCMFrameCount;
  float *pSampleData = drwav_open_file_and_read_pcm_frames_f32(
      filePath, &channels, &sampleRate, &totalPCMFrameCount, NULL);
  if (pSampleData == NULL) {
    printf("failed to load %s\n", filePath);
    exit(1);
  }

  //
  if (channels == 1)
    for (int i = 0; i < totalPCMFrameCount; i++) {
      input.push_back(pSampleData[i]);
    }
  else if (channels == 2) {
    for (int i = 0; i < totalPCMFrameCount; i++) {
      input.push_back((pSampleData[2 * i] + pSampleData[2 * i + 1]) / 2);
    }
  } else {
    printf("can't handle %d channels\n", channels);
    exit(1);
  }

  drwav_free(pSampleData, NULL);
}

bool compare(freq_amp const& i, freq_amp const& j) {
  return (i.second > j.second);
}
// //sort descending by amplitude here... i think
// from http://www.cplusplus.com/reference/algorithm/sort/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

struct Peak {
  double magnitude, frequency;
};

struct Grain {
  int begin;  // index into the vector<float> of the sound file
  int end;    // index into the vector<float> of the sound file

  float peakToPeak;
  float rms;
  float zcr;
  float centroid;
  float f0;

  // add more features?
};

std::vector<float> input;
std::vector<Grain> grain;

const int BLOCK_SIZE = 512;
const int SAMPLE_RATE = 48000;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    printf("granular-synth <.wav>");
    exit(1);
  }

  load(input, argv[1]);
  printf("input size is %ld\n", input.size());
  fflush(stdout);

  // this is how to get a command line argument
  if (argc > 2) {
    int N = std::stoi(argv[2]);
  }

  int clipSize = 2048;
  int hopSize = 1024;
  int fftSize = 8192;

  for (int n = 0; n < input.size() - clipSize; n += hopSize) {
    // XXX build a grain; add it to the vector of grains
    // - peak-to-peak amplitude
    // - root-mean-squared
    // - zero-crossing rate

    Grain grain; // one grain = 2048 samples of audio

    std::vector<double> clip(clipSize, 0.0);
    for (int i = 0; i < clip.size(); i++) {
      //          input      *       Hann window
      clip[i] = input[n + i] * (1 - cos(2 * M_PI * i / clip.size())) / 2;
    }

    // before zero-padding for fft, let's calculate zero crossing rate of the signal (otherwise, value would be off from the zero padding)
    float zcr = 0.0;
    for (int i = 0; i < clip.size(); i++) {
      if (clip[i] == 0.0) { zcr += 1.0; }
    }
    zcr = zcr / clip.size();
    grain.zcr = zcr;

    CArray data;
    data.resize(fftSize);
    for (int i = 0; i < clip.size(); i++) {
      data[i] = clip[i];
    }
    for (int i = clip.size(); i < fftSize; i++) {
      data[i] = 0.0;
    }

    // XXX the size of data really must be a power of two!
    //
    fft(data);

    // XXX here is where you might compute the spectral centroid
    //
    float spectral_centroid = 0.0; // = weighted mean of frequencies present in the signal
                                    // determined using a Fourier transform, with magnitudes used as weights
                                    // reference: https://en.wikipedia.org/wiki/Spectral_centroid
    float accum_denominator = 0.0; // sum of weighted frequency
    float accum_numerator = 0.0; // sum of weighted frequency * center frequency
    for (int i = 0; i < data.size() / 2; i++) { // for loop is the summation from n = 0 to n = N - 1. N is data.size() / 2.
      float cf = 1.0 * SAMPLE_RATE * i / data.size(); // center frequency of the bin
      float magnitude = abs(data[i]) / (clip.size() / 2); // weight (magnitude)

      accum_denominator += magnitude;
      accum_numerator += (cf * magnitude);
    }
    spectral_centroid = accum_numerator / accum_denominator; // calculate
    grain.centroid = spectral_centroid; // set spectral centroid attribute

    // rms amplitude calculation
    float rms = 0.0;
    float zcr = 0.0;
    for (int i = 0; i < data.size(); i++) { // taking all amplitudes here
      rms += std::pow(abs(data[i]) / (clip.size() / 2), 2); // take each amplitude val, square it, and sum across bins
    }
    rms = std::sqrt(rms / data.size()); 
    grain.rms = rms;

    // peak to peak calculation, as well as general peak calculation
    std::vector<Peak> peak;
    float trough = abs(data[0]) / (clip.size() / 2); // initialize with first peak amplitude value
    for (int i = 1; i < data.size() / 2; i++) {
      float amp = abs(data[i]) / (clip.size() / 2); 
      if (amp < trough) { trough = amp; } // find the lowest amplitude value

      // only accept maxima
      if (abs(data[i - 1]) < abs(data[i]))
        if (abs(data[i + 1]) < abs(data[i]))
          peak.push_back({abs(data[i]) / (clip.size() / 2),
                          1.0 * SAMPLE_RATE * i / data.size()});
    }

    std::sort(peak.begin(), peak.end(), [](Peak const &a, Peak const &b) {
      return a.magnitude > b.magnitude;
    });

    peak.resize(10);  // throw away the extras
    grain.peakToPeak = peak[0].magnitude - trough; 

    // XXX here's where you might estimate f0
    
  }

  for (int i = 0; i < grain.size() ; i++) {
    std::cout << grain[i].peakToPeak << ", " 
              << grain[i].rms << ", " 
              << grain[i].zcr << ", " 
              << grain[i].centroid << ", " 
              << grain[i].f0 << std:: endl;
  }
  return 0;
}
