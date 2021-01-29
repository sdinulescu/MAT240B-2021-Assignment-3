// Assignment 3
//
//
// -- Karl Yerkes / 2021-01-23 / MAT240B
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>  // std::sort
#include <algorithm>  // std::sort, std::min
#include <cmath>      // ::cos()
#include <complex>
#include <iostream>
#include <valarray>
#include <vector>

#include "al/app/al_App.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

const int BLOCK_SIZE = 512;
const int SAMPLE_RATE = 48000;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

using namespace al;

struct MyApp : App {
  Parameter background{"background", "", 0.0, "", 0.0f, 1.0f};
  Parameter db{"db", "", -60.0, "", -60.0f, 0.0f};
  ControlGUI gui;

  std::vector<float> input;
  std::vector<Grain> grains;
  std::vector<float> playback;


  MyApp(int argc, char *argv[]) {
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
        double d_ = input[n + i] * (1 - cos(2 * M_PI * i / clip.size())) / 2;
        clip[i] = d_;
        playback.push_back(d_);
      }

      // before zero-padding for fft, let's calculate zero crossing rate of the signal (otherwise, value would be off from the zero padding)
      float zcr = 0.0;
      for (int i = 0; i < clip.size(); i++) {
        if (clip[i] == 0.0) { zcr += 1.0; }
      }
      zcr = zcr / clip.size();
      //std::cout << "zcr: " << zcr << std::endl;
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
      //std::cout << "fft done" << std::endl;

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
      //std::cout << "sc: " << spectral_centroid << std::endl;
      grain.centroid = spectral_centroid; // set spectral centroid attribute

      // rms amplitude calculation
      float rms = 0.0;
      for (int i = 0; i < data.size(); i++) { // taking all amplitudes here
        rms += std::pow(abs(data[i]) / (clip.size() / 2), 2); // take each amplitude val, square it, and sum across bins
      }
      rms = std::sqrt(rms / data.size()); 
      //std::cout << "rms: " << rms << std::endl;
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
      float ptp = peak[0].magnitude - trough;
      //std::cout << "ptp: " << ptp << std::endl;
      grain.peakToPeak = ptp; 

      // XXX here's where you might estimate f0
      std::vector<float> differences;
      std::vector<int> frequencies; 
      for (int i = 0; i < peak.size(); i++) { // for the 10 biggest peaks
        // compute the differences between all frequency values
        for (int j = i; j < peak.size(); j++) {
          float diff = 0.0;
          diff = abs(peak[i].frequency - peak[j].frequency);
          if (differences.size() == 0) {
            differences.push_back(diff); frequencies.push_back(1);
          } else {
            auto iterator = std::find(differences.begin(), differences.end(), diff);
            if ( iterator != differences.end() ) { // difference value found in the vector
              int index = iterator - differences.begin();
              frequencies[index]++;
            } else { differences.push_back(diff); frequencies.push_back(1); }
          }
          // if the difference value has not yet been added to our vector, push it back. else, update the frequency.
          // if (differences.size() == 0) { differences.push_back(diff); frequencies.push_back(1); } 
          // else {
          //   for (int k = 0; k < differences.size(); k++) {
          //     if (differences[k] == diff) { frequencies[k] = frequencies[k] + 1; }
          //     else { differences.push_back(diff); frequencies.push_back(1); }
          //   }
          // }
        }
      }
      // std::cout << "differences found: ";
      // for (int i = 0; i < differences.size(); i++) {
      //   std::cout << differences[i] << ", ";
      // }
      // std::cout << std::endl;
      // find the most common difference
      int index = 0;
      for (int i = 0; i < frequencies.size() - 1; i++) {
        //get the biggest frequency
        if (frequencies[i + 1] > frequencies[i]) { index = i + 1; }
      }
      //std::cout << "index = " << index << " f0 = " << differences[index] << std::endl;
      grain.f0 = differences[index];

      grains.push_back(grain);
    }
  }

  void onCreate() override {
    gui << background;
    gui << db;
    gui.init();
    navControl().active(false);
  }

  void onAnimate(double dt) override {
    //
  }

  void onDraw(Graphics &g) override {
    g.clear(background);
    gui.draw(g);
  }

  int callback_index = 0;
  void onSound(AudioIOData &io) override {
    while (io()) {
      float f = 0;

      // XXX
      // f += ?
      if (!playback.empty()) { f+= playback[callback_index]; }

      if (f > 0.5) { f = 0.5; }

      //f *= dbtoa(db.get());
      io.out(0) = f;
      io.out(1) = f;

      if (callback_index < playback.size()) { callback_index++; } //sound loop index
    }
  }
  
  // read-only, won't overwrite
  // bool compareP2P(Grain const& i, Grain const& j) { return i.peakToPeak > j.peakToPeak; }
  // bool compareRMS(Grain const& i, Grain const& j) { return i.rms > j.rms; }
  // bool compareZCR(Grain const& i, Grain const& j) { return i.zcr > j.zcr; }
  // bool compareSC(Grain const& i, Grain const& j) { return i.centroid > j.centroid; }
  // bool compareF0(Grain const& i, Grain const& j) { return i.f0 > j.f0; }

  void arrangeGrains(std::vector<Grain> const& g, char* const& s) {
    callback_index = 0;
    //push back the specified grain's feature into a vector. this will be used to resynthesize the sound for playback.
    for (int i = 0; i < g.size(); i++) {
      playback.clear();
      if (strcmp(s, "p2p")) { playback[i] = g[i].peakToPeak; }
      else if (strcmp(s, "rms")) { playback[i] = g[i].rms; }
      else if (strcmp(s, "zcr")) { playback[i] = g[i].zcr; }
      else if (strcmp(s, "centroid")) { playback[i] = g[i].centroid; }
      else if (strcmp(s, "f0")) { playback[i] = g[i].f0; }
    }
    std::sort(playback.begin(), playback.end()); // sort the values
    // for (int i = 0; i < playback.size(); i++) {
    //   std::cout << playback[i] << std::endl;
    // }
  }

  bool onKeyDown(const Keyboard &k) override {
    int midiNote = asciiToMIDI(k.key());
    
     
    return true;
  }

  bool onKeyUp(const Keyboard &k) override {
    int midiNote = asciiToMIDI(k.key());


    // respond to user action to re-order the grains
    if (k.key() == '1') { // peak-to-peak
      arrangeGrains(grains, "p2p");
    } else if (k.key() == '2') { // rms
      arrangeGrains(grains, "rms");  
    } else if (k.key() == '3') { // zcr
      arrangeGrains(grains, "zcr");
    } else if (k.key() == '4') { // spectral centroid
      arrangeGrains(grains, "centroid");
    } else if (k.key() == '5') { // f0
      arrangeGrains(grains, "f0");
    } else if (k.key() == ' ') {
      callback_index = 0;
    }

    return true;
  }
};

int main(int argc, char *argv[]) {
  MyApp app(argc, argv);
  app.configureAudio(SAMPLE_RATE, BLOCK_SIZE, 2, 1);
  app.start();
  return 0;
}
