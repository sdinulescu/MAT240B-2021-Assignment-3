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
  std::vector<Grain> grain;

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

      Grain grain;

      std::vector<double> clip(clipSize, 0.0);
      for (int i = 0; i < clip.size(); i++) {
        //          input      *       Hann window
        clip[i] = input[n + i] * (1 - cos(2 * M_PI * i / clip.size())) / 2;
      }

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

      std::vector<Peak> peak;
      for (int i = 1; i < data.size() / 2; i++) {
        // only accept maxima
        //
        if (abs(data[i - 1]) < abs(data[i]))
          if (abs(data[i + 1]) < abs(data[i]))
            peak.push_back({abs(data[i]) / (clip.size() / 2),
                            1.0 * SAMPLE_RATE * i / data.size()});
      }

      std::sort(peak.begin(), peak.end(), [](Peak const &a, Peak const &b) {
        return a.magnitude > b.magnitude;
      });

      peak.resize(10);  // throw away the extras

      // XXX here's where you might estimate f0

      this->grain.emplace_back(grain);
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

  void onSound(AudioIOData &io) override {
    while (io()) {
      float f = 0;

      // XXX
      // f += ?
      //

      f *= dbtoa(db.get());
      io.out(0) = f;
      io.out(1) = f;
    }
  }

  bool onKeyDown(const Keyboard &k) override {
    int midiNote = asciiToMIDI(k.key());

    // respond to user action to re-order the grains
    //
    return true;
  }

  bool onKeyUp(const Keyboard &k) override {
    int midiNote = asciiToMIDI(k.key());
    return true;
  }
};

int main(int argc, char *argv[]) {
  MyApp app(argc, argv);
  app.configureAudio(SAMPLE_RATE, BLOCK_SIZE, 2, 1);
  app.start();
  return 0;
}
