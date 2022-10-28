#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "types.h"
#include "piano.h"
#include "sndfile.h"

using namespace std;

double Piano::freqTable[NUM_NOTES];

void Piano::fillFrequencyTable() {
    double NOTE_UP_SCALAR = pow(2.0, 1.0 / 12.0);
    double A = 6.875;    // A
    A *= NOTE_UP_SCALAR;    // A#
    A *= NOTE_UP_SCALAR;    // B
    A *= NOTE_UP_SCALAR;    // C, frequency of midi note 0
    for (int i = 0; (i < NUM_NOTES); i++)    // 128 midi notes
    {
        freqTable[i] = A;
        A *= NOTE_UP_SCALAR;
    }
}

long Piano::go(double *out, int m) {
    long n = 0;
    for (int i = 0; i < m; i++) {
        if (sample > this->samples)
            break;
        n++;
        double vString = 0.0;
        for (int k = 0; k < nStrings; k++) {
            vString += string[k]->input_velocity();
        }
        double hLoad = hammer->load(t, vString / nStrings);
        double load = 0;
        for (int k = 0; k < nStrings; k++) {
            load += (2 * Z * string[k]->go_hammer(hLoad / (2 * Z))) / (Z * nStrings + Zb);
        }
        double output = 0.0;
        for (int k = 0; k < nStrings; k++) {
            output += string[k]->go_soundboard(load);
        }

        output = soundboard->reverb(output);

        //output += filter(output,&shaping1);
        //output = filter(output,&shaping2);
        //output += filter(output,&shaping3);
        t += dt;
        sample++;
        out[i] = output * 100.0;
    }
    return n;
}

Piano::Piano(int note, double Fs, double v0, long samples) {
    this->Fs = Fs;
    this->v0 = v0;
    this->samples = samples;
    sample = 0;
    t = 0.0;
    dt = 1.0 / Fs;

    this->f = freqTable[note];

    double f0 = 27.5;
    double rho = 7850.0;
    double p = 2.0 + 1.0 * log(f / f0) / log(4192 / f0);
    double m = 0.06 - 0.058 * pow(log(f / f0) / log(4192 / f0), 0.1);
    double K = 40 / pow(0.7e-3, p);
//    double L = 1.4 - 1.32 * log(f / f0) / log(4192.0 / f0);
    double L = 0.04 + 1.4 / (1 + exp(-3.4 + 1.4 * log(f / f0)));
    double r = 0.002 * pow(1 + 0.6 * log(f / f0), -1.4);
    double rhoL = PI * r * r * rho;
    double T = (2 * L * f) * (2 * L * f) * rhoL;
    Z = sqrt(T * rhoL);
    Zb = 4000.0;
    Zh = 0;
    double E = 200e9;
//    double fLong = sqrt(E / rho) / (2.0 * L);


    double rCore = (r < 0.0006) ? r : 0.0006;
    double B = pow(M_PI, 3) * E * pow(rCore, 4) / (4.0 * L * L * T);
    double hp = 1.0 / 7.0;

    printf("f = %g, r = %g mm, L = %g, T = %g, hammer = %g, Z = %g, K = %g, B = %g\n", f, 1000 * r, L, T, hp, Z, K, B);

    if (note < 31)
        nStrings = 1;
    else if (note < 41)
        nStrings = 2;
    else
        nStrings = 3;

    double c1b = 20.0;
    double c3b = 20.0;
    for (int k = 0; k < nStrings; k++) {
        string[k] = new dwgs(f * TUNE[k], Fs, hp, B, Z, Zb + (nStrings - 1) * Z, Zh);
    }

    double a = -1.0 / 4.0;
    double mix = 1;
    double alpha = 0.1e-4 * log(f / f0) / log(4192 / f0);
    soundboard = new Reverb(c1b, c3b, a, mix, Fs);
    hammer = new Hammer(f, Fs, m, K, p, Z, alpha, v0);

    biquad(500.0, Fs, 10, notch, &shaping1);
    biquad(200.0, Fs, 1.0, high, &shaping2);
    biquad(800.0, Fs, 1.0, low, &shaping3);

}

Piano::~Piano() {
    for (int k = 0; k < nStrings; k++) {
        delete string[k];
    }
    delete hammer;
    delete soundboard;
}

void usage() {
    printf("usage: piano <time(s)> <note[21...108](A4=69)> <velocity[0-10](m/s)> [-b](binary file output to wave.out instead of wav output to out.wav)\n");
    exit(-1);
}

int main() {
    double tTot = 5.0;
    int note = 69;
    double v0 = 5.0;

    Piano::fillFrequencyTable();
    int Fs = 44100;

    long samples = (int) (tTot * Fs);

    PcmWriter *pcmWriter;
    pcmWriter = new PcmWriter((char *) "D:/a.out/piano/out1.wav", samples, Fs, 1);

    auto *piano = new Piano(note, (double) Fs, v0, samples);

    double buf[128];
    while (piano->go(buf, 128)) {
        pcmWriter->write(buf, 128);
    }
    delete pcmWriter;
}
