#ifndef PIANO_H
#define PIANO_H

double TUNE[3] = {1, 1.0003, 0.9996};
#define NUM_NOTES 128

#include "filter.h"
#include "dwgs.h"
#include "reverb.h"
#include "hammer.h"
#include "types.h"

class Piano {
public:
    Piano(int note, double Fs, double v0, long samples);

    ~Piano();

    long go(double *out, int m);

//    double getResampleRatio() = 0;

    static void fillFrequencyTable();

    static double freqTable[NUM_NOTES];

    long samples;
    long sample;
    double t;
    double dt;
    double Z;
    double Zb;
    double Zh;
    double f;

    int nStrings;
    dwgs *string[3]{};
    Hammer *hammer;
    Reverb *soundboard;
    Filter shaping1{};
    Filter shaping2{};
    Filter shaping3{};
};

#endif
