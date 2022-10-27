#ifndef REVERB_H
#define REVERB_H

class Reverb {
public:
    Reverb(double c1, double c3, double a, double mix, double Fs);

    double reverb(double in);

    double mix;
    Delay d[8]{};
    double A[8][8]{};
    double o[8]{};
    double b[8]{};
    double c[8]{};
    Filter decay[8]{};
};


#endif
