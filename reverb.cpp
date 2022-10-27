#include "filter.h"
#include "reverb.h"

Reverb::Reverb(double c1, double c3, double a, double mix, double Fs) {
    this->mix = mix;
    int lengths[8] = {37, 87, 181, 271, 359, 592, 687, 721};
    double aa[8] = {a, 1+a, a, a, a, a, a, a};

    for (int k = 0; k < 8; k++) {
        init_delay(&(d[k]), lengths[k]);
        o[k] = 0;
        b[k] = 1;
        c[k] = (k % 2 == 0) ? 1.0 / 8.0 : -1.0 / 8.0;

        loss(Fs / lengths[k], c1, c3, &(decay[k]));
    }
    for (int j = 0; j < 8; j++)
        for (int k = 0; k < 8; k++)
            A[j][k] = aa[(8 + (k - j)) % 8];
}

double Reverb::reverb(double in) {
    double i[8];

    for (int j = 0; j < 8; j++) {
        i[j] = b[j] * in;
        for (int k = 0; k < 8; k++) {
            i[j] += A[j][k] * o[k];
        }
    }

    double out = 0;
    for (int j = 0; j < 8; j++) {
        o[j] = filter(delay(i[j], &(d[j])), &(decay[j]));
        out += c[j] * o[j] * .5;
    }

    return mix * out + (1.0 - mix) * in;
}

