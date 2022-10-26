#include "hammer.h"
#include <cmath>

Hammer::Hammer(double f, double Fs, double m, double K, double p, double Z, double alpha, double v0) {
    this->S = 3;
    this->p = p;
    this->K = K;
    this->mi = 1.0 / m;
    this->alpha = alpha;
    this->Z2i = 1.0 / (2.0 * Z);
    this->dt = 1.0 / (Fs * S);
    this->dti = 1.0 / dt;
    this->x = 0;
    this->v = v0;
    this->a = 0.0;
    this->F = 0.0;
    this->upPrev = 0;
}

Hammer::~Hammer() = default;

double Hammer::load(double t, double vin) {
    for (int k = 0; k < S; k++) {
        double up;
        up = (x > 0) ? pow(x, p) : 0;
        double dup_dt = (up - upPrev) * dti;
        double v1;
        double x1;
        for (int i = 0; i < 3; i++) {
            F = K * (up + alpha * dup_dt);
            if (F < 0) F = 0;
            a = -F * mi;
            v1 = v + a * dt;
            x1 = x + (v1 - (vin + F * Z2i)) * dt;
            double upnew = (x1 > 0) ? pow(x1, p) : 0;
            double dupdtnew = (upnew - upPrev) / (2 * dt);
            double change = dupdtnew - dup_dt;
            dup_dt = dup_dt + (double) 0.5 * change;
        }
        upPrev = up;
        v = v1;
        x = x1;
    }

    return F;
}
