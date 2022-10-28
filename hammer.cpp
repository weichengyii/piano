#include "hammer.h"
#include <cmath>

Hammer::Hammer(double f, double Fs, double m, double K, double p, double Z, double alpha, double v0) {
    this->S = 10;
    this->p = p;
    this->F0 = K;
    this->m = m;
    this->alpha = alpha;
    this->Z = Z;
    this->dt = 1.0 / (Fs * S);
    this->x = 0;
    this->v = v0;
    this->a = 0.0;
    this->F = 0.0;
    this->upPrev = 0;
}

Hammer::~Hammer() = default;

double Hammer::F_out(double t, double F_in) {
    for (int k = 0; k < S; k++) {
        double up;
        up = (x > 0) ? pow(x, p) : 0;
        double dup_dt = (up - upPrev) / dt;
        double v1;
        double x1;
        for (int i = 0; i < 300; i++) {
            F = F0 * (up + alpha * dup_dt);
            if (F < 0) F = 0;
            a = - F / m;
            v1 = v + a * dt;
            x1 = x + (v1 - (F_in + F / (2 * Z))) * dt;
            double up_new = (x1 > 0) ? pow(x1, p) : 0;
            double du_pdt_new = (up_new - upPrev) / (2 * dt);
            double change = du_pdt_new - dup_dt;
            dup_dt = dup_dt + (double) 0.5 * change;
        }
        upPrev = up;
        v = v1;
        x = x1;
    }

    return F;
}
