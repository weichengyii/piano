#include "filter.h"
#include <cmath>
#include <cstring>
#include "types.h"

long choose(long n, long k) {
    long divisor = 1;
    long multiplier = n;
    long answer = 1;
    k = min(k, n - k);
    while (divisor <= k) {
        answer = (answer * multiplier) / divisor;
        multiplier--;
        divisor++;
    }
    return answer;
}

double Db(double B, double f, int M) {
    double C1, C2, k1, k2, k3;
    if (M == 4) {
        C1 = .069618;
        C2 = 2.0427;
        k1 = -.00050469;
        k2 = -.0064264;
        k3 = -2.8743;
    } else {
        C1 = .071089;
        C2 = 2.1074;
        k1 = -.0026580;
        k2 = -.014811;
        k3 = -2.9018;
    }

    double logB = log(B);
    double kd = exp(k1 * logB * logB + k2 * logB + k3);
    double Cd = exp(C1 * logB + C2);
    double halfstep = pow(2.0, 1.0 / 12.0);
    double Ikey = log(f * halfstep / 27.5) / log(halfstep);
    double D = exp(Cd - Ikey * kd);
    return D;
}


void complex_divide(double Hn[2], double Hd[2], double H[2]) {
    double magn = sqrt(Hn[0] * Hn[0] + Hn[1] * Hn[1]);
    double argn = atan2(Hn[1], Hn[0]);
    double magd = sqrt(Hd[0] * Hd[0] + Hd[1] * Hd[1]);
    double argd = atan2(Hd[1], Hd[0]);
    double mag = magn / magd;
    double arg = argn - argd;
    H[0] = mag * cos(arg);
    H[1] = mag * sin(arg);
}

double groupDelay(Filter *c, double f, double Fs) {
    double df = 5;
    double f2 = f + df;
    double f1 = f - df;
    double omega2 = 2 * PI / (Fs / f2); // FS / f2 =  n2
    double omega1 = 2 * PI / (Fs / f1); // FS / f1 = n1
    return (omega2 * phaseDelay(c, f2, Fs) - omega1 * phaseDelay(c, f1, Fs)) / (omega2 - omega1);
}

double phaseDelay(Filter *c, double f, double Fs) {
    double Hn[2];
    double Hd[2];
    double H[2];

    Hn[0] = 0.0;
    Hn[1] = 0.0;
    Hd[0] = 0.0;
    Hd[1] = 0.0;

    double omega = 2 * PI / (Fs / f); // FS / f = n
    int N = c->n;
    for (int k = 0; k <= N; k++) {
        Hn[0] += cos(k * omega) * c->b[k];
        Hn[1] += sin(k * omega) * c->b[k];
    }
    for (int k = 0; k <= N; k++) {
        Hd[0] += cos(k * omega) * c->a[k];
        Hd[1] += sin(k * omega) * c->a[k];
    }
    complex_divide(Hn, Hd, H);
    double arg = atan2(H[1], H[0]);
    if (arg < 0) arg = arg + 2 * PI;

    return arg / omega;
}

void differentiator(Filter *c) {
    c->x = new double[2];
    c->y = new double[2];
    c->a = new double[2];
    c->b = new double[2];
    memset(c->x, 0, 2 * sizeof(double));
    memset(c->y, 0, 2 * sizeof(double));

    c->a[0] = 1;
    c->a[1] = 0;
    c->b[0] = 1;
    c->b[1] = -1;

    c->n = 1;
}


void resonator(double f, double Fs, double tau, Filter *c) {
    c->x = new double[3];
    c->y = new double[3];
    c->a = new double[3];
    c->b = new double[3];
    memset(c->x, 0, 3 * sizeof(double));
    memset(c->y, 0, 3 * sizeof(double));

    double rp = exp(-1 / (tau * Fs));
    double omega = 2 * PI * f / Fs;
    c->a[0] = 1;
    c->a[1] = -2 * rp * cos(omega);
    c->a[2] = rp * rp;
    c->b[0] = 0;
    c->b[1] = sin(omega);
    c->b[2] = 0;
    c->n = 2;
}

void loss(double f0, double c1, double c3, Filter *c) {
    c->x = new double[2];
    c->y = new double[2];
    c->a = new double[2];
    c->b = new double[2];
    memset(c->x, 0, 2 * sizeof(double));
    memset(c->y, 0, 2 * sizeof(double));

    double g = 1.0 - c1 / f0;
    double b = 4.0 * c3 + f0;
    double a1 = (-b + sqrt(b * b - 16.0 * c3 * c3)) / (4.0 * c3);
    c->b[0] = g * (1 + a1);
    c->b[1] = 0;
    c->a[0] = 1;
    c->a[1] = a1;

    c->n = 1;
}

void biquad(double f0, double fs, double Q, int type, Filter *c) {
    c->x = new double[3];
    c->y = new double[3];
    c->a = new double[3];
    c->b = new double[3];
    memset(c->x, 0, 3 * sizeof(double));
    memset(c->y, 0, 3 * sizeof(double));

    double a = 1 / (2 * tan(PI * f0 / fs));
    double a2 = a * a;
    double aoQ = a / Q;
    double d = (4 * a2 + 2 * aoQ + 1);

    c->a[0] = 1;
    c->a[1] = -(8 * a2 - 2) / d;
    c->a[2] = (4 * a2 - 2 * aoQ + 1) / d;

    switch (type) {
        case pass:
            c->b[0] = 2 * aoQ / d;
            c->b[1] = 0;
            c->b[2] = -2 * aoQ / d;
            break;
        case low:
            c->b[0] = 1 / d;
            c->b[1] = 2 / d;
            c->b[2] = 1 / d;
            break;
        case high:
            c->b[0] = 4 * a2 / d;
            c->b[1] = -8 * a2 / d;
            c->b[2] = 4 * a2 / d;
            break;
        case notch:
            c->b[0] = (1 + 4 * a2) / d;
            c->b[1] = (2 - 8 * a2) / d;
            c->b[2] = (1 + 4 * a2) / d;
            break;
    }

    c->n = 2;
}

void thirian(double D, int N, Filter *c) {
    c->x = new double[N + 1];
    c->y = new double[N + 1];
    c->a = new double[N + 1];
    c->b = new double[N + 1];
    memset(c->x, 0, (N + 1) * sizeof(double));
    memset(c->y, 0, (N + 1) * sizeof(double));

    for (int k = 0; k <= N; k++) {
        double ak = (double) choose((long) N, (long) k);
        if (k % 2 == 1)
            ak = -ak;
        for (int n = 0; n <= N; n++) {
            ak *= ((double) D - (double) (N - n));
            ak /= ((double) D - (double) (N - k - n));
        }
        c->a[k] = (double) ak;
        c->b[N - k] = (double) ak;
    }
    c->n = N;
}

void thirianDispersion(double B, double f, int M, Filter *c) {
    int N = 2;
    double D;
    D = Db(B, f, M);

    if (D <= 1.0) {
        c->x = new double[N + 1];
        c->y = new double[N + 1];
        c->a = new double[N + 1];
        c->b = new double[N + 1];
        memset(c->x, 0, (N + 1) * sizeof(double));
        memset(c->y, 0, (N + 1) * sizeof(double));
        c->a[0] = 1;
        c->a[1] = 0;
        c->a[2] = 0;
        c->b[0] = 1;
        c->b[1] = 0;
        c->b[2] = 0;
    } else {
        thirian(D, N, c);
    }
}

void destroy_filter(Filter *c) {
    delete c->a;
    delete c->b;
    delete c->x;
    delete c->y;
}

double filter(double in, Filter *c) {
    double *a = c->a;
    double *b = c->b;
    int n = c->n;
    double *x = c->x + n;
    double *y = c->y + n;
    double *xend = c->x;

    while (x != xend) {
        *(x) = *(x - 1);
        *(y) = *(y - 1);
        x--;
        y--;
    }
    *x = in;
    double out = *(b) * *(x);
    b++;
    a++;
    x++;
    y++;
    xend = x + c->n;

    while (x != xend) {
        out += *(b) * *(x);
        out -= *(a) * *(y);
        b++;
        a++;
        x++;
        y++;
    }
    *(c->y) = out;

    return out;
}


//void init_delay(Delay *c, int di) {
//    // turn size into a mask for quick modding
//    c->size = 2 * di;
//    int p = 0;
//    while (c->size) {
//        c->size /= 2;
//        p++;
//    }
//    c->size = 1;
//    while (p) {
//        c->size *= 2;
//        p--;
//    }
//
//    c->x = new double[c->size];
//    c->y = new double[c->size];
//    memset(c->x, 0, c->size * sizeof(double));
//    memset(c->y, 0, c->size * sizeof(double));
//
//    c->cursor = 0;
//    c->d1 = c->size - di;
//}

//double probe_delay(Delay *c, int pos) {
//    return c->y[(c->cursor - pos + c->size) % (c->size)];
//}
//
//void destroy_delay(Delay *c) {
//    delete c->x;
//    delete c->y;
//}

//double delay(double in, Delay *c) {
//    int cursor = c->cursor;
//    int d1 = c->d1;
//    double y0 = c->x[d1];
//    c->y[cursor] = y0;
//    c->x[cursor] = in;
//    c->d1++;
////    c->d1 &= c->mask;
//    if (c->d1 == c->size) c->d1 = 0;
//    c->cursor++;
////    c->cursor &= c->mask;
//    if (c->cursor == c->size) c->cursor = 0;
//
//    return y0;
//}

void init_delay(Delay *c, int di) {
    // turn size into a mask for quick modding
    c->size = di + 1;

    c->x = new double[c->size];
    memset(c->x, 0, c->size * sizeof(double));
    c->cursor = 0;
}

double probe_delay(Delay *c, int pos) {
    return c->x[(c->cursor - pos + c->size) % (c->size)];
}

void destroy_delay(Delay *c) {
    delete c->x;
}

double delay(double in, Delay *c) {
    c->x[c->cursor] = in;
//    c->cursor++;
//    if (c->cursor == c->size) c->cursor = 0;
    c->cursor = (c->cursor + 1) % c->size;
    return c->x[c->cursor];
}