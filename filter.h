#ifndef FILTER_H
#define FILTER_H


struct Filter {
    double *x, *y;
    double *a, *b;
    int n;
};

struct Delay {
    int d1;
    int size;
    int mask;
    int cursor;
    double *x;
    double *y;
};

enum biqaudType {
    pass = 0,
    low,
    high,
    notch
};

long choose(long n, long k);

double probe_delay(Delay *c, int pos);

double groupDelay(Filter *c, double f, double Fs);

double phaseDelay(Filter *c, double f, double Fs);

void thirian(double D, int N, Filter *c);

void thirianDispersion(double B, double f, int M, Filter *c);

void resonator(double f, double Fs, double tau, Filter *c);

void differentiator(Filter *c);

void loss(double f0, double c1, double c3, Filter *c);

void biquad(double f0, double fs, double Q, int type, Filter *c);

void destroy_filter(Filter *c);

double filter(double in, Filter *c);

void init_delay(Delay *c, int di);

void destroy_delay(Delay *c);

double delay(double in, Delay *c);


#endif
