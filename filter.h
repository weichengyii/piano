#ifndef FILTER_H
#define FILTER_H


struct Filter {
    double *x, *y;
    double *a, *b;
    int n;
};

struct Delay {
    int size;
    int cursor;
    double *x;

    void init(int di);
    double delay(double in);
    [[nodiscard]] double probe(int pos) const;
    void destroy() const;
};

enum biqaudType {
    pass = 0,
    low,
    high,
    notch
};

long choose(long n, long k);

//double probe_delay(Delay *c, int pos);

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



#endif
