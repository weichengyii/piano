#ifndef HAMMER_H
#define HAMMER_H
#include <cstdio>

class Hammer {
public:
  Hammer(double f, double Fs, double m, double K, double p, double Z, double alpha, double v0);
  ~Hammer();
  
  double load(double t, double vin);

  double dt;
  double dti;
  double x;
  double v;
  double a;

  int S;
  double mi;
  double K;
  double p;
  double F;
  double upPrev;
  double alpha;
  double Z2i;
};

#endif
