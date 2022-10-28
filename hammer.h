#ifndef HAMMER_H
#define HAMMER_H
#include <cstdio>

class Hammer {
public:
  Hammer(double f, double Fs, double m, double K, double p, double Z, double alpha, double v0);
  ~Hammer();
  
  double F_out(double t, double F_in);

  double dt;
  double x;
  double v;
  double a;
  double m;

  int S;
  double F0;
  double p;
  double F;
  double upPrev;
  double alpha;
  double Z;
};

#endif
