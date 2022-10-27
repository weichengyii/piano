#ifndef DWGS_H
#define DWGS_H

#include "filter.h"

class dwgNode {
public:
    explicit dwgNode(double z);

    double z;
    double load;
    double a[2]{};
};

class dwgs;

class dwg {
public:
    dwg(double z, int del1, int del2, int commute, dwgs *parent);

    ~dwg();

    void init();

    void update() const;

    void doDelay();

    void doLoad();

    void connectLeft(dwgNode *node);

    void connectRight(dwgNode *node);

    void connectLeft(dwgNode *node, int polarity);

    void connectRight(dwgNode *node, int polarity);

    int del1;
    int del2;
    int nl;
    int nr;
    int pl[2]{};
    int pr[2]{};
    dwgNode *cl[2]{};
    dwgNode *cr[2]{};
    dwgNode *l, *r;
    double loadL{}, loadR{};
    double alphaThisL{};
    double alphaThisR{};
    double alphaL[2]{};
    double alphaR[2]{};

    Delay d[2]{};
    dwgs *parent;
    int commute;
};

class dwgs {
public:
    dwgs(double f, double Fs, double inPos, double c1, double c3, double B, double Z, double Zb, double Zh);

    ~dwgs();

    double input_velocity();

    double go_hammer(double load);

    double go_soundboard(double load);

    Filter dispersion[4]{};
    Filter lowPass{};
    Filter fracDelay{};

    int M;
    dwg *d[4]{};
};

#endif
