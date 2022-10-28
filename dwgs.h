#ifndef DWGS_H
#define DWGS_H

#include "filter.h"

class dwgNode {
public:
    explicit dwgNode(double z);

    double Z;
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

    Delay delayLine[2]{};
    dwgs *parent;
    int commute;
};

class dwgs {
public:
    dwgs(double f, double Fs, double inPos, double B, double Z, double Zb, double Zh);

    ~dwgs();

    [[nodiscard]] double input_velocity() const;

    [[nodiscard]] double go_hammer(double load) const;

    [[nodiscard]] double go_soundboard(double load) const;

    Filter dispersion[4]{};
    Filter lowPass{};
    Filter fracDelay{};

    int M;
    dwg *leftString, *rightString;
    dwg *bridge, *hammer;
};

#endif
