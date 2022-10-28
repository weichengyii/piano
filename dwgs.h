#ifndef DWGS_H
#define DWGS_H

#include "filter.h"
#include <vector>

class dwgNode {
public:
    explicit dwgNode(double z);

    double z;
    double load;
    double a[2]{};

    double alpha{};
};

class dwgs;

class dwg {
public:
    dwg(double Z, int del1, int del2, int commute, dwgs *parent);

    ~dwg();

    void init();

    void update() const;

    void doDelay();

    void doLoad();

    void connectLeft(dwgNode *node);

    void connectRight(dwgNode *node);

    int nl;
    int nr;
    dwgNode *cl[2]{};
    dwgNode *cr[2]{};
    dwgNode *lNode, *rNode;
    double lLoad{}, rLoad{};
    double alphaThisL{};
    double alphaThisR{};
    double alphaL[2]{};
    double alphaR[2]{};

    Delay lower{};
    Delay upper{};
    dwgs *parent;
    int commute;
};

class dwgs {
public:
    dwgs(double f, double Fs, double inPos, double B, double Z, double Zb, double Zh);

    ~dwgs();

    [[nodiscard]] double input_velocity() const;

    [[nodiscard]] double go_hammer(double v) const;

    [[nodiscard]] double go_soundboard(double load) const;

    Filter dispersion[4]{};
    Filter lowPass{};
    Filter fracDelay{};

    int M;
    dwg *dwgLeftString, *dwgRightString;
    dwg *dwgBridge, *dwgHammer;
};

#endif
