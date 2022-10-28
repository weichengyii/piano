#include "dwgs.h"
#include <cmath>
#include <cstdio>

dwgs::dwgs(double f, double Fs, double inPos, double B, double Z, double Zb, double Zh) {
    // delTot = 2n
    // n = L / dx = L / (a*dt) = L / (2Lf*dt) = 1/(2f*dt) = (1/dT)/2f = Fs/(2f)
    // delTot = 2 * n = Fs / f
    double delTot = Fs / f;
    int del1 = (int) (inPos * 0.5 * delTot);
    if (del1 < 2) del1 = 1;

    if (f > 400) {
        M = 1;
        thirianDispersion(B, f, M, &(dispersion[0]));
    } else {
        M = 4;
        for (int m = 0; m < M; m++)
            thirianDispersion(B, f, M, &(dispersion[m]));
    }
    double dispersionDelay = M * groupDelay(&(dispersion[0]), f, Fs);

    double c1 = 0.25;
    double c3 = 5.85;
    loss(f, c1, c3, &lowPass);
    double lowPassDelay = groupDelay(&lowPass, f, Fs);

    int del2 = (int) (0.5 * (delTot - 2.0 * del1) - dispersionDelay);
    int del3 = (int) (0.5 * (delTot - 2.0 * del1) - lowPassDelay - 5.0);
    if (del2 < 2) del2 = 1;
    if (del3 < 2) del3 = 1;

    double D = (delTot - (double) (del1 + del1 + del2 + del3 + dispersionDelay + lowPassDelay));
    thirian(D, (int) D, &fracDelay);
    double tuningDelay = groupDelay(&fracDelay, f, Fs);

    printf("total delay = %g/%g, left del = %d/%d, right del = %d/%d,"
           " dispersion delay = %g, low pass delay = %g, fractional delay = %g/%g\n",
           del1 + del1 + del2 + del3 + dispersionDelay + lowPassDelay + tuningDelay,
           delTot, del1, del1, del2, del3, dispersionDelay, lowPassDelay, tuningDelay, D);

    dwgLeftString = new dwg(Z, del1, del1, 0, this);
    dwgRightString = new dwg(Z, del2, del3, 1, this);
    dwgBridge = new dwg(Zb, 0, 0, 0, this);
    dwgHammer = new dwg(Zh, 0, 0, 0, this);

    dwgLeftString->connectRight(dwgRightString->nodeL);
    dwgRightString->connectLeft(dwgLeftString->nodeR);
    dwgRightString->connectRight(dwgBridge->nodeL);
    dwgBridge->connectLeft(dwgRightString->nodeR);

    dwgLeftString->connectRight(dwgHammer->nodeL);
    dwgRightString->connectLeft(dwgHammer->nodeL);
    dwgHammer->connectLeft(dwgLeftString->nodeR);
    dwgHammer->connectLeft(dwgRightString->nodeL);

    dwgLeftString->init();
    dwgRightString->init();
    dwgBridge->init();
    dwgHammer->init();
}

dwgs::~dwgs() {
    delete dwgLeftString;
    delete dwgRightString;
    delete dwgBridge;
    delete dwgHammer;
}

dwgNode::dwgNode(double z) {
    a[0] = 0;
    a[1] = 0;
    this->Z = z;
    this->load = 0;
}

dwg::dwg(double z, int del1, int del2, int commute, dwgs *parent) {
    this->parent = parent;

    if (del1 > 1) {
        init_delay(&(delayLine[0]), del1 - 1);
    }
    if (del2 > 1) {
        init_delay(&(delayLine[1]), del2 - 1);
    }

    this->del1 = del1;
    this->del2 = del2;
    nl = 0;
    nr = 0;
    nodeL = new dwgNode(z);
    nodeR = new dwgNode(z);
    this->commute = commute;
}

void dwg::init() {
    double zTot;

    zTot = nodeL->Z;
    for (int k = 0; k < nl; k++) {
        zTot += cl[k]->Z;
    }
    alphaThisL = 2.0 * nodeL->Z / zTot;
    for (int k = 0; k < nl; k++) {
        alphaL[k] = 2.0 * cl[k]->Z / zTot;
    }

    zTot = nodeR->Z;
    for (int k = 0; k < nr; k++) {
        zTot += cr[k]->Z;
    }
    alphaThisR = 2.0 * nodeR->Z / zTot;
    for (int k = 0; k < nr; k++) {
        alphaR[k] = 2.0 * cr[k]->Z / zTot;
    }

}

dwg::~dwg() {
    delete nodeL;
    delete nodeR;
}

void dwg::connectLeft(dwgNode *node, int polarity) {
    cl[nl] = node;
    pl[nl++] = polarity;
}

void dwg::connectRight(dwgNode *node, int polarity) {
    cr[nr] = node;
    pr[nr++] = polarity;
}

void dwg::connectLeft(dwgNode *node) {
    connectLeft(node, 0);
}

void dwg::connectRight(dwgNode *node) {
    connectRight(node, 0);
}

void dwg::doDelay() {
    double dar;
    if (del1 == 1)
        dar = nodeR->a[0];
    else
        dar = delay(nodeR->a[0], &(delayLine[0]));

    double dal;
    if (del2 == 1)
        dal = nodeL->a[1];
    else
        dal = delay(nodeL->a[1], &(delayLine[1]));

    nodeL->a[0] = dar;
    nodeR->a[1] = dal;
}

void dwg::doLoad() {
    if (nl == 0)
        loadL = 0;
    else {
        loadL = alphaThisL * nodeL->a[0];
        for (int k = 0; k < nl; k++) {
            int polarity = pl[k] ? 0 : 1;
            loadL += cl[k]->load;
            loadL += alphaL[k] * cl[k]->a[polarity];
        }
    }

    if (nr == 0)
        loadR = 0;
    else {
        loadR = alphaThisR * nodeR->a[1];
        for (int k = 0; k < nr; k++) {
            int polarity = pr[k] ? 1 : 0;
            loadR += cr[k]->load;
            loadR += alphaR[k] * cr[k]->a[polarity];
        }
    }
}

void dwg::update() const {
    double a = (loadL - nodeL->a[0]);
    if (commute) {
        for (int m = 0; m < parent->M; m++) {
            a = filter(a, &(parent->dispersion[m]));
        }
    }
    nodeL->a[1] = a;

    a = (loadR - nodeR->a[1]);
    if (commute) {
        a = filter(a, &(parent->lowPass));
        a = filter(a, &(parent->fracDelay));
    }
    nodeR->a[0] = a;
}


double dwgs::input_velocity() const {
    return dwgRightString->nodeL->a[0] + dwgLeftString->nodeR->a[1];
}

double dwgs::go_hammer(double load) const {
    dwgHammer->nodeL->load = load;

    dwgLeftString->doDelay();
    dwgRightString->doDelay();

    return dwgRightString->nodeR->a[1];
}

double dwgs::go_soundboard(double load) const {
    dwgBridge->nodeL->load = load;

    dwgLeftString->doLoad();
    dwgRightString->doLoad();
    dwgBridge->doLoad();

    dwgLeftString->update();
    dwgRightString->update();
    dwgBridge->update();

    return dwgBridge->nodeL->a[1];
}
