#include "dwgs.h"
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

    dwgLeftString->connectRight(dwgRightString->lNode);
    dwgRightString->connectLeft(dwgLeftString->rNode);
    dwgRightString->connectRight(dwgBridge->lNode);
    dwgBridge->connectLeft(dwgRightString->rNode);

    dwgLeftString->connectRight(dwgHammer->lNode);
    dwgRightString->connectLeft(dwgHammer->lNode);
    dwgHammer->connectLeft(dwgLeftString->rNode);
    dwgHammer->connectLeft(dwgRightString->lNode);

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
    this->z = z;
    this->load = 0;
}

dwg::dwg(double Z, int del1, int del2, int commute, dwgs *parent) {
    this->parent = parent;

    if (del1 > 1) {
        lower.init(del1 - 1);
    }
    if (del2 > 1) {
        upper.init(del2 - 1);
    }

    nl = 0;
    nr = 0;
    lNode = new dwgNode(Z);
    rNode = new dwgNode(Z);
    this->commute = commute;
}

void dwg::init() { // init alpha
    double zTot;

    zTot = lNode->z;
    for (int k = 0; k < nl; k++) {
        zTot += cl[k]->z;
    }
    alphaThisL = 2.0 * lNode->z / zTot;
    for (int k = 0; k < nl; k++) {
        alphaL[k] = 2.0 * cl[k]->z / zTot;
    }

    zTot = rNode->z;
    for (int k = 0; k < nr; k++) {
        zTot += cr[k]->z;
    }
    alphaThisR = 2.0 * rNode->z / zTot;
    for (int k = 0; k < nr; k++) {
        alphaR[k] = 2.0 * cr[k]->z / zTot;
    }

}

dwg::~dwg() {
    delete lNode;
    delete rNode;
}

void dwg::connectLeft(dwgNode *node) {
    cl[nl] = node;
    nl++;
}

void dwg::connectRight(dwgNode *node) {
    cr[nr] = node;
    nr++;
}

void dwg::doDelay() {
    double dar = lower.delay(rNode->a[0]);

    double dal = upper.delay(lNode->a[1]);

    lNode->a[0] = dar;
    rNode->a[1] = dal;
}

void dwg::doLoad() {
    if (nl == 0)
        lLoad = 0;
    else {
        lLoad = alphaThisL * lNode->a[0];
        for (int k = 0; k < nl; k++) {
            lLoad += cl[k]->load;
            lLoad += alphaL[k] * cl[k]->a[1];
        }
    }

    if (nr == 0)
        rLoad = 0;
    else {
        rLoad = alphaThisR * rNode->a[1];
        for (int k = 0; k < nr; k++) {
            rLoad += cr[k]->load;
            rLoad += alphaR[k] * cr[k]->a[0];
        }
    }
}

void dwg::update() const {
    double a = lLoad - lNode->a[0];
    if (commute) {
        for (int m = 0; m < parent->M; m++) {
            a = filter(a, &(parent->dispersion[m]));
        }
    }
    lNode->a[1] = a;

    a = rLoad - rNode->a[1];
    if (commute) {
        a = filter(a, &(parent->lowPass));
        a = filter(a, &(parent->fracDelay));
    }
    rNode->a[0] = a;
}


double dwgs::input_velocity() const {
    return dwgRightString->lNode->a[0] + dwgLeftString->rNode->a[1];
}

double dwgs::go_hammer(double v) const {
    dwgHammer->lNode->load = v;

    dwgLeftString->doDelay();
    dwgRightString->doDelay();

    return dwgRightString->rNode->a[1];
}

double dwgs::go_soundboard(double load) const {
    dwgBridge->lNode->load = load;

    dwgLeftString->doLoad();
    dwgRightString->doLoad();
    dwgBridge->doLoad();

    dwgLeftString->update();
    dwgRightString->update();
    dwgBridge->update();

    return dwgBridge->lNode->a[1];
}
