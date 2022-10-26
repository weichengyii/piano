#include "dwgs.h"
#include <cmath>
#include <cstdio>

dwgs::dwgs(double f, double Fs, double inPos, double c1, double c3, double B, double Z, double Zb, double Zh) {
    double delTot = Fs / f;
    int del1 = (int) (inPos * 0.5 * delTot);
    if (del1 < 2)
        del1 = 1;

    if (f > 400) {
        M = 1;
        thirianDispersion(B, f, M, &(dispersion[0]));
    } else {
        M = 4;
        for (int m = 0; m < M; m++)
            thirianDispersion(B, f, M, &(dispersion[m]));
    }
    double dispersionDelay = M * groupDelay(&(dispersion[0]), f, Fs);
    loss(f, c1, c3, &lowPass);
    double lowPassDelay = groupDelay(&lowPass, f, Fs);

    int del2 = (int) (0.5 * (delTot - 2.0 * del1) - dispersionDelay);
    int del3 = (int) (0.5 * (delTot - 2.0 * del1) - lowPassDelay - 5.0);
    if (del2 < 2)
        del2 = 1;
    if (del3 < 2)
        del3 = 1;

    double D = (delTot - (double) (del1 + del1 + del2 + del3 + dispersionDelay + lowPassDelay));
    thirian(D, (int) D, &fracDelay);
    double tuningDelay = groupDelay(&fracDelay, f, Fs);

    printf("total delay = %g/%g, left del = %d/%d, right del = %d/%d,"
           " dispersion delay = %g, low pass delay = %g, fractional delay = %g/%g\n",
           del1 + del1 + del2 + del3 + dispersionDelay + lowPassDelay + tuningDelay,
           delTot, del1, del1, del2, del3, dispersionDelay, lowPassDelay, tuningDelay, D);

    d[0] = new dwg(Z, del1, del1, 0, this);
    d[1] = new dwg(Z, del2, del3, 1, this);
    d[2] = new dwg(Zb, 0, 0, 0, this);
    d[3] = new dwg(Zh, 0, 0, 0, this);

    d[0]->connectRight(d[1]->l);
    d[1]->connectLeft(d[0]->r);
    d[1]->connectRight(d[2]->l);
    d[2]->connectLeft(d[1]->r);

    d[0]->connectRight(d[3]->l);
    d[1]->connectLeft(d[3]->l);
    d[3]->connectLeft(d[0]->r);
    d[3]->connectLeft(d[1]->l);

    d[0]->init();
    d[1]->init();
    d[2]->init();
    d[3]->init();
}

dwgs::~dwgs() {
    for (auto &k: d) {
        delete k;
    }
}

dwgNode::dwgNode(double z) {
    a[0] = 0;
    a[1] = 0;
    this->z = z;
    this->load = 0;
}

dwg::dwg(double z, int del1, int del2, int commute, dwgs *parent) {
    this->parent = parent;

    if (del1 > 1) {
        init_delay(&(d[0]), del1 - 1);
    }
    if (del2 > 1) {
        init_delay(&(d[1]), del2 - 1);
    }

    this->del1 = del1;
    this->del2 = del2;
    nl = 0;
    nr = 0;
    l = new dwgNode(z);
    r = new dwgNode(z);
    this->commute = commute;
}

void dwg::init() {
    double zTot;

    zTot = l->z;
    for (int k = 0; k < nl; k++) {
        zTot += cl[k]->z;
    }
    alphaThisL = 2.0 * l->z / zTot;
    for (int k = 0; k < nl; k++) {
        alphaL[k] = 2.0 * cl[k]->z / zTot;
    }

    zTot = r->z;
    for (int k = 0; k < nr; k++) {
        zTot += cr[k]->z;
    }
    alphaThisR = 2.0 * r->z / zTot;
    for (int k = 0; k < nr; k++) {
        alphaR[k] = 2.0 * cr[k]->z / zTot;
    }

}

dwg::~dwg() {
    delete l;
    delete r;
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
        dar = r->a[0];
    else
        dar = delay(r->a[0], &(d[0]));

    double dal;
    if (del2 == 1)
        dal = l->a[1];
    else
        dal = delay(l->a[1], &(d[1]));

    l->a[0] = dar;
    r->a[1] = dal;
}

void dwg::doLoad() {
    if (nl == 0)
        loadL = 0;
    else {
        loadL = alphaThisL * l->a[0];
        for (int k = 0; k < nl; k++) {
            int polarity = pl[k] ? 0 : 1;
            loadL += cl[k]->load;
            loadL += alphaL[k] * cl[k]->a[polarity];
        }
    }

    if (nr == 0)
        loadR = 0;
    else {
        loadR = alphaThisR * r->a[1];
        for (int k = 0; k < nr; k++) {
            int polarity = pr[k] ? 1 : 0;
            loadR += cr[k]->load;
            loadR += alphaR[k] * cr[k]->a[polarity];
        }
    }
}

void dwg::update() const {
    double a = (loadL - l->a[0]);
    if (commute) {
        for (int m = 0; m < parent->M; m++) {
            a = filter(a, &(parent->dispersion[m]));
        }
    }
    l->a[1] = a;

    a = (loadR - r->a[1]);
    if (commute) {
        a = filter(a, &(parent->lowPass));
        a = filter(a, &(parent->fracDelay));
    }
    r->a[0] = a;
}


double dwgs::input_velocity() {
    return d[1]->l->a[0] + d[0]->r->a[1];
}

double dwgs::go_hammer(double load) {
    d[3]->l->load = load;
    for (int k = 0; k < 2; k++) {
        d[k]->doDelay();
    }

    return d[1]->r->a[1];
}

double dwgs::go_soundboard(double load) {
    d[2]->l->load = load;
    for (int k = 0; k < 3; k++) {
        d[k]->doLoad();
    }

    for (int k = 0; k < 3; k++) {
        d[k]->update();
    }

    return d[2]->l->a[1];
}
