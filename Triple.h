//
// Created by 周昊 on 2023/3/8.
//

#ifndef NNGRAM_TRIPLE_H
#define NNGRAM_TRIPLE_H
//#include "./NG/NGramSolver.h"

class Triple {
public:
    long time;
    int X;
    int Y;

    Triple () {}

    Triple (long time, int x, int y): time(time), X(x), Y(y) {}

    Triple (long time, double x, double y): time(time), X((int)(x)), Y((int)(y)) {}

    static double getDis(const Triple* p1, const Triple* p2) {
        return sqrt(pow(p1->X - p2->X, 2) + pow(p1->Y - p2->Y, 2) + pow((p1->time - p2->time) / 60.0, 2));
    }

    static double getPhysicDis(const Triple* p1, int x, int y) {
        return sqrt(pow(p1->X - x, 2) + pow(p1->Y - y, 2));
    }

    static double getPhysicDis(const Triple* p1, const Triple* p2) {
        return sqrt(pow(p1->X - p2->X, 2) + pow(p1->Y - p2->Y, 2));
    }

};



#endif //NNGRAM_TRIPLE_H
