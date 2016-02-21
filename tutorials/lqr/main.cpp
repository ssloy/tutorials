#include <iostream>
#include "OpenNL_psm.h"

int main() {
    const int N = 60;
    const double xN = 2.3;
    const double x0 = .5;
    const double hard_penalty = 100.;

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, N*2);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);

    nlBegin(NL_ROW); // x0 = x0
    nlCoefficient(0, 1);
    nlRightHandSide(x0);
    nlScaleRow(hard_penalty);
    nlEnd(NL_ROW);

    nlBegin(NL_ROW); // xN = xN
    nlCoefficient((N-1)*2, 1);
    nlRightHandSide(xN);
    nlScaleRow(hard_penalty);
    nlEnd(NL_ROW);

    nlBegin(NL_ROW); // uN = 0, for convenience, normally uN is not defined
    nlCoefficient((N-1)*2+1, 1);
    nlScaleRow(hard_penalty);
    nlEnd(NL_ROW);

    for (int i=0; i<N-1; i++) {
        nlBegin(NL_ROW); // x{N+1} = xN + uN
        nlCoefficient((i+1)*2  , -1);
        nlCoefficient((i  )*2  ,  1);
        nlCoefficient((i  )*2+1,  1);
        nlScaleRow(hard_penalty);
        nlEnd(NL_ROW);
    }

    for (int i=0; i<N; i++) {
        nlBegin(NL_ROW); // ui = 0, soft
        nlCoefficient(i*2+1, 1);
        nlEnd(NL_ROW);
    }

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    nlSolve();

    for (int i=0; i<N; i++) {
        std::cout << nlGetVariable(i*2) << " " <<  nlGetVariable(i*2+1) << std::endl;
    }

    nlDeleteContext(nlGetCurrent());
    return 0;
}

