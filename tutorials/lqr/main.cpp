#include <iostream>
#include <vector>
#include "OpenNL_psm.h"

int main() {
    const int N = 60;
    const double x0 = 9.1;
    const double v0 = -2.5;
    const double hard_penalty = 100.;
    const double rho = 16.;

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, N*3);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);

    nlBegin(NL_ROW);
    nlCoefficient(0, 1); // x0 = 3.1
    nlRightHandSide(x0);
    nlScaleRow(hard_penalty);
    nlEnd(NL_ROW);

    nlBegin(NL_ROW);
    nlCoefficient(1, 1); // v0 = .5
    nlRightHandSide(v0);
    nlScaleRow(hard_penalty);
    nlEnd(NL_ROW);

    nlBegin(NL_ROW);
    nlCoefficient((N-1)*3, 1); // xN = 0
    nlScaleRow(hard_penalty);
    nlEnd(NL_ROW);

    nlBegin(NL_ROW);
    nlCoefficient((N-1)*3+1, 1); // vN = 0
    nlScaleRow(hard_penalty);
    nlEnd(NL_ROW);

    nlBegin(NL_ROW); // uN = 0, for convenience, normally uN is not defined
    nlCoefficient((N-1)*3+2, 1);
    nlScaleRow(hard_penalty);
    nlEnd(NL_ROW);

    for (int i=0; i<N-1; i++) {
        nlBegin(NL_ROW); // x{N+1} = xN + vN
        nlCoefficient((i+1)*3  , -1);
        nlCoefficient((i  )*3  ,  1);
        nlCoefficient((i  )*3+1,  1);
        nlScaleRow(hard_penalty);
        nlEnd(NL_ROW);

        nlBegin(NL_ROW); // v{N+1} = vN + uN
        nlCoefficient((i+1)*3+1, -1);
        nlCoefficient((i  )*3+1,  1);
        nlCoefficient((i  )*3+2,  1);
        nlScaleRow(hard_penalty);
        nlEnd(NL_ROW);
    }

    for (int i=0; i<N; i++) {
        nlBegin(NL_ROW); // xi = 0, soft
        nlCoefficient(i*3, 1);
        nlEnd(NL_ROW);

        nlBegin(NL_ROW); // vi = 0, soft
        nlCoefficient(i*3+1, 1);
        nlEnd(NL_ROW);

        nlBegin(NL_ROW); // ui = 0, soft
        nlCoefficient(i*3+2, 1);
        nlScaleRow(rho);
        nlEnd(NL_ROW);
    }

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    nlSolve();

    std::vector<double> solution;
    for (int i=0; i<3*N; i++) {
        solution.push_back(nlGetVariable(i));
    }

    nlDeleteContext(nlGetCurrent());

     double a = -0.0513868, b = -0.347324;

    std::cerr << a << " " << b << std::endl;

    double xi = x0;
    double vi = v0;
    for (int i=0; i<N; i++) {
        double ui = xi*a + vi*b;
        xi = xi + vi;
        vi = vi + ui;
        std::cout << (ui-solution[i*3+2]) << std::endl;
    }

    return 0;
}

