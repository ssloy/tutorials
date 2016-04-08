#include <iostream>
#include <vector>
#include "OpenNL_psm.h"

int main() {
    const int N = 1000; // 2 seconds

    const double x0 = 24.5; // 24.5 cm from the goal
    const double v0 = 0;   // zero initial speed, cm/sec

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, N*4);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);

    nlBegin(NL_SYSTEM);

    nlLockVariable(0);
    nlSetVariable(0, x0);

    nlLockVariable(1);
    nlSetVariable(1, v0);

    nlLockVariable(3);
    nlSetVariable(3, x0);

    nlLockVariable((N-1)*4);
    nlSetVariable((N-1)*4, 0);

    nlLockVariable((N-1)*4+1);
    nlSetVariable((N-1)*4+1, 0);

    nlLockVariable((N-1)*4+2);
    nlSetVariable((N-1)*4+2, 0);

    nlBegin(NL_MATRIX);

    for (int i=0; i<N-1; i++) {
        nlRowScaling(1000);
        nlBegin(NL_ROW); // x{i+1} = xi + vi
        nlCoefficient((i+1)*4  ,   -1);
        nlCoefficient((i  )*4  ,    1);
        nlCoefficient((i  )*4+1, .002);
        nlEnd(NL_ROW);

        nlRowScaling(1000);
        nlBegin(NL_ROW); // v{i+1} = a vi + b ui
        nlCoefficient((i+1)*4+1, -1);
        nlCoefficient((i  )*4+1, 0.974);
        nlCoefficient((i  )*4+2, 0.111808); // remember that x and v are measured in cm and cm/s
        nlEnd(NL_ROW);

        nlRowScaling(1000);
        nlBegin(NL_ROW); // error integral
        nlCoefficient((i+1)*4+3, -1);
        nlCoefficient((i  )*4+3,  1);
        nlCoefficient((i  )*4  , .002);
        nlEnd(NL_ROW);
    }

    for (int i=0; i<N; i++) {
        nlRowScaling(5);
        nlBegin(NL_ROW); // xi = 0, soft
        nlCoefficient(i*4, 1);
        nlEnd(NL_ROW);

        nlRowScaling(1);
        nlBegin(NL_ROW); // vi = 0, soft
        nlCoefficient(i*4+1, 1);
        nlEnd(NL_ROW);

        nlRowScaling(4);
        nlBegin(NL_ROW); // ui = 0, soft
        nlCoefficient(i*4+2, 1);
        nlEnd(NL_ROW);

        nlRowScaling(.1);
        nlBegin(NL_ROW); // zi = 0, soft
        nlCoefficient(i*4+3, 1);
        nlEnd(NL_ROW);
    }

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    nlSolve();

    std::vector<double> solution;
    for (int i=0; i<4*N; i++) {
        solution.push_back(nlGetVariable(i));
    }

    for (int i=0; i<N; i++) {
        std::cout << i*2./1000. << " " << solution[i*4] << " " << solution[i*4+1] << " " << solution[i*4+2]  << " " << solution[i*4+3] << std::endl;
    }


    nlDeleteContext(nlGetCurrent());

    // Control curves being found, it is time to find static gain matrix

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, 3);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);

    for (int i=0; i<N; i++) {
        nlBegin(NL_ROW);
        nlCoefficient(0, solution[i*4  ]/100); // bring position and speed back to m and m/s
        nlCoefficient(1, solution[i*4+1]/100);
        nlCoefficient(2, solution[i*4+3]);
        nlRightHandSide(solution[i*4+2]);
        nlEnd(NL_ROW);
    }

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    nlSolve();
    double a = nlGetVariable(0);
    double b = nlGetVariable(1);
    double c = nlGetVariable(2);

    std::cerr << a << " " << b << " " << c << std::endl;

    nlDeleteContext(nlGetCurrent());

    return 0;
}

