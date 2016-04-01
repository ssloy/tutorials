#include <iostream>
#include <vector>
#include "OpenNL_psm.h"

int main() {
    const int N = 2500;
    const double x0 = 500;
    const double v0 = 180;
    const double hard_penalty = 100.;

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, N*3);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
            nlSolverParameteri(NL_PRECONDITIONER, NL_PRECOND_SSOR);

    nlBegin(NL_SYSTEM);

    nlLockVariable(0);
    nlSetVariable(0, x0);

    nlLockVariable(1);
    nlSetVariable(1, v0);

    nlLockVariable((N-1)*3);
    nlSetVariable((N-1)*3, 0);

    nlLockVariable((N-1)*3+1);
    nlSetVariable((N-1)*3+1, 0);

    nlLockVariable((N-1)*3+2);
    nlSetVariable((N-1)*3+2, 0);

    nlBegin(NL_MATRIX);

    for (int i=0; i<N-1; i++) {
        nlRowScaling(hard_penalty);
        nlBegin(NL_ROW); // x{N+1} = xN + vN
        nlCoefficient((i+1)*3  , -1);
        nlCoefficient((i  )*3  ,  1);
        nlCoefficient((i  )*3+1,  .01);
        nlEnd(NL_ROW);

        nlRowScaling(hard_penalty);
        nlBegin(NL_ROW); // v{N+1} = vN + uN
        nlCoefficient((i+1)*3+1, -1);
        nlCoefficient((i  )*3+1,  .98);
        nlCoefficient((i  )*3+2,  .157);
        nlEnd(NL_ROW);
    }

    for (int i=0; i<N; i++) {
        nlRowScaling(1);
        nlBegin(NL_ROW); // xi = 0, soft
        nlCoefficient(i*3, 1);
//        nlScaleRow(1);
        nlEnd(NL_ROW);

        nlRowScaling(1);
        nlBegin(NL_ROW); // vi = 0, soft
        nlCoefficient(i*3+1, 1);
        nlEnd(NL_ROW);

        nlRowScaling(10);
        nlBegin(NL_ROW); // ui = 0, soft
        nlCoefficient(i*3+2, 1);
        nlEnd(NL_ROW);
    }

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    nlSolve();

    std::vector<double> solution;
    for (int i=0; i<3*N; i++) {
        solution.push_back(nlGetVariable(i));
    }



/*
    for (int i=0; i<N; i++) {
        for (int j=0; j<3; j++) {
            std::cout << solution[i*3+j] << " ";
        }
        std::cout << std::endl;
    }
*/



    nlDeleteContext(nlGetCurrent());


    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, 2);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);

    for (int i=0; i<N; i++) {
        nlBegin(NL_ROW);
        nlCoefficient(0, solution[i*3  ]);
        nlCoefficient(1, solution[i*3+1]);
        nlRightHandSide(solution[i*3+2]);
        nlEnd(NL_ROW);
    }

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    nlSolve();
    double a = nlGetVariable(0);
    double b = nlGetVariable(1);

    nlDeleteContext(nlGetCurrent());


// -0.0360376 -0.0439987

    std::cerr << a << " " << b << std::endl;

    for (int i=0; i<N; i++) {
        double ui = solution[i*3]*a + solution[i*3+1]*b;
        std::cout << solution[i*3+2] << " " << ui<< std::endl;
    }


    return 0;
}

