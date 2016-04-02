#include <iostream>
#include <vector>
#include "OpenNL_psm.h"

int main() {
    const int N = 2500; // 5 seconds

    const double x0 = 245; // 245 mm from the goal
    const double v0 = 0;   // zero initial speed

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, N*3);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);

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
        nlRowScaling(1000);
        nlBegin(NL_ROW); // x{i+1} = xi + vi
        nlCoefficient((i+1)*3  ,  -1);
        nlCoefficient((i  )*3  ,   1);
        // well, on microcontroller side the position will be measured in encoder ticks and not in mm
        // the cart moves for 1 mm per 100 encoder ticks
        // for computational stability reasons the position in this program is given in mm and not in ticks,
        // therefore the equation i am solving here is x{i+1} = xi + .01 vi,
        // since the speed is still in ticks / dt
        nlCoefficient((i  )*3+1, .01);
        nlEnd(NL_ROW);

        nlRowScaling(1000);
        nlBegin(NL_ROW); // v{i+1} = .97 vi + .218 ui
        nlCoefficient((i+1)*3+1, -1);
        nlCoefficient((i  )*3+1,  .97);
        nlCoefficient((i  )*3+2,  .218);
        nlEnd(NL_ROW);
    }

    for (int i=0; i<N; i++) {
        nlRowScaling(1);
        nlBegin(NL_ROW); // xi = 0, soft
        nlCoefficient(i*3, 1);
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

    for (int i=0; i<N; i++) {
        std::cout << i*2 << " " << solution[i*3] << " " << solution[i*3+1]/2*10 << " " << solution[i*3+2] << std::endl;
    }


    nlDeleteContext(nlGetCurrent());

    // Control curves being found, it is time to find static gain matrix

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

    std::cerr << a << " " << b << std::endl;

    nlDeleteContext(nlGetCurrent());

    return 0;
}

