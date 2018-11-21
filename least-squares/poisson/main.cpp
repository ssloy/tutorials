#include <vector>
#include <iostream>
#include <cmath>
#include "geometry.h"
#include "model.h"
#include "OpenNL_psm.h"

float f(float x, float y) {
    return sin(x*M_PI)*sin(y*M_PI);
}

int main(int argc, char** argv) {
    if (argc<2) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }

    Model m(argv[1]);
    std::vector<bool> boundary_verts(m.nverts(), false);
    for (int i=0; i<m.nhalfedges(); i++) {
        if (m.opp(i)<0) {
            boundary_verts[m.from(i)] = true;
            boundary_verts[m.to  (i)] = true;
        }
    }

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, m.nverts());
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);

    for (int v=0; v<m.nverts(); v++) {
        if (boundary_verts[v]) {
            nlBegin(NL_ROW);
            nlCoefficient(v, 10);
            nlEnd(NL_ROW);
        }
    }

    for (int hi=0; hi<m.nhalfedges(); hi++) {
        if (m.opp(hi)<0 || m.opp(hi)<hi) continue;
        Vec3f v1 = m.point(m.from(hi));
        Vec3f v2 = m.point(m.to  (hi));
        nlBegin(NL_ROW);
        nlCoefficient(m.to  (hi),  1);
        nlCoefficient(m.from(hi), -1);
        nlRightHandSide(f(v2.x, v2.y) - f(v1.x, v1.y));
        nlEnd(NL_ROW);
    }

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    nlSolve();
    for (int i=0; i<m.nverts(); i++) {
        m.point(i).z = nlGetVariable(i);
    }
    std::cout << m << std::endl;
    return 0;
}

