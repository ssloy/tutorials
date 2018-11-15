#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <limits>
#include <iostream>

#include "geometry.h"

std::vector<Vec3f> verts;
std::vector<std::vector<int> > faces;
std::vector<std::vector<int> > vvadj;

// fills vvadj (vertex-to-vertex) adjacency array
void reconstruct_adjacency() {
    vvadj = std::vector<std::vector<int> >(verts.size());
    for (int i=0; i<(int)faces.size(); i++) {
        for (int k=0; k<3; k++) {
            int v1 = faces[i][k];
            int v2 = faces[i][(k+1)%3];
            vvadj[v1].push_back(v2);
            vvadj[v2].push_back(v1);
        }
    }
    for (int i=0; i<(int)verts.size(); i++) { // remove duplicates
        std::sort(vvadj[i].begin(), vvadj[i].end());
        vvadj[i].erase(unique(vvadj[i].begin(), vvadj[i].end()), vvadj[i].end());
    }
}

// fills verts and faces arrays, supposes .obj file to have "f " entries without slashes
void load_obj(const char *filename) {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v[i];
            verts.push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            std::vector<int> f;
            int idx;
            iss >> trash;
            while (iss >> idx) {
                idx--; // in wavefront obj all indices start at 1, not zero
                f.push_back(idx);
            }
            faces.push_back(f);
        }
    }
    std::cerr << "# v# " << verts.size() << " f# "  << faces.size() << std::endl;
    reconstruct_adjacency();
}

void print_obj() {
    for (int i=0; i<(int)verts.size(); i++) {
        std::cout << "v " << verts[i] << std::endl;
    }
    for (int i=0; i<(int)faces.size(); i++) {
        std::cout << "f ";
        for (int k=0; k<3; k++) {
            std::cout << (faces[i][k]+1) << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    load_obj("face.obj");

    std::vector<Vec3f> curvature(verts.size(), Vec3f(0,0,0));
    for (int i=0; i<(int)verts.size(); i++) {
        for (int j=0; j<(int)vvadj[i].size(); j++) {
            curvature[i] = curvature[i] - verts[vvadj[i][j]];
        }
        curvature[i] = verts[i] + curvature[i] / (float)vvadj[i].size();
    }

    for (int it=0; it<100; it++) {
        for (int i=0; i<(int)verts.size(); i++) {
            Vec3f bary(0,0,0);
            for (int j=0; j<(int)vvadj[i].size(); j++) {
                bary = bary + verts[vvadj[i][j]];
            }
            bary = bary / (float)vvadj[i].size();
            verts[i] = bary + curvature[i]*2.1; // play with the coeff here
        }
    }

    print_obj();
    return 0;
}

