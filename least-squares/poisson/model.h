#ifndef __MODEL_H__
#define __MODEL_H__
#include <vector>
#include <string>
#include "geometry.h"

class Model {
private:
    std::vector<Vec3f> verts;
    std::vector<Vec3i> faces;
    std::vector<std::vector<int> > v2h;    // vertex to halfedge incidency
    std::vector<int> opposites;            // halfedges
    void compute_opposites();             
public:                                   
    Model(const char *filename);          
                                          
    int nverts();                          // number of vertices
    int nfaces();                          // number of triangles
    int nhalfedges();                      // number of halfedges = nfaces()*3
                                          
    Vec3f &point(int i);                   // coordinates of the vertex i
    int vert(int fi, int li);              // index of the vertex for the triangle fi and local index li
    void get_bbox(Vec3f &min, Vec3f &max); // bounding box for all the vertices, including isolated ones

    int first_halfedge(int vid);           // retrieve the index of the first halfedge in the star for the given vertex
    int from(int hid);                     // starting vertex for the halfedge
    int to(int hid);                       // end vertex for the halfedge
    int opp(int hid);                      // opposite halfedge or -1 if it does not exist
    int next(int hid);                     // next halfedge in the triangle
    int prev(int hid);                     // previous halfedge in the triangle, prev(hid) = next(next(hid))
};

std::ostream& operator<<(std::ostream& out, Model &m);

#endif //__MODEL_H__

