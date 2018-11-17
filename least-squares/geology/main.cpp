#include <vector>
#include <iostream>
#include "geometry.h"
#include "model.h"
#include "tgaimage.h"
#include "OpenNL_psm.h"

const int width  = 1280;
const int height = 1280;
const TGAColor white  (255, 255, 255, 255);
const TGAColor green  (  0, 255,   0, 255);
const TGAColor red    (255,   0,   0, 255);
const TGAColor blue   (  0,   0, 255, 255);
const TGAColor yellow (255, 255,   0, 255);
const TGAColor magenta(255,   0, 255, 255);
const TGAColor colors[] = {green, blue, yellow, magenta};
const int ncolors = sizeof(colors)/sizeof(TGAColor);


bool is_edge_present(Vec3f v1, Vec3f v2, Model &m) {
    for (int j=0; j<m.nhalfedges(); j++) {
        Vec3f u1 = m.point(m.from(j));
        Vec3f u2 = m.point(m.to  (j));
        float threshold = 1e-3;
        if (((u1-v1).norm()<threshold && (u2-v2).norm()<threshold) || ((u2-v1).norm()<threshold && (u1-v2).norm()<threshold)) {
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv) {
    if (argc<2) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }

    std::vector<Model> m;
    for (int i=1; i<argc; i++) {
        m.push_back(Model(argv[i]));
    }
    TGAImage frame(width, height, TGAImage::RGB);

    for (int i=0; i<m[0].nhalfedges(); i++) {
        Vec3f v[2] = { m[0].point(m[0].from(i)), m[0].point(m[0].to(i)) };
        TGAColor color = white;

        bool fault = m.size()>=2 && is_edge_present(v[0], v[1], m[1]);

        int horizon = -1;
        for (int j=2; j<(int)m.size(); j++) {
            if (is_edge_present(v[0], v[1], m[j])) {
                horizon = j-2;
                break;
            }
        }

        if (fault) color = red;
        if (horizon>=0) color = colors[horizon%ncolors];

        Vec2i s[2];
        for (int j=0; j<2; j++) {
            s[j] = Vec2i(v[j].x*width, v[j].y*height);
        }

        frame.line(s[0], s[1], color);
    }

    frame.write_tga_file("framebuffer.tga");
    return 0;

#if 0
    {
        std::vector<Model> m;
        for (int i=1; i<argc; i++) {
            m.push_back(Model(argv[i]));
        }
        Vec3f min, max;
        m[0].get_bbox(min, max);
        float maxside = std::max(max.x-min.x, std::max(max.y-min.y, max.z-min.z));
        for (int i=0; i<(int)m.size(); i++) {
            std::cout << "= " << argv[1+i] << std::endl;
            for (int v=0; v<m[i].nverts(); v++) {
                Vec3f p = m[i].point(v);
                p  = (p - min)/maxside;
                m[i].point(v) = Vec3f(p.x, p.z, p.y);
            }
            std::cout << m[i];
        }
    }
    return 0;
#endif
}

