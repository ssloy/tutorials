#include <vector>
#include <iostream>

#include "model.h"
#include "tgaimage.h"

const int width  = 1280;
const int height = 1280;
const TGAColor white(255, 255, 255, 255);
const TGAColor green(  0, 255,   0, 255);
const TGAColor red  (255,   0,   0, 255);


#include "geometry.h"


/*
*/

// well, a line is a line, no tricks here
void line(Vec2i a, Vec2i b, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(a.x-b.x)<std::abs(a.y-b.y)) {
        std::swap(a.x, a.y);
        std::swap(b.x, b.y);
        steep = true;
    }
    if (a.x>b.x) {
        std::swap(a, b);
    }
    for (int x=a.x; x<=b.x; x++) {
        float t = (x-a.x)/(float)(b.x-a.x);
        int y = a.y*(1.-t) + b.y*t;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}


int main(int argc, char** argv) {
    if (argc<2) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }
    Model model(argv[1]);
//    std::cout << model;

//    Model faults("chevron-side-faults-maxx.obj");
//    Model model("chevron/chevron-side-maxx.obj");

    return 0;
    TGAImage frame(width, height, TGAImage::RGB);

    Vec3f min, max;
    model.get_bbox(min, max);
    float maxside = std::max(max.x-min.x, std::max(max.y-min.y, max.z-min.z));

    for (int i=0; i<model.nhalfedges(); i++) {
        Vec3f v[2] = { model.point(model.from(i)), model.point(model.to(i)) };
        TGAColor color = white;

/*
        bool flag[2] = {false, false};
        for (int j=0;j<2; j++) {
            for (int k=0; k<horizons.nverts(); k++) {
                Vec3f u = horizons.point(k);
                flag[j] |= (u-v[j])*(u-v[j]) < 1e-2;
            }
        }
        if (flag[0] && flag[1]) color = green;

/*
        bool flag2[2] = {false, false};
        for (int j=0;j<2; j++) {
            for (int k=0; k<faults.nverts(); k++) {
                Vec3f u = faults.point(k);
                flag2[j] |= (u-v[j])*(u-v[j]) < 1e-2;
            }
        }
        if (flag2[0] && flag2[1]) color = red;
*/

        Vec2i s[2];
        for (int j=0; j<2; j++) {
            v[j]  = (v[j] - min)/maxside;
            s[j] = Vec2i(v[j].y*width, v[j].z*height);
        }

        line(s[0], s[1], frame, color);
    }

    frame.write_tga_file("framebuffer.tga");

//    print_obj();
    return 0;
}

