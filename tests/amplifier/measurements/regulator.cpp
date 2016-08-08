// ~$ lscpu | grep Endian
// Byte Order:            Little Endian

#include <cassert>
#include <iostream>
#include <iomanip>      // std::setprecision

#include <fstream>
#include <vector>


unsigned char reverse(unsigned char b) {
    b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
    b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
    b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
    return b;
}


int main () {
    std::ifstream is("analyzer.log");

    std::vector<std::vector<unsigned char> > data;
    data.reserve(200000);
    std::vector<long> time;
    time.reserve(200000);

    long nsamples = 0, nraising = 0, nbytes = 0;
    unsigned char c,p;

    unsigned char current = 0;
    unsigned char   bitno = 0;
    std::vector<unsigned char> packet;

    while (is.read((char *)&c, 1)) {
        if (c&1 && !(p&1)) {
            if (packet.size()==0 && bitno==0) {
                time.push_back(nsamples);
            }
            current |= (((c&4)>>2) << bitno);
            bitno++;
            nraising++;
            if (8==bitno) {
                packet.push_back(reverse(current));
                if (packet.size()==10) {
                    data.push_back(packet);
                    packet.resize(0);
                }
                bitno = 0;
                current = 0;
                nbytes++;
            }

        }
        p = c;
        nsamples++;
    }

    is.close();

    for (int i=0; i<(int)data.size(); i++) {
        unsigned short _g = ((unsigned short)data[i][0] << 8) | ((unsigned short)data[i][1]);
        unsigned short _y = ((unsigned short)data[i][2] << 8) | ((unsigned short)data[i][3]);
        unsigned int   _x = ((unsigned int)data[i][4] << 24) | ((unsigned int)data[i][5] << 16) | ((unsigned int)data[i][6] << 8) | ((unsigned int)data[i][7]);
        unsigned short _p = ((unsigned short)data[i][8] << 8) | ((unsigned short)data[i][9]);
        short g = *(short *)(void *)&_g;
        short y = *(short *)(void *)&_y;
        int x = *(int *)(void *)&_x;
        short p = *(short *)(void *)&_p;

//        std::cerr << time[i]/24000000. << "\t"<< g*5./255. << "\t" << y*2.5/255. << std::endl;
        std::cerr << std::fixed << std::setprecision(9) << time[i]/24000000. << "\t"<< g/1000. << "\t" << y/1000. << "\t" << x/1000000. << "\t" << (p/255.*24.) << std::endl;
    }


//    std::cerr << nsamples << " " << nraising << " " << nbytes << std::endl;
    assert(time.size()==data.size());
    assert(nsamples==200000000L);
    assert(nraising % (8*10)==0);
    assert(nraising == nbytes*8);

    return 0;
}

