/**
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include <fstream>
#include <iomanip>
#include <iostream>
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/symmetry.h"
#include "channelflow/utilfuncs.h"
#include "cfbasics/mathdefs.h"
#include "modules/ddc/macros.h"
using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        string purpose(
            "Read Waleffe's dataset and save it");
        
        
        ArgList args(argc, argv, purpose);
        const bool padded = args.getbool("-p", "--padded", false, "set padding modes to zero");
        const string symmstr = args.getstr("-symms", "--symmetries", "", "file of symmetries to satisfy");
        
        const string EQname = args.getstr(1, "<fieldname>", "equilibria file");
        
        args.check();
        args.save("./");

        const int Nx = 32;
        const int Ny = 35; 
        const int Nz = 32;
    
        const Real Lx = 2.0*pi/1.14;//2*pi/1.14
        const Real Lz = 4.0*pi/5.0;//4*pi/5
        const Real ymin = -1;
        const Real ymax = 1;
        
        args.check();
        args.save("./");

        SymmetryList s;
        if (symmstr.length() > 0) {
            s = SymmetryList(symmstr);
            cout << "Restricting random field to invariant subspace generated by symmetries" << endl;
            cout << s << endl;
        }

        FlowField u(Nx, Ny, Nz, 3, Lx, Lz, ymin, ymax);
        u.makePhysical();
        ///////////////////////////////////////////////////
        // https://people.math.wisc.edu/~fwaleffe/ECS/RRC-data.html
        // http://channelflow.org/dokuwiki/doku.php?id=database#fnt__2
        // to run: ./build/modules/ddc/tools/ddc_readEQ
        ifstream inputFileU(EQname+".asc");
        if (!inputFileU) {
            std::cerr << "Error opening file!" << std::endl;
        }
        for (int nx=0; nx<Nx; ++nx)
            for (int ny=0; ny<Ny; ++ny)
                for (int nz=0; nz<Nz; ++nz)
                    for (int i=0; i<3; ++i){
                        inputFileU >> u(nx,ny,nz,i);
                        // cout << u(nx,ny,nz,i) <<endl;
                    }
        inputFileU.close();
        u.makeSpectral();
        Real magn = L2Norm(u);
        for (int i = 0; i < s.length(); ++i)
            u += s[i](u);

        u *= magn / L2Norm(u);
        u.setPadded(padded);
        u.save(EQname);

        #ifdef P5
        FlowField temp(Nx, Ny, Nz, 1, Lx, Lz, ymin, ymax);
        temp.addPerturbations(temp.kxmaxDealiased(), temp.kzmaxDealiased(), 1.0, 0.5);

        for (int i = 0; i < s.length(); ++i)
            temp += s[i](temp);

        temp *= 0.01 / L2Norm(temp);
        temp.setPadded(padded);
        temp.save("t");
        #endif

        #ifdef P6
        FlowField salt(Nx, Ny, Nz, 1, Lx, Lz, ymin, ymax);
        salt.addPerturbations(salt.kxmaxDealiased(), salt.kzmaxDealiased(), 1.0, 0.5);

        for (int i = 0; i < s.length(); ++i)
            salt += s[i](salt);

        salt *= 0.01 / L2Norm(salt);
        salt.setPadded(padded);
        salt.save("s");
        #endif 
    }
    cfMPI_Finalize();
}
