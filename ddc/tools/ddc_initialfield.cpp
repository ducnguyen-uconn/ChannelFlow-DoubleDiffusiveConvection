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
#include "modules/ddc/macros.h"
#include "modules/ddc/addPerturbations.h"
using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        string purpose(
            "Construct initial fields for DNS or Nsolver");

        ArgList args(argc, argv, purpose);
        const string modelabel = args.getstr("-m", "--modelabel", "nor", "problem/mode of initial conditions");
        const bool padded = args.getbool("-p", "--padded", false, "set padding modes to zero");
        const int Nx = args.getint("-Nx", "--Nx", "# x gridpoints");
        const int Ny = args.getint("-Ny", "--Ny", "# y gridpoints");
        const int Nz = args.getint("-Nz", "--Nz", "# z gridpoints");
        const Real alpha = args.getreal("-a", "--alpha", 0, "Lx = 2 pi/alpha");
        const Real gamma = args.getreal("-g", "--gamma", 0, "Lz = 2 pi/gamma");
        const Real lx = (alpha == 0.0) ? args.getreal("-lx", "--lx", 0.0, "Lx = 2 pi lx") : 1 / alpha;
        const Real lz = (gamma == 0.0) ? args.getreal("-lz", "--lz", 0.0, "Lz = 2 pi lz") : 1 / gamma;
        const Real Lx = (lx == 0.0) ? args.getreal("-Lx", "--Lx", "streamwise (x) box length") : 2 * pi * lx;
        const Real Lz = (lz == 0.0) ? args.getreal("-Lz", "--Lz", "spanwise   (z) box length") : 2 * pi * lz;
        const Real ymin = args.getreal("-ymin", "--ymin", -0.5, "lower wall height (y is wallnormal) ");
        const Real ymax = args.getreal("-ymax", "--ymax", +0.5, "upper wall height (y is wallnormal) ");
        const int seed = args.getint("-sd", "--seed", 1, "seed for random number generator");
        const Real smooth = args.getreal("-s", "--smoothness", 0.5, "smoothness of field, 0 < s < 1");
        const Real magn = args.getreal("-m", "--magnitude", 0.05, "magnitude  of field, 0 < m < 1");
        const bool meanfl = args.getflag("-mf", "--meanflow", "perturb the mean");

        const string symmstr = args.getstr("-symms", "--symmetries", "", "file of symmetries to satisfy");

        
        #if defined(P6)
        const string uname = args.getstr(3, "<fieldname>", "output velocity file");
        const string tname = args.getstr(2, "<fieldname>", "output temperature file");
        const string sname = args.getstr(1, "<fieldname>", "output salinity file");
        #elif defined(P5)
        const string uname = args.getstr(2, "<fieldname>", "output velocity file");
        const string tname = args.getstr(1, "<fieldname>", "output temperature file");
        #else
        const string uname = args.getstr(1, "<fieldname>", "output velocity file");
        #endif
        
        args.check();
        args.save("./");

        srand48(seed);

        FlowField u(Nx, Ny, Nz, 3, Lx, Lz, ymin, ymax);
        FlowField temp(Nx, Ny, Nz, 1, Lx, Lz, ymin, ymax);
        FlowField salt(Nx, Ny, Nz, 1, Lx, Lz, ymin, ymax);
        
        // Perturb velocity field
        cout << "Perturbing velocity, temperature, and salinity fields of " << modelabel << " mode ... " << flush;
        if(modelabel=="yang2021jfm_case3"){// define initial contions of Yang2021JFM's simulation
            addRandomPerturbations(u,1e-3);
            #ifdef P5
            addSinusoidalPerturbations(temp,-0.05,6.0);
            #endif
            #ifdef P6
            addSinusoidalPerturbations(salt,-0.05,6.0);
            #endif
        } else if(modelabel=="yang2021jfm_case5"){// define initial contions of Yang2021JFM's simulation
            addRandomPerturbations(u,1e-3);
            #ifdef P5
            addSinusoidalPerturbations(temp,-0.025,8.0);
            #endif
            #ifdef P6
            addSinusoidalPerturbations(salt,-0.025,8.0);
            #endif
        } else if (modelabel=="eaves2016jfm"){
            addRandomPerturbations(u,1e-4);
            #ifdef P5
            addSinusoidalPerturbations(temp,-0.05,6.0);
            #endif
            #ifdef P6
            addSinusoidalPerturbations(salt,-0.05,6.0);
            #endif
        } else if (modelabel=="sin"){
            addSinusoidalPerturbations(u,-2*magn,1.0);
            #ifdef P5
            addSinusoidalPerturbations(temp,-magn,1.0);
            #endif
            #ifdef P6
            addSinusoidalPerturbations(salt,-magn,1.0);
            #endif
        } else{// random mode
            addRandomPerturbations(u,1e-3);
            #ifdef P5
            addRandomPerturbations(temp,1e-3);
            #endif
            #ifdef P6
            addRandomPerturbations(salt,1e-3);
            #endif
        }
        cout << "done" << endl;
        u.setPadded(padded);
        u.save(uname);
        #ifdef P5
        temp.setPadded(padded);
        temp.save(tname);
        #endif
        #ifdef P6
        salt.setPadded(padded);
        salt.save(sname);
        #endif
    }
    cfMPI_Finalize();
}
