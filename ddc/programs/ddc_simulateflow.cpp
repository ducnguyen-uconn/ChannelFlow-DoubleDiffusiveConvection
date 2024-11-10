/**
 * 
 * Original author: Duc Nguyen
 */

#include <iomanip>
#include <iostream>


#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/chebyshev.h"
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/poissonsolver.h"
#include "channelflow/symmetry.h"
#include "channelflow/tausolver.h"
#include "channelflow/utilfuncs.h"
#include "modules/ddc/ddc.h"
#include "modules/ddc/addPerturbations.h"
using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        WriteProcessInfo(argc, argv);
        string purpose(
            "integrate wall-bounded double-diffusive convection "
            "from a given initial condition and save fields to disk.");

        ArgList args(argc, argv, purpose);

        DDCFlags flags(args);
        

        args.section("Program options");
        const string modelabel = args.getstr("-m", "--modelabel", "nor", "problem/mode of initial conditions");

        const string outdir = args.getpath("-o", "--outdir", "data/", "output directory");
        const string ulabel = args.getstr("-ul", "--ulabel", "u", "output velocity field prefix");
        const string tlabel = args.getstr("-tl", "--tlabel", "t", "output temperature field prefix");
        const string slabel = args.getstr("-sl", "--slabel", "s", "output salinity field prefix");

        const int Nx_ = args.getint("-Nx", "--Nx", 32, "# x gridpoints");
        const int Ny_ = args.getint("-Ny", "--Ny", 31, "# y gridpoints");
        const int Nz_ = args.getint("-Nz", "--Nz", 32, "# z gridpoints");
        const Real Lx_ = args.getreal("-Lx", "--Lx", 1, "streamwise (x) box length");
        const Real Lz_ = args.getreal("-Lz", "--Lz", 1, "spanwise   (z) box length");
        const Real ymin_ = args.getreal("-ymin", "--ymin", 0.0, "lower wall height (y is wallnormal) ");
        const Real ymax_ = args.getreal("-ymax", "--ymax", 1.0, "upper wall height (y is wallnormal) ");
       
        args.check();
        args.save("./");
        mkdir(outdir);
        args.save(outdir);
        flags.save(outdir);
        

        TimeStep dt(flags);
        fftw_loadwisdom();

        cout << "Parameters: " << endl;
        cout << "Ra = " << flags.Ra << endl;
        cout << "Le = " << flags.Le << endl;
        cout << "Pr = " << flags.Pr << endl;
        cout << "Rrho = " << flags.Rrho << endl;
        if(P7!=0.0)cout << "Rsep = " << flags.Rsep << endl;
        cout << "Ua = " << flags.ulowerwall << endl;
        cout << "Wa = " << flags.wlowerwall << endl;
        cout << "Ub = " << flags.uupperwall << endl;
        cout << "Wb = " << flags.wupperwall << endl;
        cout << "Ta = " << flags.tlowerwall << endl;
        cout << "Tb = " << flags.tupperwall << endl;
        cout << "Sa = " << flags.slowerwall << endl;
        cout << "Sb = " << flags.supperwall << endl;

        // Construct data fields: 3d velocity and 1d pressure
        cout << "Building velocity, temperature, salinity and pressure fields..." << flush;
        vector<FlowField> fields = {
            FlowField(Nx_, Ny_, Nz_, 3, Lx_, Lz_, ymin_, ymax_), // velocity
            FlowField(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ymin_, ymax_), // temperature
            FlowField(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ymin_, ymax_), // salinity
            FlowField(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ymin_, ymax_)};// pressure
        cout << "done" << endl;

        // Perturb velocity field
        cout << "Perturbing velocity, temperature, and salinity fields ... " << flush;
        if(modelabel=="yang2021jfm"){// define initial contions of Yang2021JFM's simulation
            addRandomPerturbations(fields[0],1e-3);
            addSinusoidalPerturbations(fields[1],-0.05,6.0);
            addSinusoidalPerturbations(fields[2],-0.05,6.0);
        }else{// mormal mode
            addRandomPerturbations(fields[0],1e-3);
            addRandomPerturbations(fields[1],1e-3);
            addRandomPerturbations(fields[2],1e-3);
        }
        cout << "done" << endl;

        // Construct Navier-Stoke integrator, set integration method
        cout << "Building DDC-DNS..." << flush;
        DDC ddc(fields, flags);
        cout << "done" << endl;

        int count=0;
        // Real cfl = ddc.CFL(fields[0]);
        for (Real t = flags.t0; t <= flags.T; t += dt.dT()) {
            cout << "         t == " << t << endl;
            cout << "       CFL == " << ddc.CFL(fields[0]) << endl;
            cout << " L2Norm(u) == " << L2Norm(fields[0]) << endl;
            cout << "divNorm(u) == " << divNorm(fields[0]) << endl;
            // cout << "      dPdx == " << ddc.dPdx() << endl;
            // cout << "     Ubulk == " << ddc.Ubulk() << endl;

            // Write velocity and modified pressure fields to disk
            if (P7!=0.0){
                fields[0].save(outdir + ulabel + i2s(int(count)));//<<--- save only fluctuations
                fields[1].save(outdir + tlabel + i2s(int(count)));
                fields[2].save(outdir + slabel + i2s(int(count)));
            }else{
                FlowField u_tot = totalVelocity(fields[0], flags); u_tot.save(outdir + ulabel + i2s(int(count)));//<<--- save total fields
                FlowField temp_tot = totalTemperature(fields[1], flags); temp_tot.save(outdir + tlabel + i2s(int(count)));
                FlowField salt_tot = totalSalinity(fields[2], flags); salt_tot.save(outdir + slabel + i2s(int(count)));
            }
            count+=1;

            

            // Take n steps of length dt
            ddc.advance(fields, dt.n());

            if (dt.variable() &&
                dt.adjust(ddc.CFL(fields[0])))  // TODO: dt.variable()==true is checked twice here, remove it.
                ddc.reset_dt(dt);
            cout << endl;
        }
    }
    cfMPI_Finalize();
}
