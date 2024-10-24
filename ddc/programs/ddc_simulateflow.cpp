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
        const int Nx_ = args.getint("-Nx", "--Nx", 32, "# x gridpoints");
        const int Ny_ = args.getint("-Ny", "--Ny", 31, "# y gridpoints");
        const int Nz_ = args.getint("-Nz", "--Nz", 32, "# z gridpoints");
        const Real Lx_ = args.getreal("-Lx", "--Lx", 1, "streamwise (x) box length");
        const Real Lz_ = args.getreal("-Lz", "--Lz", 1, "spanwise   (z) box length");
        const Real ymin_ = args.getreal("-ymin", "--ymin", 0.0, "lower wall height (y is wallnormal) ");
        const Real ymax_ = args.getreal("-ymax", "--ymax", 1.0, "upper wall height (y is wallnormal) ");
       
        args.check();
        args.save("./");
        

        TimeStep dt(flags);

        fftw_loadwisdom();
       
        

        cout << "Parameters: " << endl;
        cout << "Ra = " << flags.Ra << endl;
        cout << "Le = " << flags.Le << endl;
        cout << "Pr = " << flags.Pr << endl;
        cout << "Rrho = " << flags.Rrho << endl;
        cout << "Ua = " << flags.ulowerwall << endl;
        cout << "Wa = " << flags.wlowerwall << endl;
        cout << "Ub = " << flags.uupperwall << endl;
        cout << "Wb = " << flags.wupperwall << endl;
        cout << "Ta = " << flags.tlowerwall << endl;
        cout << "Tb = " << flags.tupperwall << endl;
        cout << "Sa = " << flags.slowerwall << endl;
        cout << "Sb = " << flags.supperwall << endl;

        // Construct data fields: 3d velocity and 1d pressure
        cout << "building velocity, temperature, salinity and pressure fields..." << flush;
        vector<FlowField> fields = {
            FlowField(Nx_, Ny_, Nz_, 3, Lx_, Lz_, ymin_, ymax_), // velocity
            FlowField(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ymin_, ymax_), // temperature
            FlowField(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ymin_, ymax_), // salinity
            FlowField(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ymin_, ymax_)};// pressure
        cout << "done" << endl;


        // Define size and smoothness of initial disturbance
        Real spectralDecay = 0.5;
        Real magnitude = 0.05;
        int kxmax = 3;
        int kzmax = 3;
        // Perturb velocity field
        cout << "Perturbing velocity field..." << flush;
        fields[0].addPerturbations(kxmax, kzmax, 1.0, spectralDecay);
        fields[0] *= magnitude / L2Norm(fields[0]);
        fields[1].addPerturbations(kxmax, kzmax, 1.0, spectralDecay);
        fields[1] *= magnitude / L2Norm(fields[1]);
        fields[2].addPerturbations(kxmax, kzmax, 1.0, spectralDecay);
        fields[2] *= magnitude / L2Norm(fields[2]);
        cout << "done" << endl;

        // Construct Navier-Stoke integrator, set integration method
        cout << "Building DDC-DNS..." << flush;
        DDC ddc(fields, flags);
        cout << "done" << endl;

        mkdir("data");
        Real cfl = ddc.CFL(fields[0]);
        for (Real t = flags.t0; t <= flags.T; t += dt.dT()) {
            cout << "         t == " << t << endl;
            cout << "       CFL == " << cfl << endl;
            cout << " L2Norm(u) == " << L2Norm(fields[0]) << endl;
            cout << "divNorm(u) == " << divNorm(fields[0]) << endl;
            // cout << "      dPdx == " << ddc.dPdx() << endl;
            // cout << "     Ubulk == " << ddc.Ubulk() << endl;

            // Write velocity and modified pressure fields to disk
           // fields[0].save("data/u" + i2s(int(t)));//<<--- save only fluctuations
            // fields[1].save("data/t" + i2s(int(t)));
            // fields[2].save("data/s" + i2s(int(t)));

            FlowField u_tot = totalVelocity(fields[0], flags); u_tot.save("data/u" + i2s(int(t)));//<<--- save total fields
            FlowField temp_tot = totalTemperature(fields[1], flags); temp_tot.save("data/t" + i2s(int(t)));
            FlowField salt_tot = totalSalinity(fields[2], flags); salt_tot.save("data/s" + i2s(int(t)));

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
