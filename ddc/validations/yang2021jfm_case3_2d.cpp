/**
 * Layering and vertical transport in sheared double-diffusive convection in the diffusive regime
    Y Yang, R Verzicco, D Lohse, CP Caulfield. Journal of Fluid Mechanics, 2022
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
        cout << "================================================================\n";
        cout << "This program integrates a plane Thermolhaline-shear flow from a random\n";
        cout << "initial condition. Velocity fields are saved at intervals dT=1.0\n";
        cout << "in a data/ directory, in channelflow's binary data file\n";
        cout << "format." << endl << endl;

        // Define gridsize
        const int Nx = 384;
        const int Ny = 385;
        const int Nz = 10;

        // Define box size
        const Real Lx = 2.0;
        const Real a = 0.0;
        const Real b = 1.0;
        const Real Lz = 0.02;

        // Define flow parameters
        const Real Pr = 10.0;
        const Real Rrho = 2.0;
        const Real Ra = 1e4;
        const Real Le = 100.0; // similar to \tau=1/Le=0.01
        

        

        fftw_loadwisdom();
        // Define DNS parameters
        DDCFlags flags;
        
        flags.Pr = Pr;
        flags.Rrho = Rrho;
        flags.Ra = Ra;
        flags.Le = Le;

        
        flags.uupperwall = 0.0;
        flags.ulowerwall = 0.0;
        flags.wupperwall = 0.0;
        flags.wlowerwall = 0.0;
        flags.tupperwall = 0.0;
        flags.tlowerwall = 1.0;
        flags.supperwall = 0.0;
        flags.slowerwall = 1.0;
        
        // flags.timestepping = SBDF3;
        // flags.initstepping = CNRK2;
        // flags.dealiasing = DealiasXZ;
        // flags.constraint = PressureGradient;
        // flags.taucorrection = true;
        flags.nonlinearity = Convection;
        flags.t0 = 0;
        flags.dt = 0.01;
        flags.T = 30;
        flags.dtmin = 1e-3;

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
        

        TimeStep dt(flags);

        

        // Construct data fields: 3d velocity and 1d pressure
        cout << "building velocity, temperature, salinity and pressure fields..." << flush;
        vector<FlowField> fields = {
            FlowField(Nx, Ny, Nz, 3, Lx, Lz, a, b), // velocity
            FlowField(Nx, Ny, Nz, 1, Lx, Lz, a, b), // temperature
            FlowField(Nx, Ny, Nz, 1, Lx, Lz, a, b), // salinity
            FlowField(Nx, Ny, Nz, 1, Lx, Lz, a, b)};// pressure
        cout << "done" << endl;


        // Define size and smoothness of initial disturbance
        Real spectralDecay = 0.5;
        Real magnitude = 0.05;
        int kxmax = 4;
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

        mkdir("yang2021jfm_case3_2d");
        Real cfl = ddc.CFL(fields[0]);
        for (Real t = flags.t0; t <= flags.T; t += dt.dT()) {
            cout << "         t == " << t << endl;
            cout << "       CFL == " << cfl << endl;
            cout << " L2Norm(u) == " << L2Norm(fields[0]) << endl;
            cout << "divNorm(u) == " << divNorm(fields[0]) << endl;
            // cout << "      dPdx == " << ddc.dPdx() << endl;
            // cout << "     Ubulk == " << ddc.Ubulk() << endl;

            // Write velocity and modified pressure fields to disk
            // fields[0].save("yang2021jfm_case3_2d/u" + i2s(int(t)));//<<--- save only fluctuations
            // fields[1].save("yang2021jfm_case3_2d/t" + i2s(int(t)));
            // fields[2].save("yang2021jfm_case3_2d/s" + i2s(int(t)));

            FlowField u_tot = totalVelocity(fields[0], flags); u_tot.save("yang2021jfm_case3_2d/u" + i2s(int(t)));//<<--- save total fields
            FlowField temp_tot = totalTemperature(fields[1], flags); temp_tot.save("yang2021jfm_case3_2d/t" + i2s(int(t)));
            FlowField salt_tot = totalSalinity(fields[2], flags); salt_tot.save("yang2021jfm_case3_2d/s" + i2s(int(t)));

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
