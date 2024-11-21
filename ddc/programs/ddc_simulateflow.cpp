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
#include "modules/ddc/ddcdsi.h"
#include "modules/ddc/ddc.h"
#include "modules/ddc/addPerturbations.h"
#include "modules/ddc/boundaryCondition.h"
#include "modules/ddc/turbulenceStatistics.h"
using namespace std;
using namespace chflow;

string printdiagnostics(FlowField& u, const DNS& dns, Real t, const TimeStep& dt, Real nu, Real umin, bool vardt,
                        bool pl2norm, bool pchnorm, bool pdissip, bool pshear, bool pdiverge, bool pUbulk, bool pubulk,
                        bool pdPdx, bool pcfl);

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

        const bool pcfl = args.getflag("-cfl", "--cfl", "print CFL number each dT");
        const bool pl2norm = args.getflag("-l2", "--l2norm", "print L2Norm(u) each dT");
        const bool pchnorm = args.getbool("-ch", "--chebyNorm", true, "print chebyNorm(u) each dT");
        const bool pdissip = args.getflag("-D", "--dissipation", "print dissipation each dT");
        const bool pshear = args.getflag("-I", "--input", "print wall shear power input each dT");
        const bool pdiverge = args.getflag("-dv", "--divergence", "print divergence each dT");
        const bool pubulk = args.getflag("-u", "--ubulk", "print ubulk each dT");
        const bool pUbulk = args.getflag("-Up", "--Ubulk-print", "print Ubulk each dT");
        const bool pdPdx = args.getflag("-p", "--pressure", "print pressure gradient each dT");
        const Real umin = args.getreal("-u", "--umin", 0.0, "stop if chebyNorm(u) < umin");

        const Real ecfmin = args.getreal("-e", "--ecfmin", 0.0, "stop if Ecf(u) < ecfmin");

        const int Nx_ = args.getint("-Nx", "--Nx", 32, "# x gridpoints");
        const int Ny_ = args.getint("-Ny", "--Ny", 31, "# y gridpoints");
        const int Nz_ = args.getint("-Nz", "--Nz", 32, "# z gridpoints");
        const Real Lx_ = args.getreal("-Lx", "--Lx", 1, "streamwise (x) box length");
        const Real Lz_ = args.getreal("-Lz", "--Lz", 1, "spanwise   (z) box length");
        const Real ymin_ = args.getreal("-ymin", "--ymin", 0.0, "lower wall height (y is wallnormal) ");
        const Real ymax_ = args.getreal("-ymax", "--ymax", 1.0, "upper wall height (y is wallnormal) ");

        const string outdir_profiles = args.getpath("-op", "--outdir_profiles", "profiles/", "output directory of mean profiles");
        const bool savetot = args.getflag("-savetot", "--savetotfields", "save total fields");
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
        cout << "Ri = " << flags.Ri << endl;
        cout << "Rrho = " << flags.Rrho << endl;
        #ifdef P7
        cout << "Rsep = " << flags.Rsep << endl;
        #endif
        cout << "Ua = " << flags.ulowerwall << endl;
        cout << "Wa = " << flags.wlowerwall << endl;
        cout << "Ub = " << flags.uupperwall << endl;
        cout << "Wb = " << flags.wupperwall << endl;
        #ifdef P5
        cout << "Ta = " << flags.tlowerwall << endl;
        cout << "Tb = " << flags.tupperwall << endl;
        #endif
        #ifdef P6
        cout << "Sa = " << flags.slowerwall << endl;
        cout << "Sb = " << flags.supperwall << endl;
        #endif

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
        if(modelabel=="yang2021jfm_case3"){// define initial contions of Yang2021JFM's simulation
            addRandomPerturbations(fields[0],1e-3);
            #ifdef P5
            addSinusoidalPerturbations(fields[1],-0.05,6.0);
            #endif
            #ifdef P6
            addSinusoidalPerturbations(fields[2],-0.05,6.0);
            #endif
        } else if(modelabel=="yang2021jfm_case5"){// define initial contions of Yang2021JFM's simulation
            addRandomPerturbations(fields[0],1e-3);
            #ifdef P5
            addSinusoidalPerturbations(fields[1],-0.025,8.0);
            #endif
            #ifdef P6
            addSinusoidalPerturbations(fields[2],-0.025,8.0);
            #endif
        } else if (modelabel=="eaves2016jfm"){
            addRandomPerturbations(fields[0],1e-4);
            #ifdef P5
            addSinusoidalPerturbations(fields[1],-0.05,6.0);
            #endif
            #ifdef P6
            addSinusoidalPerturbations(fields[2],-0.05,6.0);
            #endif
        } else{// mormal mode
            // addRandomPerturbations(fields[0],1e-3);
            // #ifdef P5
            // addRandomPerturbations(fields[1],1e-3);
            // #endif
            // #ifdef P6
            // addRandomPerturbations(fields[2],1e-3);
            // #endif

            addSinusoidalPerturbations(fields[0],-0.2,1.0);
            #ifdef P5
            addSinusoidalPerturbations(fields[1],-0.1,1.0);
            #endif
            #ifdef P6
            addSinusoidalPerturbations(fields[2],-0.1,1.0);
            #endif
        }
        cout << "done" << endl;

        // Construct Navier-Stoke integrator, set integration method
        cout << "Building DDC-DNS..." << flush;
        DDC ddc(fields, flags);
        cout << "done" << endl;

        ios::openmode openflag = (flags.t0 > 0) ? ios::app : ios::out;
        ofstream eout, x0out;
        openfile(eout, outdir + "energy.asc", openflag);
        eout << ddcfieldstatsheader_t("t", flags) << endl;
        
        #ifdef P6
        #ifdef SAVESTATS
        DDCTurbStats meanProfiles(Ny_);
        mkdir(outdir_profiles);
        #endif
        #endif

        int count=0;
        for (Real t = flags.t0; t <= flags.T; t += dt.dT()) {
            // cout << "         t == " << t << endl;
            // cout << "       CFL == " << ddc.CFL(fields[0]) << endl;
            // cout << " L2Norm(u) == " << L2Norm(fields[0]) << endl;
            // cout << "divNorm(u) == " << divNorm(fields[0]) << endl;
            // cout << "      dPdx == " << ddc.dPdx() << endl;
            // cout << "     Ubulk == " << ddc.Ubulk() << endl;

            string s;
            s = printdiagnostics(fields[0], ddc, t, dt, flags.nu, umin, dt.variable(), pl2norm, pchnorm, pdissip,
                                 pshear, pdiverge, pUbulk, pubulk, pdPdx, pcfl);
            if (ecfmin > 0 && Ecf(fields[0]) < ecfmin) {
                cferror("Ecf < ecfmin == " + r2s(ecfmin) + ", exiting");
            }

            cout << s;
            s = ddcfieldstats_t(fields[0], fields[1], fields[2], t, flags);
            eout << s << endl;
            #ifdef P6
            #ifdef SAVESTATS
            meanProfiles.addSnapshot(fields, flags);
            // save horizontially averaged fields
            meanProfiles.saveTurbStats(outdir_profiles + "meanprofile" + i2s(int(count)), fields);
            #endif
            #endif
            
            // Write velocity and modified pressure fields to disk
            if(!savetot){
                fields[0].save(outdir + ulabel + i2s(int(count)));//<<--- save only fluctuations
                #ifdef P5
                fields[1].save(outdir + tlabel + i2s(int(count)));
                #endif
                #ifdef P6
                fields[2].save(outdir + slabel + i2s(int(count)));
                #endif
            }else{
                FlowField u_tot = totalVelocity(fields[0], flags); u_tot.save(outdir + ulabel + i2s(int(count)));//<<--- save total fields
                #ifdef P5
                FlowField temp_tot = totalTemperature(fields[1], flags); temp_tot.save(outdir + tlabel + i2s(int(count)));
                #endif
                #ifdef P6
                FlowField salt_tot = totalSalinity(fields[2], flags); salt_tot.save(outdir + slabel + i2s(int(count)));
                #endif
            }
            count+=1;

            // Take n steps of length dt
            #ifndef FREESLIP
            ddc.advance(fields, dt.n());
            #endif
            #ifdef FREESLIP // this is not really used
            for (int step = 0; step < dt.n(); ++step) {
                cout << "." << flush;
                freeslipBC(fields);// add free-slip boundary conditions
                ddc.advance(fields, 1);
            }  // End of time stepping loop
            #endif

            if (dt.variable() &&
                dt.adjust(ddc.CFL(fields[0])))  // TODO: dt.variable()==true is checked twice here, remove it.
                ddc.reset_dt(dt);
            cout << endl;
        }
    }
    cfMPI_Finalize();
}


string printdiagnostics(FlowField& u, const DNS& dns, Real t, const TimeStep& dt, Real nu, Real umin, bool vardt,
                        bool pl2norm, bool pchnorm, bool pdissip, bool pshear, bool pdiverge, bool pUbulk, bool pubulk,
                        bool pdPdx, bool pcfl) {
    // Printing diagnostics
    stringstream sout;
    sout << "           t == " << t << endl;
    if (vardt)
        sout << "          dt == " << Real(dt) << endl;
    if (pl2norm)
        sout << "   L2Norm(u) == " << L2Norm(u) << endl;

    if (pchnorm || umin != 0.0) {
        Real chnorm = chebyNorm(u);
        sout << "chebyNorm(u) == " << chnorm << endl;
        if (chnorm < umin) {
            cout << "Exiting: chebyNorm(u) < umin." << endl;
            exit(0);
        }
    }
    Real h = 0.5 * (u.b() - u.a());
    u += dns.Ubase();
    if (pl2norm)
        sout << "   energy(u+U) == " << 0.5 * L2Norm(u) << endl;
    if (pdissip)
        sout << "   dissip(u+U) == " << dissipation(u) << endl;
    if (pshear)
        sout << "wallshear(u+U) == " << abs(wallshearLower(u)) + abs(wallshearUpper(u)) << endl;
    if (pdiverge)
        sout << "  divNorm(u+U) == " << divNorm(u) << endl;
    if (pUbulk)
        sout << "mean u+U Ubulk == " << dns.Ubulk() << endl;
    u -= dns.Ubase();
    if (u.taskid() == u.task_coeff(0, 0)) {
        if (pubulk)
            sout << "         ubulk == " << Re(u.profile(0, 0, 0)).mean() << endl;
    }
    if (pdPdx)
        sout << "          dPdx == " << dns.dPdx() << endl;
    if (pl2norm)
        sout << "     L2Norm(u) == " << L2Norm(u) << endl;
    if (pl2norm)
        sout << "   L2Norm3d(u) == " << L2Norm3d(u) << endl;

    Real cfl = dns.CFL(u);
    if (u.taskid() == u.task_coeff(0, 0)) {
        ChebyCoeff U = dns.Ubase();
        ChebyCoeff W = dns.Wbase();

        U.makeSpectral();
        U += Re(u.profile(0, 0, 0));
        Real Ucenter = U.eval(0.5 * (u.a() + u.b()));
        Real Uwall = pythag(0.5 * (U.eval_b() - U.eval_a()), 0.5 * (W.eval_b() - W.eval_a()));
        Real Umean = U.mean();
        // sout << "        1/nu == " << 1 / nu << endl;
        sout << "  Uwall h == " << Uwall * h << endl;
        sout << "  Ubulk h == " << dns.Ubulk() * h << endl;
        sout << "  Umean h == " << Umean << " * " << h << endl;
        sout << "  Umean h == " << Umean * h << endl;
        sout << " Uparab h == " << 1.5 * dns.Ubulk() * h << endl;
        sout << "Ucenter h == " << Ucenter * h << endl;
    }
    sout << "         CFL == " << cfl << endl;
    return sout.str();
}