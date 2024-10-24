/**
 * Dynamical System Interface for the DDC module
 *
 * Original author: Duc Nguyen
 */


#include "modules/ddc/ddcdsi.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "cfbasics/mathdefs.h"
#include "channelflow/diffops.h"
// #include "viscoelastic/veutils.h"
#include "modules/ddc/ddc.h"

using namespace std;

namespace chflow {

/*utility functions*/

std::vector<Real> ddcstats(const FlowField& u, const FlowField& temp, const FlowField& salt, const DDCFlags flags) {
    double l2n = L2Norm(u);
    if (std::isnan(l2n)) {
        cferror("L2Norm(u) is nan");
    }

    FlowField u_tot = totalVelocity(u, flags);
    FlowField temp_tot = totalTemperature(temp, flags);
    FlowField salt_tot = totalSalinity(salt, flags);

    std::vector<Real> stats;
    stats.push_back(L2Norm(u));
    stats.push_back(heatcontent(temp_tot, flags));// averaged temp
    stats.push_back(saltcontent(salt_tot, flags));// averaged salt
    stats.push_back(L2Norm(u_tot));
    stats.push_back(L2Norm(temp));
    stats.push_back(L2Norm(temp_tot));
    stats.push_back(L2Norm(salt));
    stats.push_back(L2Norm(salt_tot));
    stats.push_back(L2Norm3d(u));
    stats.push_back(Ecf(u));
    stats.push_back(wallshear(u_tot));
    return stats;
}

string ddcfieldstatsheader(const DDCFlags flags) {
    stringstream header;
    header << setw(14) << "L2(u')" << setw(10) << "<T>(y=" << flags.ystats << ")"  // change position with L2(T')
           << setw(14) << "L2(u)" << setw(14) << "L2(T)" << setw(14) << "e3d" << setw(14) << "ecf" << setw(14)
           << "ubulk" << setw(14) << "wbulk" << setw(14) << "wallshear" << setw(14) << "buoyPowIn" << setw(14)
           << "totalDiss" << setw(14) << "heatinflux" << setw(14) << "L2(T')" << setw(10) << "Nu(y=" << flags.ystats
           << ")";
    return header.str();
}

string ddcfieldstatsheader_t(const string tname, const DDCFlags flags) {
    stringstream header;
    header << setw(6) << "#(" << tname << ")" << ddcfieldstatsheader(flags);
    return header.str();
}

string ddcfieldstats(const FlowField& u, const FlowField& temp, const FlowField& salt, const DDCFlags flags) {
    std::vector<Real> stats = ddcstats(u, temp, salt, flags);
    // Return string
    stringstream s;
    for (uint i = 0; i < stats.size(); i++) {
        s << setw(14) << stats[i];
    }
    return s.str();
}

string ddcfieldstats_t(const FlowField& u, const FlowField& temp, const FlowField& salt, const Real t, const DDCFlags flags) {
    std::vector<Real> stats = ddcstats(u, temp, salt, flags);
    // Return string
    stringstream s;
    s << setw(8) << t;
    for (uint i = 0; i < stats.size(); i++) {
        s << setw(14) << stats[i];
    }
    return s.str();
}

FlowField totalVelocity(const FlowField& velo, const DDCFlags flags) {
    // copy
    FlowField u(velo);
    FlowField tmp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());

    // get base flow
    DDC ddc({u, tmp, tmp, tmp}, flags);
    ChebyCoeff Ubase = ddc.Ubase();
    ChebyCoeff Wbase = ddc.Wbase();

    // add base flow (should be identical to code in function temperatureNL in DDE
    for (int ny = 0; ny < u.Ny(); ++ny) {
        if (u.taskid() == u.task_coeff(0, 0)) {
            u.cmplx(0, ny, 0, 0) += Complex(Ubase(ny), 0.0);
            u.cmplx(0, ny, 0, 2) += Complex(Wbase(ny), 0.0);
        }
    }
    if (u.taskid() == u.task_coeff(0, 0)) {
        u.cmplx(0, 0, 0, 1) -= Complex(flags.Vsuck, 0.);
    }

    return u;
}

FlowField totalTemperature(const FlowField& temp, const DDCFlags flags) {
    // copy
    FlowField T(temp);
    FlowField tmp(T.Nx(), T.Ny(), T.Nz(), 3, T.Lx(), T.Lz(), T.a(), T.b(), T.cfmpi());

    // get base flow
    DDC ddc({tmp, T, T, T}, flags);
    ChebyCoeff Tbase = ddc.Tbase();

    // add base flow (should be identical to code in function temperatureNL in DDE
    for (int ny = 0; ny < T.Ny(); ++ny) {
        if (T.taskid() == T.task_coeff(0, 0))
            T.cmplx(0, ny, 0, 0) += Complex(Tbase(ny), 0.0);
    }

    return T;
}

FlowField totalSalinity(const FlowField& salt, const DDCFlags flags) {
    // copy
    FlowField S(salt);
    FlowField tmp(S.Nx(), S.Ny(), S.Nz(), 3, S.Lx(), S.Lz(), S.a(), S.b(), S.cfmpi());

    // get base flow
    DDC ddc({tmp, S, S, S}, flags);
    ChebyCoeff Sbase = ddc.Sbase();

    // add base flow (should be identical to code in function temperatureNL in DDE
    for (int ny = 0; ny < S.Ny(); ++ny) {
        if (S.taskid() == S.task_coeff(0, 0))
            S.cmplx(0, ny, 0, 0) += Complex(Sbase(ny), 0.0);
    }

    return S;
}


Real heatcontent(const FlowField& ttot, const DDCFlags flags) {
    assert(ttot.ystate() == Spectral);
    int N = 100;
    Real dy = (flags.ystats - ttot.a()) / (N - 1);
    Real avt = 0;  // average temperature
    if (ttot.taskid() == ttot.task_coeff(0, 0)) {
        ChebyCoeff tprof = Re(ttot.profile(0, 0, 0));
        for (int i = 0; i < N; ++i) {
            Real y = ttot.a() + i * dy;
            avt += tprof.eval(y);
        }
        avt *= 1.0 / N;
    }
#ifdef HAVE_MPI
    MPI_Bcast(&avt, 1, MPI_DOUBLE, ttot.task_coeff(0, 0), ttot.cfmpi()->comm_world);
#endif
    return avt;
}
Real saltcontent(const FlowField& stot, const DDCFlags flags) {
    assert(stot.ystate() == Spectral);
    int N = 100;
    Real dy = (flags.ystats - stot.a()) / (N - 1);
    Real avt = 0;  // average temperature
    if (stot.taskid() == stot.task_coeff(0, 0)) {
        ChebyCoeff sprof = Re(stot.profile(0, 0, 0));
        for (int i = 0; i < N; ++i) {
            Real y = stot.a() + i * dy;
            avt += sprof.eval(y);
        }
        avt *= 1.0 / N;
    }
#ifdef HAVE_MPI
    MPI_Bcast(&avt, 1, MPI_DOUBLE, stot.task_coeff(0, 0), stot.cfmpi()->comm_world);
#endif
    return avt;
}


/* Begin of ddcDSI class*/

ddcDSI::ddcDSI() {}

ddcDSI::ddcDSI(DDCFlags& ddcflags, FieldSymmetry sigma, PoincareCondition* h, TimeStep dt, bool Tsearch, bool xrelative,
               bool zrelative, bool Tnormalize, Real Unormalize, const FlowField& u, const FlowField& temp, const FlowField& salt, ostream* os)
    : cfDSI(ddcflags, sigma, h, dt, Tsearch, xrelative, zrelative, Tnormalize, Unormalize, u, os),
      ddcflags_(ddcflags) {}

Eigen::VectorXd ddcDSI::eval(const Eigen::VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);

    Real T;
    extractVectorDDC(x, u, temp, salt, sigma_, T);

    FlowField Gu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gtemp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gsalt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    G(u, temp, salt, T, h_, sigma_, Gu, Gtemp, Gsalt, ddcflags_, dt_, Tnormalize_, Unormalize_, fcount_, CFL_, *os_);
    Eigen::VectorXd Gx(Eigen::VectorXd::Zero(x.rows()));
    //   Galpha *= 1./vednsflags_.b_para;
    field2vector(Gu, Gtemp, Gsalt, Gx);  // This does not change the size of Gx and automatically leaves the last entries zero

    return Gx;
}

Eigen::VectorXd ddcDSI::eval(const Eigen::VectorXd& x0, const Eigen::VectorXd& x1, bool symopt) {
    FlowField u0(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField u1(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp0(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp1(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt0(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt1(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    Real T0, T1;
    FieldSymmetry sigma0, sigma1;
    extractVectorDDC(x0, u0, temp0, salt0, sigma0, T0);
    extractVectorDDC(x1, u1, temp1, salt1, sigma1, T1);

    FlowField Gu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gtemp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gsalt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);

    f(u0, temp0, salt0, T0, h_, Gu, Gtemp, Gsalt, ddcflags_, dt_, fcount_, CFL_, *os_);
    if (symopt) {
        Gu *= sigma0;
        if (sigma0.sy() == -1) {
            // wall-normal mirroring in velocity requires sign change in temperature
            FieldSymmetry inv(-1);
            sigma0 *= inv;
        }
        Gtemp *= sigma0;
        Gsalt *= sigma0;
    }
    Gu -= u1;
    Gtemp -= temp1;
    Gsalt -= salt1;

    // normalize
    if (Tnormalize_) {
        Gu *= 1.0 / T0;
        Gtemp *= 1.0 / T0;
        Gsalt *= 1.0 / T0;
    }
    if (Unormalize_ != 0.0) {
        Real funorm = L2Norm3d(Gu);
        Gu *= 1. / sqrt(abs(funorm * (Unormalize_ - funorm)));
        // u should stay off zero, so normalize with u for now - temp should also stay away from zero
        Gtemp *= 1. / sqrt(abs(funorm * (Unormalize_ - funorm)));
        Gsalt *= 1. / sqrt(abs(funorm * (Unormalize_ - funorm)));
    }

    Eigen::VectorXd Gx(Eigen::VectorXd::Zero(x0.rows()));
    field2vector(Gu, Gtemp, Gsalt, Gx);  // This does not change the size of Gx and automatically leaves the last entries zero

    return Gx;
}

void ddcDSI::save(const Eigen::VectorXd& x, const string filebase, const string outdir, const bool fieldsonly) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);

    u.save(outdir + "u" + filebase);
    temp.save(outdir + "t" + filebase);
    salt.save(outdir + "s" + filebase);

    if (!fieldsonly) {
        string fs = ddcfieldstats(u, temp, salt, ddcflags_);
        if (u.taskid() == 0) {
            if (xrelative_ || zrelative_ || !sigma.isIdentity())
                sigma.save(outdir + "sigma" + filebase);
            if (Tsearch_)
                chflow::save(T, outdir + "T" + filebase);
            // sigma.save (outdir+"sigmaconverge.asc", ios::app);
            ofstream fout((outdir + "fieldconverge.asc").c_str(), ios::app);
            long pos = fout.tellp();
            if (pos == 0)
                fout << ddcfieldstatsheader() << endl;
            fout << fs << endl;
            fout.close();
            ddcflags_.save(outdir);
        }
    }
}

string ddcDSI::stats(const Eigen::VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);
    return ddcfieldstats_t(u, temp, salt, mu_, ddcflags_);
}

pair<string, string> ddcDSI::stats_minmax(const Eigen::VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gu(u);
    FlowField Gtemp(temp);
    FlowField Gsalt(salt);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);

    std::vector<Real> stats = ddcstats(u, temp, salt, ddcflags_);
    std::vector<Real> minstats(stats);
    std::vector<Real> maxstats(stats);

    // quick hack to avoid new interface or creating simple f() for DDC
    TimeStep dt = TimeStep(ddcflags_.dt, 0, 1, 1, 0, 0, false);
    int fcount = 0;
    PoincareCondition* h = 0;
    Real CFL = 0.0;
    std::ostream muted_os(0);
    Real timep = T / 100.0;

    *os_ << "Using flag -orbOut: Calculate minmax-statistics of periodic orbit." << endl;
    for (int t = 0; t < 100; t++) {
        f(u, temp, salt, timep, h, Gu, Gtemp, Gsalt, ddcflags_, dt, fcount, CFL, muted_os);
        stats = ddcstats(Gu, Gtemp, Gsalt, ddcflags_);
        for (uint i = 0; i < stats.size(); i++) {
            minstats[i] = (minstats[i] < stats[i]) ? minstats[i] : stats[i];
            maxstats[i] = (maxstats[i] > stats[i]) ? maxstats[i] : stats[i];
        }
        u = Gu;
        temp = Gtemp;
        salt = Gsalt;
    }
    // Return string
    stringstream smin;
    stringstream smax;
    smin << setw(8) << mu_;
    smax << setw(8) << mu_;
    for (uint i = 0; i < stats.size(); i++) {
        smin << setw(14) << minstats[i];
        smax << setw(14) << maxstats[i];
    }

    pair<string, string> minmax;
    minmax = make_pair(smin.str(), smax.str());
    return minmax;
}

string ddcDSI::statsHeader() { return ddcfieldstatsheader_t(ddc_cPar2s(ddc_cPar_), ddcflags_); }

/// after finding new solution fix phases
void ddcDSI::phaseShift(Eigen::VectorXd& x) {
    if (xphasehack_ || zphasehack_) {
        FlowField unew(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
        FlowField tnew(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
        FlowField snew(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
        FieldSymmetry sigma;
        Real T;
        extractVectorDDC(x, unew, tnew, snew, sigma, T);
        // vector2field (x,unew);
        const int phasehackcoord = 0;  // Those values were fixed in continuesoln anyway
        const parity phasehackparity = Odd;
        const Real phasehackguess = 0.0;

        if (zphasehack_) {
            FieldSymmetry tau = zfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
            cout << "fixing z phase of potential solution with phase shift tau == " << tau << endl;
            unew *= tau;
            tnew *= tau;
            snew *= tau;
        }
        if (xphasehack_) {
            FieldSymmetry tau = xfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
            cout << "fixing x phase of potential solution with phase shift tau == " << tau << endl;
            unew *= tau;
            tnew *= tau;
            snew *= tau;
        }
        if (uUbasehack_) {
            cout << "fixing u+Ubase decomposition so that <du/dy> = 0 at walls (i.e. Ubase balances mean pressure "
                    "gradient))"
                 << endl;
            Real ubulk = Re(unew.profile(0, 0, 0)).mean();
            if (abs(ubulk) < 1e-15)
                ubulk = 0.0;

            ChebyCoeff Ubase = laminarProfile(ddcflags_.nu, ddcflags_.constraint, ddcflags_.dPdx,
                                              ddcflags_.Ubulk - ubulk, ddcflags_.Vsuck, unew.a(), unew.b(),
                                              ddcflags_.ulowerwall, ddcflags_.uupperwall, unew.Ny());

            fixuUbasehack(unew, Ubase);
        }
        makeVectorDDC(unew, tnew, snew, sigma, T, x);
    }
}

void ddcDSI::phaseShift(Eigen::MatrixXd& y) {
    if (xphasehack_ || zphasehack_) {
        FlowField unew(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
        FlowField tnew(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
        FlowField snew(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
        Eigen::VectorXd yvec;
        FieldSymmetry sigma;
        Real T;

        const int phasehackcoord = 0;  // Those values were fixed in continuesoln anyway
        const parity phasehackparity = Odd;
        const Real phasehackguess = 0.0;

        FieldSymmetry taux(0.0, 0.0);
        FieldSymmetry tauz(0.0, 0.0);

        extractVectorDDC(y.col(0), unew, tnew, snew, sigma, T);

        if (xphasehack_) {
            taux = xfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
            cout << "fixing x phase of potential solution with phase shift tau == " << taux << endl;
        }
        if (zphasehack_) {
            tauz = zfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
            cout << "fixing z phase of potential solution with phase shift tau == " << tauz << endl;
        }

        for (int i = 0; i < y.cols(); i++) {
            extractVectorDDC(y.col(i), unew, tnew, snew, sigma, T);
            unew *= taux;
            tnew *= taux;
            snew *= taux;
            unew *= tauz;
            tnew *= tauz;
            snew *= tauz;
            makeVectorDDC(unew, tnew, snew, sigma, T, yvec);
            y.col(i) = yvec;
        }
    }
}

Real ddcDSI::extractT(const Eigen::VectorXd& x) {  // inefficient hack
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);
    return T;
}

Real ddcDSI::extractXshift(const Eigen::VectorXd& x) {  // inefficient hack
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);
    return sigma.ax();
}

Real ddcDSI::extractZshift(const Eigen::VectorXd& x) {  // inefficient hack
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);
    return sigma.az();
}

void ddcDSI::makeVectorDDC(const FlowField& u, const FlowField& temp, const FlowField& salt, const FieldSymmetry& sigma, const Real T,
                           Eigen::VectorXd& x) {
    if (u.Nd() != 3)
        cferror("ddcDSI::makeVector(): u.Nd() = " + i2s(u.Nd()) + " != 3");
    if (temp.Nd() != 1)
        cferror("ddcDSI::makeVector(): temp.Nd() = " + i2s(temp.Nd()) + " != 1");
    if (salt.Nd() != 1)
        cferror("ddcDSI::makeVector(): salt.Nd() = " + i2s(salt.Nd()) + " != 1");
    int taskid = u.taskid();

    int uunk = field2vector_size(u, temp, salt);                   // # of variables for u and alpha unknonwn
    const int Tunk = (Tsearch_ && taskid == 0) ? uunk : -1;  // index for T unknown
    const int xunk = (xrelative_ && taskid == 0) ? uunk + Tsearch_ : -1;
    const int zunk = (zrelative_ && taskid == 0) ? uunk + Tsearch_ + xrelative_ : -1;
    int Nunk = (taskid == 0) ? uunk + Tsearch_ + xrelative_ + zrelative_ : uunk;
    if (x.rows() < Nunk)
        x.resize(Nunk);
    field2vector(u, temp, salt, x);
    if (taskid == 0) {
        if (Tsearch_)
            x(Tunk) = T;
        if (xrelative_)
            x(xunk) = sigma.ax();
        if (zrelative_)
            x(zunk) = sigma.az();
    }
}

void ddcDSI::extractVectorDDC(const Eigen::VectorXd& x, FlowField& u, FlowField& temp, FlowField& salt, FieldSymmetry& sigma, Real& T) {
    int uunk = field2vector_size(u, temp, salt);  // number of components in x that corresond to u and alpha
    vector2field(x, u, temp, salt);
    const int Tunk = uunk + Tsearch_ - 1;
    const int xunk = uunk + Tsearch_ + xrelative_ - 1;
    const int zunk = uunk + Tsearch_ + xrelative_ + zrelative_ - 1;
    Real ax = 0;
    Real az = 0;
    if (u.taskid() == 0) {
        T = Tsearch_ ? x(Tunk) : Tinit_;
        ax = xrelative_ ? x(xunk) : axinit_;
        az = zrelative_ ? x(zunk) : azinit_;
    }
#ifdef HAVE_MPI
    MPI_Bcast(&T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&az, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    sigma = FieldSymmetry(sigma_.sx(), sigma_.sy(), sigma_.sz(), ax, az, sigma_.s());
}

Eigen::VectorXd ddcDSI::xdiff(const Eigen::VectorXd& a) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(a, u, temp, salt);
    Eigen::VectorXd dadx(a.size());
    dadx.setZero();
    u = chflow::xdiff(u);
    temp = chflow::xdiff(temp);
    salt = chflow::xdiff(salt);
    field2vector(u, temp, salt, dadx);
    dadx *= 1. / L2Norm(dadx);
    return dadx;
}

Eigen::VectorXd ddcDSI::zdiff(const Eigen::VectorXd& a) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(a, u, temp, salt);
    Eigen::VectorXd dadz(a.size());
    dadz.setZero();
    u = chflow::zdiff(u);
    temp = chflow::zdiff(temp);
    salt = chflow::zdiff(salt);
    field2vector(u, temp, salt, dadz);
    dadz *= 1. / L2Norm(dadz);
    return dadz;
}

Eigen::VectorXd ddcDSI::tdiff(const Eigen::VectorXd& a, Real epsDt) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    // quick hack to avoid new interface or creating simple f() for DDC
    TimeStep dt = TimeStep(epsDt, 0, 1, 1, 0, 0, false);
    int fcount = 0;
    PoincareCondition* h = 0;
    Real CFL = 0.0;
    FlowField edudtf(u);
    FlowField edtempdtf(temp);
    FlowField edsaltdtf(salt);
    std::ostream muted_os(0);
    //   vector2field (a, u, temp);
    extractVectorDDC(a, u, temp, salt, sigma, T);
    // use existing f() instead of simple
    f(u, temp, salt, epsDt, h, edudtf, edtempdtf, edsaltdtf, ddcflags_, dt, fcount, CFL, muted_os);
    //   f (temp, 1,epsDt, edtempdtf, dnsflags_, *os_);
    edudtf -= u;
    edtempdtf -= temp;
    edsaltdtf -= salt;
    Eigen::VectorXd dadt(a.size());
    field2vector(edudtf, edtempdtf, edsaltdtf, dadt);
    dadt *= 1. / L2Norm(dadt);
    return dadt;
}

void ddcDSI::updateMu(Real mu) {
    DSI::updateMu(mu);
    if (ddc_cPar_ == ddc_continuationParameter::none) {
        cfDSI::updateMu(mu);
    }else if (ddc_cPar_ == ddc_continuationParameter::Lx) {
        Lx_ = mu;
    } else if (ddc_cPar_ == ddc_continuationParameter::Lz) {
        Lz_ = mu;
    }else {
        throw invalid_argument("ddcDSI::updateMu(): continuation parameter is unknown");
    }
}

void ddcDSI::chooseMuDDC(string muName) {
    ddc_continuationParameter ddc_cPar = s2ddc_cPar(muName);

    if (ddc_cPar == ddc_continuationParameter::none)
        ddcDSI::chooseMu(muName);
    else
        chooseMuDDC(ddc_cPar);
}

void ddcDSI::chooseMuDDC(ddc_continuationParameter mu) {
    switch (mu) {
        case ddc_continuationParameter::Lx:
            updateMu(Lx_);
            break;
        case ddc_continuationParameter::Lz:
            updateMu(Lz_);
            break;
        case ddc_continuationParameter::none:
            throw invalid_argument(
                "ddcDSI::chooseMu(): continuation parameter is none, we should not reach this point");
        default:
            throw invalid_argument("ddcDSI::chooseMu(): continuation parameter is unknown");
    }
}

ddc_continuationParameter ddcDSI::s2ddc_cPar(string muname) {
    std::transform(muname.begin(), muname.end(), muname.begin(), ::tolower);  // why is the string made lower case?
    if (muname == "lx")
        return ddc_continuationParameter::Lx;
    else if (muname == "lz")
        return ddc_continuationParameter::Lz;
    else
        // cout << "ddcDSI::s2ddc_cPar(): ddc_continuation parameter '"+muname+"' is unknown, defaults to 'none'" <<
        // endl;
        return ddc_continuationParameter::none;
}

string ddcDSI::printMu() { return ddc_cPar2s(ddc_cPar_); }

string ddcDSI::ddc_cPar2s(ddc_continuationParameter ddc_cPar) {
    if (ddc_cPar == ddc_continuationParameter::none)
        return cfDSI::cPar2s(cPar_);
    else if (ddc_cPar == ddc_continuationParameter::Lx)
        return "Lx";
    else if (ddc_cPar == ddc_continuationParameter::Lz)
        return "Lz";
    else
        throw invalid_argument("ddcDSI::ddc_cPar2s(): continuation parameter is not convertible to string");
}

void ddcDSI::saveParameters(string searchdir) {
    // cfDSI::saveParameters (searchdir);
    ddcflags_.save(searchdir);
}

void ddcDSI::saveEigenvec(const Eigen::VectorXd& ev, const string label, const string outdir) {
    FlowField efu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField eft(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efs(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(ev, efu, eft, efs);
    efu *= 1.0 / L2Norm(efu);
    eft *= 1.0 / L2Norm(eft);
    efs *= 1.0 / L2Norm(efs);
    efu.save(outdir + "efu" + label);
    eft.save(outdir + "eft" + label);
    efs.save(outdir + "efs" + label);
}

void ddcDSI::saveEigenvec(const Eigen::VectorXd& evA, const Eigen::VectorXd& evB, const string label1,
                          const string label2, const string outdir) {
    FlowField efAu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efBu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efAt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efBt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efAs(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efBs(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(evA, efAu, efAt, efAs);
    vector2field(evB, efBu, efBt, efBs);
    Real cu = 1.0 / sqrt(L2Norm2(efAu) + L2Norm2(efBu));
    Real ct = 1.0 / sqrt(L2Norm2(efAt) + L2Norm2(efBt));
    Real cs = 1.0 / sqrt(L2Norm2(efAs) + L2Norm2(efBs));
    efAu *= cu;
    efBu *= cu;
    efAt *= ct;
    efBt *= ct;
    efAs *= cs;
    efBs *= cs;
    efAu.save(outdir + "efu" + label1);
    efBu.save(outdir + "efu" + label2);
    efAt.save(outdir + "eft" + label1);
    efBt.save(outdir + "eft" + label2);
    efAs.save(outdir + "efs" + label1);
    efBs.save(outdir + "efs" + label2);
}

/* OUTSIDE CLASS */

// G(x) = G(u,sigma) = (sigma f^T(u) - u) for orbits
void G(const FlowField& u, const FlowField& temp, const FlowField& salt, Real& T, PoincareCondition* h, const FieldSymmetry& sigma,
       FlowField& Gu, FlowField& Gtemp, FlowField& Gsalt, const DDCFlags& ddcflags, const TimeStep& dt, bool Tnormalize, Real Unormalize,
       int& fcount, Real& CFL, ostream& os) {
    f(u, temp, salt, T, h, Gu, Gtemp, Gsalt, ddcflags, dt, fcount, CFL, os);
    Real funorm = L2Norm3d(Gu);
    Gu *= sigma;
    Gu -= u;
    if (sigma.sy() == -1) {
        // wall-normal mirroring in velocity requires sign change in temperature and salinity
        FieldSymmetry tsigma(sigma);
        FieldSymmetry ssigma(sigma);
        FieldSymmetry inv(-1);
        tsigma *= inv;
        ssigma *= inv;
        Gtemp *= tsigma;
        Gsalt *= ssigma;
    } else {
        Gtemp *= sigma;
        Gsalt *= sigma;
    }
    Gtemp -= temp;
    Gsalt -= salt;

    if (Tnormalize) {
        Gu *= 1.0 / T;
        Gtemp *= 1.0 / T;
        Gsalt *= 1.0 / T;
    }
    if (Unormalize != 0.0) {
        Gu *= 1. / sqrt(abs(funorm * (Unormalize - funorm)));
        // u should stay off zero, so normalize with u for now - temp should also stay away from zero
        Gtemp *= 1. / sqrt(abs(funorm * (Unormalize - funorm)));
        Gsalt *= 1. / sqrt(abs(funorm * (Unormalize - funorm)));
    }
}

void f(const FlowField& u, const FlowField& temp, const FlowField& salt, Real& T, PoincareCondition* h, FlowField& f_u, FlowField& f_temp, FlowField& f_salt,
       const DDCFlags& ddcflags_, const TimeStep& dt_, int& fcount, Real& CFL, ostream& os) {
    if (!isfinite(L2Norm(u))) {
        os << "error in f: u is not finite. exiting." << endl;
        exit(1);
    }
    DDCFlags flags(ddcflags_);
    flags.logstream = &os;
    TimeStep dt(dt_);
    vector<FlowField> fields = {u, temp, salt, FlowField(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi())};

    //   f_u = u;
    //   f_temp = temp;
    //   f_salt = salt;
    // No Poincare section, just integration to time T
    if (h == 0) {
        if (T < 0) {
            os << "f: negative integration time T == " << T << endl
               << "returning f(u,T) == (1+abs(T))*u" << endl
               << "returning f(temp,T) == (1+abs(T))*temp" << endl
               << "returning f(salt,T) == (1+abs(T))*salt" << endl;
            fields[0] *= 1 + abs(T);
            fields[1] *= 1 + abs(T);
            fields[2] *= 1 + abs(T);
            return;
        }
        // Special case #1: no time steps
        if (T == 0) {
            os << "f: T==0, no integration, returning u, temp, and salt" << endl;
            return;
        }
        dt.adjust_for_T(T, false);
        flags.dt = dt;
        // Adjust dt for CFL if necessary
        DDC ddc(fields, flags);
        ddc.advance(fields, 1);
        if (dt.variable()) {
            dt.adjust(ddc.CFL(fields[0]), false);
            ddc.reset_dt(dt);
        }
        //  t == current time in integration
        //  T == total integration time
        // dT == CFL-check interval
        // dt == DNS time-step
        //  N == T/dT,  n == dT/dt;
        //  T == N dT, dT == n dt
        //  t == s dT (s is loop index)

        os << "f^T: " << flush;
        for (int s = 1; s <= dt.N(); ++s) {
            Real t = s * dt.dT();
            CFL = ddc.CFL(fields[0]);
            if (s % 10 == 0)
                os << iround(t) << flush;
            else if (s % 2 == 0) {
                if (CFL < dt.CFLmin())
                    os << '<' << flush;
                else if (CFL > dt.CFLmax())
                    os << '>' << flush;
                else
                    os << '.' << flush;
            }
            ddc.advance(fields, dt.n());
            if (dt.variable() && dt.adjust(CFL, false))
                ddc.reset_dt(dt);
        }

    }
    // Poincare section computation: return Poincare crossing nearest to t=T, with Tmin < t < Tmax.
    else {
        cout << "Poincare sectioning not yet implemented (markd as experimental)." << endl;
        exit(1);
        /*    // Adjust dt for CFL if necessary
            DNSPoincare dns (f_u, h, flags);
            if (dt.variable()) {
              dns.advance (f_u, p, 1);
              dt.adjust (dns.CFL());
              dns.reset_dt (dt,u);
              f_u = u;
            }
            // Collect all Poincare crossings between Tfudgemin and Tfudgemax
            // If we don't find one in that range, go as far as Tlastchance
            Real dTfudge = 1.05;
            Real Tfudgemin = lesser (0.90*T, T - dTfudge*dt.dT());
            Real Tfudgemax = Greater (1.02*T, T + dTfudge*dt.dT());
            Real Tlastchance = 10*T;

            vector<FlowField> ucross;
            vector<Real> tcross;
            vector<int>  scross;
            int s=0;
            int crosssign = 0; // look for crossings in either direction

            os << "f^t: " << flush;

            for (Real t=0; t<=Tlastchance; t += dt.dT(), ++s) {

              CFL = dns.CFL();

              if (s % 10 == 0)   os << iround (t);
              else if (s % 2 == 0) {
                if (CFL > dt.CFLmax())  os << '>';
                else if (CFL < dt.CFLmin())  os << '<';
                else  os << '.';
                os << flush;
              }

              // Collect any Poincare crossings
              bool crossed = dns.advanceToSection (f_u, p, dt.n(), crosssign, Tfudgemin);
              if (crossed && t >= Tfudgemin) {
                ucross.push_back (dns.ucrossing());
                tcross.push_back (dns.tcrossing());
                scross.push_back (dns.scrossing());
              }

              // If we've found at least one crossing within the fudge range, stop.
              // Otherwise continue trying until Tlastchance
              if (ucross.size() > 0 && t >= Tfudgemax)
                break;
            }

            if (ucross.size() <1) {

              os << "\nError in f(u, T, f_u, flags, dt, fcount, CFL, os) :\n";
              os << "the integration did not reach the Poincare section.\n";
              os << "Returning laminar solution and a b.s. value for the crossing time.\n";
              os << "I hope you can do something useful with them." << endl;
              f_u.setToZero();
              T = dns.time();
              ++fcount;
              return;
            }
            os << "  " << flush;

            // Now select the crossing that is closest to the estimated crossing time
            FlowField ubest = ucross[0];
            Real  Tbest = tcross[0];
            int   sbest = scross[0];
            int   nbest = 0;

            for (uint n=1; n<ucross.size(); ++n) {
              if (abs (tcross[n]-T) < abs (Tbest-T)) {
                ubest = ucross[n];
                Tbest = tcross[n];
                sbest = scross[n];
                nbest = n;
              }
            }
            os << nbest << (sbest==1 ? '+' : '-') << " at t== " << Tbest << flush;

            T = Tbest;
            f_u = ubest;

            // Now check if there are any crossings of opposite sign close by.
            // This signals near-tangency to Poincare section, which'll mess up
            // the search. Just print warning and let user intervene manually.
            Real l2distubestucross = 0;
            for (uint n=0; n<ucross.size(); ++n) {
              l2distubestucross = L2Dist (ubest, ucross[n]);
              if ( (u.taskid() == 0) && (scross[n] != sbest)) {
                os << "\nWARNING : There is a nearby Poincare crossing of opposite sign," << endl;
                os << "signalling near-tangency to section. You should probably switch " << endl;
                os << "to another Poincare crossing." << endl;
                os << "(ubest, unear) signs == " << sbest << ", " << scross[n] << endl;
                os << "(ubest, unear) times == " << Tbest << ", " << tcross[n] << endl;
                os << "(ubest, unear) dist  == " << l2distubestucross << endl;
              }
            }*/
    }

    if (!isfinite(L2Norm(f_u))) {
        os << "error in f: f(u,t) is not finite. exiting." << endl;
        exit(1);
    }
    if (!isfinite(L2Norm(f_temp))) {
        os << "error in f: f(temp,t) is not finite. exiting." << endl;
        exit(1);
    }
    if (!isfinite(L2Norm(f_salt))) {
        os << "error in f: f(salt,t) is not finite. exiting." << endl;
        exit(1);
    }

    ++fcount;
    f_u = fields[0];
    f_temp = fields[1];
    f_salt = fields[2];
    return;
}

}  // namespace chflow