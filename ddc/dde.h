/**
 * System class of double-diffusive equations for standard channel flows
 *
 * Original author: Duc Nguyen
 */

#ifndef DDE_H
#define DDE_H

#include "channelflow/diffops.h"
#include "channelflow/flowfield.h"
#include "channelflow/nse.h"
#include "channelflow/tausolver.h"
#include "modules/ddc/ddcflags.h"

namespace chflow {

// nonlinear term of NSE plus the linear coupling term to the temperature equation
void momentumNL(const FlowField& u, const FlowField& T, const FlowField& S, 
                ChebyCoeff Ubase, ChebyCoeff Wbase,
                FlowField& f,FlowField& tmp, DDCFlags flags);

// nonlinear term of heat equation plus the linear coupling term to the momentum equation
void temperatureNL(const FlowField& u, const FlowField& T, 
                   ChebyCoeff Ubase, ChebyCoeff Wbase, ChebyCoeff Tbase,
                   FlowField& f, FlowField& tmp, DDCFlags flags);

// nonlinear term of salt equation plus the linear coupling term to the momentum equation
void salinityNL(const FlowField& u, const FlowField& S, 
                ChebyCoeff Ubase, ChebyCoeff Wbase, ChebyCoeff Sbase,
                FlowField& f, FlowField& tmp, DDCFlags flags);

class DDE : public NSE {
   public:
    
    DDE(const std::vector<FlowField>& fields, const DDCFlags& flags);
    DDE(const std::vector<FlowField>& fields, const std::vector<ChebyCoeff>& base, const DDCFlags& flags);
    virtual ~DDE();

    void nonlinear(const std::vector<FlowField>& infields, std::vector<FlowField>& outfields) override;
    void linear(const std::vector<FlowField>& infields, std::vector<FlowField>& outfields) override;

    // calls a tausolver for each Fourier mode
    void solve(std::vector<FlowField>& outfields, const std::vector<FlowField>& infields, const int i = 0) override;

    // redefines the tausolver objects with new time-stepping constant (allocates memory for tausolver at first use)
    void reset_lambda(const std::vector<Real> lambda_t) override;

    // vector of RHS is smaller than of fields because of missing pressure equation
    std::vector<FlowField> createRHS(const std::vector<FlowField>& fields) const override;

    // returns vector of symmetries confining the vector of fields to a subspace
    std::vector<cfarray<FieldSymmetry>> createSymmVec() const override;

    
    const ChebyCoeff& Ubase() const override; // (y-dependent)
    const ChebyCoeff& Wbase() const override; // constant
    const ChebyCoeff& Tbase() const; // constant
    const ChebyCoeff& Sbase() const; // constant

   protected:
    HelmholtzSolver*** heatsolver_;  // 3d cfarray of tausolvers, indexed by [i][mx][mz] for substep, Fourier Mode x,z
    HelmholtzSolver*** saltsolver_;

    DDCFlags flags_;  // User-defined integration parameters
    

    // additional base solution profiles
    ChebyCoeff Tbase_;    // temperature base profile (physical)
    ChebyCoeff Tbaseyy_;  // 2. deriv. of temperature base profile
    ChebyCoeff Sbase_;    // salinity base profile (physical)
    ChebyCoeff Sbaseyy_;  // 2. deriv. of salinity base profile

    ChebyCoeff Pbasey_;  //  wall normal pressure gradient (y-dependent)

    // constant terms
    ComplexChebyCoeff Cu_;  // constant
    ComplexChebyCoeff Cw_;
    ComplexChebyCoeff Ct_;
    ComplexChebyCoeff Cs_;
    bool nonzCu_;
    bool nonzCw_;
    bool nonzCt_;
    bool nonzCs_;

    ComplexChebyCoeff Tk_;
    ComplexChebyCoeff Rtk_;

    ComplexChebyCoeff Sk_;
    ComplexChebyCoeff Rsk_;

   private:
    void createDDCBaseFlow();
    void initDDCConstraint(const FlowField& u);  // method called only at construction
    void createConstants();

    bool baseflow_;
    bool constraint_;
};

// Construct laminar flow profile for given flow parameters.
// [a,b]   == y position of [lower, upper] walls
// ChebyCoeff ShearVelocityProfile(Real a, Real b, int Ny, DDCFlags flags);
ChebyCoeff laminarVelocityProfile(Real dPdx, Real Ubulk, Real Ua, Real Ub, Real a, Real b, int Ny,
                                  DDCFlags flags);
ChebyCoeff linearTemperatureProfile(Real a, Real b, int Ny, DDCFlags flags);
ChebyCoeff linearSalinityProfile(Real a, Real b, int Ny, DDCFlags flags);
ChebyCoeff PressureGradientY(Real a, Real b, int Ny, DDCFlags flags);

}  // namespace chflow
#endif