/**
 * This class creates initial perturbations
 *
 * Original author: Duc Nguyen
 */

#ifndef ADD_PERTURB_H
#define ADD_PERTURB_H

#include "channelflow/diffops.h"
#include "channelflow/flowfield.h"

namespace chflow {
void addSinusoidalPerturbations(FlowField& f, Real mag, Real k) {
    
    Vector y = f.ygridpts();
    
    if (mag == 0.0)
        return;

    if (f.Nd() == 3) {
        f.makePhysical();
        srand(time(0));
        // Add a div-free perturbation to the base flow.
        lint Nz = f.Nz();
        lint nxlocmin = f.nxlocmin();
        lint nxlocmax = f.nxlocmin() + f.Nxloc();
        lint nylocmin = f.nylocmin();
        lint nylocmax = f.nylocmax();
        for (lint ny = nylocmin; ny < nylocmax; ++ny)
            for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
                for (lint nz = 0; nz < Nz; ++nz) {
                    f(nx, ny, nz, 0) = mag*sin(2*pi*k*y[ny]);
                }
        f.makeSpectral();
    } else if (f.Nd() == 1) {
        f.makePhysical();
        lint Nz = f.Nz();
        lint nxlocmin = f.nxlocmin();
        lint nxlocmax = f.nxlocmin() + f.Nxloc();
        lint nylocmin = f.nylocmin();
        lint nylocmax = f.nylocmax();
        for (lint ny = nylocmin; ny < nylocmax; ++ny)
            for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
                for (lint nz = 0; nz < Nz; ++nz) {
                    f(nx, ny, nz, 0) = mag*sin(2*pi*k*y[ny]);
                }
        f.makeSpectral();
    } else {
        exit(1);
    }
    
}
void addRandomPerturbations(FlowField& f, Real mag) {
    
    if (mag == 0.0)
        return;

    if (f.Nd() == 3) {
        f.makePhysical();
        srand(time(0));
        // Add a div-free perturbation to the base flow.
        lint Nz = f.Nz();
        lint nxlocmin = f.nxlocmin();
        lint nxlocmax = f.nxlocmin() + f.Nxloc();
        lint nylocmin = f.nylocmin();
        lint nylocmax = f.nylocmax();
        for (lint ny = nylocmin; ny < nylocmax; ++ny)
            for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
                for (lint nz = 0; nz < Nz; ++nz) {
                f(nx, ny, nz, 0) = mag*randomReal(0.0,1.0);
                f(nx, ny, nz, 1) = mag*randomReal(0.0,1.0);
                f(nx, ny, nz, 2) = mag*randomReal(0.0,1.0);
            }
        f.makeSpectral();
    } else if(f.Nd() == 1){
        f.makePhysical();
        srand(time(0));
        // Add a div-free perturbation to the base flow.
        lint Nz = f.Nz();
        lint nxlocmin = f.nxlocmin();
        lint nxlocmax = f.nxlocmin() + f.Nxloc();
        lint nylocmin = f.nylocmin();
        lint nylocmax = f.nylocmax();
        for (lint ny = nylocmin; ny < nylocmax; ++ny)
            for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
                for (lint nz = 0; nz < Nz; ++nz) {
                f(nx, ny, nz, 0) = mag*randomReal(0.0,1.0);
            }
        f.makeSpectral();
    } else {
        exit(1);
    }
    
}
}  // namespace chflow
#endif