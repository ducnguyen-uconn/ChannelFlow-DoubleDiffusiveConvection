/**
 * free-slip boundary condition on walls
 *
 * Original author: Duc Nguyen
 */

#ifndef FREESLIPBC_H
#define FREESLIPBC_H

#include "channelflow/diffops.h"
#include "channelflow/flowfield.h"
#include "modules/ddc/macros.h"
#ifdef FREESLIP
namespace chflow {
void freeslipBC(std::vector<FlowField>& fields) {
    FlowField u = fields[0];
    
    if (u.Nd() != 3) {
        cferror("freeslipBC(u): u.Nd() = " + i2s(u.Nd()) + " != 3");
        exit(1);
    }
    
    u.makePhysical();
    #ifdef P5
    FlowField temp = fields[1];
    temp.makePhysical();
    #endif
    #ifdef P6
    FlowField salt = fields[2];
    salt.makePhysical();
    #endif

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();
    
    // top
    if(nylocmin == 0){
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx){
            for (lint nz = 0; nz < Nz; ++nz) {
                u(nx, nylocmin, nz, 0) = u(nx, nylocmin+1, nz, 0); // du/dy=0
                u(nx, nylocmin, nz, 1) = 0; // No-penetration
                u(nx, nylocmin, nz, 2) = u(nx, nylocmin+1, nz, 2); // dw/dy=0
                #ifdef P5
                temp(nx, nylocmin, nz, 0) = 0; // constant temperature
                #endif
                #ifdef P6
                salt(nx, nylocmin, nz, 0) = 0; // constant salinity
                #endif
            }
        }
    }

    // bottom
    if(nylocmax == u.Ny()-1){
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx){
            for (lint nz = 0; nz < Nz; ++nz) {
                u(nx, nylocmax, nz, 0) = u(nx, nylocmax-1, nz, 0); // du/dy=0
                u(nx, nylocmax, nz, 1) = 0; // No-penetration
                u(nx, nylocmax, nz, 2) = u(nx, nylocmax-1, nz, 2); // dw/dy=0
                #ifdef P5
                temp(nx, nylocmax, nz, 0) = 0; // constant temperature
                #endif
                #ifdef P6
                salt(nx, nylocmax, nz, 0) = 0; // constant salinity
                #endif
            }
        }
    }

    u.makeSpectral();
    fields[0] = u; // update velocity
    #ifdef P5
    temp.makeSpectral();
    fields[1] = temp; // update velocity
    #endif
    #ifdef P6
    salt.makeSpectral();
    fields[2] = salt; // update velocity
    #endif
}
}  // namespace chflow
#endif
#endif