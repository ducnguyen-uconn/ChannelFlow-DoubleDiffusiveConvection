/**
 * turbulence statistics for turbulent flow
 *
 * Original author: Duc Nguyen
 */

#ifndef TUR_STATS_H
#define TUR_STATS_H

#include "channelflow/turbstats.h"
#include "channelflow/flowfield.h"
#include "modules/ddc/macros.h"
#ifdef P6
#ifdef SAVESTATS
namespace chflow {
    class DDCTurbStats : public TurbStats {
        public:
        DDCTurbStats(const int Ny): 
            temp_m(Ny),
            salt_m(Ny),
            temp_grad(Ny),
            salt_grad(Ny),
            temp_flux(Ny),
            salt_flux(Ny) { };
        void saveTurbStats(const string& filebase, vector<FlowField>& fieldsn){
            string filename(filebase);
            filename += string(".csv");
            ofstream os(filename.c_str());
            os << setprecision(REAL_DIGITS);
            char s = ',';
            int Ny = temp_m.length();
            Vector y = fieldsn[0].ygridpts();
            
            os  << "ypoints" << s
                << "temp_m" << s << "salt_m" << s 
                << "temp_grad" << s << "salt_grad" << s 
                << "temp_flux" << s << "salt_flux" << '\n' ;
            for (int ny = 0; ny < Ny; ++ny) {
                os  << y[ny] << s
                    << temp_m[ny] << s << salt_m[ny] << s
                    << temp_grad[ny] << s << salt_grad[ny] << s
                    << temp_flux[ny] << s << salt_flux[ny] << '\n';
            }
        };
        void addSnapshot(vector<FlowField>& fieldsn, const DDCFlags flags){
            FlowField u = totalVelocity(fieldsn[0], flags);
            FlowField t = totalTemperature(fieldsn[1], flags);
            FlowField dtdy = ydiff(t);
            FlowField s = totalSalinity(fieldsn[2], flags);
            FlowField dsdy = ydiff(s);
            
            u.makePhysical();
            t.makePhysical();
            dtdy.makePhysical();
            s.makePhysical();
            dsdy.makePhysical();
            lint Nx = u.Nx();
            lint Ny = u.Ny();
            lint Nz = u.Nz();
            lint nxlocmin = u.nxlocmin();
            lint nxlocmax = u.nxlocmin() + u.Nxloc();
            lint nylocmin = u.nylocmin();
            lint nylocmax = u.nylocmax();
            for (lint ny = 0; ny < Ny; ++ny){
                Vector sum(6);
                if (nylocmin <= ny && ny < nylocmax){
                    for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
                        for (lint nz = 0; nz < Nz; ++nz) {
                            sum[0] += t(nx,ny,nz,0);
                            sum[1] += s(nx,ny,nz,0);
                            sum[2] += dtdy(nx,ny,nz,0);
                            sum[3] += dsdy(nx,ny,nz,0);
                            sum[4] += u(nx,ny,nz,1)*t(nx,ny,nz,0);
                            sum[5] += u(nx,ny,nz,1)*s(nx,ny,nz,0);
                        }
                }
                // Sum up results from all processes    
                #ifdef HAVE_MPI
                Vector tmp(sum) ;
                MPI_Allreduce(&tmp[0], &sum[0], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                MPI_Allreduce(&tmp[1], &sum[1], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                MPI_Allreduce(&tmp[2], &sum[2], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                MPI_Allreduce(&tmp[3], &sum[3], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                MPI_Allreduce(&tmp[4], &sum[4], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                MPI_Allreduce(&tmp[5], &sum[5], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                #endif
                temp_m[ny] = sum[0]/(Nx*Nz);
                salt_m[ny] = sum[1]/(Nx*Nz);
                temp_grad[ny] = sum[2]/(Nx*Nz);
                salt_grad[ny] = sum[3]/(Nx*Nz);
                temp_flux[ny] = sum[4]/(Nx*Nz);
                salt_flux[ny] = sum[5]/(Nx*Nz);
            }
        };

        private:
        Vector temp_m;// horizontally averaged profile of temperature [<T>h]
        Vector salt_m;// horizontally averaged profile of salinity [<S>h]
        // horizontally averaged profile of density [<S-Rrho*T>h/(1-Rrho)]

        Vector temp_grad;// horizontally averaged profile of temperature gradient [<dT/dy>h]
        Vector salt_grad;// horizontally averaged profile of salinity gradient [<dS/dy>h]

        Vector temp_flux;// horizontally averaged profile of temperature's convective flux [<v*T>h]
        Vector salt_flux;// horizontally averaged profile of salinity's convective flux [<v*S>h]
    
        // horizontally averaged profiles of turbulent diffusivity [<v*T>h/<dT/dy>h] and [<v*S>h/<dS/dy>h]
    };
}  // namespace chflow
#endif
#endif
#endif