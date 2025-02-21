"""

To run, restart, and plot using e.g. 4 processes:
    $ mpiexec -n 8 python3 ./ddc/validations/yang2021jfm_case3_2d.py
    $ mpiexec -n 8 python3 ./ddc/validations/yang2021jfm_case3_2d.py --restart
    $ mpiexec -n 8 python3 _.py _/*.h5
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import math
import h5py
import logging
logger = logging.getLogger(__name__)


# Allow restarting via command line
restart = (len(sys.argv) > 1 and sys.argv[1] == '--restart')

dealias = 3/2  
pi = np.pi
Pr, Ra, Rp, Le = 10.0, 1e6, 2.0, 100.0 # Yang2021JFM case 3 
Lx, Lz, Nx, Nz = 2.0, 1.0, 384, 384
stop_sim_time = 700 + 300*restart # Stopping criteria
# Bases
coords = d3.CartesianCoordinates('x','z')
dist = d3.Distributor(coords, dtype=np.float64)
# define the coordinate system
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
zbasis = d3.ChebyshevT(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)
# define fields
p = dist.Field(name='p', bases=(xbasis,zbasis)) # pressure
u = dist.VectorField(coords, name='u', bases=(xbasis,zbasis)) # velocity
sa = dist.Field(name='sa', bases=(xbasis,zbasis)) # salinity
te = dist.Field(name='te', bases=(xbasis,zbasis)) # temperature
baru = dist.Field(bases=(zbasis))

# Substitutions
x, z = dist.local_grids(xbasis, zbasis) # get coordinate arrays in horizontal and vertical directions
ex, ez = coords.unit_vector_fields(dist) # get unit vectors in horizontal and vertical directions
# define vertical velocity component
w = u @ ez

# create constant sub-field for incompressible flow condition's equation
tau_p = dist.Field(name='tau_p') 
tau_t1 = dist.Field(name='tau_t1', bases=xbasis)
tau_t2 = dist.Field(name='tau_t2', bases=xbasis)
tau_s1 = dist.Field(name='tau_s1', bases=xbasis)
tau_s2 = dist.Field(name='tau_s2', bases=xbasis)
tau_u1 = dist.VectorField(coords, name='tau_u1', bases=xbasis)
tau_u2 = dist.VectorField(coords, name='tau_u2', bases=xbasis)
# because this term is only a contant added to the equation, we don't need to instantiate it for bases system


lift_basis = zbasis.derivative_basis(1)
lift = lambda A: d3.Lift(A, lift_basis, -1)
grad_te = d3.grad(te)  + ez*lift(tau_t1) # First-order reduction
grad_sa = d3.grad(sa)  + ez*lift(tau_s1) # First-order reduction
grad_u = d3.grad(u)  + ez*lift(tau_u1) # First-order reduction
# First-order form: "lap(f)" becomes "div(grad_f)"
lap_u = d3.div(grad_u)
lap_te = d3.div(grad_te)
lap_sa = d3.div(grad_sa)
# First-order form: "div(A)" becomes "trace(grad_A)"

dx = lambda A: d3.Differentiate(A, coords['x']) 
dz = lambda A: d3.Differentiate(A, coords['z']) 

baru['g'] = 0.5*(z-0.5)

# Problem
problem = d3.IVP([p, tau_p, u, te, sa,tau_t1, tau_t2,tau_s1, tau_s2, tau_u1, tau_u2], namespace=locals())
problem.add_equation("trace(grad_u) + tau_p = 0")
problem.add_equation("integ(p) = 0") # Pressure gauge
problem.add_equation("dt(u) + baru*dx(u) + w*dz(baru)*ex + grad(p) - Pr*lap_u - Pr*Ra*(te-sa/Rp)*ez + lift(tau_u2) = - u@grad_u")
problem.add_equation("dt(te) + baru*dx(te) - lap_te - w  + lift(tau_t2)= - u@grad_te")
problem.add_equation("dt(sa) + baru*dx(sa) - (1./Le)*lap_sa - Rp*w  + lift(tau_s2)= - u@grad_sa")

problem.add_equation("u(z=0) = 0")
problem.add_equation("u(z=Lz) =  0")
problem.add_equation("te(z=0) = 0")
problem.add_equation("te(z=Lz) =  0")
problem.add_equation("sa(z=0) = 0")
problem.add_equation("sa(z=Lz) =  0")

# timestepper = d3.RK443 # 3rd-order 4-stage DIRK+ERK scheme [Ascher 1997 sec 2.8] https://doi-org.ezproxy.lib.uconn.edu/10.1016/S0168-9274(97)00056-1
timestepper = d3.RK222

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Initial conditions
if not restart:
    p.fill_random('g', seed=42, distribution='normal', scale=0.7*1e-4) # Random noise
    u.fill_random('g', seed=42, distribution='normal', scale=0.7*1e-4) # Random noise
    te['g'] = -0.05*np.sin(2.0*pi*6.0*z)
    sa['g'] = -0.05*np.sin(2.0*pi*6.0*z)
    file_handler_mode = 'overwrite'
    initial_timestep = 0.02
    max_timestep = 0.02
else:
    write, initial_timestep = solver.load_state('checkpoints/checkpoints_s1.h5')
    max_timestep = 0.02
    file_handler_mode = 'append'
    logger.info('Imported last-step data successfully')

# store data for analysis later
dataset = solver.evaluator.add_file_handler('snapshots', sim_dt=10.0, max_writes=1000, mode=file_handler_mode)
dataset.add_task(te, name='temperature')
dataset.add_task(sa, name='salinity')
# dataset.add_task(p, name='pressure')
# dataset.add_task(u@ex, name='velocity_u')
# dataset.add_task(u@ez, name='velocity_w')
# dataset.add_task(-d3.div(d3.skew(u)), name='vorticity')
# store data to restart later
checkpoints = solver.evaluator.add_file_handler('checkpoints', sim_dt=100, max_writes=1, mode=file_handler_mode)
checkpoints.add_tasks(solver.state)

# CFL
CFL = d3.CFL(solver, initial_timestep, cadence=10,max_dt=max_timestep,
             safety=0.2, threshold=0.1,
             max_change=1.5, min_change=0.5
             )
CFL.add_velocity(u)

xg = xbasis.global_grid(dist, scale=dealias)
zg = zbasis.global_grid(dist, scale=dealias)
Tg = te.allgather_data('g')

# Main loop
oldtime = 0.
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        solver.step(timestep)  
        if (solver.iteration-1) % 100 == 0:
            logger.info('Completed iteration {}, time={:.3f}, dt={:.10f}'.format(solver.iteration, solver.sim_time, timestep))        
            ########################### <--- plot instantaneous temperature distribution
            Tg = te.allgather_data('g')
            Sg = sa.allgather_data('g')
            if dist.comm.rank == 0:
                # plot temperature distribution
                plt.figure(figsize=(8,3))
                plt.pcolormesh(xg.ravel(),zg.ravel(),Sg.transpose(),cmap='jet')
                plt.colorbar() 
                plt.xlabel(r'$x$')
                plt.ylabel(r'$z$')
                plt.title("t = {:.3f}".format(solver.sim_time))
                plt.savefig('snapshots/S_time={:.3f}.png'.format(solver.sim_time), bbox_inches='tight')
                plt.close()
                # plot horizontaly averaged density profiles
                # meanT = np.mean(Tg,axis=0,keepdims=True)
                # meanS = np.mean(Sg,axis=0,keepdims=True)
                # meanDensity = (meanS-meanT)+(1-Rp)*np.linspace(0,1.,int(Nz*dealias))
                # plt.figure(figsize=(3,3))
                # plt.plot(np.copy(meanDensity[0]),np.linspace(0,1.,int(Nz*dealias)))
                # plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)# Hide xticks
                # plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
                # plt.title("t = {:.3f}".format(solver.sim_time))
                # plt.savefig('snapshots/meanDensity_time={:.3f}.png'.format(solver.sim_time), bbox_inches='tight')
                # plt.close()
            ###########################
            if math.isnan(np.max(Tg)):
                logger.error('NaN values')
                break
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()