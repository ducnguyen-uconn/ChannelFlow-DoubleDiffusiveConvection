A exact executable command likes this:
```bash
mpiexec -n <ncpu> ./<exename> <option1> <option2> ...
```

|Option  | Default   | Description |
|:------------------------|:----------|:------------------------------------------------------------------|
|`-Nx <value>`| $200$| Number of points along x-direction |
|`-Ny <value>`| $101$ | Number of points along y-direction, notes that $N_y\%2=1$ |
|`-Nz <value>`| $10$ | Number of points along z-direction |
|`-Lx <value>`| $2$| Streamwise length |
|`-Lz <value>`| $0.05$ | Spanwise length |
|`-Pr <value>` | $10$ | Prandtl number $Pr=\frac{\nu}{\kappa_T}$|
|`-Ra <value>`| $10^3$ | Thermal Rayleigh number $Ra_T=\frac{g\alpha\Delta T H^3}{\nu\kappa_T}$|
|`-Le <value>`| $100$ | Lewis number $Le=\frac{\kappa_T}{\kappa_S}$ |
|`-Rr <value>`| $2$ | Density stability ratio $R_\rho=\frac{\alpha \Delta T}{\beta \Delta S}$ |
|`-Ua <value>`| $0$ | X-velocity at lower wall, U(y=a) |
|`-Ub <value>`| $0$ | X-velocity at upper wall, U(y=b) |
|`-Wa <value>`| $0$ | Z-velocity at lower wall, W(y=a) |
|`-Wb <value>`| $0$ | Z-velocity at upper wall, W(y=b) |
|`-Ta <value>`| $0$ | Temperature at lower wall, T(y=a) |
|`-Tb <value>`| $1$ | Temperature at upper wall, T(y=b) |
|`-Sa <value>`| $0$ | Salinity at lower wall, S(y=a) |
|`-Sb <value>`| $1$ | Salinity at upper wall, S(y=b) |
|`-T0 <value>`| $0$ | Start time of DNS |
|`-T <value>`| $20$ | Final time of DNS |
|`-dt <value>`| $0.03125$ | Timestep |
|`-dT <value>`| $1$ | Save interval |
|`-nl <value>`| "rot" | Method of calculating  nonlinearity, one of [rot\|conv\|div\|skew\|alt\|linear] |


Examples:
```bash
mpiexec -n 16 ./ddc_simulateflow -Pr 10 -Ra 1000 -Le 100 -Rr 2 -dt 0.02 -dT 1 -T 100 \
-Nx 200 -Ny 81 -Nz 10 -Lx 2 -Lz 0.02 -bf "laminar" -nl "conv"
```