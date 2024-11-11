<h2 style="text-align:center;">Double-Diffusive Convection: An Extensional Module of Channelflow</h2>

<div style="text-align: center">
    <figure>
        <img src="images/finger2d.mp4" alt="Salinity behavior in 2D fingering mode" width="50%"/>
        <figcaption>Salinity behavior of 2D finger regime</figcaption>
    </figure>
    ![Video](<video width="640" height="480" controls>
    <source src="images/finger2d.mp4" type="video/mp4">
    Salinity behavior of 2D finger regime.
    </video>)
</div>


<div style="text-align: center">
    <figure>
        <img src="images/s_yang2021jfm_case3_2d_noslip.gif" alt="Salinity" width="50%" align="left">
        <img src="images/t_yang2021jfm_case3_2d_noslip.gif" alt="Image 2" width="50%" align="right">
        <figcaption>Salinity and temperature behaviors of 2D diffusive regime in channel flow configuration.</figcaption>
    </figure>
</div>


ChannelFlow-DoubleDiffusiveConvection is an extensional module of Channelflow 2.0 for wall-bounded double-component problems like Double-Diffusive Convection and Binary Fluid Convection. To use this code, pls read following instruction to [install Channelflow](INSTALL.md) (this also contains setup on a [HPC](HPCsetup.md))

After knowing how to install the standard Channelflow, you can clone them to your local machine and add present `ddc` module inside by
```bash
git clone https://github.com/epfl-ecps/channelflow.git

cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
mkdir -p ./channelflow/modules/
rm -rf ./channelflow/modules/ddc
cp -r ./ddc ./channelflow/modules/ddc
``` 
And build them, for example
```bash
mkdir -p build
cd build
cmake ../channelflow -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DWITH_DDC=ON -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/user/local/ -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm -Wno-unused-variable " -DWITH_SHARED=OFF -DWITH_HDF5CXX=OFF
make -j16
```

### Building DNS
First of all, you need to define governing equations of problem. In this code, we offer nondimensional governing equations, which have form:

$`\begin{align}
    \frac{\partial \boldsymbol{u}}{\partial t} + \boldsymbol{u}_{\text{tot}} \cdot \nabla \boldsymbol{u}_{\text{tot}} &= -\nabla p + p_1\nabla^2\boldsymbol{u}+ p_1p_2(p_3\theta - p_4s) \mathbf{j},\\
    \frac{\partial \theta}{\partial t} + \boldsymbol{u}_{\text{tot}} \cdot \nabla \theta_{\text{tot}} &= p_5\nabla^2 \theta,\\
    \frac{\partial s}{\partial t} + \boldsymbol{u}_{\text{tot}} \cdot \nabla s_{\text{tot}} &= p_6 \nabla^2 s + p_7\nabla^2 \theta,\\
    \nabla \cdot \boldsymbol{u} &= 0,
\end{align}`$

where $\boldsymbol{u}$, $\theta$, and $s$ are fluctuations of the velocity, temperature, and a third field. If you can not see equations, let read [pdf](README.pdf) file instead. The third field $s$ may be the salinity in double-diffusive convection (Radko [2013](https://doi.org/10.1017/CBO9781139034173)) or the convective mass flux in binary fluid convection (Mercader [2013](https://doi.org/10.1017/jfm.2013.77)). The subscript `tot` indicates the total value of fields, which is defined as sum of base flow and fluctuation of each field. Also, another suggestion of governing equations is introduced for DDC in channel flow (Yang [2021](https://doi.org/10.1017/jfm.2021.1091)) with wall's boundary velocity normalized into unit velocity, $U_0=1$. Because this code offers two options DDC and BFC, so first you need to define the problem you want to use in this code. This is perfomed by modifying controling parameters ($p_i$) via a header file `ddc/macros.h`.


| Parameters | Double-diffusive convection  | Moving-wall bounded Double-diffusive convection   | Binary fluid convection  | Description                                                       |
|:------------------------|:--------|:----------|:----------|:------------------------------------------------------------------|
| $p_1$ | $Pr_T$  | $\sqrt{\frac{Pr_T}{Ra_T}}$ | $Pr$ | $Pr_T=\frac{\nu}{\kappa_T}$ is Prandtl number |
| $p_2$ | $Ra_T$  | $\sqrt{\frac{Ra_T}{Pr_T}}$ | $Ra_T$ | $Ra_T=\frac{g\alpha \Delta_T H^3}{\nu\kappa_T}$ is Thermal Rayleigh number |
| $p_3$ | $1$  | $1$ | $1+R_{sep}$ | $R_{sep}$ is Separation ratio  |
| $p_4$ | $\frac{1}{R_\rho}$ | $\frac{1}{R_\rho}$ | $R_{sep}$ | $R_\rho=\frac{\beta\Delta_S}{\alpha\Delta_T}$ is Density stability ratio  |
| $p_5$ | $1$  | $\frac{1}{\sqrt{Pr_T Ra_T}}$ | $1$ | |
| $p_6$ | $\frac{1}{Le}$  | $\frac{1}{Le\sqrt{Pr_T Ra_T}}$ | $\frac{1}{Le}$ | $Le=\frac{\kappa_T}{\kappa_S}$ is Lewis number |
| $p_7$ | $0$  | $0$    | $1$  |  |

A correct defination looks like this
```cpp
// Example 1: Stationary-wall bounded double-diffusive convection [Yang2016PNAS]
// Boundary velocity is not normalized into unit velocity, don't know exact U0
// Only use this form for stationary walls
#define P1 Pr 
#define P2 Ra
#define P3 1.0
#define P4 1.0/Rrho
#define P5 1.0
#define P6 1.0/Le
#define P7 0.0

// Example 2: Moving-wall bounded double-diffusive convection [Yang2021JFM]
// Boundary velocity is normalized into unit velocity, U0=1.0
// For example, Ua=-0.5 Ub=0.5 or Ua=0 Ub=1
// #define P1 sqrt(Pr/Ra) 
// #define P2 (1.0/sqrt(Pr/Ra))
// #define P3 1.0
// #define P4 (1.0/Rrho)
// #define P5 (1.0/sqrt(Pr*Ra))
// #define P6 (1.0/(Le*sqrt(Pr*Ra)))
// #define P7 0.0

// Example 3: Binary fluid convection [Mercader2013JFM]
// #define P1 Pr 
// #define P2 Ra
// #define P3 (1.0+Rsep)
// #define P4 Rsep
// #define P5 1.0
// #define P6 1.0/Le
// #define P7 1.0
```


### Running executable files
A exact executable command likes this:
```bash
mpiexec -n <ncpu> ./<exename> <option1> <option2> ...
```
with controling parameters

|Option  | Default   | Description |
|:------------------------|:----------|:------------------------------------------------------------------|
|`-Nx <value>`| $200$| Number of points along x-direction |
|`-Ny <value>`| $101$ | Number of points along y-direction, notes that $N_y$ is odd number |
|`-Nz <value>`| $6$ | Number of points along z-direction, minimal number is 6 ( can be used for 2D setup) |
|`-Lx <value>`| $2$| Streamwise length |
|`-Lz <value>`| $0.004$ | Spanwise length, default setup is a small length representing 2D domain |
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
mpiexec -n 16 ./ddc_simulateflow -Pr 10 -Ra 1000 -Le 100 -Rr 2 -dt 0.02 -dT 1 -T 100 -Nx 200 -Ny 81 -Nz 10 -Lx 2 -Lz 0.02 -nl "conv"
```