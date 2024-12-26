#ifndef MACROS_H
#define MACROS_H
/**
 * @brief this file defines governing equations of double-component systems like double-diffusive convection.
 * 
 * A system of double components has form:
 * Momentum equation (variable U, P):           dt(U) + div(grad(U)) = - grad(P) + p1*lap(U) + p2(p3*T-p4*S)*ey
 * First component's equation (variable T):     dt(T) + div(grad(T)) =             p5*lap(T)
 * Second component's equation (variable S):    dt(S) + div(grad(S)) =             p6*lap(S) + p7*lap(T)
 * where U, = Ubase + u, is total velocity.
 * 
 */

// Example 1: stationary-wall bounded double-diffusive convection [Yang2016PNAS]: p1 = Pr, p2 = Ra, p3 = 1, p4 = 1/Rrho, p5 = 1, p6 = 1/Le, p7 = 0
// #define P1 Pr 
// #define P2 Pr*Ra
// #define P3 1.0
// #define P4 1.0/Rrho
// #define P5 1.0
// #define P6 1.0/Le


// Example 2: Moving-wall bounded double-diffusive convection [Yang2021JFM]
// wall's boundary velocity is normalized into unit velocity, so we can set up U0=1.0
// for example, Ua=-0.5 Ub=0.5
// #define P1 sqrt(Pr/Ra) 
// #define P2 Ri/(Rrho-1)
// #define P3 1.0
// #define P4 Rrho
// #define P5 1.0/sqrt(Pr*Ra)
// #define P6 1.0/(Le*sqrt(Pr*Ra))


// Example 3: Binary fluid convection [Mercader2013JFM]: p1 = Pr, p2 = Ra, p3 = 1+Rsep, p4 = 1, p5 = 1, p6 = 1.0/Le, p7 = 1
// #define P1 Pr 
// #define P2 Pr*Ra
// #define P3 (1.0+Rsep)
// #define P4 Rsep
// #define P5 1.0
// #define P6 1.0/Le
// #define P7 1.0


// Example 4: Couette flow 
// #define P1 1.0/Rey


/* Example 5: Stratified plane Couette flow [Langham2019JFM] */
// #define P1 1.0/Rey
// #define P2 (-1.0)
// #define P3 Ri
// #define P5 (1.0/(Rey*Pr))

// Example 6: RBC [Zheng JFM 2024]
#define P1 sqrt(Pr/Ra) 
#define P2 1.0
#define P3 1.0
#define P5 1.0/sqrt(Pr*Ra)

// #define FREESLIP
// #define SAVESTATS
// #define EAVES2016JFM
#endif