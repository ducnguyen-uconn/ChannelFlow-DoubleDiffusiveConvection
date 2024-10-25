#ifndef MACROS_H
#define MACROS_H
/**
 * @brief this file define governing equations of double-component systems like double-diffusive convection.
 * 
 * A system of double components has form:
 * Momentum equation (variable U, P):           dt(U) + div(grad(Utot)) = - grad(P) + p1*lap(U) + p1*p2(p3*T-p4*S)*ey
 * First component's equation (variable T):     dt(T) + div(grad(Ttot)) =             p5*lap(T)
 * Second component's equation (variable S):    dt(S) + div(grad(Stot)) =             p6*lap(S) + p7*lap(T)
 * where Utot, = Ubase + U, is total velocity.
 * 
 * Example 1: for double-diffusive convection [Yang2016PNAS]: p1 = Pr, p2 = Ra, p3 = 1, p4 = 1/Rrho, p5 = 1, p6 = 1/Le, p7 = 0
 * Example 2: for binary fluid convection [Mercader2013JFM]: p1 = Pr, p2 = Ra, p3 = 1+Rsep, p4 = 1, p5 = 1, p6 = tau = Le, p7 = 1
 * 
 */

// Example 1: for double-diffusive convection [Yang2016PNAS]: p1 = Pr, p2 = Ra, p3 = 1, p4 = 1/Rrho, p5 = 1, p6 = 1/Le, p7 = 0
#define P1 Pr 
#define P2 Ra
#define P3 1.0
#define P4 1.0/Rrho
#define P5 1.0
#define P6 1.0/Le
#define P7 0.0

// Example 2: for binary fluid convection [Mercader2013JFM]: p1 = Pr, p2 = Ra, p3 = 1+Rsep, p4 = 1, p5 = 1, p6 = Le, p7 = 1
// #define P1 Pr 
// #define P2 Ra
// #define P3 (1.0+Rsep)
// #define P4 Rsep
// #define P5 1.0
// #define P6 Le
// #define P7 1.0

#endif