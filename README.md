# RSPA-2024-0558
Data for RSPA manuscript "Analytic Representation of the Polarization Tensor: Application to the Plasmon Resonance of Nanoparticles"

This is a repository of the codes and figures produced for this manuscript.

## Example 1: corrugated disk
### code: 
spectral.f90
### usage:
./spectral  <rad> <omega> <Ndpt>  <Nmom> <idir> 

      - rad/omega: geometric parameters for r= 1+rad*cos(omega*theta)
      - Ndpt : Number of boundary points per unit length
      - Nmom: Number of moments (even integer)
      - idir: set to i=idir for alpha[i,i]
      
### example:
./spectral .05 16 200 10 1
#### goal: 
- find the coefficients of the expansion $s [1- \nu_1 s + \nu_2 s^2 + \cdots]$ of the polarization $\alpha_{1,1}$ and $\alpha_{2,2}$
#### output file: 
- numoments.dat
#### post-processing: 
- run "maple spect.mw" find the PAs [m-1/m] and [m/m] from the moments $\{\nu[k]\} k=0,1,...,M$ of the variable $s = (eps1/eps2 -1)$, then plot $\text{Im}([m-1/m])$ or $[m/m]$ vs $z$, with $s= -1/(z+1/2)$, $z =x+iy_0$, $-0.5<x< 0.5$, $y_0=0.01$;


## Example 2: lens geometry

### code: 
lens.f90
### usage: 

./lens  rad psi Npan Nsub Nleg Nmom
  
  - rad: radius lens
  - psi: half-angle of lens (in degrees)
  - Npan: number of panels
  - Nsub: number of subpanels
  - Nleg: order of Gauss-Legendre quadrature on each subpanel (8, or 16)
  - Nmom: Number of moments (even integer)
  
#### example:
./lens 1.0 120 16 32 16 10
#### goal: 
find the coefficients of the expansion $2\lambda [1+ \mu_1 \lambda + \mu_2 \lambda^2 + \cdots]$ of the polarization $\alpha_{1,1}$ and $\alpha_{2,2}$
#### output file: 
moments.dat
#### post-processing: 
- run "maple polarization.mw" to find a Padé approximant of $\text{Im}[\alpha_{1,1}(z)]$ and $\text{Im}[\alpha_{2,2}(z)]$ with $z= -1/(2\lambda) = x+ iy_0$,  $y_0 \ll 1$
- run "maple diffpade.mw" to find a 1st-order differential Padé approximants 
