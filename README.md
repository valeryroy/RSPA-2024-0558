# RSPA-2024-0558
Data for RSPA manuscript "Analytic Representation of the Polarization Tensor: Application to the Plasmon Resonance of Nanoparticles"

This is a repository of the codes and figures produced for this manuscript.

## Example 1: lense geometry

### code: 

lens.f90

#### usage: 

./lens  rad psi Npan Nsub Nleg Nmom
  
  rad: radius lens
  
  psi: half-angle of lens (in degrees)
  
  Npan: number of panels
  
  Nsub: number of subpanels
  
  Nleg: order of Gauss-Legendre quadrature on each subpanel (8, or 16)
  
  Nmom: Number of moments (even integer)
  
#### example:
./lens 1.0 120 16 32 16 10

#### goal: 
find the coefficients of the expansion $2\lambda [1+ \mu_1 \lambda + \mu_2 lambda^2 + \cdots]$ of the polarization $\alpha_{1,1}$ and $\alpha_{2,2}$

#### output file: 
moments.dat
#### post-processing: 
run "maple polarization.mw" to find Im[alpha[1,1](z)] and Im[alpha[2,2](z)} with z= -1/(2*lambda) = x+ i*y0,  y0<<1
