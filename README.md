# cooper-verner

This is a simple implementation of the Cooper-Verner algorithm for numerically solving initial-value ODE systems.

The algorithm is a fixed-step explicit 8th-order member of the Runge-Kutta family of ODE solvers that uses 11 function evaluations per step.

Provided with the method is a simple 4th-order polynomial interpolant for dense output.

The method has been sourced from: 
Cooper, G. J., and Verner, J. H. “Some Explicit Runge-Kutta Methods of High Order.” SIAM Journal on Numerical Analysis, vol. 9, no. 3, September 1972, pp. 389–405. JSTOR, http://www.jstor.org/stable/2156139.
