# MHD-Riemann-solvers 
Simulate the behaviour of compressible fluid in the presence or absence of magnetic field. Contain features like
1) Positive preserving schemes for density and pressure
2) Constrained transport method to ensure solenoidal magnetic field
3) Compact third-order slope limiter (TVD)
4) OpenMP implementation for parallel processing
5) Iteration by column-major order for better cache utilization and improved performance

Types of problem to solve:
1) Shock Tube
2) Rotor problem
3) Blast wave

Example of MHD rotor simulation: 
https://www.youtube.com/watch?v=Ac8E4WVf0NM&ab_channel=Spaceduck496
