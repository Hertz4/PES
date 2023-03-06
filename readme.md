<h1>Robust Analytic Continuation of Green's functions</h1>
<h2>Projection, Estimation and Semidefinite relaxation (PES)</h2>

Reference: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.107.075151

Requirement of code:

1. cvx: http://cvxr.com/cvx/download/
2. Chebfun: https://www.chebfun.org/download/

Our code could deal with both clean and noisy Matsubara data, both scalar and matrix data. See the references for detailed choices.

There are three folders. 

Folder 1: clean (noiseless) case: ES method for noiseless fermionic Green's function. 

       matlab>>> test_noiseless
        
Folder 2: noisy case: PES method for noisy fermionic Green's function. 

       matlab>>> test_noisy1
       matlab>>> test_noisy2
       matlab>>> test_noisy3
Folder 3: bosonic: ES/PES method for Bosonic case.

       matlab>>> test_bosonic
