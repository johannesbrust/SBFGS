# SBFGS
Structured-BFGS Optimization Algorithms

Algorithms for problems where the nonlinear objective function
has the form
                                                                       
                            minimize {f(x) = u(x) + k(x)}                                   

with x in R^n and where the grad(u), grad(k) and Hess(k) are known,
but Hess(u) is estimated via a structured compact approximation
 
Details of the methods are in `[1]`


## Organization
Files are organized in 6 folders:

    ALG/            (algorithms from [2])      
    ALG_COMPACT/    (algorithms from [1])    
    EXAMPLE/        (example)    
    EXTERN_SOLVERS/ (third party software)    
    INTERFACE/      (interfaces to external software)    
    MISC/           (structured objective functions, etc.)    
    REL_EXPERIMENTS/(experiments from [1] for release)

  Note: The main algorithms from the article are in `CSBP3_D.m` (L-S-BFGS-P)
  and in `CSBM3_SP_MF_SCALE.m` (L-S-BFGS-M)

## Example
To try-out an example navigate inside folder "EXAMPLE/".
From within the folder you can run the example

	>> example

The outcomes of running example_1 look like:
```
>> example

Running Problem: P1 
 Running: CSBMSV3_MF_SCALE on    n=1000 
Iter       obj          Resid        ||d||         alpha       #LS UseZoom Reg     PriR
0	0.000000e+00	  3.38e+00	
1	-2.562288e+01	  1.98e+01	3.38e+00	4.95e-02	1	0	0	-
2	-6.977861e+01	  1.80e+01	1.59e-01	1.00e+00	1	0	0	-
3	-2.028126e+02	  1.69e+01	9.95e-01	1.00e+00	1	0	0	-
4	-2.264990e+02	  1.14e+01	5.03e-01	3.48e-01	1	0	0	-
5	-2.361670e+02	  6.01e+00	2.18e-02	1.00e+00	1	0	0	-
6	-2.428693e+02	  2.86e+00	4.56e-02	1.00e+00	1	0	0	-
7	-2.436621e+02	  3.24e+00	2.56e-02	1.00e+00	1	0	0	-
8	-2.442691e+02	  1.76e+00	3.29e-03	1.00e+00	1	0	0	-
9	-2.449016e+02	  1.59e+00	6.74e-03	1.00e+00	1	0	0	-
10	-2.457273e+02	  2.61e+00	1.21e-02	1.00e+00	1	0	0	-
11	-2.473509e+02	  3.43e+00	3.30e-02	1.00e+00	1	0	0	-
12	-2.479824e+02	  4.11e+00	4.18e-02	5.07e-01	1	0	0	-
13	-2.488017e+02	  1.25e+00	2.62e-02	1.00e+00	1	0	0	-
14	-2.489275e+02	  6.57e-01	3.16e-03	1.00e+00	1	0	0	-
15	-2.489737e+02	  4.86e-01	2.41e-03	1.00e+00	1	0	0	-
16	-2.490154e+02	  5.25e-01	2.16e-03	1.00e+00	1	0	0	-
17	-2.490686e+02	  4.16e-01	4.68e-03	1.00e+00	1	0	0	-
18	-2.490787e+02	  1.98e-01	3.53e-03	3.27e-01	1	0	0	-
19	-2.490820e+02	  5.85e-02	8.73e-04	1.00e+00	1	0	0	-
20	-2.490825e+02	  2.83e-02	2.62e-04	1.00e+00	1	0	0	-
.             .             .           .           .       .   .     .       .
.             .             .           .           .       .   .     .       .
65	-2.490868e+02	  4.85e-05	3.17e-07	1.00e+00	1	0	0	-
66	-2.490868e+02	  2.86e-05	2.66e-07	1.00e+00	1	0	0	-
67	-2.490868e+02	  1.66e-05	1.44e-07	1.00e+00	1	0	0	-
68	-2.490868e+02	  1.09e-05	4.83e-08	1.00e+00	1	0	0	-
69	-2.490868e+02	  1.08e-05	6.15e-08	1.00e+00	1	0	0	-
70	-2.490868e+02	  1.20e-05	1.50e-07	1.00e+00	1	0	0	-
71	-2.490868e+02	  1.12e-05	7.13e-08	1.00e+00	1	0	0	-
72	-2.490868e+02	  1.57e-06	3.71e-08	1.00e+00	1	0	0	-
Finished: CSBMSV3_MF_SCALE in time=1.37
```

## References
[1] J.J.Brust, Z.Di, S.Leyffer and C.G.Petra, Compact representations of structured BFGS matrices,
 Comput. Optim. App. 80, 55-88, (2021)

[2] C.G. Petra, N.Y. Chiang, M.Anitescu, A Structured Quasi-Newton Algorithm for Optimizing with Incomplete 
Hessian Information, SIAM J. Optim, 29(2), (2019)
