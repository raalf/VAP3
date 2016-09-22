# VAP
The MATLAB version of the higher-order potential flow code FreeWake.

Branches:
- Master
    - Working version of VAP2.0
    - Optimized to make use of MATLAB's features

- VAP1.0
    - Contains files as copied from the C++ version of FreeWake into MATLAB
	- Provides identical results to the C++ version of FreeWake, but runs slow


TO-DO:

- Remove structures
- Vectorize functions

We can implement some of the low-level functions such as vortex sheet induction from HVAP which are written well. 

Areas to focus on:

- DVE creation
    - generateDVEs_v2
    - generateSingfct
    - panelBoundary
    - generateWakeDVEs

- D-Matrix
    - fcnAssembleWingD
    - fcnDVEKinCond
    - fcnNewWakeVorticity
    - fcnSolveDWakeMatrix

