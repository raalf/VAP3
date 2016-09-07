# VAP
The MATLAB version of the higher-order potential flow code FreeWake.


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

