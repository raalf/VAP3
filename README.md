# VAP
The MATLAB version of the higher-order potential flow code FreeWake.

Branches:
- Master
    - Working version of VAP2.0
    - Optimized to make use of MATLAB's features

- VAP1.0
    - Contains files as copied from the C++ version of FreeWake into MATLAB
	- Provides similar results to the C++ version of FreeWake, but runs slow
    - Known issues with VAP1.0:
        - If relaxing the wake, the very first row of elements emitted (timestep 0) are not connected to the second row. 
        - Splits in the wings are not handled correctly.
