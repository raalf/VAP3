# VAPXML
This markdown file holds infomation and records of VAPXML, development of a new input format for future versions of VAP. Similar to previous input files, parameters of the flight vehicle geometry are stored in ASCII text file but in XML 1.0 language. This will provide freedom for future versions of VAP and enable analysis of complex vechicle geometry and configuration. 


# VAPXML File Format Specifications
## Legend
- (\*)  Multiple Entry
- [FW] Parameter exists in FreeWake
- [VAP2] Parameter introduced to VAP2
- [VAP3] Parameter introduced to VAP3


## Condition
- `condition`


## Vehicle
- `vehicle*`
	- `metrics`
		- `area` Reference area of the vehicle (m)
		- `span` Reference span of the vehicle (m)
		- `cmac` Mean aerodynamic chord of the vehicle (m)
	- `wing*`
		- `panel*`
			- `symmetry` Symmetry edge (0, 1 or 2)
			- `N` Number of spanwise elements
			- `M` Number of chordwise elements
			- `section*`
				- `X` X position of section leading edge (m)
				- `Y` Y position of section leading edge (m)
				- `Z` Z position of section leading edge (m)
				- `chord` Chord length of section (m)
				- `epsilon` Twist angle of section (deg)
				- `airfoil`

				
---

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
