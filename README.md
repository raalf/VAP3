# VAPXML
This markdown file holds infomation and records of VAPXML, development of a new input format for future versions of VAP. Similar to previous input files, parameters of the flight vehicle geometry are stored in ASCII text file but in XML 1.0 language. This will provide freedom for future versions of VAP and enable analysis of complex vechicle geometry and configuration. 


# VAPXML File Format Specifications
## Legend
- (\*)  Multiple Entry
- [FW] Parameter exists in FreeWake
- [VAP2] Parameter introduced to VAP2
- [VAP3] Parameter introduced to VAP3

## Notes
### Default Unit
All numerical input would have a default unit assigned. (e.g. meters for the reference span). This can be setup with a `unit` tag under the numerical parameter such that specific units can be converted to the default unit in *future* VAP releases.

For example:
`<span unit='feet'>5</span>` 

would yield the same result as
`<span>1.524</span>`


## VAP Settings
- `VAP`
	- `flagRELAX` (True | False)
	- `flagSTEADY` (True | False)
	- `flagTRI` (True | False)
	- `valMAXTIME` Maximum number of time steps
	- `valMINTIME` Minimum number of time steps
	- `valDELTIME` Inteval of each time step (second)
	- `valDELTAE` Convergence delta-span effic. (0 if only timestepping)



## Condition
- `condition`
	- `valDENSITY` Fluid density (kg/m^3)
	- `valKINV` Kinematic viscosity (m^2/s)


## Vehicle
- `vehicle*`
	- `alpha` Overall angle of attack of the entire vehicle (deg)
	- `beta` Overall angle of sideslip of the entire vehicle (deg)
	- `roll` Overall angle of roll of the entire vehicle (deg)
	- `pitch` Overall angle of pitch of the entire vehicle (deg)
	- `yaw` Overall angle of yaw (Heading) of the entire vehicle (deg)
	- `radius` Turning radius of the vehicle (zero = straight) (m)
	- `vinf` Freestream velocity of the vehicle (m/s)
	- `x` X position of the entire vehicle (Position of vehicle's origin in global reference frame)
	- `y` Y position of the entire vehicle 
	- `z` Z position of the entire vehicle 

	
	- `wing*`
		- `incidence` Angle of incidence of the root (deg)
		- `trimable` Defines a surface that can be trimmed (True | False)
		- `area` Reference area of the wing (m^2)
		- `span` Reference span of the wing (m)
		- `cmac` Mean aerodynamic chord of the wing (m)
		- `panel*`
			- `symmetry` Symmetry edge (0, 1 or 2)
			- `N` Number of spanwise elements
			- `M` Number of chordwise elements
			- `section*`
				- `X` X position of section leading edge with respect to vehicle origin(m)
				- `Y` Y position of section leading edge with respect to vehicle origin(m)
				- `Z` Z position of section leading edge with respect to vehicle origin(m)
				- `chord` Chord length of section (m)
				- `epsilon` Twist angle of section (deg)
				- `airfoil`

	- `rotor*`
		- `RPM` Rotational speed of rotor (RPM) 
		- `dia` Rotor diameter (m) 
		- `xhub` X location of center of rotor with respect to vehicle origin (m)
		- `yhub` Y location of center of rotor with respect to vehicle origin (m)
		- `zhub` Z location of center of rotor with respect to vehicle origin (m)
		- `axisx` Axis of rotation unit vector X-component with respect to vehicle origin (m)
		- `axisy` Axis of rotation unit vector Y-component with respect to vehicle origin (m)
		- `axisz` Axis of rotation unit vector Z-component with respect to vehicle origin (m)
		- `blades` Number of blades
		- `panel*`
			- `N` Number of spanwise elements
			- `M` Number of chordwise elements
			- `section*`
				- `X` X position of section leading edge with respect to `hub` and XY plane (m)
				- `Y` Y position of section leading edge with respect to `hub` and XY plane (m)
				- `Z` Z position of section leading edge with respect to `hub` and XY plane (m)
				- `chord` Chord length of section (m)
				- `epsilon` Twist angle of section (deg)
				- `airfoil`
	- `fuselage*`
	- `viscous`

				
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
