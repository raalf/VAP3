# VAP 3.5
The MATLAB version of the higher-order potential flow code FreeWake.

This markdown file holds information and records of VAP 3, development of a new input format for future versions of VAP. Similar to previous input files, parameters of the flight vehicle geometry are stored in ASCII text file but in XML 1.0 language. This will provide freedom for future versions of VAP and enable analysis of complex vehicle geometry and configuration. 


# VAP 3.5 File Format Specifications
## Legend
- (\*)  Multiple Entry
- [FW] Parameter exists in FreeWake
- [VAP2] Parameter introduced to VAP2
- [VAP3] Parameter introduced to VAP3

## Input Structure
### Setting
- `Setting`
	- `relax` (True | False) Toggle wake relaxation or fixed wake
	- `steady` (True | False)
	- `maxtime` Maximum number of time steps
	- `delta_time` Inteval of each time step (second)
	- `start_force` Timestep number to begin force calculations
	- `stiff_wing` (1 | 0) Toggle on (1) or off (0) structure model
	- `fixed_lift` (1 | 0) Toggle fixed lift where lift = weight and freestream is calcualted
	- `gust_mode` Select gust mode, 0 - No gust, 1 - Sine wave gust, 2 - 1-cosine gust, 3 - Sharp edge gust

### Condition
- `condition`
	- `density` Fluid density (kg/m^3)
	- `kin_viscosity` Kinematic viscosity (m^2/s)
	- `gust_amplitude`
	- `gust_length`
	- `gust_start`

### Vehicle
- `vehicle*`
	- `global_x` X position of the entire vehicle (Position of vehicle's origin in global reference frame)
	- `global_y` Y position of the entire vehicle 
	- `global_z` Z position of the entire vehicle
	- `weight` Weight of vehicle (N)
	- `interference_drag` Percent of interference drag (often due to fuselage)
	- `speed` Freestream velocity of the vehicle (m/s)
	- `alpha` Overall angle of attack of the entire vehicle (deg)
	- `beta` Overall angle of sideslip of the entire vehicle (deg)
	- `roll` Overall angle of roll of the entire vehicle (deg)
	- `fpa` FPA (Flight Path Angle) = Pitch - AOA (deg)
	- `track` Track heading of the entire vehicle (deg)
	- `radius` Turning radius of the vehicle (zero = straight) (m) 
	- `ref_area` Reference area of the wing (m^2)
	- `ref_span` Reference span of the wing (m)
	- `ref_cmac` Mean aerodynamic chord of the wing (m)
	- `wing*`
		- `symmetry` (True | False) 
		- `incidence` Angle of incidence of the root (deg)
		- `trimable` Defines a surface that can be trimmed (True | False)
		- `triangular_elements` Uses triangular elements, use with caution (True | False)
		- `chordwise_elements` Number of chordwise elements
		- `vehicle_x` X reference position of wing
		- `vehicle_y` Y reference position of wing
		- `vehicle_z` Z reference position of wing
		- `panel*`
			- `spanwise_elements` Number of spanwise elements in this panel
			- `strip_airfoil` Airfoil name for this panel
			- `section*`
				- `wing_x` X position of section leading edge with respect to vehicle_x (m)
				- `wing_y` Y position of section leading edge with respect to vehicle_y (m)
				- `wing_z` Z position of section leading edge with respect to vehicle_z (m)
				- `chord` Chord length of section (m)
				- `twist` Twist angle of section (deg)
	- `rotor*`
		- `rpm` Rotational speed of rotor (RPM) 
		- `collective` Collective blade pitch (deg)
		- `ref_diam` Rotor diameter (m) 
		- `rotation_direction` Rotational direction, CCW - Counterclockwise, CW - clockwise
		- `veh_x_hub` X location of center of rotor with respect to vehicle origin (m)
		- `veh_y_hub` Y location of center of rotor with respect to vehicle origin (m)
		- `veh_z_hub` Z location of center of rotor with respect to vehicle origin (m)
		- `veh_x_axis` Axis of rotation unit vector X-component with respect to vehicle origin (m)
		- `veh_y_axis` Axis of rotation unit vector Y-component with respect to vehicle origin (m)
		- `veh_z_axis` Axis of rotation unit vector Z-component with respect to vehicle origin (m)
		- `blades` Number of blades
		-  `chordwise_elements` Number of chordwise elements
		- `panel*`
			- `strip_airfoil` Airfoil name for this panel
			- `spanwise_elements` Number of spanwise elements
			- `section*`
				- `rotor_x` X position of section leading edge with respect to `hub` and XY plane (m)
				- `rotor_y` Y position of section leading edge with respect to `hub` and XY plane (m)
				- `rotor_z` Z position of section leading edge with respect to `hub` and XY plane (m)
				- `chord` Chord length of section (m)
				- `twist` Twist angle of section (deg)

## Notes
### Default Unit
All numerical input would have a default unit assigned. (e.g. meters for the reference span). This can be setup with a `unit` tag under the numerical parameter such that specific units can be converted to the default unit in ***future*** VAP releases.

For example:
`<span unit='feet'>5</span>` 

would yield the same result as
`<span>1.524</span>`
	
	
---