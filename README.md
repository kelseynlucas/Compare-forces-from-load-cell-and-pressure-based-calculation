# Compare-forces-from-load-cell-and-pressure-based-calculation
Scripts used to compare forces and torques measured by a force-torque sensor on a flapping foil apparatus, and those calculated based on simultaneously-collected PIV velocity data.  Velocity data were used to calculate pressure fields with queen2 (http://dabirilab.com/software/) and then forces and torques (https://github.com/kelseynlucas/Pressure-based-force-calculation-for-foils).

Compare_Flapper_Pcode.py - the Python script used to perform cross-correlation and RMSE calculations. Reads in data files containing the forces from a flapping foil apparatus, and forces derived from pressure-based calculations (https://github.com/kelseynlucas/Pressure-based-force-calculation-for-foils).

For preparing data for comparisons:

downsample_flapper_data.py codes read in data from the flapping foil apparatus, and write out every 10th data point (because equal sampling was required for cross-correlation, and pressure-based calculations were performed at 100Hz, compared to the flapping foil apparatus' 1000Hz).
