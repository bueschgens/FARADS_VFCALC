# FARADS_VFCALC

## FARADS - Fast and Accurate calculation of RADiation heat transfer between arbitrary three-dimensional Surfaces

FARADS provides a framework for calculating view factors of arbitrary three-dimensional geometries. With these view factors the additional application of the net radiation method for calculating exchanged heat fluxes can be used.

other modules of FARADS:
- FARADS_GEOM
- FARADS_PLOT
- FARADS_VFCALC
- FARADS_QRAD

This module provides functions and types for:
- calculating view factors
- checking for blocking between element pairs (ray triangle interection based on [1])

This is core of the FARADS framework.

[1] Möller, T.; Trumbore, B.: "Fast, Minimum Storage Ray/Triangle Intersection", Journal of Graphics Tools, vol. 2(1) (1997), 21-28
