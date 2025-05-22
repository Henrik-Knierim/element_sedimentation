# helium_rain
An implementation of helium rain into MESA.

## Quickstart
Copy `helium_rain.f90` and `Schötter+2018.csv` into the `src` directory your MESA project. Then, add the functions and parameters of the `run_star_extras.f90` file in this repo to the `run_star_extras.f90` file of your MESA project.

Next, add the compiling info of `helium_rain.f90` to your makefile. An example of how to do this can be found in `makefile` in this repo.

## Details
The code computes an array $dY$, that is the change in the helium mass fraction due to sedimentation. Once obtained, the mass abundances are changed according to
$$dX_i = \frac{- dY}{1-Y} X_i$$
where $X_i$ is the mass fraction of the $i$-th chemical species and $Y$ is the helium mass fraction.