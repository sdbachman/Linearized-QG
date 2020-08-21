# Linearized-QG
A 20-layer linearized QG model

This Dedalus script creates a 512 x 512 20-layer quasigeostrophic model, where the nonlinear advection, dissipation, and bottom drag are omitted.  This was used to output growth rates and eddy flux angles in Bachman (2020, Fluids).

This script writes output every 100 time steps.  The dynamics are mainly governed by the planetary vorticity gradient, "beta", whose value is set by the Python controller script, "controller.py". The horizontal resolution is set via the choice of the layerwise deformation radius, "Ld".

This repository comes with the Python controller script and a few MATLAB scripts used to produce the plots in the aforementioned paper.
