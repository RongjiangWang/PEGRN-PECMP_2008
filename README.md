FORTRAN code for calculating co- and post-seismic deformation in multi-layered poroelastic half-space based on the poroelastic dislocation theory. The code structure is similar as the viscoelastic variant PSGRN-PSCMP.

Highlights:

(1) orthonormal propagator algorithm for numerical stability

(2) finite fault model

(3) azimuthal-equidistant projection between spherical and local cartesian coordinates

(4) output of complete geophysical observables (time-dependent displacement, stress/strain, tilt, block rotation, pore pressure and Darcy flux)

For Windows user, the executable file is provided under folder "WindowsEXE". Linux user may compile the source codes with "gfortran" via a single command like, e.g.,

~>cd .../PEGRN2020/SourceCode

~>gfortran -o pEgrn2020 *.f -O3

to get the excutable code pegrn2020.

After start the executable code, the program ask for an input file in the ASCII format. An example input file is provided under folder "InputFile". You may change the input data included in this file for your own applications.

References

Wang, R. and H.-J. Kümpel (2003). Poroelasticity: Efficient modelling of strongly coupled, slow defor-mation processes in a multilayered half-space, Geophysics, 68(2), 705-717

Wang, R., F. Lorenzo-Martin and F. Roth (2006). PSGRN/PSCMP - a new code for calculating co- and post-seismic deformation, geoid and gravity changes based on the viscoelastic-gravitational dislocation theory, Computers and Geosciences, 32, 527-541. DOI:10.1016/j.cageo.2005.08.006
