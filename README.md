# IntensityOnlyImagingwMUSIC
Collection of Matlab codes used to compute the examples shown in the paper, "Intensity-only inverse scattering using MUSIC" by A. D. Kim and C. Tsogka (2019).

The codes are currently set to compute the 3 ellipsoids problem.

1. ForwardProblem3D.m computes the solution of the forward problem using the Method of Fundamental Solutions for the suite of experiments corresponding to different incident fields. Because this code computes the solution of the linear system directly using "backslash", it can be slow when considering many scattering obstacles.

This main driver code requires the following supporting codes.

* Ellipsoid.m -- This code computes a random distribution of points on an ellipsoid and their corresponding unit, outward normals used by the Method of Fundamental Solutions.

 * ComputeMFSCoeffs.m -- This code assembles and solves the linear system of equations for the Method of Fundamental Solutions expansion coefficients.
 
 * ComputeMeasurements.m -- This code uses the Method of Fundamental Solutions expansion coefficients to compute the scattered fields on the measurement plane.
 
The result of running this code is a MAT file containing the data to be used by the inverse scattering code.

2. MUSICImaging3D.m uses the results from ForwardProblem3D.m to compute the imaging procedure given in Section 3 of "Intensity-only inverse scattering using MUSIC" by A. D. Kim and C. Tsogka and plots the results similarly to what is shown in Section 5 of that paper.

This main driver code requires the following supporting code.

* ellipsoid_fit.m -- This code computes a surface of the ellipsoidal scattering objects to be reconstructed. This surface is plotted along with the isosurface of the imaging result to compare results.

Also included in this repository is the sample MAT file called "ForwardData.mat" to run "MUSICImaging3D.m" without having to run "ForwardProblem3D.m".
