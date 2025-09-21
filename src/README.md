# Source Code (`src/`)

This directory contains the MATLAB functions that implement the reconstruction of a 3D intensity grid in reciprocal space.  
Before running these functions, add this folder to your MATLAB path:

addpath(genpath('src'));

Description of the functions:

Main Function
	IntensityGrid_plusSymmetry_parallel_Lorentz_V3.m
	Purpose: Reconstructs a 3D intensity grid in reciprocal space.

	Features:
		Applies detector corrections (solid angle, polarization, Lorentz factor).
		Incorporates symmetrization.
		Handles hot-pixel masking.
		Supports parallel execution with parfor.

	Outputs:
		hgrid, kgrid, lgrid – reciprocal-space voxel grid edges
		IG – intensity grid
		sigmaIG – voxel-wise intensity uncertainty
		VG – total voxel volume calculated as the sum of all the pixels contributing into the voxel 
		CG – number of contributing pixels per voxel

Supporting Functions
	prepare_voxel_grid.m
		Creates the 3D grid in h, k, l space.
		Defines voxel sizes and grid boundaries.

	LoadMATImage.m
		Loads raw detector images saved in .mat format.
		Applies the hot-pixel mask.
		Returns intensity and geometry metadata.

	adjust_geometry.m
		Adjusts detector geometry parameters (P0) using instrument settings.

	Pixel_to_XYZ_TAU.m
		Converts detector pixel coordinates (X, Y) into scattering vectors in the lab frame.
		Calculates additional geometry (e.g., oblique angles, rotation matrices).

	Lorentz.m
		Computes the Lorentz correction factor and voxel volume contributions.