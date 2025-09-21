# 3D Reciprocal-Space Intensity Reconstruction

This repository contains MATLAB code for reconstructing a three-dimensional intensity grid in reciprocal space from diffraction images.
The method includes detector corrections, symmetry operations, and Lorentz-factor weighting, and is designed to reconstruct the intensity distribution in reciprocal space.

IMPORTANT: This code is provided for general information only. At present, it is configured for specific diffractometer geometries (ESRF ID28 and the four-circle X-ray diffractometer at Tel Aviv University).

The code assumes that the raw experimental data are stored as .mat files. We provide helper scripts to convert CBF files to MAT files; this conversion must be performed before running the reconstruction.
The CBF→MAT conversion relies on cbfread.m (distributed by Oliver Bunk and available online).

---

## Main Script

The central routine is:

**`IntensityGrid_plusSymmetry_parallel_Lorentz_V3.m`**  
Located in: `src/`

**Purpose:**  
Reconstructs the 3D intensity distribution in reciprocal space by combining diffraction images with the instrument geometry.

**Inputs:**  
FolderName : Path to the folder containing the experimental data (MATLAB .mat files).
UB (3x3 double): Crystal orientation matrix. Each column gives the reciprocal-lattice basis vectors expressed in the laboratory coordinate system. 
P0 (1x6 double): Geometry of the experiment 
   P0(1) sample-to-detector distance (mm)
   P0(2:3) X and Y coordinates (pixels) of the direct beam at zero detector angle
   P0(4) detector pixel size (mm); for a PILATUS detector this is 0.172 mm
   P0(5) detector rotation angle; set to 0 to read from the image header
   P0(6) wavelength (Å)   

hL, Nh: range and number of intervals along a*
kL, Nk: range and number of intervals along b*
lL, Nl: range and number of intervals along c* 

SymMatrix (3x3xNsym double): Stack of symmetry operators. For the m-3m point group, use a 3×3×48 array, it is available for loading in variables folder.


**Outputs:**  
- `hgrid, kgrid, lgrid` — voxel-grid centers in reciprocal space  
- `IG` — voxel-wise average intensity grid 
- `sigmaIG` — voxel-wise intensity uncertainty  
- `VG` — voxel volumes obtained as the sum of all collected pixels falling into the voxel  
- `CG` — number of contributing pixels per voxel  

---

**Minimal usage example**

addpath('src');

FolderName = 'path/to/mat/files';

UB = [ ... % 3x3 orientation matrix ];

P0 = [244, 516, 527,  0.172, 0, 0.6968];

hL = [-2, 2];  Nh = 201;
kL = [-2, 2];  Nk = 201;
lL = [-2, 2];  Nl = 201;

load('variables/SymMatrix.mat', 'SymMatrix'); % 3x3x48

[hgrid,kgrid,lgrid,IG,sigmaIG,VG,CG] = ...
    IntensityGrid_plusSymmetry_parallel_Lorentz_V3( ...
        FolderName, UB, P0, hL, Nh, kL, Nk, lL, Nl, SymMatrix);



