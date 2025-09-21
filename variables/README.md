# Variables

This folder contains auxiliary `.mat` files used by the reconstruction routines.  
These files provide symmetry operations and detector-specific masks for hot pixels.

---

## Files

### `SymMatrix.mat`
- **Contents:** A 3×3×48 array.
- **Description:** Each 3x3 slice is a symmetry matrix of the *m-3m* cubic point group (48 operations in total).
- **Usage:** Supplied as the `SymMatrix` argument when calling  
  `IntensityGrid_plusSymmetry_parallel_Lorentz_V3`.

### `HotPixels_ID28.mat`
- **Contents:** List of detector pixel positions corresponding to *hot pixels*.
- **Description:** Used to mask out bad pixels during intensity reconstruction.  
  The provided file is specific to the **ID28 beamline**.
- **Customization:**  
  - You can create your own hot-pixel mask by listing the `(X, Y)` pixel positions for your instrument.  
  - Save it in a `.mat` file with the same variable structure and replace this file.
  - If your detector does not have any hot pixels, please save this variable as an empty matrix. 

---

## Notes
- Both files are small but essential for reproducible analysis.  
- Place them on the MATLAB path or adjust your scripts to load them directly from this folder.  