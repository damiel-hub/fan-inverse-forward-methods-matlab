# Inverse Computational Morphology of Debris and Alluvial Fans

The forward method, developed by [Chen and Capart (2022)](https://doi.org/10.1016/j.cageo.2022.105228), simulates alluvial fan topography using a fast algorithm based on the visibility polygon. This method requires initial terrain data, apex location, and the elevation-distance relationship (elevation profile) as inputs to create a detailed topographic model of the alluvial fan.

This project focuses on the inverse approach, which extracts the elevation-distance relationship from fan topography using the visibility polygon. This extracted profile can then be used as input for the forward method to accurately reconstruct the fan topography.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Example Commands](#example-commands)
- [More Functions and Examples](#more-functions-and-examples)
- [How to Cite](#how-to-cite)
- [Reference](#reference)


## Installation

### Prerequisites

Before you begin, ensure you have the following:

- **MATLAB** (version R2021a or later recommended)
- **Required Toolbox**: Mapping Toolbox, Parallel Computing Toolbox


## Usage

### Basic Usage

To use this project, follow these steps:

1. Add the 'functions' folder to your MATLAB path.
2. Navigate to the main directory of the project.
3. Run the demo script by entering the following command in the MATLAB Command Window:

    ```matlab
    run('demo_inverse_forward_methods.m')
    ```

The demo script includes:
- Loading example GeoTiff data files: `topo_post_event` (alluvial fan terrain), `shape_fan_boundary` (fan boundary polygon shapefile), and `topo_pre_event` (initial terrain).
- Clipping the alluvial fan topography using the boundary polygon.
- Mapping shortest path distances within and along the boundary.
- Applying median filtering and quadratic fitting to the elevation-distance relationship.
- Using the forward method to reconstruct the fan topography with the extracted profile.

## Example Commands

**Clip the Fan Topography:**

```matlab
[xMesh_crop, yMesh_crop, zMesh_post_crop] = clipGeoTiff(topo_post_event, shape_fan_boundary);
```
- Input: 
  - `topo_post_event`: Path to the GeoTIFF file containing the alluvial fan topography.
  - `shape_fan_boundary`: Path to the shapefile defining the boundary of the fan.
- Output:
  - `xMesh_post_crop`: X mesh coordinates of the clipped alluvial fan topography.
  - `yMesh_post_crop`: Y mesh coordinates of the clipped alluvial fan topography.
  - `zMesh_post_crop`: Z mesh coordinates (elevation) of the clipped alluvial fan topography.



**Map the Shortest Path Distance Within the Boundary:**

```matlab
sMap = shortest_path_distance_within_boundary(xMesh_crop, yMesh_crop, zMesh_post_crop, 1);
```

- **Inputs:**
  - `xMesh_crop`, `yMesh_crop`, `zMesh_post_crop`: Clipped mesh coordinates from the previous step.
  - `0` or `1`: Specifies whether to plot the shortest path results. Use `0` to skip plotting and `1` to display the plot.
- **Outputs:**
  - `sMap`: A 2D matrix of shortest path distances within the boundary of the alluvial fan.

**Map the Shortest Path Distance Along the Boundary:**

```matlab
xysBoundary = shortest_path_distance_along_boundary(xMesh_crop, yMesh_crop, zMesh_post_crop);
```

- **Inputs:**
  - `xMesh_crop`, `yMesh_crop`, `zMesh_post_crop`: Clipped mesh coordinates from the previous step.

- **Outputs:**
  - `xysBoundary`: An `n × 3` matrix with X and Y coordinates of the boundary in the first two columns, and shortest path distances in the third column.


**Median Filtering and Quadratic Fitting:**

```matlab
% With median filter applied
fitting_s_z_within_boundary = process_s_z_relationship(sMap, zMesh_post_crop, bin_size, ds, outlength, 1);

% Without median filter
fitting_s_z_along_boundary = process_s_z_relationship(xysBoundary(:,3), zBoundary, bin_size, ds, outlength, 1, 'medianFilter', 0);
```

- **Inputs:**
  - `sMap` and `xysBoundary(:,3)`: 2D matrix or 1D vector representing shortest path distances.
  - `zMesh_post_crop` and `zBoundary`: Elevation values of the fan topography correspoding to the previous parameter.
  - `bin_size`: Size of bins for filtering useless when without median filter.
  - `ds`: Sampling distance.
  - `outlength`: Length of the straight-line fitting beyond the fan toe and apex.
  - `0` or `1`: Controls plotting. Use `0` to skip the plot, `1` to display it.
  - `'medianFilter'`: Controls median filtering. Use `0` to disable, `1` or omit to enable.

- **Outputs:**
  - `fitting_s_z_within_boundary` and `fitting_s_z_along_boundary`: Processed elevation-distance relationships fitted with a quadratic function.



**Reconstruct the Fan Using the Forward Method:**

```matlab
[zTopo_sim, heightAG_Volume_All] = reconstruct_fan_surface(xMesh, yMesh, zMesh_pre, xApex, yApex, fanSimVolume, guessHeightAboveGround_top, guessHeightAboveGround_bottom, fitting_s_z_within_boundary, "fanBoundarySHP", shape_fan_boundary);
```

- **Inputs:**
  - `xMesh`, `yMesh`, `zMesh_pre`: Mesh grid coordinates (`xMesh`, `yMesh`) and initial elevation data (`zMesh_pre`).
  - `xApex`, `yApex`: Coordinates of the fan's apex point.
  - `fanSimVolume`: Expected volume of the simulated fan.
  - `guessHeightAboveGround_top`: Initial estimate of the apex height above ground for the upper boundary.
  - `guessHeightAboveGround_bottom`: Initial estimate of the apex height above ground for the lower boundary.
  - `fitting_s_z_within_boundary`: 2D array representing the distance-elevation relationships obtained in the previous step.
  - `"fanBoundarySHP"`: (string) A name-value pair where `"fanBoundarySHP"` specifies the boundary shapefile `shape_fan_boundary` for volume calculation. Default is `nan`.
  - `"tol"`: (scalar) Tolerance for the difference between the expected volume ($V_e$) and the simulated volume ($V_s$). The simulaiton stops when the condition $\text{abs}(V_s - V_e) / V_e \leq \text{tol}$ is met. Default is `0.03`.
  - `"debug"`: (boolean) If set to `1`, saves two guessed fan topographies as GeoTIFF files. Default is `0`.
  - `"epgs_code"`: (scalar) The EPSG code to use when saving GeoTIFF files. This is only applicable if `"debug"` is set to `1`. Default is `3826`.

- **Outputs:**
  - `zTopo`: Reconstructed topography of the alluvial fan.
  - `heightAG_Volume_All`: 2D array containing the apex height above ground and the corresponding calculated fan volume for each iteration.


## More Functions and Examples

**Fan topography simulation (Forward method)**

```matlab
[zTopo,kTopoAll,xyzkApexAll,xyzVisPolygon,xyVisPolygonAll,thetaMesh] = FanTopo(xMesh,yMesh,zMesh,xApexM,yApexM,zApexM,options)
```
- **Inputs:**
  - `xMesh`, `yMesh`, `zMesh`: Mesh grid coordinates (`xMesh`, `yMesh`) and the initial elevation data (`zMesh`).

  - `xApexM`, `yApexM`, `zApexM`: Vectors containing the coordinates (`xApexM`, `yApexM`) and elevation (`zApexM`) of the apexes for multiple fans, allowing for the simultaneous simulation of multiple fan surfaces.



  - **`options`**: A structure containing optional parameters that customize the fan morphology and simulation process. The fields in `options` depend on the chosen `caseName` and can include the following:

    - **`caseName`**: (string) Specifies the type of fan morphology to generate. Possible values include:
      - **`'cone'`**: Generates a conical fan shape.
        - **`tanAlphaM`**: (vector) Defines the slope angles (tangents) for each apex, which determine the steepness of the fan.
      - **`'concave'`**: Generates a concave-shaped fan.
        - **`tanAlphaM`**: (vector) Defines the slope angles (tangents) for each apex.
        - **`KM`**: (vector) Concavity factors for each apex, which control the curvature of the fan.
      - **`'infinite'`**: Generates a fan with an infinite slope.
        - **`tanAlphaM`**: (vector) Defines the slope angles (tangents) for each apex.
        - **`KM`**: (vector) Concavity factors for each apex.
        - **`tanInfiniteM`**: (vector) Slope values for cases where the tangent approaches infinity, useful for modeling extreme or vertical slopes.
      - **`'myProfile'`**: Uses custom profiles for fan generation.
        - **`dz_interpM`**: (cell array) Contains interpolation values for elevation, used to create spline-based morphologies for the fan.

    - **`dispflag`**: (scalar) A flag to control the display of the generated topography. Set to `1` to visualize the topography, or `0` to skip the visualization. Default is set as `0`.

    - **`saveVisPolygon`**: (scalar) A flag to control whether the visibility polygons should be saved. Set to `1` to save the polygons, or `0` to skip saving them. Default is set as `0`.

    This structure allows for flexible fan generation, enabling the simulation of various fan types by adjusting parameters specific to each morphology type.
- **Outputs:**
  - `zTopo`: 2D matrix of final fan topography (elevation after aggradation).
  - `kTopoAll`: 2D matrix with indices of the apex dominating each mesh grid point.
  - `xyzkApexAll`: Matrix of apex coordinates and indices (including child apexes).
  - `xyzVisPolygon`: Cell array of 3D coordinates (`x`, `y`, `z`) for visibility polygons. Only generated if `saveVisPolygon` is set to `1`.
  - `xyVisPolygonAll`: Matrix of `x` and `y` coordinates for all visibility polygons. Only generated if `saveVisPolygon` is set to `1`.
  - `thetaMesh`: 2D matrix showing the angle of each point relative to the apex(es).




## How to Cite

If you use this project in your research or work, please give us credit by citing it with the following BibTeX entries:

### Inverse method: extracting the distance-elevation relationship
```
Chiu, Y.H. & Capart, H. (In process). *Inverse Computational Morphology of Debris and Alluvial Fans*. Retrieved from [https://github.com/damiel-hub/geometric-alluvial-fan].
```
### Forward method: simulating fan topography

```bibtex
@article{chenComputationalMorphologyDebris2022,
title = {Computational morphology of debris and alluvial fans on irregular terrain using the visibility polygon},
journal = {Computers & Geosciences},
volume = {169},
pages = {105228},
year = {2022},
issn = {0098-3004},
doi = {https://doi.org/10.1016/j.cageo.2022.105228},
url = {https://www.sciencedirect.com/science/article/pii/S0098300422001777},
author = {Tzu-Yin Kasha Chen and Hervé Capart},
keywords = {Morphology, Debris flow, Alluvial fan, Surface of revolution, Eikonal equation, Visibility polygon}
}
```


## Reference

- Chen, T.Y.K. & Capart, H. (2022). [Computational morphology of debris and alluvial fans on irregular terrain using the visibility polygon](https://doi.org/10.1016/j.cageo.2022.105228). *Computers & Geosciences*, *169*, 105228.


