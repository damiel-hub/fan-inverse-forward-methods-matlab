# Inverse Computational Morphology of Confined Debris and Alluvial Fans Using the Visibility Polygon

The forward method, developed by [Chen and Capart (2022)](https://doi.org/10.1016/j.cageo.2022.105228), simulates alluvial fan topography using a fast algorithm based on the visibility polygon. This method requires initial terrain data, apex location, and the elevation-distance relationship (elevation profile) as inputs to create a detailed topographic model of the alluvial fan.

This project focuses on the inverse approach, which extracts the elevation-distance relationship from fan topography using the visibility polygon. This extracted profile can then be used as input for the forward method to accurately reconstruct the fan topography.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Example Commands](#example-commands)
- [Cite This](#cite-this)
- [Reference](#reference)


## Installation

### Prerequisites

Before you begin, ensure you have the following:

- **MATLAB** (version R2021a or later recommended)
- **Required Toolbox**: Mapping Toolbox

### Installing Additional Dependencies

1. Download the `inpoly` function from [GitHub](https://github.com/dengwirda/inpoly).
2. Extract the ZIP file and place the `inpoly-master` folder in your project directory.
3. Add the `inpoly-master` folder to your MATLAB path using the following command:

    ```matlab
    addpath(genpath('path_to_your_project/inpoly-master'))
    ```

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
- Loading example data files: `"topo_fan.tif"` (alluvial fan terrain), `"fan_extent.shp"` (fan boundary polygon shapefile), and `"topo_initial.tif"` (initial terrain).
- Clipping the alluvial fan topography using the boundary polygon.
- Mapping shortest path distances within and along the boundary.
- Applying median filtering and quadratic fitting to the elevation-distance relationship.
- Using the forward method to reconstruct the fan topography with the extracted profile.

## Example Commands

**Clip the Fan Topography:**

```matlab
[xMesh_post_crop, yMesh_post_crop, zMesh_post_crop] = clipGeoTiff('data/topo_fan.tif', 'data/shape/fan_extent.shp');
```
- Input: 
  - `'data/topo_fan.tif'`: Path to the GeoTIFF file containing the alluvial fan topography.
  - `'data/shape/fan_extent.shp'`: Path to the shapefile defining the boundary of the fan.
- Output:
  - `xMesh_post_crop`: X mesh coordinates of the clipped alluvial fan topography.
  - `yMesh_post_crop`: Y mesh coordinates of the clipped alluvial fan topography.
  - `zMesh_post_crop`: Z mesh coordinates (elevation) of the clipped alluvial fan topography.



**Map the Shortest Path Distance Within the Boundary:**

```matlab
sMap = shortest_path_distance_within_boundary(xMesh_post_crop, yMesh_post_crop, zMesh_post_crop, 0);
```

- **Inputs:**
  - `xMesh_post_crop`, `yMesh_post_crop`, `zMesh_post_crop`: Clipped mesh coordinates from the previous step.
  - `0` or `1`: Specifies whether to plot the shortest path results. Use `0` to skip plotting and `1` to display the plot.
- **Outputs:**
  - `sMap`: A 2D matrix of shortest path distances within the boundary of the alluvial fan.

**Map the Shortest Path Distance Along the Boundary:**

```matlab
xysBoundary = shortest_path_distance_along_boundary(xMesh_post_crop, yMesh_post_crop, zMesh_post_crop);
```

- **Inputs:**
  - `xMesh_post_crop`, `yMesh_post_crop`, `zMesh_post_crop`: Clipped mesh coordinates from the previous step.

- **Outputs:**
  - `xysBoundary`: An `n × 3` matrix with X and Y coordinates of the boundary in the first two columns, and shortest path distances in the third column.

**Median Filtering and Quadratic Fitting:**

```matlab
fitting_s_z_within_boundary = process_s_z_relationship(sMap, zMesh_post_crop, bin_size, ds, outlength, 1);
```

- **Inputs:**
  - `sMap`: A 2D matrix of shortest path distances within the boundary.
  - `zMesh_post_crop`: Elevation values of the clipped fan topography.
  - `bin_size`: Size of the bins used for filtering.
  - `ds`: Sampling distance.
  - `outlength`: The length of the straight line fitting beyond the fan toe and apex.
  - `0` or `1`: Indicates whether to plot the median filter and fitting result. Use `0` to skip plotting and `1` to show the plot.

- **Outputs:**
  - `fitting_s_z_within_boundary`: Processed elevation-distance relationship fitted with a quadratic function.

**Reconstruct the Fan Using the Forward Method:**

```matlab
[zTopo, ~, ~, ~, ~, ~] = FanTopo_slope_bd(xMesh_pre, yMesh_pre, zMesh_pre, xApex, yApex, zApex, 'caseName', 'myProfile', 'dz_interpM', {fitting_s_z_within_boundary});
```

- **Inputs:**
  - `xMesh_pre`, `yMesh_pre`, `zMesh_pre`: Initial mesh coordinates and elevation data of the fan before reconstruction.
  - `xApex`, `yApex`, `zApex`: Coordinates of the fan’s apex.
  - `'caseName'`: Identifier for the case, such as `'cone'`, `'concave'`, `'infinite'`, or `'myProfile'`.
  - `'dz_interpM'`: A cell array containing the distance-elevation relationships `{fitting_s_z_within_boundary, ...}` for different fan apexes.

- **Outputs:**
  - `zTopo`: Reconstructed topography of the alluvial fan.


## How to Cite

If you use this project in your research or work, please give us credit by citing it with the following BibTeX entries:

### Inverse method: extracting the distance-elevation relationship
```
Chiu, Y.H. & Capart, H. (In process). *Inverse Computational Morphology of Confined Debris and Alluvial Fans Using the Visibility Polygon*. Retrieved from [https://github.com/damiel-hub/geometric-alluvial-fan].
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


