# Inverse computational morphology of confined debris and alluvial fans using the visibility polygon

The forward method, developed by [Chen and Capart (2022)](https://doi.org/10.1016/j.cageo.2022.105228), simulates alluvial fan topography using a fast algorithm based on the visibility polygon. This approach requires initial terrain data, apex location, and the elevation-distance relationship (elevation profile) as input. By utilizing these inputs, the method generates a detailed topographic model of the alluvial fan.

In contrast, this project focuses on the inverse method, which aims to extract the elevation-distance relationship from fan topography using the visibility polygon. The extracted profile serves as the input data for the forward method, enabling accurate reconstruction of fan topography.


## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
- [Reference](#reference)

## Installation

### Prerequisites

Before you begin, ensure you have the following:

- **MATLAB** (version R2021a or later recommended)
- **Required Toolboxes**: Verify that you have the following MATLAB toolboxes installed:
  - Mapping Toolbox

### Installing Additional Dependencies
You also need to download ['inpoly'](https://github.com/dengwirda/inpoly) function from GitHub and add the 'inpoly-master' folder to your MATLAB path. 



## Usage

### Basic Usage

To use this project, follow the steps below:

1. Ensure that 'functions' folder is added to your MATLAB path.
2. Navigate to the main directory of the project.
3. Run the demo script by entering the following command in the MATLAB Command Window:

```matlab
run('demo_inverse_forward_methods.m')
```

The demo script demonstrates the basic 


## Reference