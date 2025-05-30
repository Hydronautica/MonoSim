# MonoSim User Guide

## Getting Started

### Prerequisites

- MATLAB R2019b or newer
- Signal Processing Toolbox (for spectral analysis)
- Optional: Parallel Computing Toolbox (for faster simulations)

### Installation

1. Clone or download the MonoSim repository
2. Add the MonoSim directory and all subdirectories to your MATLAB path:

```matlab
addpath(genpath('/path/to/monoSim'));
```

## Running a Simulation

### Basic Usage

To run a simulation with default parameters:

```matlab
% Navigate to the MonoSim directory
cd /path/to/monoSim

% Run the simulator
monopileSimulator
```

### Customizing Parameters

To customize the simulation parameters:

1. Open `initializeParameters.m`
2. Modify the parameters as needed
3. Run `monopileSimulator.m`

## Key Parameters

### Structural Parameters

- `params.L`: Total length of the monopile [m]
- `params.D`: Diameter of the monopile [m]
- `params.t`: Wall thickness of the monopile [m]
- `params.E`: Young's modulus [Pa]
- `params.rho`: Material density [kg/m³]

### Mesh Parameters

- `params.nElem`: Number of finite elements
- `params.embedDepth`: Depth of embedment into the soil [m]

### Time Integration Parameters

- `params.T`: Total simulation time [s]
- `params.dt`: Time step size [s]
- `params.betaN`: Newmark-β parameter
- `params.gammaN`: Newmark-γ parameter

### Environmental Parameters

#### Wave Parameters

- `params.waveHeight`: Wave height [m]
- `params.wavePeriod`: Wave period [s]
- `params.waterDepth`: Water depth [m]
- `params.irregular`: Boolean flag for irregular waves
- `params.secondOrder`: Boolean flag for second-order wave theory

#### Wind Parameters

- `params.windSpeed`: Reference wind speed [m/s]
- `params.windRefHeight`: Reference height for wind speed [m]
- `params.windExponent`: Power law exponent for wind profile

#### Soil Parameters

- `params.soilType`: Type of soil model
- `params.soilLayers`: Definition of soil layers

## Output and Visualization

### Available Plots

The simulator automatically generates several plots:

1. Tip displacement time history
2. Bending stress distribution
3. Nodal acceleration time history
4. Surface elevation (for wave loads)

### Animations

Set `params.video = true` to generate animations:

1. 2D animation of the monopile deflection
2. 3D animation showing forces and deflection

### Accessing Results

After running the simulation, the following variables are available in the workspace:

- `U`: Displacement history [nDOF × nSteps]
- `V`: Velocity history [nDOF × nSteps]
- `Aacc`: Acceleration history [nDOF × nSteps]
- `FHyd_history`: Hydrodynamic force history [nDOF × nSteps]
- `FAero_history`: Aerodynamic force history [nDOF × nSteps]
- `Fsoil_history`: Soil reaction force history [nDOF × nSteps]

## Advanced Usage

### Modifying Force Models

To implement custom force models:

1. Create a new function in the appropriate subdirectory (aero, hydro, or soil)
2. Modify the corresponding function call in `monopileSimulator.m`

### Adding New Visualization Methods

To add custom visualization:

1. Create a new plotting function in the FEM directory
2. Add a call to your function at the end of `monopileSimulator.m`

## Troubleshooting

### Common Issues

1. **Numerical Instability**
   - Try reducing the time step size (`params.dt`)
   - Adjust the Newmark parameters (`params.betaN` and `params.gammaN`)

2. **Memory Issues**
   - Reduce the simulation time or increase the time step
   - Use a coarser mesh (reduce `params.nElem`)

3. **Slow Execution**
   - Enable parallel computing if available
   - Simplify the force models or reduce the complexity of the simulation