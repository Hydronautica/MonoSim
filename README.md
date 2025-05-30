# MonoSim

A MATLAB-based simulator for offshore monopile structures subjected to environmental loads.

## Overview

MonoSim is a comprehensive simulator for analyzing the dynamic response of monopile foundations used in offshore wind turbines. The simulator integrates multiple physical domains:

- Aerodynamic loading on the wind turbine structure
- Hydrodynamic loading from waves and currents
- Soil-structure interaction at the foundation level
- Structural dynamics using finite element methods

The simulation uses the Newmark-β time integration method to solve the equations of motion and provides various visualization tools for analyzing the results.

## Project Structure

```
monoSim/
├── monopileSimulator.m       # Main simulation script
├── initializeParameters.m    # Parameter initialization
├── aero/                     # Aerodynamic modules
├── hydro/                    # Hydrodynamic modules
├── soil/                     # Soil dynamics modules
├── FEM/                      # Finite element modules
├── Integration/              # Time integration modules
└── docs/                     # Documentation
```

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/MonoSim.git
   ```

2. Add the MonoSim directory and all subdirectories to your MATLAB path:
   ```matlab
   addpath(genpath('/path/to/MonoSim'));
   ```

## Usage

To run a simulation with default parameters:

```matlab
% Navigate to the MonoSim directory
cd /path/to/MonoSim

% Run the simulator
monopileSimulator
```

To customize the simulation parameters:

1. Open `initializeParameters.m`
2. Modify the parameters as needed
3. Run `monopileSimulator.m`

## Documentation

Detailed documentation is available in the `docs` directory:

- [User Guide](docs/user-guide.md)
- [Aerodynamics](docs/aerodynamics.md)
- [Hydrodynamics](docs/hydrodynamics.md)
- [Soil Dynamics](docs/soil-dynamics.md)
- [Finite Element Method](docs/finite-element-method.md)
- [Time Integration](docs/time-integration.md)

## Requirements

- MATLAB R2019b or newer
- Signal Processing Toolbox (for spectral analysis)
- Optional: Parallel Computing Toolbox (for faster simulations)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use this code in your research, please cite it as:

```
@software{MonoSim,
  author = {Jacob Fontaine},
  title = {MonoSim: A MATLAB-based simulator for offshore monopile structures},
  url = {https://github.com/yourusername/MonoSim},
  year = {2023},
}
```