# MonoSim Documentation

This documentation provides an overview of the MonoSim project, a MATLAB-based simulator for offshore monopile structures subjected to environmental loads.

## Table of Contents

1. [Overview](#overview)
2. [Aerodynamics](#aerodynamics)
3. [Hydrodynamics](#hydrodynamics)
4. [Soil Dynamics](#soil-dynamics)
5. [Finite Element Method](#finite-element-method)

## Overview

MonoSim is a comprehensive simulator for analyzing the dynamic response of monopile foundations used in offshore wind turbines. The simulator integrates multiple physical domains:

- Aerodynamic loading on the wind turbine structure
- Hydrodynamic loading from waves and currents
- Soil-structure interaction at the foundation level
- Structural dynamics using finite element methods

The simulation uses the Newmark-β time integration method to solve the equations of motion and provides various visualization tools for analyzing the results.

### Project Structure

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

### How to Use

To run the simulation:

1. Ensure all dependencies are installed
2. Modify parameters in `initializeParameters.m` as needed
3. Run `monopileSimulator.m` in MATLAB

## Next Steps

Explore the specific modules in the documentation sections below for detailed information on the theoretical background and implementation details.