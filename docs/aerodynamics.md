# Aerodynamics

## Theory

The aerodynamic module in MonoSim calculates the wind loads acting on the monopile structure. The implementation includes:

### Wind Load Models

The simulator uses several approaches to model wind loads:

1. **Steady Wind Loads**: Constant wind forces applied to the structure based on a specified wind profile.

2. **Dynamic Wind Loads**: Time-varying wind forces that can include turbulence effects.

3. **Guy Wire Forces**: Additional forces from guy wires if they are present in the structure.

### Wind Profile

The wind velocity profile typically follows a power law or logarithmic profile with height:

```
V(z) = V_ref * (z/z_ref)^α
```

Where:
- V(z) is the wind velocity at height z
- V_ref is the reference wind velocity at reference height z_ref
- α is the power law exponent (typically 0.12 for open sea)

### Aerodynamic Force Calculation

The aerodynamic force on a structural element is calculated using:

```
F = 0.5 * ρ_air * Cd * A * V^2
```

Where:
- ρ_air is the air density
- Cd is the drag coefficient
- A is the projected area of the element
- V is the relative velocity between the wind and the structure

## Implementation

The aerodynamic forces are implemented in the following files:

- `computeAerodynamicForces.m`: Calculates time-varying wind forces on the structure
- `computeConstantForces.m`: Applies constant forces (including steady wind loads)
- `computeGuyWireForces.m`: Calculates forces from guy wires if present

## Validation

*This section will contain validation studies comparing the aerodynamic model results with experimental data or other numerical models.*

### Planned Validation Studies

1. Comparison with wind tunnel tests
2. Validation against field measurements
3. Benchmarking against industry-standard tools