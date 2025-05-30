# Soil Dynamics

## Theory

The soil dynamics module in MonoSim models the interaction between the monopile foundation and the surrounding soil. The implementation includes:

### p-y Curves

The lateral soil resistance is modeled using p-y curves, which represent the relationship between lateral soil resistance (p) and lateral displacement (y). The simulator implements:

1. **Static p-y Curves**: Based on standard design codes (e.g., API, DNV).

2. **Dynamic p-y Curves**: Account for cyclic loading effects.

3. **Hysteretic Behavior**: Models energy dissipation in the soil during cyclic loading.

### Mathematical Formulation

The p-y relationship is typically expressed as:

```
p = A * p_u * tanh(k * y / (A * p_u))
```

Where:
- p is the lateral soil resistance per unit length
- p_u is the ultimate lateral resistance
- A is an empirical factor
- k is the initial modulus of subgrade reaction
- y is the lateral displacement

### Hysteretic Models

For cyclic loading, hysteretic models are implemented to capture the nonlinear behavior of soil:

1. **Simple Hysteretic Model**: Uses rules for loading, unloading, and reloading paths.

2. **Memory-Based Model**: Tracks the loading history to determine the appropriate path.

## Implementation

The soil dynamics are implemented in the following files:

- `computeSoilForces.m`: Main function for calculating soil reaction forces
- `initializeSoilStates.m`: Sets up initial soil state variables
- `pySimpleHysteretic.m`: Implements the p-y hysteretic model for lateral resistance
- `mrSimpleHysteretic.m`: Implements the moment-rotation hysteretic model for rotational resistance

## Validation

*This section will contain validation studies comparing the soil model results with experimental data or other numerical models.*

### Planned Validation Studies

1. Comparison with centrifuge tests
2. Validation against full-scale field measurements
3. Benchmarking against finite element soil models
4. Verification against analytical solutions for simple cases
5. Comparison with industry-standard design approaches