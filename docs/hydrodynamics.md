# Hydrodynamics

## Theory

The hydrodynamic module in MonoSim calculates the wave and current loads acting on the monopile structure. The implementation includes:

### Wave Theories

The simulator supports multiple wave theories:

1. **Linear (Airy) Wave Theory**: First-order approximation for small amplitude waves.

2. **Second-Order Wave Theory**: Includes second-order effects for more accurate representation of nonlinear waves.

3. **Irregular Waves**: Representation of realistic sea states using wave spectra (e.g., JONSWAP, Pierson-Moskowitz).

### Morison Equation

The hydrodynamic forces on slender cylindrical structures are calculated using the Morison equation:

```
F = F_D + F_M = 0.5 * ρ * C_D * D * |u| * u + C_M * ρ * (π * D²/4) * du/dt
```

Where:
- F_D is the drag force component
- F_M is the inertia force component
- ρ is the water density
- C_D is the drag coefficient
- C_M is the inertia coefficient
- D is the cylinder diameter
- u is the water particle velocity
- du/dt is the water particle acceleration

### Second-Order Wave Elevation

For more accurate representation of wave kinematics, second-order wave elevation is calculated using:

```
η₂(t) = Σ Σ B⁺(ωᵢ,ωⱼ) * cos[(ωᵢ+ωⱼ)t + (φᵢ+φⱼ)] + B⁻(ωᵢ,ωⱼ) * cos[(ωᵢ-ωⱼ)t + (φᵢ-φⱼ)]
```

Where:
- B⁺ and B⁻ are the second-order transfer functions
- ωᵢ and ωⱼ are wave frequencies
- φᵢ and φⱼ are wave phases

## Implementation

The hydrodynamic forces are implemented in the following files:

- `computeMorisonForces.m`: Calculates wave and current forces using the Morison equation
- `secondOrderElevation.m`: Computes second-order wave elevation
- `computeBplus.m` and `computeBminus.m`: Calculate second-order transfer functions
- `plotSurfaceElevation.m`: Visualizes the wave surface elevation

## Validation

*This section will contain validation studies comparing the hydrodynamic model results with experimental data or other numerical models.*

### Planned Validation Studies

1. Comparison with wave tank experiments
2. Validation against field measurements
3. Benchmarking against computational fluid dynamics (CFD) simulations
4. Verification against analytical solutions for simple cases