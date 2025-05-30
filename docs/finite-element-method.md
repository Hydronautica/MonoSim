# Finite Element Method

## Theory

The Finite Element Method (FEM) module in MonoSim models the structural dynamics of the monopile foundation. The implementation includes:

### Beam Elements

The monopile structure is discretized using beam elements with:

1. **Degrees of Freedom**: Each node has translational and rotational degrees of freedom.

2. **Element Formulation**: Euler-Bernoulli or Timoshenko beam theory depending on the slenderness ratio.

### Dynamic Equation

The structural dynamics are governed by the equation of motion:

```
M * ü + C * u̇ + K * u = F(t)
```

Where:
- M is the mass matrix
- C is the damping matrix
- K is the stiffness matrix
- u, u̇, and ü are the displacement, velocity, and acceleration vectors
- F(t) is the external force vector

### Damping Model

Rayleigh damping is used to model structural damping:

```
C = α * M + β * K
```

Where:
- α and β are the Rayleigh damping coefficients

### Natural Frequencies

The natural frequencies of the structure are calculated by solving the eigenvalue problem:

```
(K - ω² * M) * φ = 0
```

Where:
- ω is the natural frequency
- φ is the mode shape vector

## Implementation

The FEM components are implemented in the following files:

- `assembleMeshAndMatrices.m`: Creates the mesh and assembles global matrices
- `calculateNaturalFrequencies.m`: Computes natural frequencies and mode shapes
- `plotNodalAcceleration.m`: Visualizes nodal accelerations
- `plotTipDisplacementAndStress.m`: Visualizes tip displacement and bending stress
- `generate2DAnimation.m`: Creates 2D animation of the structure's response
- `generate3DAnimation.m`: Creates 3D animation of the structure's response

## Validation

*This section will contain validation studies comparing the FEM model results with analytical solutions or other numerical models.*

### Planned Validation Studies

1. Comparison with analytical solutions for simple beam problems
2. Validation against commercial FEM software
3. Verification of natural frequencies and mode shapes
4. Benchmarking against experimental data when available