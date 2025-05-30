# DDPG Reinforcement Learning Controller for Monopile Guy Wire Control

This directory contains a complete Deep Deterministic Policy Gradient (DDPG) reinforcement learning implementation for controlling guy wire cables to minimize tower top displacement in offshore monopile structures.

## Overview

The DDPG controller uses LSTM neural networks to learn optimal control policies for two guy wire cables based on:
- **Observations**: Tower top displacement and wave elevation
- **Actions**: Control displacements for guy wires 1 and 2 (±0.5m range)
- **Objective**: Minimize tower top displacement while minimizing control effort

## Files Description

### Core DDPG Implementation
- `trainDDPGController.m` - Main training script for the DDPG agent
- `testDDPGController.m` - Testing and performance evaluation script
- `createDDPGSimulinkModel.m` - Creates Simulink model with integrated DDPG controller

### Network Architecture
- **Actor Network**: LSTM (64 units) → FC(32) → FC(16) → FC(2) → Tanh → Scaling
- **Critic Network**: 
  - Observation path: LSTM (64 units) → FC(32)
  - Action path: FC(32)
  - Combined: Addition → FC(64) → FC(32) → FC(1)
- **Sequence Length**: 100 time steps
- **Sample Time**: 0.01 seconds

## Quick Start Guide

### 1. Training the DDPG Agent

```matlab
% Navigate to the Simulink directory
cd('/Users/jacobfontaine/Documents/Graduate Projects/PhD/monoSim/Simulink')

% Start training (this will take several hours)
trainDDPGController();
```

**Training Parameters:**
- Episodes: 1000
- Steps per episode: 2000
- Experience buffer: 1M transitions
- Mini-batch size: 128
- Learning rates: 1e-4 (critic), 1e-5 (actor)

### 2. Testing the Trained Agent

```matlab
% Test the trained agent
testDDPGController();

% This will generate:
% - Performance plots
% - Statistical analysis
% - ddpg_test_results.mat file
```

### 3. Using in Simulink

```matlab
% Create Simulink model with DDPG controller
createDDPGSimulinkModel();

% Open the created model
open_system('MonopileDDPGControl');

% Run simulation
sim('MonopileDDPGControl');
```

## Training Details

### Environment Setup
- **State Space**: 2D continuous (displacement, wave elevation)
- **Action Space**: 2D continuous (guy wire controls, ±0.5m)
- **Reward Function**: `r = -|displacement| - 0.1*(|u1| + |u2|)`
- **Episode Length**: 2000 steps (20 seconds)
- **Termination**: Excessive displacement (>5m) or max steps

### DDPG Configuration
- **Target Network Update**: Soft update (τ = 1e-3)
- **Exploration Noise**: Ornstein-Uhlenbeck process
- **Noise Variance**: 0.1 with decay rate 1e-5
- **Discount Factor**: 0.99 (default)

### LSTM Network Benefits
- **Memory**: Captures temporal dependencies in wave-structure interaction
- **Sequence Processing**: Uses 100-step history for informed decisions
- **Compact Size**: 64 hidden units for computational efficiency

## Performance Metrics

The system tracks several key performance indicators:

1. **Displacement Metrics**:
   - Maximum absolute displacement
   - RMS displacement
   - Displacement distribution

2. **Control Metrics**:
   - Control effort (RMS)
   - Control action smoothness
   - Guy wire force levels

3. **Learning Metrics**:
   - Episode rewards
   - Training convergence
   - Policy stability

## Expected Results

A well-trained DDPG agent should achieve:
- **Displacement Reduction**: 30-50% compared to uncontrolled case
- **Smooth Control**: Minimal high-frequency control actions
- **Stability**: Consistent performance across different wave conditions
- **Efficiency**: Low computational overhead for real-time operation

## Troubleshooting

### Common Issues

1. **Training Not Converging**:
   - Reduce learning rates
   - Increase exploration noise
   - Check reward function scaling

2. **Agent File Not Found**:
   - Ensure `trainedDDPGAgent.mat` exists
   - Run training script first
   - Check file path

3. **Simulink Integration Issues**:
   - Verify MATLAB Function block syntax
   - Check signal dimensions
   - Ensure proper initialization

4. **Memory Issues**:
   - Reduce experience buffer size
   - Decrease mini-batch size
   - Use smaller network architectures

### Performance Tuning

1. **Improve Learning Speed**:
   - Increase learning rates (carefully)
   - Use prioritized experience replay
   - Implement curriculum learning

2. **Better Control Performance**:
   - Tune reward function weights
   - Adjust action limits
   - Modify network architecture

3. **Reduce Computational Load**:
   - Decrease sequence length
   - Use smaller LSTM units
   - Implement model compression

## Advanced Usage

### Custom Wave Conditions

Modify the wave generation in `testDDPGController.m`:

```matlab
function wave_elev = generateWaveElevation(t)
    % Your custom wave model here
    H_s = 7.0; % Increase significant wave height
    T_p = 12.0; % Change peak period
    % ... additional wave components
end
```

### Reward Function Modification

Adjust the reward function in `MonopileStepFcn`:

```matlab
% Example: Add velocity penalty
reward = -abs(tip_displacement) - 0.05*abs(tip_velocity) - 0.1*(abs(deltaL_1) + abs(deltaL_2));
```

### Network Architecture Changes

Modify `createActorLSTMNetwork` or `createCriticLSTMNetwork` for different architectures:

```matlab
% Example: Add dropout for regularization
lstmLayer(hiddenSize, 'Name', 'lstm1', 'OutputMode', 'last')
dropoutLayer(0.2, 'Name', 'dropout1')
```

## Integration with Existing Code

The DDPG system integrates seamlessly with the existing monopile simulation:
- Uses `monopileSimulinkFunction.m` as the environment
- Maintains compatibility with all force models
- Preserves full dynamics and physics

## Future Enhancements

1. **Multi-Agent Systems**: Control multiple monopiles simultaneously
2. **Adaptive Control**: Online learning and adaptation
3. **Robust Control**: Uncertainty quantification and robust policies
4. **Hardware Integration**: Real-time deployment on embedded systems

## References

- Lillicrap, T. P., et al. "Continuous control with deep reinforcement learning." ICLR 2016.
- Hochreiter, S., & Schmidhuber, J. "Long short-term memory." Neural computation 9.8 (1997): 1735-1780.
- MATLAB Reinforcement Learning Toolbox Documentation

---

**Note**: Training times can vary significantly based on hardware. GPU acceleration is recommended for faster training. The system has been tested with MATLAB R2023a and Reinforcement Learning Toolbox.