# Backstepping

Backstepping control implementations for personalized 6-MP chemotherapy dosing.

**Location:** `~/UofT/Research/Leukemia/Code/Backstepping`

## Quick Start

```matlab
cd 'Backstepping'
BacksteppingMultiStart
```

## Control Variants

| Variant | Description | States |
|---------|-------------|--------|
| Adaptive | Parameter estimation with projection | 2-state |
| State Feedback | Full state feedback | 8-state |
| Output Feedback | Observer-based (ANC only) | 8-state |
| LastFourStates | PD subsystem only | 4-state |

## Key Scripts

| Script | Purpose |
|--------|---------|
| `BacksteppingMultiStart.m` | Multi-patient simulation |
| `BacksteppingPlots.m` | Control law, time domain, phase plots |
| `BacksteppingStatistics.m` | Settling time, safety compliance |
| `adaptive_backstepping_projection.m` | Robust adaptive control |

## Control Methods

### Adaptive Backstepping (2-state)
```
Model: u -> P -> N
Parameters: theta_1=k_tr*Base^gamma, theta_2=k_tr*slope*Base^gamma, theta_3=k_tr
Gains: g_1=0.00001, g_2=0.0001, g_3=0.0001
```

### State Feedback (8-state)
```
Error: z_i = x_i - alpha_i (i=1...8)
Lyapunov: V = 1/2 * sum(C_i * z_i^2)
Gains: C_1=k_ma/2, C_2-8 = [0.01, 0.01, 0.1, 10, 10, 10, 10]
Dose saturation: [0, 250] mg
```

### Output Feedback
```
Observer gain: place(A', C', desired_poles)'
Measurement: Only ANC (x_8) required
Design: Linearization, pole placement, observability check
```

## Model Structure

**Full Jost Model (8 states):**
```
gut6MP -> blood6MP -> blood6TGN -> prol -> tr1 -> tr2 -> tr3 -> circ
```

## Performance Summary

| Metric | Value |
|--------|-------|
| Target ANC | 1.5 G/L |
| Median settling time | < 100 days |
| Safety range | 0.5-2.5 G/L |
| Patients tested | 215 |
| Initial conditions | 5 per patient |
| Generated plots | 351 |

## Results

- Multi-patient MAT files
- Figures: ControlLaw, SimulationResults, PhasePlot

## Dependencies

- Symbolic Math Toolbox (control law derivation)
- Control System Toolbox (observer design)
- Simulink (.slx models)
- Parallel Computing Toolbox

## Key Findings

- Adaptive with projection prevents parameter drift
- Output feedback clinically practical (single ANC measurement)
- Robust convergence from 5 initial conditions per patient
