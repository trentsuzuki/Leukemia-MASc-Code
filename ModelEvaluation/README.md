# ModelEvaluation

MATLAB-based framework for model order reduction by fitting lower-order PK dynamics to higher-order clinical model responses.

**Location:** `~/UofT/Research/Leukemia/Code/ModelEvaluation`

## Quick Start

```matlab
pk_matching
```

## Key Features

- Model order reduction: 3rd order -> 1st order approximation
- Error metrics: MAE, RMSE for model fit quality
- Pharmacokinetic matching analysis
- Friberg model step response comparison
- Parameter optimization for simplified dynamics

## Evaluation Methods

| Method | Description |
|--------|-------------|
| Model Reduction | Fit 1st/2nd order models to 3rd order Jost dynamics |
| Parameter Fitting | Optimize F/V and ke for simplified PK |
| Error Metrics | MAE, RMSE for approximation quality |
| Visual Comparison | Step response overlay plots |

## Error Metrics

```
Absolute Error: |y_3rd_order - y_1st_order|
MAE: mean(|y_3rd_order - y_1st_order|)
RMSE: sqrt(mean((y_3rd_order - y_1st_order)^2))
```

## Key Scripts

| Script | Purpose |
|--------|---------|
| `pk_matching.mlx` | Fit 1st order PK to 3rd order Jost dynamics |
| `evaluate_model.m` | Parallel simulation framework |
| `plot_nlmixr_data.m` | Predicted vs measured comparison |

## Model Reduction

**Original Jost Model (3rd order):**
```
A = [-k_a, 0, 0; k_a, -k_20, 0; 0, FM3kme, -CL6tgn]
```

**Simplified Model (1st order):**
```
sys_1st_order = ss(k_e, F_over_V, 1, 0)
```

**Optimization Goal:** Minimize MAE/RMSE between step responses

## Input Data

| File | Size | Description |
|------|------|-------------|
| Treatment Data | 2.9 MB | 55,306 rows, nlmixr2 output with IPRED |
| Patient Parameters | 7.6 KB | 118 rows, BASE, KTR, SLOPE, GAMMA |

## Dependencies

- Control System Toolbox
- Statistics and Machine Learning Toolbox
- MATLAB Live Editor (.mlx files)

## Integration

- **Leukemia-PKPD:** nlmixr2 parameter estimates
- **NMPC/Backstepping:** Simplified dynamics for faster control design
