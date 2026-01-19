# Leukemia Project MASc Code

This repository contains a collection of subfolders relating to the leukemia treatment project.

| **Module** | **Purpose** | **Key Method** |
|--------|---------|------------|
| `friberg-control` | Control validation using Friberg PK model | State/output feedback, 7-day dosing |
| `NMPC` | Nonlinear MPC dose optimization | GA + Joint UKF for state/parameter estimation |
| `ModelEvaluation` | Model order reduction: 3rd order -> 1st order PK approximation | MAE/RMSE, step response matching |
| `Backstepping` | Backstepping control variants | Adaptive, state-feedback, output-feedback |
| `Robust` | H-infinity/H2 robust control for uncertain params | LMI synthesis over parameter ranges |
| `Stability Check` | Friberg model stability analysis | Eigenvalue analysis for patient cohorts |
| `Jost 2020 Data Analysis` | Clinical data analysis (116 patients) | ANC trajectories, dose-response analysis |
| `JostModelWithParams.m` | Core 8-state ODE model | gut6MP -> blood6MP -> blood6TGN -> prol -> tr1-tr3 -> circ |
