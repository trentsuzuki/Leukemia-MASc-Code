# Jost 2020 Data Analysis

Analysis of patient data from Jost et al. 2020 paper on childhood ALL maintenance therapy.

**Location:** `~/UofT/Research/Leukemia/Code/Jost 2020 Data Analysis`

## Quick Start

```matlab
cd 'Jost 2020 Data Analysis'
Jost2020Analysis
treatment_outcome_figures
```

## Key Features

- MATLAB statistical analysis and visualization
- 116 individual patient outcome figures
- Includes complete supplementary materials from paper (Data Sheets 1-4)
- Treatment protocol analysis with safety monitoring

## Data Description

| Dataset | Records | Description |
|---------|---------|-------------|
| Data Sheet 1 | 116 patients | ID, Base, ktr, slope, gamma parameters |
| Data Sheet 2 | 55,306 | ANC measurements, dosing events, BSA |

## Key Scripts

| Script | Purpose |
|--------|---------|
| `Jost2020Analysis.m` | ANC statistics, BMI calculation, dose analysis |
| `treatment_outcome_figures.m` | 116 individual patient outcome plots |
| `SteveChanDataProcessing.m` | NONMEM format conversion |
| `Jost2020DataAnalysis.mlx` | Interactive exploratory analysis |

## Clinical Protocol

- **Drug:** 6-Mercaptopurine (6-MP), 50 mg/m² nominal dose
- **Protocol:** AIEOP-BFM 2009
- **Safe ANC range:** 0.5-2.0 G/L
- **Dose adjustment:** Lower dose when ANC < 0.5 G/L

## Dependencies

- MATLAB with Statistics and Machine Learning Toolbox
- Excel readtable capability

## Outputs

- Population ANC distribution plots
- Individual treatment trajectories (116 PNG files)
- Safety analysis summaries
