# Friberg Control

This folder contains scripts to validate control designs from [Model-Based Control Strategies for Personalized Leukemia Treatment](https://utoronto.scholaris.ca/items/5a67ad1e-bc95-4891-a07e-9771bdd3fe97).

### Important Scripts

- **friberg_model.m**: Dynamics for the Friberg model, as described in the thesis.
- **stability_regions.m**: Generates figure showing identified stability regions for the Friberg model.
- **state_fb_allPatients.m**: State feedback simulations and plots for patients from Jost 2020.
- **state_fb_7days_allPatients.m**: State feedback simulations + plots, with doses changing *only every 7 days*.
- **output_fb_allPatients.m**: Output feedback simulations and plots for patients from Jost 2020.
- **output_fb_7days_allPatients.m**: Output feedback simulations + plots, with doses changing *only every 7 days*.
