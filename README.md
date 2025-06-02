# Doubly Robust Weighted and Centered Least Squares (DR-WCLS)

This repository contains implementation and simulation code for the Doubly Robust Weighted and Centered Least Squares (DR-WCLS) method.

## Repository Structure

### Main Simulation Code
The core simulation code for the primary experiments in the paper is located in the `Main Simulation` folder:

- **`sim-omit.R`**  
  The main execution file that runs the primary simulation study. Execute this script to reproduce the main results.

- **Supporting Files**  
  All other files in this directory provide necessary functions, utilities, and data generation routines required to execute `sim-omit.R`.

### Appendix Simulations
The `Appendices` folder contains supplemental simulation code:
- **Appendix D.2**: Simulation code for large time horizons
- **Appendix D.3**: Simulation code for large time horizons with a clipped bandit treatment policy

## Usage Instructions

### Running Main Simulation
1. Navigate to the `Main Simulation` directory
2. Execute the primary script in R:
   ```R
   setwd("~/Main Simulation")
   source("sim-omit.R")