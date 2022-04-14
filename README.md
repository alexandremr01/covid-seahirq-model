# COVID-ITA-TASK

Run environment with `docker compose run main`.

Install the code with 

```bash
mkdir build
cd build
cmake ..
make
```


# Parameters

`phi` - Ratio between hospitalization and infected per age, assumed proportional to the IFR of the age group.
`tau` - Proportionality constant between IFR and phi.
`zeta` - Ratio between individuals that died without and with hospitalization, assumed constant for all age groups.
`theta` - Ratio between infected and hospitalized mortality. 

# Code Reference

`/src` contains all the C++ code, divided into 3 modules

## SEAHIRQ Lib

This is the core domain:

- `/src/seahirq_lib/scenario.cpp` - Orchestrates a model, scenario parameters and a solver 
- `/src/seahirq_lib/models.cpp` - Model differential equations
- `/src/seahirq_lib/solvers.cpp` - RungeKutta solvers for 1 and 2 dimensional arrays

## Python Interface

Files that will generate the python interface for the model.

- `/src/seahirq_lib/bind.cpp` - Declares the interface between Python and C code

The remaining functions are in majority handlers.

## Commands

Files that will generate the executable interface for the model.

- `/src/cmd/spatial_covid0d_estrat.cpp` - Main function
- `/src/seahirq_lib/real_data_reader.cpp` - Implements I/O
