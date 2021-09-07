# COVID-ITA-TASK

Run environment with `docker-compose run main`.

Install the code with `cd build`, `premake4 gmake` and finally `make`.


# Parameters

`phi` - Ratio between hospitalization and infected per age, assumed proportional to the IFR of the age group.
`tau` - Proportionality constant between IFR and phi.
`zeta` - Ratio between individuals that died without and with hospitalization, assumed constant for all age groups.
`theta` - Ratio between infected and hospitalized mortality. 