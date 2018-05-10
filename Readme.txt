A bucket model based on Folland and Parker, 1995. But is expanded to take in
hourly resolved data to simulate the diurnal cycle of changes in bucket measurements
by Duo Chan.

Currenly, the model takes ERA-interim 1985-2014 climatology t2m, d2m, and w10m to run.
An alternative choice is to take ICOADS 1950-1990 climatology.

Insolation is from ERA-interim 1985-2014 climatology.
Cloud field is from ICOADS 1950-1990 climatology.
Initial SST is from 1982-2014 OI-SST.

All of these input enviromental fields are regridded onto 5X5 degrees.

BKT: bucket
MD: model
STP: step
DVR: driver
PREP: prepare driver
RUN: runscript

You need m_map toolbox to run this model, which can be found here:
https://www.eoas.ubc.ca/~rich/map.html

Also, smooth3CD,  conv2nan,
