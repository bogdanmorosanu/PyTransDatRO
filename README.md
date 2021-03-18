# PyTransDatRO

## Description

This is a Python library which can be used for the coordinate transformation between ETRS89 and Stereo70 coordinate reference systems (CRSs).

It implments the official algorithm published by the Romanian *National Agency for Cadastre and Land Registration*, *Department of Geodesy and Cartography*, *Geodetic Service*.
For more details about the transformation please visit the [wiki](https://github.com/bogdanmorosanu/PyTransDatRO/wiki) section.

## Requirements

Python 3.6 or later.

## Intallation

    $ pip install pytransdatro

## Usage

This example does a transformation 
```python
import pytransdatro

# Input/Output units of measure:
#   - ETRS89: radians
#   - Stereo70: meters

t = pytransdatro.TransRO()

# Stereo70 -> ETRS89
# (n, e[, h]) -> (lat, lon[, h])
lat, lon = t.st70_to_etrs89(693771.731, 310723.518)
lat, lon, h = t.st70_to_etrs89(693771.731, 310723.518, 122.714)

# ETRS89 -> Stereo70
# (lat, lon[, h]) -> (n, e[, h])
n, e = t.etrs89_to_st70(0.832795488097716, 0.3922724455407368)
n, e, h = t.etrs89_to_st70(0.832795488097716, 0.3922724455407368, 162.0874)
```
## Features

The library is optimized for transformation of coordinates inside the same grid cell, that is, some parameters used in computation are cached for multiple calls against an instance of the class *TransRO*. This boosts the performance when transforming multiple coordinates covered by the same grid cell. The 2D transformation (no elevation/height) uses a grid cell of 11 km size. In the case of a 3D transformation (elevation/height inluced) an additional grid of 0.0333333 degrees is used.

