A geospatial workflow for wildfire-risk assessment in Deliblato Sands, Serbia.

The project combines a static fire-susceptibility layer with dynamic environmental modifiers from reanalysis, weather forecast, and climate projection datasets.

## Components

- **Static susceptibility**
  - vegetation flammability
  - road influence
  - buildings
  - population density

- **Previous measurements**
  - ERA5-Land hourly reanalysis
  - ERA5 cloud cover
  - pixel-based dynamic modifier

- **Short-term forecast**
  - ECMWF near-real-time forecast
  - 15-day forecast window
  - weather variables averaged over the protected area

- **Long-term climate**
  - NASA GDDP-CMIP6
  - SSP245 and SSP585 scenarios
  - 2031–2041 monthly climate tendency modifier

## Links

- [Interactive Earth Engine app](https://cirroco.projects.earthengine.app/view/firerisk)
- [Previous measurements notebook](https://colab.research.google.com/github/minumesaros/DeliblatoSandsWildfire/blob/main/notebooks/past_fire_risk.ipynb)
- [Short-term forecast notebook](https://colab.research.google.com/github/minumesaros/DeliblatoSandsWildfire/blob/main/notebooks/short_term_forecast.ipynb)
- [Long-term climate notebook](https://colab.research.google.com/github/minumesaros/DeliblatoSandsWildfire/blob/main/notebooks/long_term_climate.ipynb)

## Main data sources

- Google Earth Engine
- ERA5-Land hourly reanalysis
- ECMWF IFS near-real-time forecast
- NASA GDDP-CMIP6 climate projections
- Static raster layers prepared for Deliblato Sands
