## Repository Details
This repository contains all the code to reproduce the results in our paper, "[Validation of the Bond et al. (2010) SDSS-derived kinematic models for the Milky Way's disk and halo stars with Gaia Data Release 3 proper motion and radial velocity data](https://arxiv.org/abs/2406.03541)".

We validate the Bond et. al. (2010) kinematic models for the Milky Way's disk and halo stars with Gaia Data Release 3 data. Bond et al. constructed models for stellar velocity distributions using stellar radial velocities measured by the Sloan Digital Sky Survey (SDSS) and stellar proper motions derived from SDSS and the Palomar Observatory Sky Survey astrometric measurements. These models describe velocity distributions as functions of position in the Galaxy, with separate models for disk and halo stars that were labeled using SDSS photometric and spectroscopic metallicity measurements. We find that the Bond et al. model predictions are in good agreement with recent measurements of stellar radial velocities and proper motions by the Gaia survey. In particular, the model accurately predicts the skewed non-Gaussian distribution of rotational velocity for disk stars and its vertical gradient, as well as the dispersions for all three velocity components. Additionally, the spatial invariance of velocity ellipsoid for halo stars when expressed in spherical coordinates is also confirmed by Gaia data at galacto-centric radial distances of up to 15 kpc.

All our results from this paper can be reproduced with this notebook. The data needed to run the code can be found in https://doi.org/10.5281/zenodo.11483943

## Quickstart
1. Clone this repository: ```git clone https://github.com/sidchaini/MWKinematicsFGKM.git```
2. Download data: ```wget -P data/ -i data/download.txt```
3. Install required packages ```pip install -r requirements.txt```
4. Open and run the notebook ```[validationBond+2010_GaiaDR3.ipynb](validationBond+2010_GaiaDR3.ipynb)```

## Authors
Bruno Domínguez, Siddharth Chaini, Karlo Mrakovčić, Brandon Sallee, and Željko Ivezić



