# Introduction
This python programs takes Van Allen Probes (mageis and rept) particle data, and fits PADs with a sin^n function.
Currently a static folder tree. When setting up, name location of data. (ex: /user/youruser/data)
Then folder structure is, e.g. RBSP/rept/L2/A, RBSP/mageis/L3/B, RBSP/ephemeris/A

Writes to e.g. RBSP/rept/L4PAI/A

Code will currently download ephemeris files from https://rbsp-ect.newmexicoconsortium.org
if they don't exist in the folder structure, but other files need to be downloaded.



# Install commands
sudo python3 -m pip install -e .
python3 -m sin_PADs config
sudo python3 setup.py install


import sin_PADs
import datetime
from sin_PADs import create_fits


sin_PADs.create_fits.load_data('A','2012-10-26','rept',rewrite='yes')
#takes sc_id, date, instrument, and rewrite cdf?

[![DOI](https://zenodo.org/badge/618530805.svg)](https://zenodo.org/badge/latestdoi/618530805)

