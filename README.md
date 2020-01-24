Analysis of current meter data collected by the Monterey Bay Aquarium Research Institute (MBARI) at Station M, located on the abyssal plain 220 km offshore of central California. Current meter data were collected near the seabed at ~4000 m depth from October 2014–October 2018. Data, methods and results are described in the publication:

Connolly, T. P., P. R. McGill, R. G. Henthorn, D. A. Burrier, C. Michaud, Near-bottom currents at Station M in the abyssal Northeast Pacific, *submitted to Deep Sea Research II*

### Data

#### Current meter data

Station M current meter data are publicly available on Zenodo:
http://dx.doi.org/10.5281/zenodo.3612575

#### Wind data

Winds from the NCEP North American Regional Reanalysis (NARR) product are used in [rover_compare_wind.ipynb](rover_compare_wind.ipynb), available at:
https://www.esrl.noaa.gov/psd/thredds/catalog/Datasets/NARR/monolevel/catalog.html

The analysis uses the 10 m wind velocity in the `uwnd.10.*.nc` and `vwnd.10.*.nc` files from 2014-2018. These can be stored locally, or the code can be modified to use the ESRL OpenDAP server.

### Using the data in the analysis

Download all current meter data files to one common base directory and extract the zip files as three separate sub-directories. Modify the `data_dir` variable in [datapath.py](datapath.py) to point to the location of this base directory. If using the NCEP NARR winds, modify the `ncep_dir` variable to point to the location of the locally-stored .nc files (they should all be in one directory).

### Required packages

The Python notebooks and modules in this repository make use of the following packages:

* xarray (http://xarray.pydata.org)

* Python tools for wavelet analysis (https://github.com/regeirk/pycwt)

* GSW toolbox (https://github.com/TEOS-10/GSW-Python)

* UTide (https://github.com/wesleybowman/UTide)

* LMFIT (https://lmfit.github.io/lmfit-py/)

* physoce (https://github.com/physoce/physoce-py)

The conda environment can be recreated with the `environment.yml` file, following the instructions at https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

### Analysis

#### Analysis of benthic rover data

[rover_analysis.ipynb](rover_analysis.ipynb)
Figure 1 - Example time series
Figure 2 - Running average time series
Figure 3 - Seasonal time series
Figure 4 - Benthic rover spectrum

[rover_compare_wind.ipynb](rover_compare_wind.ipynb)
Figure 5 - Wind forcing and wavelets
* First run [rover_analysis.ipynb](rover_analysis.ipynb) to create `rover_processed_2015_2018.csv`
* Requires NCEP NARR output (see above)

#### Analysis of ADCP data

[ADCP_spectrum_complexpca.ipynb](ADCP_spectrum_complexpca.ipynb)
Figure 6 - ADCP rotary spectra
Figure 10 - complex PCA

#### Analysis of World Ocean Atlas data

[stratification_woa.ipynb](stratification_woa.ipynb)
Figure 7 - Stratification (World Ocean Atlas)

#### Bottom boundary layer - tidal analysis and idealized modeling

[ADCP_utide_analysis.ipynb](ADCP_utide_analysis.ipynb)
Run UTide analysis on ADCP data and save results in a NetCDF file

[ADCP_utide_results_theory_M2.ipynb](ADCP_utide_results_theory_M2.ipynb)
Figure 8 - Plot results of harmonic analysis and modeling for M2 constituent
* First run [ADCP_utide_analysis.ipynb](ADCP_utide_analysis.ipynb)

[ADCP_utide_results_theory_K1.ipynb](ADCP_utide_results_theory_K1.ipynb)
Figure 9 - Plot results of harmonic analysis and modeling for K1 constituent
* First run [ADCP_utide_analysis.ipynb](ADCP_utide_analysis.ipynb)

[ADCP_utide_testcases_soulsby.ipynb](ADCP_utide_testcases_soulsby.ipynb)
Shows results of analytical model, which can be compared with examples in Soulsby (1983)

#### Bottom boundary layer - time dependent numerical model

[time_dependent_model_loggrid_multiple_zo.ipynb](time_dependent_model_loggrid_multiple_zo.ipynb)
Run time-dependent numerical model for different values of $z_o$

[time_dependent_model_analysis.ipynb](time_dependent_model_analysis.ipynb)
Figure 11 - Plot RMSE for different $$ values modeled friction velocity for one month
* First run [time_dependent_model_loggrid_multiple_zo.ipynb](time_dependent_model_loggrid_multiple_zo.ipynb)
