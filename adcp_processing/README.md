### Python code for processing RDI ADCP data.

* [ADCP.py](ADCP.py) contains general tools for creating an xarray dataset from a text file created by RDI tools.

* [process_StationM_ADCP.py](process_StationM_ADCP.py) is script for processing Station M data (including low pass filtering and tidal analysis) and saving as a NetCDF file. The tidal analysis is performed with UTIde using eight constituents.