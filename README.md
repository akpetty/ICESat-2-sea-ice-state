# New estimates of chord length and lead fraction with ICESat-2
Contact: Alek Petty / alek.a.petty@nasa.gov / www.alekpetty.com
Code contributors: Alek Petty and Marco Bagnardi

Scripts and notebooks for deriving new information regarding the polar sea ice state from ICESat-2, used in Petty et al., (2021).   

Petty, A. A., Bagnardi, M., Kurtz, N., Tilling, R., Fons, S., Armitage, T., et al. (2021). Assessment of ICESat‐2 sea ice surface classification with Sentinel‐2 imagery: implications for freeboard and new estimates of lead and floe geometry. Earth and Space Science, 8, e2020EA001491. https://doi.org/10.1029/2020EA001491'.
 
#### Versions

v1.0: This initial repository was used to compare ATL07 with S-2 and to process ICESat-2 lead fraction and chord length estimates for the Arctic and Southern Ocean as shown in Petty et al., (2021)

### Getting Started

If you are familiar with using conda to manage your Python environment you should be able to just import the conda environment included in this repository as:
```
conda env create -f environment.yml
```
then simply activate this environment to work in (to avoid using you own Python environment which will likely have dependency issues)
```
source activate pyicesat2
```
Alternatively, the environment file simply lists the Python modules needed if you want to go it alone.


### Repo layout

The ```/scripts/``` folder contains the primary scripts to process ATL07 data. The scripts contain basic descriptions. SHared functions are stored in utils.py.

The ```/data/``` folder to store the raw and gridded data output (the data are not included as they are too large for GitHub). Data for the v1.0 release are archived on Zenodo: 10.5281/zenodo.4553552

The ```/figures/``` folder contains...figures!

The ```/notebooks/``` folder will eventually be populated with some Jupyter Notebooks 


