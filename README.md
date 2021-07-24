SMARPs
======
We derived Space-weather MDI Active Region Patches, or SMARPs, from maps of the solar surface magnetic field taken by the Michelson Doppler Imager (MDI) aboard the ESA/NASA Solar and Heliospheric Observatory (SoHO). These data include maps that track every solar active region observed by MDI, along with keywords that describe the physical characteristics of each active region. These data are stored in a [publicly-available, web-accessible pSQL database](http://jsoc.stanford.edu/ajax/lookdata.html) at Stanford University.

We designed the SMARP data for use in concert with another data product called the Space-Weather HMI Active Region Patches (SHARPs, [Bobra et al. 2014](https://doi.org/10.1007/s11207-014-0529-3)), derived from photospheric magnetic field data taken by the Helioseismic and Magnetic Imager instrument aboard the NASA Solar Dynamics Observatory. Combined, the SMARP and SHARP databases provide a continuous, seamless set of active region data from 1996 until the present day. 

Users can access these data with a [SunPy](https://sunpy.org/) affiliated package called [drms](https://drms.readthedocs.io/en/stable/). If you use `drms` in your research, please cite [The SunPy Community et al. 2020](https://dx.doi.org/10.3847/1538-4357/ab4f7a) and [Glogowski et al. 2019](https://joss.theoj.org/papers/10.21105/joss.01614).

### Contents

This repository contains two folders. The `example_gallery` folder contains several notebooks and functions designed to help users understand the SMARP data and how to use the SMARP and SHARP data together. The `paper` folder contains notebooks that reproduce the figures and analysis in the SMARP paper (Bobra et al. 2021). To use these notebooks together with all the requisite Python packages, create a new [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) called smarp using the provided `smarp.yml` environment file like this:

```
> conda env create -f smarp.yml
```

**Example Gallery**

* The notebook `Create_SMARP_and_SHARP_maps_with_SunPy.ipynb` is a good place to get started. This notebook queries SMARP and SHARP data using the SunPy affiliated package `drms`, generates a Pandas dataframe of all the associated keyword metadata, and creates [SunPy Maps](https://docs.sunpy.org/en/stable/code_ref/map.html) of the SMARP and SHARP image data.
* The notebook `Compare_SMARP_and_SHARP_bitmaps.ipynb` compares the activity bitmaps in the SMARP and SHARP data sets. The activity bitmaps encode the spatial distribution of the active regions -- called TARP and HARP regions, respectively -- as well as additional information about magnetic field strength and, for the TARP bitmaps, photometric features observed in the continuum intensity data.
* The notebook `Extract_SMARP_and_SHARP_maps_from_full-disk_data.ipynb` shows how to extract the SMARP and SHARP partial-disk maps in Helioprojective Cartesian coordinates directly from full-disk maps.
* The notebook `Compare_SMARP_and_SHARP_coordinates.ipynb` shows how to map a coordinate in the SHARP data to the same physical feature in the SMARP data. It also explains how optical distortion in the MDI data limits the ability for high-precision alignment between the two data series.
* The notebook `Coalign_SMARP_and_SHARP_maps.ipynb` shows how to co-align SMARP and SHARP magnetic field maps.
* The functions in `calculate_smarpkeys.py` calculate spaceweather keywords from the line-of-sight magnetic field data in the SMARP series. Sample data are included in this repository under the `files` directory.

**Paper**

* The notebooks `Figure1.ipynb`, `Figure2.ipynb`, and `Figure3.ipynb` reproduce the three figures in the SMARP paper, submitted to the *Astrophysical Journal Supplement Series*.
* The notebook `Matching_TARPs_and_HARPs.ipynb` shows how to match HARP regions with TARP regions during the overlap period when both MDI and HMI took data. The SMARP and SHARP data overlap for half a year, between 1 May 2010 and 28 October 2010.

### Citation

If you use the Space-weather MDI Active Region Patch data in your research, please consider this repository. 