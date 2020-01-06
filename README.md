# **curvit**
> create lightcurves from UVIT data.

**curvit** is an open source python package to produce lightcurves from **UVIT (Ultraviolet imaging Telescope)** data.  The events list from the **official UVIT L2 pipeline (version 5.6 onwards)** is required as in input to the package.  Other pipelines are not yet supported (*but do mail me about your requirement, we can figure something out!*). 

## Installation
Linux, OS X, and Windows:
```sh
pip install curvit --user
```

## Getting started

**curvit**'s capabilities can be best demonstrated by examples. First, we need to get events list to be provided as input (an events list is a FITS file containing events). Go to ISSDC's [AstroBrowse website](https://astrobrowse.issdc.gov.in/astro_archive/archive/Home.jsp) and download UVIT data of your interest. Here, for example, the publicly available Level2 data of FO Aqr was chosen (Observation ID: `G06_084T01_9000000710`). 

This dataset, `LEVL2AS1UVT20161005G06_084T01_9000000710.zip`, is a compressed file which needs to be extracted. Once extracted, a directory named `20161005_G06_084T01_9000000710_level2` can be found. This directory has the following structure. 

```
20161005_G06_084T01_9000000710_level2/
└── uvit
    ├── RAS_NUV
    │   ├── pipeline
    │   ├── uvt_01
    │   ├── uvt_02
    │   ├── uvt_03
    │   ├── ...
    │   ├── ...
    │   └── uvt_ci
    ├── RAS_VIS
    │   ├── pipeline
    │   ├── uvt_01
    │   ├── uvt_02
    │   ├── uvt_03
    │   ├── ...
    │   ├── ...
    │   └── uvt_ci
    ├── DISCLAIMER.txt
    ├── LEVL1AS1UVT20161005G06_084T01_9000000710_05546_V2.2_dqr.xml
    └── README.txt
```
Please read the `README.txt` for details on Level2 data and what it contains. `RAS_VIS` directory contains images that were corrected for satellite drift by using the VIS (visible) channel. For the images inside `RAS_NUV` directory, NUV (near-ultraviolet) channel was used. For most cases, the data from `RAS_VIS` would be suitable. If you have downloaded a dataset that is different than the one mentioned above, check the statistics inside `DISCLAIMER.txt` to decide what to use. 

Our directory of interest,`RAS_VIS`, has the following contents. 
```
RAS_VIS/
├── pipeline
│   └── LEVL2AS1UVT20161005G06_084_L2_DM_params.txt
├── uvt_01
│   ├── F_01
│   │   ├── AS1G06_084T01_9000000710uvtFIIPC00F1A_l2err.fits
│   │   ├── AS1G06_084T01_9000000710uvtFIIPC00F1A_l2exp.fits
│   │   ├── AS1G06_084T01_9000000710uvtFIIPC00F1A_l2img.fits
│   │   ├── AS1G06_084T01_9000000710uvtFIIPC00F1I_l2img.fits
│   │   └── AS1G06_084T01_9000000710uvtFIIPC00F1_l2ce.fits
│   ├── N_01
│   │   ├── AS1G06_084T01_9000000710uvtNIIPC00F2A_l2err.fits
│   │   ├── AS1G06_084T01_9000000710uvtNIIPC00F2A_l2exp.fits
│   │   ├── AS1G06_084T01_9000000710uvtNIIPC00F2A_l2img.fits
│   │   ├── AS1G06_084T01_9000000710uvtNIIPC00F2I_l2img.fits
│   │   └── AS1G06_084T01_9000000710uvtNIIPC00F2_l2ce.fits
│   └── V_01
│       └── AS1G06_084T01_9000000710uvtVIIIM00F2_l2dr.fits
├── ...
├── ...
```

<!--- `pipeline` contains the parameter file for UVIT L2 pipeline (UL2P) used to generate the science data products inside `RAS_VIS`.  -->

Inside the directory `uvt_01`, data are organized in separate folders, each corresponding to overlapping time-ranges in UV and VIS channels, as available in Level1 dataset (`F_01`: FUV; `N_01`: NUV; `V_01`: VIS). 

The sufixes of the FITS files has the following meaning. 

* `...A_l2img.fits`: Image file in astronomical coordinates.
* `...I_l2img.fits`: Image file in instrument coordinates.
* `...A_l2exp.fits`: Exposure map for `A_l2img.fits`.
* `...A_l2err.fits`: Error map for `A_l2img.fits`.
* `...l2ce.fits` : Corrected events list.
* `...l2dr.fits` : the Relative Aspects Series (RAS) file. 

This structure of subdirectories shall repeat for all sets - `uvt_01`, `uvt_02`, `uvt_03`, etc.

For the examples given below, we will be using FUV events list (`...l2ce.fits`) from `uvt_06` as input to **curvit**.  

### `makecurves`

The `makecurves` function of **curvit** can automatically detect sources from events list and create light curves. 

```
import curvit
curvit.makecurves(events_file = 'AS1G06_084T01_9000000710uvtFIIPC00F1_l2ce.fits', how_many = 4)
```







