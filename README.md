# PanCanadian grid processing

This repository contains the scripts and instructions used to maintain the PanCanadian dataset.  
The dataset is a collection of L3b (level-3 binned) satellite images at 4km resolution, for the following sensors / years / variables:  

    MODIS (NASA)          / 2003-present  / CHL_OCX, CHL_POLY4, CHL_GSM_GS, PAR, RRS  
    SeaWiFS (NASA)        / 1997-2010     / CHL_OCX, CHL_POLY4, CHL_GSM_GS, PAR, RRS  
    VIIRS-SNPP (NASA)     / 2012-present  / CHL_OCX, CHL_POLY4, CHL_GSM_GS, PAR, RRS  
    OLCI-A (GlobColour)   / 2016-present  / CHL1, CHL2, CHL-OC5, RRS  
    OLCI-B (GlobColour)   / 2019-present  / CHL1  

where:  

    CHL_OCX       = chlorophyll-a using global NASA OCx algorithm  
    CHL_POLY4     = chlorophyll-a using regional POLY4 algorithm (Clay et al, 2019)  
    CHL_GSM_GS    = chlorophyll-a using regional GSM_GS algorithm (Clay et al, 2019)  
    PAR           = photosynthetically active radiation  
    RRS           = remote sensing reflectance  
    CHL1          = chlorophyll-a in case 1 waters  
    CHL2          = chlorophyll-a in case 2 waters (neural network algorithm)  
    CHL-OC5       = chlorophyll-a, OC5 algorithm  


*CHL_POLY4* and *CHL_GSM_GS* are regional algorithms, so they each have two sets of coefficients that are optimized to the NWA (Northwest Atlantic) or NEP (Northeast Pacific). The NWA is very large, so for some projects (such as PhytoFit), it's further reduced to a region called "Atlantic". Regions defined below:  

| Region name | Abbreviation | Longitude bounds | Latitude bounds | Number of pixels (4km, 9km) |
| ----------- | ------------ | --------------- | ---------------- | --------------------------- |
| PanCanadian grid | PANCAN | -147, -41 | 39, 86 | 529797, 132579 |
| Northwest Atlantic | NWA | -95, -42 | 39, 82 | 295425, 73771 |
| Northeast Pacific | NEP | -140, -122 | 46, 60 | 48854, 12199 |
| Atlantic | atlantic | -71, -42 | 39, 63 | 183824, unused |
| Gulf of Saint Lawrence | GoSL | -75, -49 | 41, 53 | 68067, 16962 |


Information on default flags used in MODIS/SeaWiFS/VIIRS-SNPP L3b files can be found here:  
https://oceancolor.gsfc.nasa.gov/docs/format/  
*GlobColour OLCI files include the flags, but they are omitted from NASA MODIS/SeaWiFS/VIIRS-SNPP*  

Binning scheme explained here:  
https://oceancolor.gsfc.nasa.gov/docs/format/l3bins/  
*Note that GlobColour binned files are zero-indexed, unlike NASA binned files*  

As of 2020-11-20, the raw downloaded files are stored on hecla in `/mnt/data2/claysa`.  
After subsetting them to the panCanadian grid, they are stored in `/mnt/data3/claysa`.  

#### TO-DO LIST:  

- For CHL_GSM_GS processing, use mean or median of X closest pixels (maybe 25) for each Rrs waveband used in a single data point, instead of only using the single closest point.  



--------------------------------------------------------------------------------

## INSTRUCTIONS

Before running an R (.R), Python (.py), or Shell (.sh) script, open the script in a text editor or RStudio and adjust the variables at the top of the script.  
__IMPORTANT:__ To run Shell or Python scripts from the command line, you must change the hashbang at the top of the script to the path containing your own bash and python executables, i.e. change these lines for Shell scripts and Python scripts respectively:  

    #!/bin/bash  
    #!/home/claysa/anaconda3/bin/python  

To find the right path for your own system, try typing on the command line:  

    which bash
    which python

When updating the main PanCanadian dataset, first check the last downloaded or processed day for each sensor and variable by running **00_check_latest_file.R**.  

__Input/output file paths:__  
You can specify the base input and output paths in processing, but the subfolders in those must be in the same format:  
For example, on hecla, the base path for downloads is `/mnt/data2/claysa`, and the base path for processed files (extracted to PanCanadian grid) is `/mnt/data3/claysa`. Inside each of those folders, the files are sorted into subfolders, following this system:  

- Download subfolders: sensor / variable / year  

    - Sensor = MODIS, SeaWiFS, VIIRS-SNPP, OLCI-A, or OLCI-B  
    - Variable = CHL, PAR, RRS, or SST (for MODIS/SeaWiFS/VIIRS-SNPP), CHL1, CHL2, CHL-OC5, or RRS (for OLCI-A), and CHL1 (for OLCI-B)  
    - Year = sensor-specific  
    
- Processed subfolders: sensor / variable / region / year  

    - Sensor = MODIS, SeaWiFS, VIIRS-SNPP, OLCI-A, or OLCI-B  
    - Variable = CHL_OCX, CHL_POLY4, CHL_GSM_GS, PAR, RRS, or SST (for MODIS/VIIRS-SNPP, and SeaWiFS with the exception of SST), CHL1, CHL2, CHL-OC5, or RRS (for OLCI-A), and CHL1 (for OLCI-B)  
    - Region = PANCAN, NWA, or NEP (NWA and NEP are for CHL_POLY4 and CHL_GSM_GS only, PANCAN for everything else)  
    - Year = sensor-specific  


### STEP 1  

Download / sort / convert files to L3b at 4km-resolution using the following scripts:  

- __01a_download_daily_L3b.sh__  
- __01b_seawifs_L2_to_L3b.sh__  
- __01c_olci_file_sorting.R__  

##### MODIS/VIIRS-SNPP:  

Run **01a_download_daily_L3b.sh** to download the data for your selected sensor(s), year(s), day(s), and variable(s), and make sure the output path is correct. Note that you must have an Earthdata account to download these files (see the instructions at the top of the script). To run this script from the command line, navigate to the panCan_processing project folder on the command line and type:  

`./01a_download_daily_L3b.sh`

##### SeaWiFS:  

SeaWiFS L3b files are in 9km resolution, so you need to download the L2 (level-2) files and manually process them to L3b using the OCSSW l2bin script to get 4km resolution. Instructions:  

1. __Order the files__  
- https://oceancolor.gsfc.nasa.gov/  
- Data --> Level 1&2 Browser  
- Select SeaWiFS (MLAC), then click "Reconfigure page"  
- Select the year(s) and range of days, and fill in the PanCanadian boundaries, then click "Find swaths"  
- Click "ORDER DATA"  
- If prompted, log in with your Earthdata username and password  
- Click "DO (extract my order for me)" so that the files are limited to the PanCanadian bounds  
- Select Level 2 (OC), chlorophyll a, remote sensing reflectances, and photosynthetically available radiation  
- Click "Review order"  
- Check details, then click "Submit order"  
- When confirmation email arrives, follow the instructions and confirm by clicking the link  

2. __Download the files__  
- When files are ready to download, go to https://oceancolor.gsfc.nasa.gov/  
- Services --> Order manager  
- Download "http_manifest.txt"  
- Move "http_manifest.txt" to the location where you want to download the files (currently /mnt/data2/claysa/SeaWiFS/L2)  
- Navigate to that folder on the command line and type the following to download the files:  
`wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition --wait=0.5 --random-wait -i http_manifest.txt`  
- If the files are tarred (ending in .tar), untar them by typing this:  
`for file in *.tar; do tar -xvf $file; done`  
- If they are zipped (ending in .bz2), unzip them by typing this:  
`bunzip2 *.bz2`  
- Then sort them into their appropriate "year" subfolders.  

3. __Convert files from L2 to 4km-resolution L3b__  
- Adjust the parameters in **01c_seawifs_L2_to_L3b.sh** (select the day(s), year(s), and variable(s) to process, and make sure the input and output paths are right)  
- Run the script by navigating to the panCan_processing project folder on the command line and typing the following command:  
`./01b_seawifs_L2_to_L3b.sh`


##### GlobColour OLCI-A/B:  

1. __Order the files__  
- http://hermes.acri.fr/index.php?class=archive  
- Options to select: Global (coordinates: NORTH 86, SOUTH 38, WEST -147, EAST -41), Projection: Sinusoidal (L3b), Resolution: 4km, Date(s), Binning period: daily, Sensor type: olcia/olcib, Biochemical: CHL1/CHL2/CHL-OC5, Ocean surface optical: all RRS and NRRS  
- "Search"  
- "Order products"  
- Confirm through first email  

2. __Download the files__  
- Get download link in second email (example: `ftp://ftp.hermes.acri.fr/690462250/*`)  
- On command line, navigate to the folder where you want to download the files (currently /mnt/data2/claysa)  
- Download files with the following command (changing the example link to the actual link from the second email, and including the /* at the end so that it downloads all the files from that link): `wget --user='ftp_hermes' --password='hermes%' ftp://ftp.hermes.acri.fr/690462250/*`  

3. __Sort them into folders__  
**01b_olci_file_sorting.R**  

--------------------------------------------------------------------------------

### STEP 2  

Subset files to the PanCanadian grid.  

#### MODIS / VIIRS-SNPP / SeaWiFS

Adjust variables and run the following script:  
__02a_subset_L3b_to_panCan_grid.py__  
*(PanCanadian bin and coordinate information is contained within 02a_CAN_4km.nc)*  

#### OLCI_A/B

Adjust variables and run the following script:  
__02b_olciA_global_to_pancan.R__  
*(PanCanadian bin and coordinate information is contained within 02b_panCan_ColsRows_from_ISIN.csv)*  

--------------------------------------------------------------------------------

### STEP 3  

Create CHL files using POLY4 and GSM_GS algorithms, for MODIS/SeaWiFS/VIIRS-SNPP and regions NWA and NEP.  

- __03a_create_CHL_POLY4_files.R__  
- __03b_create_CHL_GSM_GS_files.R__  

Band ratio coefficients and GSM exponents (sensor and region-specific) are contained within the following files:  

- 03a_poly_coefs_2019.xlsx  
- 03b_gsm_exponents_2019.xlsx  

--------------------------------------------------------------------------------

### STEP 4 (optional)  

Merge daily images into a single grid where each row is a pixel and each column is a day of year, then flatten the grid and write to an *fst* file for smaller file size and faster loading. These files are formatted for easier use in certain scripts (a year of data can be loaded with only one command, rather than opening 365 netCDF files).  

__04_format_annual_fst.R__  

--------------------------------------------------------------------------------

### STEP 5 (optional)  

Make climatology images.  

__05_make_climatologies.R__  

--------------------------------------------------------------------------------

### STEP 6  

Backup images to external drive.  

__06_backup_files.sh__  

--------------------------------------------------------------------------------

## TROUBLESHOOTING

**Message:**  

*Error in Rsx_nc4_put_vara_double: NetCDF: Numeric conversion not representable*  
*Error in ncvar_put(L3b, poly4_var, vals = poly4_chl) : *  
  *C function Rsx_nc4_put_vara_double returned error*  
  
**Affected scripts:**  

- 03a_create_CHL_POLY4_files.R  
- 03b_create_CHL_GSM_GS_files.R  
  
**Multiple files could not be created**  

