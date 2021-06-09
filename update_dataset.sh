#! /bin/bash

# Stephanie.Clay@dfo-mpo.gc.ca
# 8 June 2021

# This runs a series of scripts to download new files, subset them to the PanCanadian grid,
# create POLY4, GSM_GS, and EOF chlorophyll files, create annual fst files, and update the
# files formatted for use in PhytoFit.

# Data downloaded: Daily L3b (level-3 binned)
# Sensors: MODIS, VIIRS-SNPP
# Variables: CHL (OCx, POLY4, GSM_GS, EOF), PAR, RRS

# This script will be run automatically on the second day of each month to process the 
# day from the previous month.
# NOTE: SECOND DAY, NOT FIRST - give it a day to make sure NASA has updated its dataset.
# 0 3 2 * * /home/claysa/panCan_processing/update_dataset.sh >> /home/claysa/panCan_processing/update_dataset.log 2>&1

# To do:
#       - add OLCI-A and OLCI-B to regular processing
#       - figure out a system to update SST regularly (CHL/PAR/RRS are updated the next day,
#           but SST has a lag time of ~2.5-3 months)



#*******************************************************************************
# STEP 1: SET SCRIPT INPUT VARIABLES

script_dir="/home/claysa/panCan_processing/"
processing_date=$(date +'%Y-%m-%d')

year=$(date +'%Y')

# day ranges to download chl, par, rrs
doy1=$(date -d "-1 month -1 day" +'%j')
doy2=$(date -d "-2 days" +'%j')

sensors="MODIS VIIRS-SNPP"
vars_to_download="CHL PAR RRS"

vars_to_process1="CHL_OCX PAR RRS"
regions1="PANCAN"
vars_to_process2="CHL_POLY4 CHL_GSM_GS"
regions2="NWA NEP"
vars_to_process3="CHL_EOF"
regions3="GoSL"

phytofit_vars="CHL_OCX CHL_POLY4 CHL_GSM_GS CHL_EOF"
phytofit_regions="atlantic pacific"



echo
echo "***************************************************************************"
echo "*************************"${processing_date}" DATASET UPDATE*************************"
echo "***************************************************************************"
echo



# STEP 2: DOWNLOAD FILES
echo
echo "=================DOWNLOADING FILES================="
echo

for sensor in $sensors; do
    for varname in $vars_to_download; do
        ${script_dir}01a_download_daily_L3b.sh $year $doy1 $doy2 $sensor $varname
    done;
done;



# STEP 3: SUBSET TO PANCANADIAN GRID
echo
echo "=================SUBSETTING TO PANCAN GRID================="
echo

for sensor in $sensors; do
    for varname in $vars_to_download; do
        ${script_dir}02a_subset_L3b_to_panCan_grid.py $year $doy1 $doy2 $sensor $varname
    done;
done;



# STEP 4: CREATE POLY4, GSM_GS AND EOF CHLA FILES
echo
echo "=================CREATING NEW CHL FILES================="
echo

for sensor in $sensors; do
    for varname in $vars_to_process2; do
        for region in $regions2; do
            Rscript ${script_dir}03_create_new_chl_files.R $year $doy1 $doy2 $sensor $varname $region
        done;
    done;
    Rscript ${script_dir}03_create_new_chl_files.R $year $doy1 $doy2 $sensor $vars_to_process3 $regions3
done;



# STEP 5: CREATE ANNUAL FST FILES
echo
echo "=================MERGING DAILY FILES TO ANNUAL FST================="
echo

for_phytofit=FALSE
for sensor in $sensors; do
    for varname in $vars_to_process1; do
        for region in $regions1; do
            Rscript ${script_dir}04_format_annual_fst.R $for_phytofit $year $sensor $varname $region
        done;
    done;
    for varname in $vars_to_process2; do
        for region in $regions2; do
            Rscript ${script_dir}04_format_annual_fst.R $for_phytofit $year $sensor $varname $region
        done;
    done;
    for varname in $vars_to_process3; do
        for region in $regions3; do
            Rscript ${script_dir}04_format_annual_fst.R $for_phytofit $year $sensor $varname $region
        done;
    done;
done;



# STEP 6: UPDATE PHYTOFIT FILES
echo
echo "=================UPDATING PHYTOFIT FILES================="
echo

for_phytofit=TRUE
for sensor in $sensors; do
    for varname in $phytofit_vars; do
        for region in $phytofit_regions; do
            Rscript ${script_dir}04_format_annual_fst.R $for_phytofit $year $sensor $varname $region
        done;
    done;
done;






# MORE STEPS:
# backup to cloud
# backup to git repo





