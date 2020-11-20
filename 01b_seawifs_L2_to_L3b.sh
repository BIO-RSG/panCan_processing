#!/bin/bash

# Stephanie.Clay@dfo-mpo.gc.ca
# 16 Sep 2020

#----------------------------------------------------------------------------
# VARIABLES TO CHANGE:

# Days to download
# NOTE: DAYS SHOULD NOT BE ZERO-PADDED (i.e. use 5 instead of 005)
# Possible formats:
#       days=339
#       days=$(seq 354 1 356)
#       days=$(seq 154 1 158)" 152 154 159"
#       days=$(seq 1 1 152)" "$(seq 245 1 366) --> to download days 1 to 152 and 245 to 366
days=$(seq 119 1 119)

# Years to download. See "days" variable for list formatting examples.
years=$(seq 2003 1 2003)

# Variable name (CHL, PAR, or RRS)
output_varnames="CHL"

base_input_folder="/mnt/data2/claysa/SeaWiFS/L2/"
base_output_folder="/mnt/data2/claysa/SeaWiFS/"

#----------------------------------------------------------------------------
. ~/seadas-7.5.3/config/seadas.env

for output_varname in $output_varnames; do
    
    # Get the name of the variable within the L2 file.
    if [ $output_varname = 'CHL' ]; then
        input_varname="chlor_a"
    elif [ $output_varname = 'PAR' ]; then
        input_varname="par"
    elif [ $output_varname = 'RRS' ]; then
        input_varname="Rrs_412,Rrs_443,Rrs_490,Rrs_510,Rrs_555,Rrs_670"
    fi
    
    for year in ${years}; do
        
        echo "$year..."
        
        input_folder="${base_input_folder}${year}/"
        output_folder="${base_output_folder}${output_varname}/${year}/"
        
        mkdir -p ${output_folder}
        
        # Loop through days and create a txt file containing a list of netCDF filenames based on year, day of year, sensor, and variable
        #for d in {339..365}
        for d in ${days}; do

            # Make sure "day" is 3 digits. Help source: https://stackoverflow.com/questions/3191067/in-bash-how-could-i-add-integers-with-leading-zeroes-and-maintain-a-specified-b
            day=`echo $d | awk '{printf("%03d", $1)}'`
            
            find ${input_folder} -maxdepth 1 -name "S${year}${day}*" > seawifs_L2_to_L3b_tmp.txt
            output_file="S${year}${day}.L3b_DAY_${output_varname}.nc"
            
            l2bin ifile="seawifs_L2_to_L3b_tmp.txt" ofile="${output_folder}${output_file}" eday="${year}${day}" sday="${year}${day}" l3bprod=${input_varname} prodtype="day" resolution="4" flaguse="ATMFAIL,LAND,HIGLINT,HILT,HISATZEN,STRAYLIGHT,CLDICE,COCCOLITH,HISOLZEN,LOWLW,CHLFAIL,NAVWARN,MAXAERITER,CHLWARN,ATMWARN,NAVFAIL"
            
            # NOBS AND NSCENES TOO?
            

        done;

    done;
    
done;
