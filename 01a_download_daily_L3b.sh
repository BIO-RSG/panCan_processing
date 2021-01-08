#!/bin/bash

#----------------------------------------------------------------------------
# DESCRIPTION:
#
# Download daily CHL/PAR/RRS L3bin files for MODIS and VIIRS-SNPP from NASA
# (https://oceandata.sci.gsfc.nasa.gov/)
#
# UPDATE JAN 2020: EARTHDATA USERNAME/PASSWORD REQUIRED TO DOWNLOAD NASA FILES
# (see instructions below)
# https://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=11107
# https://oceancolor.gsfc.nasa.gov/data/download_methods/
#
# ERRORS:
#   As of 31 Jan 2020, the -N (timestamping) option for wget has a bug, so
#   it's not used in this script.

#----------------------------------------------------------------------------
# INSTRUCTIONS:
#
#       1.  Get an EarthData account here if you don't already have one: https://urs.earthdata.nasa.gov/home
#
#       2.  On the command line, type these lines to create netrc and urs_cookies files containing your info, replacing USERNAME and PASSWD with your own username and password:
#           echo "machine urs.earthdata.nasa.gov login USERNAME password PASSWD" > ~/.netrc ; > ~/.urs_cookies
#           chmod  0600 ~/.netrc
#
#       3.  Adjust the variables below
#
#       4.  TO RUN THE SCRIPT:
#           On the command line, navigate to the folder containing the script (currently /home/claysa/pan_canadian_grid/) and type this:
#           ./download_daily_L3b.sh

#----------------------------------------------------------------------------
# VARIABLES TO CHANGE:

# Days to download
# NOTE: DAYS SHOULD NOT BE ZERO-PADDED (i.e. use 5 instead of 005)
# Possible formats:
#       days=339
#       days=$(seq 354 1 356)
#       days=$(seq 154 1 158)" 152 154 159"
#       days=$(seq 1 1 152)" "$(seq 245 1 366) --> to download days 1 to 152 and 245 to 366
days=$(seq 299 1 366)

# Years to download. See "days" variable for list formatting examples.
years=2020

# Sensor names, case-sensitive (MODIS or VIIRS-SNPP)
# Example format: sensors="MODIS VIIRS-SNPP"
sensors="MODIS VIIRS-SNPP"

# Variable names, case-sensitive (CHL, PAR, or RRS)
# Example format: varnames="CHL PAR RRS SST"
varnames="CHL PAR RRS SST"

# Files will be downloaded to this location, in a subfolder for sensor/variable/year
# (currently this is used to download files on hecla for pan-Canadian grid processing)
output_folder='/mnt/data2/claysa'


#----------------------------------------------------------------------------

# https://superuser.com/questions/231692/how-to-convert-from-day-of-year-and-year-to-a-date-yyyymmdd
jul () { date -d "$1-01-01 +$2 days -1 day" "+%Y%m%d"; }
# And here is a bash function for UNIX-based systems, such as macOS:
# jul () { (( $2 >=0 )) && local pre=+; date -v$pre$2d -v-1d -j -f "%Y-%m-%d" $1-01-01 +%Y%m%d; }

isleap() { 
    year=$1
    (( !(year % 4) && ( year % 100 || !(year % 400) ) )) && echo "true" || echo "false"
}

# If a txt list of download filenames exists, remove it to start over
rm -f L3b_to_download.txt

for sensor in $sensors; do
    
    echo "$sensor..."
    
    
    for varname in $varnames; do
        
        echo "$varname..."
        
        # Create variables for:
        #       --the filename's first letter, based on sensor
        #       --the filename end, slightly different for VIIRS-SNPP
        # sensor letter: MODIS=A, VIIRS=V
        if [ $varname = 'SST' ]; then
            filenameSuffix="L3b.DAY.${varname}"
            if [ $sensor = 'MODIS' ]; then
                filenamePrefix='AQUA_MODIS.'
            elif [ $sensor = 'VIIRS-SNPP' ]; then
                filenamePrefix='SNPP_VIIRS.'
            fi
        else
            if [ $sensor = 'MODIS' ]; then
                filenamePrefix='A'
                filenameSuffix="L3b_DAY_${varname}"
            elif [ $sensor = 'VIIRS-SNPP' ]; then
                filenamePrefix='V'
                filenameSuffix="L3b_DAY_SNPP_${varname}"
            fi
        fi
        
        
        for year in ${years}; do
            
            echo "$year..."
            
            # Make the output folder for this sensor/variable/year L3b dataset if it doesn't already exist
            # Default directory: data2 on hecla
            # Note the '-p' option checks if the directory exists before creating it
            current_path=${output_folder}/${sensor}/${varname}/${year}
            mkdir -p ${current_path}

            check_year=$(isleap $year)

            # Loop through days and create a txt file containing a list of netCDF filenames based on year, day of year, sensor, and variable
            #for d in {339..365}
            for d in ${days}; do
                
                # If it's not a leap year but you're on day 366, end the loop
                if [[ $check_year = false && $d = 366 ]]; then
                    break
                fi
                
                # Make sure "day" is 3 digits. Help source: https://stackoverflow.com/questions/3191067/in-bash-how-could-i-add-integers-with-leading-zeroes-and-maintain-a-specified-b
                day=`echo $d | awk '{printf("%03d", $1)}'`

                if [ $varname = 'SST' ]; then
                    newDate=$(jul $year $day)
                else
                    newDate=${year}${day}
                fi
                
                # Write the current filename to the txt list of filenames
                echo "${filenamePrefix}${newDate}.${filenameSuffix}.nc" >> L3b_to_download.txt

            done;
            
            # Download the list of files for this year, using cookies to pass login credentials
            # Files will be downloaded to the path given by the current_path variable
            wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition --base=https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/ --wait=0.5 --random-wait --directory-prefix=${current_path} -i L3b_to_download.txt

            # Remove the txt list of files to start over for the next year
            rm -f L3b_to_download.txt
            
        done;
        
    done;
    
done;
