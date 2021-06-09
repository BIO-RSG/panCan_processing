# Stephanie.Clay@dfo-mpo.gc.ca
# 07 Sep 2020

# CURRENTLY FOR 4KM RESOLUTION ONLY

# This script can create the files that are used in the PhytoFit app (see below
# to set the variable for_phytofit to TRUE if this is the case).

# Description:
# Merge daily images into a single grid where each row is a pixel and each column
# is a day of year, then flatten the grid and write to fst for smaller file size
# and faster loading. These files are formatted for easier use in certain scripts
# (a year of data can be loaded with only one command, rather than opening 365 netCDF files).

# clear the environment
rm(list=ls())
gc()

library(ncdf4)
library(fst)
library(oceancolouR)


#*******************************************************************************
# VARIABLES TO CHANGE ####

# If this is TRUE and the output file for the selected year
# exists but is incomplete, the data for the new days will be
# appended to the end of it
# NOTE: this assumes you're downloading and adding days in order
data_append <- TRUE

base_input_path <- "/mnt/data3/claysa"
# NOTE: The full input path will be of the form: base_input_path/sensor/variable/region/year/
#       and the output path will be of the form: base_input_path/sensor/variable/region/annual_fst/

phytofit_output_path <- "/home/claysa/PhytoFit/data"
# NOTE: files will be sorted into atlantic or pacific subfolders in the output path

# # are the files being created for PhytoFit?
# for_phytofit <- TRUE
# 
# if (for_phytofit) {
#     # OPTION 1: MAKE PHYTOFIT FILES
#     #       sensors <- c("MODIS", "SeaWiFS", "VIIRS-SNPP")
#     #       variables <- c("CHL_OCX", "CHL_POLY4", "CHL_GSM_GS", "CHL_EOF")
#     #       regions <- c("atlantic", "pacific")
#     years <- 2021
#     sensors <- c("MODIS", "VIIRS-SNPP")
#     variables <- c("CHL_EOF")
#     regions <- c("atlantic")
# } else {
#     # OPTION 2: Make other files
#     #     sensors <- c("MODIS", "SeaWiFS", "VIIRS-SNPP", "OLCI-A", "OLCI-B")
#     #     variables <- c("CHL_OCX", "CHL_POLY4", "CHL_GSM_GS", "CHL_EOF", PAR", "RRS",
#     #                    "CHL1", "CHL2", "CHL-OC5", "SST")
#     #     regions <- c("PANCAN", "NWA", "NEP", "GoSL")
#     years <- 2021
#     sensors <- c("MODIS", "VIIRS-SNPP")
#     variables <- c("CHL_OCX", "CHL_POLY4", "CHL_GSM_GS", "CHL_EOF", "PAR", "RRS", "SST")
#     regions <- c("PANCAN", "NWA", "NEP", "GoSL")
# }


all_args <- commandArgs(trailingOnly=TRUE)
for_phytofit <- as.logical(all_args[1])
years <- as.numeric(all_args[2])
sensors <- all_args[3]
variables <- all_args[4]
regions <- all_args[5]



#*******************************************************************************
# PHYTOFIT-SPECIFIC STUFF ####

# Input files have different lengths (some are for the whole PanCanadian grid,
# or the NWA, or NEP, or GoSL), and you need to subset them based on the
# output region you selected.

# Load lat/lon vectors of the input files to use in subsetting
data("nwa_lats_4km", package="oceancolouR")
data("nwa_lons_4km", package="oceancolouR")
data("nep_lats_4km", package="oceancolouR")
data("nep_lons_4km", package="oceancolouR")
data("gosl_lats_4km", package="oceancolouR")
data("gosl_lons_4km", package="oceancolouR")
data("pancan_lats_4km", package="oceancolouR")
data("pancan_lons_4km", package="oceancolouR")
# # for 1km gosl - remaining code is not adapted for this yet
# nc <- nc_open("data/gsl_1km.nc")
# gosl_lats_1km <- ncvar_get(nc, "lat")
# gosl_lons_1km <- ncvar_get(nc, "lon")
# nc_close(nc)

# get valid indices of input data for the output data file, based on lats/lons
phytofit_index <- function(region, variable) {
    
    # lats/lons of the input data
    if (variable=="CHL_OCX") {
        lon <- pancan_lons_4km
        lat <- pancan_lats_4km
    } else if (variable=="CHL_EOF" & region=="atlantic") {
        lon <- gosl_lons_4km
        lat <- gosl_lats_4km
    } else {
        if (region=="atlantic") {
            lon <- nwa_lons_4km
            lat <- nwa_lats_4km
        } else if (region=="pacific") {
            lon <- nep_lons_4km
            lat <- nep_lats_4km
        }
    }
    
    # boundaries of the output data
    tmp_reg <- ifelse(region=="pacific", "NEP", ifelse(variable=="CHL_EOF", "GoSL", "atlantic"))
    lon_range <- lon_bounds[[tmp_reg]]
    lat_range <- lat_bounds[[tmp_reg]]
    
    # get indices to extract from files, based on input file lat/lon
    # vectors and the lat/lon range that you want for the output files
    ok <- lon >= lon_range[1] & lon <= lon_range[2] & lat >= lat_range[1] & lat <= lat_range[2]
    
    return(ok)
    
}

# function to get phytofit variable name for output filename
phytofit_variable <- function(v) {
    ifelse(v=="CHL_OCX", "ocx",
           ifelse(v=="CHL_POLY4", "poly4",
                  ifelse(v=="CHL_GSM_GS", "gsmgs",
                         ifelse(v=="CHL_EOF", "eof",
                                ifelse(v=="CHL1", "chl1",
                                       ifelse(v=="CHL2", "chl2",
                                              ifelse(v=="CHL-OC5", "chloc5", NA)))))))
}

# function to get phytofit sensor for output filename
phytofit_sensor <- function(s) {
    ifelse(s=="MODIS", "modis",
           ifelse(s=="VIIRS-SNPP", "viirs",
                  ifelse(s=="SeaWiFS", "seawifs",
                         ifelse(s=="OLCI-A", "olcia", NA))))
}

# column name in the output dataframe
if (for_phytofit) {
    output_df_colname <- "chl"
} else {
    output_df_colname <- "var"
}


#*******************************************************************************
# MAIN CODE ####

for (region in regions) {
    
    print(region)
    
    
    for (sensor in sensors) {
        
        print(sensor)
        
        waves <- all_lambda[[sensor]]
        
        
        for (variable in variables) {
            
            print(variable)
            
            if (variable=="CHL_EOF" & region=="pacific") {next}
            
            if (for_phytofit) {
                tmp_reg <- ifelse(variable=="CHL_OCX", "PANCAN",
                                  ifelse(region=="pacific", "NEP",
                                         ifelse(variable=="CHL_EOF", "GoSL", "NWA")))
                input_dir <- file.path(base_input_path, sensor, variable, tmp_reg)
                output_dir <- file.path(phytofit_output_path, region)
            } else {
                input_dir <- file.path(base_input_path, sensor, variable, region)
                output_dir <- file.path(input_dir, "annual_fst")
            }
            
            # check if there is data for this combination
            dir_check <- list.files(input_dir)
            dir_check <- dir_check[dir_check %in% years]
            if (length(dir_check)==0) {next}
            
            dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
            
            if (variable=="RRS") {
                input_variable <- paste0("Rrs_", waves)
                output_variable <- paste0("Rrs_", waves)
            } else if (variable=="PAR") {
                input_variable <- "par"
                output_variable <- "PAR"
            } else if (startsWith(variable, "CHL")) {
                input_variable <- "chlor_a"
                output_variable <- variable
            } else if (variable=="SST") {
                input_variable <- "sst"
                output_variable <- variable
            }
            
            
            # If this is for PhytoFit, you must subset the data to either the
            # atlantic, pacific, or gosl regions, depending on user-selected
            # region and variable
            if (!for_phytofit) {
                
                tmp_num_pix <- num_pix[[region]][["4km"]]
                ssok <- rep(TRUE, tmp_num_pix)
                
            } else {
                
                # get indices to extract from files, based on input file lat/lon
                # vectors and the lat/lon range that you want for the output files
                ssok <- phytofit_index(region, variable)
                tmp_num_pix <- sum(ssok)
                
            }
            
            
            for (k in 1:length(years)) {
            
                year <- years[k]
                
                print(year)
                
                # get a list of files for this year
                files <- list.files(path = file.path(input_dir, year), pattern = '*.nc')
                
                if (length(files)==0) {next}
                
                # get the day of each existing netcdf filename
                if (sensor %in% c("OLCI-A", "OLCI-B")) {
                    file_days <- as.numeric(format(as.Date(sapply(files, substr, start=5, stop=12), format="%Y%m%d"), "%j"))
                } else {
                    file_days <- as.numeric(sapply(files, substr, start=6, stop=8))
                }
                
                # don't include last day of a leap year
                noleap_ind <- file_days < 366
                file_days <- file_days[noleap_ind]
                files <- files[noleap_ind]
                
                
                for (j in 1:length(input_variable)) {
                    
                    print(output_variable[j])
                    
                    if (for_phytofit) {
                        output_fname <- paste0(output_dir, "/", region, "_", phytofit_sensor(sensor), "_", phytofit_variable(output_variable[j]), "_", year, ".fst")
                    } else {
                        output_fname <- paste0(output_dir, "/", region, "_", sensor, "_", output_variable[j] , "_", year, ".fst")
                    }
                    
                    # check if the output file for this year already exists, and if it does, (optionally) append the new days to it
                    if (file.exists(output_fname) & data_append) {
                        
                        # load the existing file
                        pts <- read_fst(output_fname)
                        
                        # reshape temporarily
                        var_mat <- matrix(pts[,output_df_colname], nrow=tmp_num_pix)
                        
                        # get a list of days that haven't been added yet to the output file
                        day_list <- c()
                        if (ncol(var_mat) < max(file_days)) {
                            day_list <- (ncol(var_mat)+1):(max(file_days))
                        }
                        
                        var_mat <- cbind(var_mat, matrix(nrow=nrow(var_mat), ncol=length(day_list)))
                        
                    } else {
                        
                        # start at beginning of year
                        day_list <- 1:(max(file_days))
                        var_mat <- matrix(data = NA, nrow = tmp_num_pix, ncol = length(day_list))
                        
                    }
                    
                    
                    for (i in day_list) { # loop through each day
                        
                        # if the current day has an existing netcdf file associated with it,
                        # add the data to the grid, otherwise add a row of NAs
                        if (i %in% file_days) {
                            d <- nc_open(filename = file.path(input_dir, year, files[which(file_days==i)]))
                            tmp_var_mat <- ncvar_get(nc = d, varid = input_variable[j])
                            var_mat[,i] <- tmp_var_mat[ssok]
                            nc_close(d)
                        } else {
                            var_mat[,i] <- rep(NA, tmp_num_pix)
                        }
                        
                    }
                    
                    # overwrite invalid data
                    if (variable != "SST") {
                        var_mat[var_mat < 0] <- NA
                    }
                    
                    # reshape and write to fst file
                    tmp <- data.frame(as.numeric(var_mat), stringsAsFactors = FALSE)
                    colnames(tmp) <- output_df_colname
                    write_fst(tmp, path=output_fname, compress=100)
                    
                    # remove large data files before processing the next year in the loop
                    rm("var_mat", "tmp")
                    gc()
                    
                }
                
            }
            
        }
        
    }
    
}
