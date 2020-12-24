# Stephanie.Clay@dfo-mpo.gc.ca
# 07 Sep 2020

# Merge daily images into a single grid where each row is a pixel and each column
# is a day of year, then flatten the grid and write to fst for smaller file size
# and faster loading. These files are formatted for easier use in certain scripts
# (a year of data can be loaded with only one command, rather than opening 365 netCDF files).

#*******************************************************************************

# clear the environment
rm(list=ls())
gc()

library(ncdf4)
library(fst)
library(oceancolouR)

sensors <- c("MODIS", "SeaWiFS", "VIIRS-SNPP", "OLCI-A", "OLCI-B")
variables <- c("CHL_OCX", "CHL_POLY4", "CHL_GSM_GS", "PAR", "RRS", "CHL1", "CHL2", "CHL-OC5", "SST")
regions <- c("PANCAN", "NWA", "NEP")

all_years <- list("MODIS"=2003:2020,
                  "SeaWiFS"=1997:2010,
                  "VIIRS-SNPP"=2012:2020,
                  "OLCI-A"=2016:2020,
                  "OLCI-B"=2019:2020)

base_input_path <- "/mnt/data3/claysa/"


#*******************************************************************************

for (region in regions) {
    
    print(region)
    
    tmp_num_pix <- num_pix[[region]]
    
    for (sensor in sensors) {
        
        print(sensor)
        
        waves <- all_lambda[[sensor]]
        years <- all_years[[sensor]]
        
        for (variable in variables) {
            
            print(variable)
            
            if (variable=="RRS") {
                input_variable <- paste0("Rrs_", waves)
                output_variable <- paste0("Rrs_", waves)
            } else if (variable=="PAR") {
                input_variable <- "par"
                output_variable <- "PAR"
            } else if (variable=="CHL_GSM_GS") {
                input_variable <- "chl_GSM_GS"
                output_variable <- variable
            } else if (startsWith(variable, "CHL")) {
                input_variable <- "chlor_a"
                output_variable <- variable
            } else if (variable=="SST") {
                input_variable <- "sst"
                output_variable <- variable
            }
            
            input_dir <- paste0(base_input_path, sensor,"/",variable,"/",region,"/")
            output_dir <- get_dir(paste0(input_dir, "annual_fst"))
            
            
            for (k in 1:length(years)) {
            
                # get list of data available for this year
                year <- years[k]
                
                print(year)
                
                # get a list of files for this year
                files <- list.files(path = paste0(input_dir, year, "/"), pattern = '*.nc')
            
                if (length(files)==0) {next}
                
                # get the day of each existing netcdf filename
                if (sensor %in% c("OLCI-A", "OLCI-B")) {
                    file_days <- as.numeric(format(as.Date(sapply(files, substr, start=5, stop=12), format="%Y%m%d"), "%j"))
                } else {
                    file_days <- as.numeric(sapply(files, substr, start=6, stop=8))
                }
                
                for (j in 1:length(input_variable)) {
                    
                    print(output_variable[j])
                    
                    output_fname <- paste0(output_dir, "/", region, "_", sensor, "_", output_variable[j] , "_", year, ".fst")
                    
                    # start at beginning of year
                    day_list <- 1:(max(file_days))
                    var_mat <- matrix(data = NA, nrow = tmp_num_pix, ncol = length(day_list))
                    
                    for (i in day_list) { # loop through each day
                        
                        # if the current day has an existing netcdf file associated with it,
                        # add the data to the grid, otherwise add a row of NAs
                        if (i %in% file_days) {
                            d <- nc_open(filename = paste0(input_dir, year, "/", files[which(file_days==i)]))
                            tmp_var_mat <- ncvar_get(nc = d, varid = input_variable[j])
                            var_mat[,i] <- tmp_var_mat
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
                    tmp <- data.frame(var=as.numeric(var_mat), stringsAsFactors = FALSE)
                    write_fst(tmp, path=output_fname, compress=100)
                    
                    # remove large data files before processing the next year in the loop
                    rm("var_mat", "tmp")
                    gc()
                    
                }
                
            }
            
        }
        
    }
    
}
