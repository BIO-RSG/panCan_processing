# Stephanie.Clay@dfo-mpo.gc.ca
# 15 Oct 2020

# On hecla, GlobColour OLCI files are downloaded to /mnt/data2/claysa/ using the
# wget command. They then need to be sorted by sensor, year, and variable using
# the code below, creating subdirectories where necessary.

years <- 2020

sensors <- c("OLCI-A", "OLCI-B")

all_variables <- list("OLCI-A"=c("CHL1", "CHL2", "CHL-OC5", "RRS"),
                      "OLCI-B"=c("CHL1"))

base_path <- "/mnt/data2/claysa"


#*******************************************************************************

# get a list of newly downloaded files
all_files <- list.files(base_path)
# remove directories
file_ind <- sapply(all_files, endsWith, suffix=".nc")
all_files <- all_files[file_ind]

setwd(base_path)

for (sensor in sensors) {
    
    # get the list of variables specific to this sensor
    variables <- all_variables[[sensor]]
    
    # make sure sensor subfolder exists
    sensor_path <- oceancolouR::get_dir(file.path(base_path, sensor))
    
    # reduce file list by sensor
    sensor_files <- all_files[grep(ifelse(sensor=="OLCI-A","OLA_","OLB_"), all_files)]
    
    for (variable in variables) {
        
        # make sure variable subfolder exists
        variable_path <- oceancolouR::get_dir(file.path(sensor_path, variable))
        
        # reduce file list by variable
        variable_files <- sensor_files[grep(variable, sensor_files)]
        
        for (year in years) {
            
            # make sure year subfolder exists
            output_dir <- oceancolouR::get_dir(file.path(variable_path, year))
            
            # reduce file list by year
            year_files <- variable_files[grep(paste0("_", year), variable_files)]
            
            # move file from base_path to year subfolder
            for (year_file in year_files) {
                file.rename(from = year_file, to = file.path(output_dir, year_file))
            }
            
        }
        
    }
    
    
    # FOR CHL1, SEPARATE GSM AND AV
    for (year in years) {
        tmp_path <- file.path(sensor_path, "CHL1", year)
        gsm_files <- list.files(tmp_path, pattern="GSM")
        gsm_file_ind <- sapply(gsm_files, endsWith, suffix=".nc")
        gsm_files <- gsm_files[gsm_file_ind]
        gsm_path <- oceancolouR::get_dir(file.path(tmp_path, "GSM"))
        for (gsm_file in gsm_files) {
            file.rename(from = file.path(tmp_path, gsm_file), to = file.path(gsm_path, gsm_file))
        }
    }
    
}

