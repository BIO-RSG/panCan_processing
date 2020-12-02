library(openxlsx)
library(ncdf4)
library(raster)
library(dplyr)
library(oceancolouR)
library(stringr)


#===============================================================================
# DESCRIPTION
#===============================================================================
#
# Add a layer to a NASA L3b PanCan image: POLY4 chl-a, using coefficients
# computed in Clay et al 2019
#
# Here, POLY4 requires:
#       Rrs from a L3b image
#       optimal coefficients
#       sensor name (MODIS, SeaWiFS, or VIIRS)
#       vector of standard wavelengths for that sensor (options below)
#
#
#===============================================================================
# VARIABLES
#===============================================================================

# These are both case-sensitive: use only the options listed
sensor <- "MODIS" # MODIS, SeaWiFS, or VIIRS-SNPP
region <- "NWA" # NWA or NEP

years <- 2003:2020

days <- 1:366

path <- paste0("/mnt/data3/claysa/", sensor, "/")

# Filename containing optimal POLY4 chla coefficients.
coef_file <- "03a_poly_coefs_2019.xlsx"


#*******************************************************************************

# Get coefficients, sensor and region-dependent
coefs <- read.xlsx(coef_file, sheet=paste0(region, "_", sensor))
coefs <- coefs$POLY4

# Get wavelengths, sensor-dependent
all_wvs <- list("MODIS"=list("blue"=c(443,488), "green"=547),
                "SeaWiFS"=list("blue"=c(443,490,510), "green"=555),
                "VIIRS-SNPP"=list("blue"=c(443,486), "green"=551))
wvs <- all_wvs[names(all_wvs)==sensor][[sensor]]


#********************************
# Get indices of lat/lon so the data can be restricted by region (NWA or NEP)

if (region=="NWA") {
    lon_range <- c(-95, -42)
    lat_range <- c(39, 82)
} else if (region=="NEP") {
    lon_range <- c(-140, -122)
    lat_range <- c(46, 60)
}

# get lats/lons for panCanadian grid at 4km resolution
data("pancan_lats_4km")
data("pancan_lons_4km")

# get lat/lon indices for the selected region (NWA or NEP)
ssok <- abs(pancan_lons_4km) >= abs(lon_range[2]) & abs(pancan_lons_4km) <= abs(lon_range[1]) & pancan_lats_4km >= lat_range[1] & pancan_lats_4km <= lat_range[2]

# To catch files that give the following error when writing the new data to netCDF:
#   Error in Rsx_nc4_put_vara_double: NetCDF: Numeric conversion not representable
bad_files <- data.frame(matrix(nrow=0, ncol=2))
colnames(bad_files) <- c("year", "day")


#********************************
# LOOP THROUGH YEARS

for (i in 1:length(years)) {
    
    year <- years[i]
    
    in_path_year <- paste0(path, "RRS/PANCAN/", year, "/")
    out_path_year <- paste0(path, "CHL_POLY4/", region, "/", year, "/")
    
    # create the directory, if necessary
    get_dir(out_path_year)
    
    L3b_files_year <- data.frame(files=list.files(in_path_year), stringsAsFactors = FALSE)
    
    if (nrow(L3b_files_year)==0) {next}
    
    
    #********************************
    # LOOP THROUGH DAYS
    
    for (j in 1:length(days)) {
        
        day <- days[j]
        
        # Subset list of files to this day
        L3b_files_day <- L3b_files_year %>%
            dplyr::filter(grepl(paste0(year,str_pad(day,width=3,side="left",pad="0")), files))
        
        if (nrow(L3b_files_day)==0) {next}
        
        
        #********************************
        # LOOP THROUGH FILES FOR THIS DAY (if more than one exists)
        
        for (fx in 1:length(L3b_files_day)) {
            
            L3b_name <- L3b_files_day[fx,"files"]
            
            cat(paste0("Getting L3b data from ",L3b_name,"...\n"))
            
            #********************************
            # GET L3B DATA
            
            all_rrs <- paste0("Rrs_", c(wvs$blue, wvs$green))
            in_file <- paste0(in_path_year,L3b_name)
            
            L3b <- nc_open(in_file)
            rrs = ncvar_get(L3b, all_rrs[1])[ssok]
            for (i in 2:length(all_rrs)) {
                rrs <- cbind(rrs, ncvar_get(L3b, all_rrs[i])[ssok])
            }
            nc_close(L3b)
            colnames(rrs) <- all_rrs
            
            L3b_dim <- as.integer(c(nrow(rrs), 1))
            
            
            #********************************
            # CALCULATE NEW POLY4 CHLA 
            cat("Computing POLY4 chl-a...\n")
            
            poly4_chl <- ocx(rrs, paste0("Rrs_", wvs$blue), paste0("Rrs_", wvs$green), coefs)
            attributes(poly4_chl)$dim <- L3b_dim
            
            
            #********************************
            # ADD NEW POLY4 CHLA TO OUTPUT NETCDF
            cat("Adding POLY4 chlor_a layer to L3b NetCDF...\n\n")
            
            old_file <- paste0(in_path_year, L3b_name)
            new_file <- paste0(out_path_year,
                               strsplit(L3b_name, "[.]")[[1]][1],
                               ifelse(sensor=="VIIRS-SNPP", ".L3b_DAY_SNPP_CHL_POLY4_", ".L3b_DAY_CHL_POLY4_"),
                               region, ".nc")
            
            file_creation <- try({
                
                # Copy the original L3b file to a new test file.
                file.copy(from=old_file,
                          to=new_file,
                          overwrite=F)
                
                # Add new variable(s) to the NetCDF file
                # https://rdrr.io/cran/ncdf4/man/ncvar_add.html
                L3b <- nc_open(new_file, write=T)
    
                dim_time <- L3b$dim[['time']]
                dim_bindata <- L3b$dim[['binDataDim']]
                # adjust bindata dimension to the selected area
                dim_bindata$name <- "binDataSubDim"
                dim_bindata$len <- length(poly4_chl)
                dim_bindata$vals <- 1:length(poly4_chl)
                # create output dims
                outdims <- list(dim_time, dim_bindata) # keep the dimensions in this order

                # NEW POLY4 CHL
                # Define variable
                poly4_var <- ncvar_def(name="chlor_a",
                                       units="mg m^-3",
                                       dim=outdims,
                                       missval=-32767,
                                       longname=paste0("Chlorophyll-a Concentration, ",
                                                       "POLY4 model with coefficients optimized to ",
                                                       sensor, " and ", region))
                # Add variable to output NetCDF
                L3b <- ncvar_add(L3b, poly4_var)
                # Write data to the variable
                ncvar_put(L3b, poly4_var, vals=poly4_chl)
                
                nc_close(L3b)
                
            }, silent=TRUE)
            
            
            # If the file wasn't created properly, remove it and add this to the list of bad files.
            if (class(file_creation)=="try-error") {
                file.remove(new_file)
                bad_files <- rbind(bad_files, data.frame(year=year, day=day))
            }
            
        }
        
    }
    
}


if (nrow(bad_files) > 0) {
    bad_file_fname <- paste0("poly4_bad_files_", sensor, "_", region, "_", paste0(range(years), collapse="-"), paste0(range(days), collapse="-"), ".csv")
    write.csv(bad_files, bad_file_fname, row.names=FALSE)
}
