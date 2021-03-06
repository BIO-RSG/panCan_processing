# Stephanie.Clay@dfo-mpo.gc.ca
# 18 Dec 2020

# Use 4km-resolution satellite Rrs panCanadian L3b files to create subsetted
# chlorophyll files using new chlorophyll-a algorithms.

# CHL_POLY4 requires:
#       optimal coefficients
#       vector of standard wavelengths for the selected sensor (options below)

# CHL_EOF requies:
#       vector of standard wavelengths for the selected sensor (options below)
#       training_set dataframe, created in a region (water type) / time period similar to the L3b

# CHL_GSM_GS requires:
#       optimal exponents
#       vector of standard wavelengths for the selected sensor (options below)
#       the following values specific to those wavelengths:
#           g coefficients
#           aphstar
#           aw and bbw (NOT temperature and salinity dependent)

rm(list=ls())
library(ncdf4)
library(raster)
library(dplyr)
library(oceancolouR)
library(stringr)

#*******************************************************************************
# VARIABLES TO CHANGE

# # These are case-sensitive: use only the options listed
# sensor <- "VIIRS-SNPP" # MODIS, SeaWiFS, or VIIRS-SNPP
# region <- "NEP" # NWA or NEP (for CHL_POLY4 or CHL_GSM_GS), or GoSL (for CHL_EOF)
# variable <- "CHL_GSM_GS" # CHL_POLY4, CHL_GSM_GS, or CHL_EOF
# 
# years <- 2021
# 
# days <- 97:122


all_args <- commandArgs(trailingOnly=TRUE)
years <- as.numeric(all_args[1])
days <- seq(as.numeric(all_args[2]), as.numeric(all_args[3]))
sensor <- all_args[4]
variable <- all_args[5]
region <- all_args[6]


path <- paste0("/mnt/data3/claysa/", sensor)
repo_path <- "/home/claysa/panCan_processing"

# acceptable range of calculated chl (anything outside this will be converted to NA)
# **note: this is required to avoid anomalous huge values (e.g. 1e124) that cause an error when writing to netCDF
chl_range <- c(0,100)


#**************************
# FOR CHL_GSM GS ONLY
# Filename containing optimal GSM_GS CHLA exponents.
exp_file <- file.path(repo_path, "data/gsm_exponents_2019.csv")

# g coefficients: gs or gc
#   gs: spectrally-dependent g coefficients
#   gc: constant g coefficients with g3=2 (making the model a quadratic)
gtype <- "gs"

# Number of clusters to use in parallel processing. Depends on the number of
# processors in your system. Will use min(num_cl, detectCores()-1)
num_cl <- 10

# 1st guesses for chl, adg, bbp respectively (starting points for Gauss-Newton
# nonlinear least squares with the nls() function). Defaults: c(0.01, 0.03, 0.019)
iop3 <- c(0.01, 0.03, 0.019)



#*******************************************************************************
# MAIN CODE

if (variable == "CHL_POLY4") {
    
    # Get coefficients, sensor and region-dependent
    coefs <- get_ocx_coefs(sensor=ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)),
                           region=tolower(region),
                           alg="poly4")
    
    # Get wavelengths, sensor-dependent
    all_wvs <- list("MODIS"=list("blue"=c(443,488), "green"=547),
                    "SeaWiFS"=list("blue"=c(443,490,510), "green"=555),
                    "VIIRS-SNPP"=list("blue"=c(443,486), "green"=551))
    wvs <- all_wvs[[sensor]]
    
    all_rrs <- paste0("Rrs_", c(wvs$blue, wvs$green))
    
    chl_longname <- paste0("Chlorophyll-a Concentration, ",
                           "POLY4 model with coefficients optimized to ",
                           sensor, " and ", region)
    
    
} else if (variable == "CHL_EOF") {
    
    # Get wavelengths, sensor-dependent
    wvs <- all_lambda[[sensor]]
    all_rrs <- paste0("Rrs_", wvs)
    
    # Get the training set filename
    tset_fname <- list.files(paste0(repo_path, "/data"), pattern=paste0("EOF_training_set_", region, "_", sensor))
    
    # Check if the file exists for this region and sensor
    if (length(tset_fname)==0) {
        stop(paste("No training set file for", region, sensor))
    } else if (!endsWith(tset_fname, ".csv")) {
        stop("Training set file must be in csv format")
    }
    
    eof_training_set <- read.csv(file.path(repo_path, "data", tset_fname), header=TRUE)
    
    # Make sure it's in the right format and has at least a few points
    if (!all(c(paste0("Rrs_",wvs), "chla") %in% colnames(eof_training_set))) {
        stop("Training set is missing some necessary Rrs or chla column names")
    } else if (nrow(eof_training_set) < 5) {
        stop("Training set must have at least 5 rows")
    }
    
    chl_longname <- "Chlorophyll-a Concentration, EOF model"
    
    
} else if (variable == "CHL_GSM_GS") {
    
    library(parallel)
    library(pbapply)
    library(compiler)
    
    num_cl <- min(num_cl, detectCores()-1)
    
    # Compile gsm function to speed it up.
    gsm_cmp <- cmpfun(gsm)
    
    # GET OPTIMAL EXPONENTS
    # Based on region, sensor, and choice of g coefs - constant or spectrally-dependent?
    exps <- read.csv(exp_file, stringsAsFactors = FALSE)
    exps <- as.numeric(unlist(exps[exps$region==region & exps$sensor==sensor & exps$gtype==gtype,4:6]))
    chl_exp <- exps[1]
    adg_exp <- exps[2]
    bbp_exp <- exps[3]
    
    wvs <- all_lambda[[sensor]]
    all_rrs <- paste0("Rrs_", wvs)
    
    chl_longname <- paste0("Chlorophyll-a Concentration, GSM_GS model with exponents ",
                           "and coefficients optimized to ", sensor, " and ", region)
    
}



#*******************************************************************************
# Get indices of lat/lon so the data can be restricted by region

lon_range <- lon_bounds[[region]]
lat_range <- lat_bounds[[region]]

# get lats/lons for panCanadian grid at 4km resolution
data("pancan_lats_4km")
data("pancan_lons_4km")

# get lat/lon indices for the selected region
ssok <- abs(pancan_lons_4km) >= abs(lon_range[2]) & abs(pancan_lons_4km) <= abs(lon_range[1]) & pancan_lats_4km >= lat_range[1] & pancan_lats_4km <= lat_range[2]

# To catch files that give the following error when writing the new data to netCDF:
#   Error in Rsx_nc4_put_vara_double: NetCDF: Numeric conversion not representable
bad_files <- data.frame(matrix(nrow=0, ncol=2))
colnames(bad_files) <- c("year", "day")



#********************************
# LOOP THROUGH YEARS

for (i in 1:length(years)) {
    
    year <- years[i]
    
    in_path_year <- paste0(path, "/RRS/PANCAN/", year)
    out_path_year <- file.path(path, variable, region, year)
    
    # create the output directories, if necessary
    dir.create(out_path_year, showWarnings=FALSE, recursive=TRUE)
    
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
        
        for (fx in 1:nrow(L3b_files_day)) {
            
            L3b_name <- L3b_files_day[fx,"files"]
            
            old_file <- file.path(in_path_year, L3b_name)
            new_file <- file.path(out_path_year,
                                  paste0(strsplit(L3b_name, "[.]")[[1]][1],
                                         ifelse(sensor=="VIIRS-SNPP", paste0(".L3b_DAY_SNPP_", variable, "_"), paste0(".L3b_DAY_", variable, "_")),
                                         region, ".nc"))
            
            if (file.exists(new_file)) {
                cat(paste0(new_file, " already exists, skipping to next file...\n"))
                next
            }
            
            cat(paste0("Getting L3b data from ",L3b_name,"...\n"))
            
            #*******************************************************************
            # GET L3B DATA
            
            in_file <- file.path(in_path_year,L3b_name)
            
            L3b <- nc_open(in_file)
            rrs <- list()
            for (i in 1:length(all_rrs)) {
                rrs[[i]] <- ncvar_get(L3b, all_rrs[i])[ssok]
            }
            nc_close(L3b)
            rrs <- do.call(cbind, rrs)
            colnames(rrs) <- all_rrs
            
            L3b_dim <- as.integer(c(nrow(rrs), 1))
            
            # remove rows with invalid Rrs
            nonNA_ind <- apply(is.finite(rrs), MARGIN=1, FUN=sum)==ncol(rrs)
            if (sum(nonNA_ind) > 1) {
                rrs <- rrs[nonNA_ind,]
            } else if (sum(nonNA_ind)==1) {
                rrs <- t(as.matrix(rrs[nonNA_ind,]))
            } else {
                next
            }
            
            
            #*******************************************************************
            # CALCULATE NEW CHLA 
            cat(paste0("Computing ", variable, "...\n"))
            
            # create empty vector of the appropriate length
            new_chl <- rep(NA, num_pix[[region]][["4km"]])
            
            if (variable == "CHL_POLY4") {
                
                chl_valid <- ocx(rrs, paste0("Rrs_", wvs$blue), paste0("Rrs_", wvs$green), coefs)
                
            } else if (variable == "CHL_EOF") {
                
                chl_valid <- eof_chl(rrs=as.data.frame(rrs), training_set=eof_training_set)
                
            } else if (variable == "CHL_GSM_GS") {
                
                # Convert values to below sea level.
                rrs <- rrs/(0.52 + 1.7*rrs)
                
                # Initiate cluster.
                cl <- makeCluster(num_cl)
                # Load necessary variables and libraries into cluster.
                clusterExport(cl, c('gsm_cmp','gsm_model','rrs','wvs','iop3',
                                    'chl_exp','adg_exp','bbp_exp','gtype'))
                iops <- pbapply(X=rrs,
                                MARGIN=1,
                                FUN=gsm_cmp,
                                lambda=wvs,
                                iop3=iop3,
                                adg_exp=adg_exp,
                                bbp_exp=bbp_exp,
                                chl_exp=chl_exp,
                                gtype=gtype,
                                cl=num_cl)
                # Return processing power to other operations.
                stopCluster(cl)
                
                # Get new GSM chl from output
                chl_valid <- as.numeric((t(iops))[,1])
                
            }
            
            # add calculated chlorophyll to valid indices of output vector
            new_chl[nonNA_ind] <- chl_valid
            attributes(new_chl)$dim <- L3b_dim
            
            new_chl[!is.finite(new_chl) | new_chl < chl_range[1] | new_chl > chl_range[2]] <- NA
            
            
            #*******************************************************************
            # ADD NEW CHLA TO OUTPUT NETCDF
            cat(paste0("Adding ", variable, " layer to L3b NetCDF...\n\n"))
            
            file_creation <- try({
                
                # Copy the original L3b file to a new file.
                file.copy(from=old_file, to=new_file, overwrite=FALSE)
                
                # Add new variable(s) to the NetCDF file
                # https://rdrr.io/cran/ncdf4/man/ncvar_add.html
                L3b <- nc_open(new_file, write=TRUE)
                
                dim_time <- L3b$dim[['time']]
                dim_bindata <- L3b$dim[['binDataDim']]
                # adjust bindata dimension to the selected area
                dim_bindata$name <- "binDataSubDim"
                dim_bindata$len <- length(new_chl)
                dim_bindata$vals <- 1:length(new_chl)
                # create output dims
                outdims <- list(dim_time, dim_bindata) # keep the dimensions in this order
                
                # Define new variable
                chla_var <- ncvar_def(name="chlor_a",
                                      units="mg m^-3",
                                      dim=outdims,
                                      missval=-32767,
                                      longname=chl_longname)
                
                # Add variable to output NetCDF
                L3b <- ncvar_add(L3b, chla_var)
                # Write data to the variable
                ncvar_put(L3b, chla_var, vals=new_chl)
                
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
    bad_file_fname <- paste0(variable, "_bad_files_", sensor, "_", region, "_", paste0(range(years), collapse="-"), paste0(range(days), collapse="-"), ".csv")
    write.csv(bad_files, bad_file_fname, row.names=FALSE)
}
