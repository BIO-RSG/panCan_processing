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


library(ncdf4)
library(raster)
library(dplyr)
library(oceancolouR)
library(stringr)

#*******************************************************************************
# VARIABLES TO CHANGE

# These are case-sensitive: use only the options listed
sensor <- "SeaWiFS" # MODIS, SeaWiFS, or VIIRS-SNPP
region <- "GoSL" # NWA or NEP (for CHL_POLY4 or CHL_GSM_GS), or GoSL (for CHL_EOF)
variable <- "CHL_EOF" # CHL_POLY4, CHL_GSM_GS, or CHL_EOF

years <- 1997:2010

days <- 1:366

path <- paste0("/mnt/data3/claysa/", sensor)


#**************************
# FOR CHL_POLY4 ONLY
# Filename containing optimal POLY4 chla coefficients.
coef_file <- "03a_poly4_coefs_2019.csv"


#**************************
# FOR CHL_GSM GS ONLY
# Filename containing optimal GSM_GS CHLA exponents.
exp_file <- "03b_gsm_exponents_2019.csv"

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

# # Use single centre pixel (single-rrs) or median of 3x3 box around centre pixel (median-rrs)?
# rformat <- "single-rrs"



#*******************************************************************************
# MAIN CODE

if (variable == "CHL_POLY4") {
    
    # Get coefficients, sensor and region-dependent
    coefs <- read.csv(coef_file, stringsAsFactors = FALSE)
    coefs <- as.numeric(unlist(coefs[coefs$region==region & coefs$sensor==sensor,3:7]))
    
    # Get wavelengths, sensor-dependent
    all_wvs <- list("MODIS"=list("blue"=c(443,488), "green"=547),
                    "SeaWiFS"=list("blue"=c(443,490,510), "green"=555),
                    "VIIRS-SNPP"=list("blue"=c(443,486), "green"=551))
    wvs <- all_wvs[names(all_wvs)==sensor][[sensor]]
    
    all_rrs <- paste0("Rrs_", c(wvs$blue, wvs$green))
    
    chl_longname <- paste0("Chlorophyll-a Concentration, ",
                           "POLY4 model with coefficients optimized to ",
                           sensor, " and ", region)
    
    
} else if (variable == "CHL_EOF") {
    
    # Get wavelengths, sensor-dependent
    wvs <- all_lambda[[sensor]]
    all_rrs <- paste0("Rrs_", wvs)
    
    # Get the training set filename
    tset_fname <- list.files(pattern=paste0("03c_EOF_training_set_", region, "_", sensor))
    
    # Check if the file exists for this region and sensor
    if (length(tset_fname)==0) {
        stop(paste("No training set file for", region, sensor))
    } else if (!endsWith(tset_fname, ".csv")) {
        stop("Training set file must be in csv format")
    }
    
    eof_training_set <- read.csv(tset_fname, header=TRUE)
    
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
    
    # # Function to get median Rrs of 3x3 pixel box around each pixel.
    # median_rrs <- function(i,ind,rrs,lambda) {
    #     
    #     # Get a vector of Rrs for each wavelength, where each value is the median of
    #     # the 3x3 pixel box around a point.
    #     rrs_mat <- sapply(1:length(lambda),function(j) {rrs[[j]][ind[i,"rmin"]:ind[i,"rmax"],ind[i,"cmin"]:ind[i,"cmax"]]})
    #     rrs_mat[rrs_mat < 0] <- NA
    #     
    #     # # IF RRS(41X) < 0.0005 AND RRS(4XX)/RRS(41X) > 3, SET TO NA SO THIS IS REMOVED
    #     # rrs41X <- rrs_mat[,1]
    #     # rrs4XX <- rrs_mat[,2]
    #     # #rrs_mat[(is.na(rrs41X) | (rrs41X < 0.0005 & rrs4XX/rrs41X > 3)), 1] <- NA
    #     # rrs_mat[(is.na(rrs41X) | rrs4XX/rrs41X > 4), 1] <- NA
    #     
    #     rrs_median <- apply(rrs_mat,2,median,na.rm=T)
    #     return(rrs_median)
    #     
    # }
    
    # Compile gsm function to speed it up.
    gsm_cmp <- cmpfun(gsm)
    
    #********************************
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
# Get indices of lat/lon so the data can be restricted by region (NWA or NEP)

lon_range <- lon_bounds[[region]]
lat_range <- lat_bounds[[region]]

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
        
        for (fx in 1:length(L3b_files_day)) {
            
            L3b_name <- L3b_files_day[fx,"files"]
            
            cat(paste0("Getting L3b data from ",L3b_name,"...\n"))
            
            #*******************************************************************
            # GET L3B DATA
            
            in_file <- file.path(in_path_year,L3b_name)
            
            L3b <- nc_open(in_file)
            rrs = ncvar_get(L3b, all_rrs[1])[ssok]
            for (i in 2:length(all_rrs)) {
                rrs <- cbind(rrs, ncvar_get(L3b, all_rrs[i])[ssok])
            }
            nc_close(L3b)
            colnames(rrs) <- all_rrs
            
            L3b_dim <- as.integer(c(nrow(rrs), 1))
            
            # For GSM or EOF, remove pixels where any Rrs are NA
            if (variable == "CHL_POLY4") {
                nonNA_ind <- rep(TRUE, nrow(rrs))
            } else {
                nonNA_ind <- apply(is.finite(rrs) & rrs >= 0, MARGIN=1, FUN=sum)==ncol(rrs)
                if (sum(nonNA_ind)==0) {
                    next
                } else if (sum(nonNA_ind)==1) {
                    rrs <- t(as.matrix(rrs[nonNA_ind,]))
                } else {
                    rrs <- rrs[nonNA_ind,]
                }
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
                
                # # FORMAT RRS
                # if (rformat=="median-rrs") {
                #     
                #     cat("Formatting Rrs using the median of the 3x3 pixel box around the centre pixel...\n\n")
                #     
                #     # Get indices of pixels in the 3x3 box around the centre pixel.
                #     ind <- data.frame(ind)
                #     ind$rmin <- apply(cbind(ind$row - 1,rep(0,nrow(ind))),1,max)
                #     ind$rmax <- apply(cbind(ind$row + 1,rep(l2_dim[1],nrow(ind))),1,min)
                #     ind$cmin <- apply(cbind(ind$col - 1,rep(0,nrow(ind))),1,max)
                #     ind$cmax <- apply(cbind(ind$col + 1,rep(l2_dim[2],nrow(ind))),1,min)
                #     
                #     # Use parallel processing to get matrix of median rrs
                #     cl <- makeCluster(num_cl)
                #     clusterExport(cl, c('median_rrs','rrs','wvs'))
                #     rrs_sat <- parSapply(cl,1:nrow(ind),median_rrs,ind=ind,rrs=rrs,lambda=wvs)
                #     stopCluster(cl)
                #     
                # }
                
                # Convert values to below sea level.
                rrs <- rrs/(0.52 + 1.7*rrs)
                
                num_cl <- min(num_cl, detectCores()-1)
                
                # Initiate cluster.
                cl <- makeCluster(num_cl)
                # Load necessary variables and libraries into cluster.
                clusterExport(cl, c('gsm_cmp','gsm_model','rrs','wvs','iop3', 'chl_exp','adg_exp','bbp_exp','gtype'))
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
            
            
            #*******************************************************************
            # ADD NEW CHLA TO OUTPUT NETCDF
            cat(paste0("Adding ", variable, " layer to L3b NetCDF...\n\n"))
            
            old_file <- file.path(in_path_year, L3b_name)
            new_file <- file.path(out_path_year,
                                  paste0(strsplit(L3b_name, "[.]")[[1]][1],
                                  ifelse(sensor=="VIIRS-SNPP", paste0(".L3b_DAY_SNPP_", variable, "_"), paste0(".L3b_DAY_", variable, "_")),
                                  region, ".nc"))
            
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
