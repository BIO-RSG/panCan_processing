library(openxlsx)
library(ncdf4)
library(raster)
library(parallel)
library(pbapply)
library(compiler)
library(dplyr)
library(oceancolouR)
library(stringr)

#===============================================================================
# DESCRIPTION
#===============================================================================
#
# Add a layer to a NASA L2 image: GSM chl-a, using g coefficients and optimal
# exponents computed in Clay et al 2019
#
# GSM requires:
#       Rrs from an L2 image
#       optimal exponents
#       sensor name (MODIS, SeaWiFS, or VIIRS)
#       vector of standard wavelengths for that sensor (options below)
#       the following values specific to those wavelengths:
#           g coefficients
#           aphstar
#           aw and bbw (NOT temperature and salinity dependent)
#
#
#===============================================================================
# VARIABLES
#===============================================================================

# These are both case-sensitive: use only the options listed
sensor <- "VIIRS-SNPP" # MODIS, SeaWiFS, or VIIRS-SNPP
region <- "NWA" # NWA or NEP

years <- 2012:2018

days <- 237:298

path <- paste0("/mnt/data3/claysa/", sensor, "/")

# Filename containing optimal GSM_GS CHLA exponents.
exp_file <- "03b_gsm_exponents_2019.xlsx"

# g coefficients: gs or gc
#   gs: spectrally-dependent g coefficients
#   gc: constant g coefficients with g3=2 (making the model a quadratic)
gtype <- "gs"

# Number of clusters to use in parallel processing. Depends on the number of
# processors in your system (Loki has 64 or some other huge number of processors
# and 12 seems to be the optimal number to use for the best speed).
num_cl <- 10

# 1st guesses for chl, adg, bbp respectively (starting points for Gauss-Newton
# nonlinear least squares with the nls() function). Defaults: c(0.01, 0.03, 0.019)
iop3 <- c(0.01, 0.03, 0.019)

# # Use single centre pixel (single-rrs) or median of 3x3 box around centre pixel (median-rrs)?
# rformat <- "single-rrs"


#===============================================================================
# MAIN CODE
#===============================================================================

# Wavelengths (must match sensor selection)
if (sensor=="MODIS") {
    lambda <- c(412,443,469,488,531,547,555,645,667,678)
} else if (sensor=="VIIRS-SNPP") {
    lambda <- c(410,443,486,551,671)
} else if (sensor=="SeaWiFS") {
    lambda <- c(412,443,490,510,555,670)
}

sensorname <- tolower(sensor)
if (sensorname=="viirs-snpp") {sensorname <- "viirs"}

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

exps <- read.xlsx(exp_file,sheet="comparisons",rowNames=T)

if (region=="NWA") {tmpreg <- "extNA"
} else if (region=="NEP") {tmpreg <- "west-coast"}

chl_exp <- exps[paste0(tmpreg,"_",sensorname,"_",gtype),1]
adg_exp <- exps[paste0(tmpreg,"_",sensorname,"_",gtype),2]
bbp_exp <- exps[paste0(tmpreg,"_",sensorname,"_",gtype),3]


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
    out_path_year <- paste0(path, "CHL_GSM_GS/", region, "/", year, "/")
    
    # create the directory, if necessary
    get_dir(out_path_year)
    
    L3b_files_year <- data.frame(files=list.files(in_path_year), stringsAsFactors = FALSE)
    
    if (nrow(L3b_files_year)==0) {next}
    
    
    #********************************
    # LOOP THROUGH DAYS
    
    for (j in 1:length(days)) {
        
        day <- days[j]
        daystr <- str_pad(day,width=3,side="left",pad="0")
        
        # Subset list of files to this day
        L3b_files_day <- L3b_files_year %>%
            dplyr::filter(grepl(paste0(year,daystr), files))
        
        if (nrow(L3b_files_day)==0) {next}
        
        
        #********************************
        # LOOP THROUGH FILES FOR THIS DAY (if more than one exists)
        
        for (fx in 1:length(L3b_files_day)) {
            
            L3b_name <- L3b_files_day[fx,"files"]
            
            cat(paste0("Getting L3b data from ",L3b_name,"...\n\n"))
            
            
            #********************************
            # GET L3B DATA
            
            in_file <- paste0(in_path_year,L3b_name)
            L3b <- nc_open(in_file)
            rrs <- matrix(ncvar_get(L3b,paste0("Rrs_",lambda[1]))[ssok], ncol=1)
            L3b_dim <- dim(rrs)
            for (i in 2:length(lambda)) {
                rrs <- cbind(rrs, matrix(ncvar_get(L3b,paste0("Rrs_",lambda[i]))[ssok], ncol=1))
            }
            nc_close(L3b)
            
            # Remove pixels where any Rrs are NA
            nonNA_ind <- which(apply(rrs, 1, function(i) {all(is.finite(i) & i >= 0)}))
            
            
            #********************************
            # FORMAT RRS
            
            # if (rformat=="median-rrs") {
            #     
            #     cat("Formatting Rrs using the median of the 3x3 pixel box around the centre pixel...\n\n")
            #     
            #     
            #     # Get indices of pixels in the 3x3 box around the centre pixel.
            #     ind <- data.frame(ind)
            #     ind$rmin <- apply(cbind(ind$row - 1,rep(0,nrow(ind))),1,max)
            #     ind$rmax <- apply(cbind(ind$row + 1,rep(l2_dim[1],nrow(ind))),1,min)
            #     ind$cmin <- apply(cbind(ind$col - 1,rep(0,nrow(ind))),1,max)
            #     ind$cmax <- apply(cbind(ind$col + 1,rep(l2_dim[2],nrow(ind))),1,min)
            #     
            #     # Use parallel processing to get matrix of median rrs
            #     ptm <- proc.time()
            #     cl <- makeCluster(num_cl)
            #     clusterExport(cl, c('median_rrs','rrs','lambda'))
            #     rrs_sat <- parSapply(cl,1:nrow(ind),median_rrs,ind=ind,rrs=rrs,lambda=lambda)
            #     stopCluster(cl)
            #     total_time <- proc.time() - ptm
            #     print(total_time)
            #     
            # }
            
            
            # Convert values to below sea level.
            rrs <- rrs/(0.52 + 1.7*rrs)
            
            #********************************
            # GET NEW GSM CHLA USING PARALLEL PROCESSING
            cat("\nComputing GSM inherent optical properties...\n\n")
            
            # To not use parallel, comment out "makeCluster", "clusterExport" and "stopCluster"
            # lines, remove cl from parLapply, and change parSapply to lapply.
            # Initiate cluster.
            cl <- makeCluster(num_cl)
            # Load necessary variables and libraries into cluster.
            clusterExport(cl, c('gsm_cmp','gsm_model','rrs','lambda','iop3', 'chl_exp','adg_exp','bbp_exp','gtype'))
            iops <- pbapply(X=rrs[nonNA_ind,],
                            MARGIN=1,
                            FUN=gsm_cmp,
                            lambda=lambda,
                            iop3=iop3,
                            adg_exp=adg_exp,
                            bbp_exp=bbp_exp,
                            chl_exp=chl_exp,
                            gtype=gtype,
                            cl=num_cl)
            # Return processing power to other operations.
            stopCluster(cl)
            
            iops <- t(iops)
            
            # Get new GSM chl from output
            gsm_gs_chl <- rep(NA, nrow(rrs))
            gsm_gs_chl[nonNA_ind] <- as.numeric(iops[,1])
            attributes(gsm_gs_chl)$dim <- L3b_dim
            
            # # Get adg443
            # adg443 <- rep(NA, nrow(rrs))
            # adg443[nonNA_ind] <- as.numeric(as.numeric(iops[,2]))
            # attributes(adg443)$dim <- L3b_dim
            # 
            # # Get bbp443
            # bbp443 <- rep(NA, nrow(rrs))
            # bbp443[nonNA_ind] <- as.numeric(as.numeric(iops[,3]))
            # attributes(bbp443)$dim <- L3b_dim
            # 
            # # Create aph
            # aph <- (gsm_gs_chl^chl_exp)*mean(aphstar)
            # attributes(aph)$dim <- L3b_dim
            
            
            #********************************
            # ADD NEW GSM_GS CHLA TO OUTPUT NETCDF
            cat("\nAdding GSM_GS chlor_a layer to L3b NetCDF...\n\n")
            
            # Copy the original L3b file to a new test file.
            old_file <- paste0(in_path_year, L3b_name)
            new_file <- paste0(out_path_year,
                               strsplit(L3b_name, "[.]")[[1]][1],
                               ifelse(sensor=="VIIRS-SNPP", ".L3b_DAY_SNPP_CHL_GSM_GS_", ".L3b_DAY_CHL_GSM_GS_"),
                               region, ".nc")
            
            file_creation <- try({
                    
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
                dim_bindata$len <- length(gsm_gs_chl)
                dim_bindata$vals <- 1:length(gsm_gs_chl)
                # create output dims
                outdims <- list(dim_time, dim_bindata) # keep the dimensions in this order
                NA_val <- -32767
                GSMmodel <- paste0("GSM_GS model with exponents and coefficients optimized to ",
                                   sensor, " and ", region)
                
                # NEW GSM CHL
                # Define variable
                chl_gsm_var <- ncvar_def(name="chl_GSM_GS",
                                         units="mg m^-3",
                                         dim=outdims,
                                         missval=NA_val,
                                         longname=paste("Chlorophyll-a Concentration,",GSMmodel))
                # Add variable to output NetCDF
                L3b <- ncvar_add(L3b, chl_gsm_var)
                # Write data to the variable
                ncvar_put(L3b,chl_gsm_var,vals=gsm_gs_chl)
                
                
                # # ADD ADG443
                # adg_gsm_var <- ncvar_def(name="geophysical_data/newGSM_adg443",
                #                          units="m^-1",
                #                          dim=outdims,
                #                          missval=NA_val,
                #                          longname=paste("Absorption at 443nm,",GSMmodel))
                # L3b <- ncvar_add(L3b, adg_gsm_var)
                # ncvar_put(L3b,adg_gsm_var,vals=adg443)
                # 
                # # ADD BBP443
                # bbp_gsm_var <- ncvar_def(name="geophysical_data/newGSM_bbp443",
                #                          units="m^-1",
                #                          dim=outdims,
                #                          missval=NA_val,
                #                          longname=paste("Particulate backscattering at 443nm,",GSMmodel))
                # L3b <- ncvar_add(L3b, bbp_gsm_var)
                # ncvar_put(L3b,bbp_gsm_var,vals=bbp443)
                # 
                # # ADD APH
                # aph_var <- ncvar_def(name="geophysical_data/newGSM_aph",
                #                      units="m^-1",
                #                      dim=outdims,
                #                      missval=NA_val,
                #                      longname=paste("Absorption of phytoplankton using mean aph* across spectrum,",GSMmodel))
                # L3b <- ncvar_add(L3b, aph_var)
                # ncvar_put(L3b,aph_var,vals=aph)
                
                nc_close(L3b)
                
            }, silent=TRUE)
            
            # If the file wasn't created properly, remove it and add this to the list of bad files.
            if (class(file_creation)=="try-error") {
                file.remove(new_file)
                bad_files <- rbind(bad_files, data.frame(year=year, day=day))
            }
            
            gc()
            
        }
        
    }
    
}


if (nrow(bad_files) > 0) {
    bad_file_fname <- paste0("gsm_gs_bad_files_", sensor, "_", region, "_", paste0(range(years), collapse="-"), paste0(range(days), collapse="-"), ".csv")
    write.csv(bad_files, bad_file_fname, row.names=FALSE)
}
