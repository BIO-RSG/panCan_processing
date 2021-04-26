# Stephanie.Clay@dfo-mpo.gc.ca
# 18 Dec 2020

# Use 4km-resolution satellite Rrs panCanadian L2 files to create subsetted
# chlorophyll files using new chlorophyll-a algorithms.

# CHL_POLY4 requires:
#       optimal coefficients
#       vector of standard wavelengths for the selected sensor (options below)

# CHL_EOF requies:
#       vector of standard wavelengths for the selected sensor (options below)
#       training_set dataframe, created in a region (water type) / time period similar to the L2

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
library(parallel)

#*******************************************************************************
# VARIABLES TO CHANGE

# These are case-sensitive: use only the options listed
sensor <- "MODIS" # MODIS, SeaWiFS, or VIIRS-SNPP
region <- "NWA" # NWA or NEP (for CHL_POLY4 or CHL_GSM_GS), or GoSL (for CHL_EOF)
variable <- "CHL_GSM_GS" # CHL_OCX_RSG, CHL_POLY4, CHL_GSM_GS, or CHL_EOF

years <- 2020

days <- 103

path <- "L2" # do not include last forward slash

# acceptable range of calculated chl (anything outside this will be converted to NA)
chl_range <- c(0,100)

# for parallel processing, number of cores to use
num_cl <- 10


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

# Use single centre pixel (single-rrs) or median of 3x3 box around centre pixel (median-rrs)?
rformat <- "single-rrs"



#*******************************************************************************
# MAIN CODE

median_rrs <- function(i,ind,rrs,lambda) {
    
    # Get a vector of Rrs for each wavelength, where each value is the median of
    # the 3x3 pixel box around a point.
    rrs_mat <- sapply(1:length(lambda),function(j) {rrs[[j]][ind[i,"rmin"]:ind[i,"rmax"],ind[i,"cmin"]:ind[i,"cmax"]]})
    rrs_mat[rrs_mat < 0] <- NA
    
    # # IF RRS(41X) < 0.0005 AND RRS(4XX)/RRS(41X) > 3, SET TO NA SO THIS IS REMOVED
    # rrs41X <- rrs_mat[,1]
    # rrs4XX <- rrs_mat[,2]
    # #rrs_mat[(is.na(rrs41X) | (rrs41X < 0.0005 & rrs4XX/rrs41X > 3)), 1] <- NA
    # rrs_mat[(is.na(rrs41X) | rrs4XX/rrs41X > 4), 1] <- NA
    
    rrs_median <- apply(rrs_mat,2,median,na.rm=T)
    return(rrs_median)
    
}

single_rrs <- function(i,ind,rrs,lambda) {
    
    # Get a vector of Rrs for each wavelength.
    rrs_mat <- sapply(1:length(lambda),function(j) {rrs[[j]][ind[i,"row"],ind[i,"col"]]})
    rrs_mat[rrs_mat < 0] <- NA
    
    # # IF RRS(41X) < 0.0005, SET TO NA SO THIS IS REMOVED
    # rrs41X <- rrs_mat[1]
    # rrs4XX <- rrs_mat[2]
    # #if (is.na(rrs_mat[1]) | rrs_mat[1] < 0.0005) {rrs_mat[1]  <- NA}
    # if (is.na(rrs_mat[1]) | rrs4XX/rrs41X > 4) {rrs_mat[1]  <- NA}
    
    return(rrs_mat)
    
}

if (variable == "CHL_OCX_RSG") {
    
    # Get coefficients, sensor and region-dependent
    coefs <- get_ocx_coefs(sensor = ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)),
                           region = tolower(region),
                           alg = "ocx")
    
    # Get wavelengths, sensor-dependent
    wvs <- get_ocx_lambda(sensor = ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)), use_443nm = TRUE)
    hu_bands <- get_ci_bands(ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)))
    
    all_rrs <- sort(unique(c(wvs$blue, wvs$green, paste0("Rrs_", hu_bands))))
    
    chl_longname <- "Chlorophyll-a Concentration, OCI model"
    
    
} else if (variable == "CHL_POLY4") {
    
    # Get coefficients, sensor and region-dependent
    coefs <- get_ocx_coefs(sensor = ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)),
                           region = tolower(region),
                           alg = "poly4")
    
    # Get wavelengths, sensor-dependent
    wvs <- get_ocx_lambda(sensor = ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)), use_443nm = FALSE)
    
    all_rrs <- c(wvs$blue, wvs$green)
    
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

# To catch files that give the following error when writing the new data to netCDF:
#   Error in Rsx_nc4_put_vara_double: NetCDF: Numeric conversion not representable
bad_files <- c()



#********************************
# LOOP THROUGH YEARS

for (i in 1:length(years)) {
    
    year <- years[i]
    
    in_path_year <- file.path(path, year)
    out_path_year <- file.path(path, year, variable)
    
    # create the output directories, if necessary
    dir.create(out_path_year, showWarnings=FALSE, recursive=TRUE)
    
    L2_files_year <- data.frame(files=list.files(in_path_year, pattern="L2_LAC_OC.x.nc"), stringsAsFactors = FALSE)
    
    if (nrow(L2_files_year)==0) {next}
    
    
    #********************************
    # LOOP THROUGH DAYS
    
    for (j in 1:length(days)) {
        
        day <- days[j]
        
        # Subset list of files to this day
        L2_files_day <- L2_files_year %>%
            dplyr::filter(grepl(paste0(year,str_pad(day,width=3,side="left",pad="0")), files))
        
        if (nrow(L2_files_day)==0) {next}
        
        
        #********************************
        # LOOP THROUGH FILES FOR THIS DAY (if more than one exists)
        
        for (fx in 1:nrow(L2_files_day)) {
            
            L2_name <- L2_files_day[fx,"files"]
            
            cat(paste0("Getting L2 data from ",L2_name,"...\n"))
            
            #*******************************************************************
            # GET L2 DATA
            
            in_file <- file.path(in_path_year,L2_name)
            l2 <- nc_open(in_file)
            variables <- names(l2$var)
            l2_chl <- ncvar_get(l2,"geophysical_data/chlor_a")
            rrs <- list()
            for (i in 1:length(all_rrs)) {
                rrs[[i]] <- ncvar_get(l2,paste0("geophysical_data/",all_rrs[i]))
            }
            nc_close(l2)
            
            # Get index of non-NA pixels from L2 chlor_a
            L2_dim <- dim(l2_chl)
            l2_chl[is.na(l2_chl) | (l2_chl > chl_range[2]) | (l2_chl < chl_range[1])] <- NA
            ind <- which(!is.na(l2_chl),arr.ind=TRUE)
            
            
            # COLLAPSE L2 INTO 1D VECTORS, ONE FOR EACH WAVEBAND
            # either take a single rrs, or the median of a 3x3 box
            
            if (rformat=="median-rrs") {
                # Get indices of pixels in the 3x3 box around the centre pixel.
                ind <- data.frame(ind)
                ind$rmin <- apply(cbind(ind$row - 1,rep(0,nrow(ind))),1,max)
                ind$rmax <- apply(cbind(ind$row + 1,rep(L2_dim[1],nrow(ind))),1,min)
                ind$cmin <- apply(cbind(ind$col - 1,rep(0,nrow(ind))),1,max)
                ind$cmax <- apply(cbind(ind$col + 1,rep(L2_dim[2],nrow(ind))),1,min)
                # Use parallel processing to get matrix of median rrs
                cl <- makeCluster(min(num_cl, detectCores()-1))
                clusterExport(cl, c('median_rrs','rrs','all_rrs'))
                rrs_sat <- parSapply(cl,1:nrow(ind),median_rrs,ind=ind,rrs=rrs,lambda=all_rrs)
                stopCluster(cl)
            } else if (rformat=="single-rrs") {
                cl <- makeCluster(min(num_cl, detectCores()-1))
                clusterExport(cl, c('single_rrs','rrs','all_rrs'))
                rrs_sat <- parSapply(cl,1:nrow(ind),single_rrs,ind=ind,rrs=rrs,lambda=all_rrs)
                stopCluster(cl)
            }
            
            # Transpose so that rows = records (pixels), columns = Rrs at each wavelength
            rrs <- t(rrs_sat)
            colnames(rrs) <- all_rrs
            
            
            #*******************************************************************
            # CALCULATE NEW CHLA 
            cat(paste0("Computing ", variable, "...\n"))
            
            if (variable == "CHL_OCX_RSG") {
                
                chl_oci <- oci(rrs, wvs$blues, wvs$green, coefs, use_443nm=TRUE,
                               sensor = ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)),
                               CI_coef_version = 1)
                chl_valid <- chl_oci$oci_chl
                
            } else if (variable == "CHL_POLY4") {
                
                chl_valid <- ocx(rrs, wvs$blue, wvs$green, coefs)
                
            } else if (variable == "CHL_EOF") {
                
                chl_valid <- eof_chl(rrs=as.data.frame(rrs), training_set=eof_training_set)
                
            } else if (variable == "CHL_GSM_GS") {
                
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
            new_chl_mat <- matrix(nrow=L2_dim[1],ncol=L2_dim[2])
            new_chl_mat[!is.na(l2_chl)] <- chl_valid
            attributes(new_chl_mat)$dim <- L2_dim
            
            # fix values so there are no errors
            new_chl_mat[!is.finite(new_chl_mat) | new_chl_mat < chl_range[1] | new_chl_mat > chl_range[2]] <- NA
            
            
            #*******************************************************************
            # ADD NEW CHLA TO OUTPUT NETCDF
            cat(paste0("Adding ", variable, " layer to L2 NetCDF...\n\n"))
            
            # Copy the original L2 file to a new test file.
            old_file <- file.path(in_path_year, L2_name)
            new_file <- file.path(out_path_year, paste0(substr(L2_name, 1, nchar(L2_name)-3), "_", variable, "_", rformat, ".nc"))
            
            
            file_creation <- try({
                
                # Copy the original L2 file to a new file.
                file.copy(from=old_file, to=new_file, overwrite=FALSE)
                
                # Add new variable(s) to the NetCDF file
                # https://rdrr.io/cran/ncdf4/man/ncvar_add.html
                L2 <- nc_open(new_file, write=TRUE)
                dim_numlines <- L2$dim[['number_of_lines']]
                dim_pixlines <- L2$dim[['pixels_per_line']]
                outdims <- list(dim_pixlines,dim_numlines) # keep the dimensions in this order
                # Define new variable
                chla_var <- ncvar_def(name=paste0("geophysical_data/", variable),
                                      units="mg m^-3",
                                      dim=outdims,
                                      missval=-32767,
                                      longname=chl_longname)
                # Add variable to output NetCDF
                L2 <- ncvar_add(L2, chla_var)
                # Write data to the variable
                ncvar_put(L2, chla_var, vals=new_chl_mat)
                nc_close(L2)
                
            }, silent=TRUE)
            
            
            # If the file wasn't created properly, remove it and add this to the list of bad files.
            if (class(file_creation)=="try-error") {
                file.remove(new_file)
                bad_files <- c(bad_files, old_file)
            }
            
        }
        
    }
    
}


# if (nrow(bad_files) > 0) {
#     bad_file_fname <- paste0(variable, "_L2_bad_files_", sensor, "_", region, "_", paste0(range(years), collapse="-"), paste0(range(days), collapse="-"), ".csv")
#     write.csv(matrix(bad_files, ncol=1), bad_file_fname, row.names=FALSE)
# }


stop()

# testing differences between NASA OCI and RSG OCI
# they should be the same, otherwise I'm doing something different than NASA

library(ggplot2)
library(patchwork)

x <- l2_chl[is.finite(l2_chl)]
chl_oci <- oci(rrs, wvs$blues, wvs$green, coefs, use_443nm=TRUE,
               sensor = ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)),
               CI_coef_version = 1)
chl_valid <- chl_oci$oci_chl
y <- chl_valid
diff <- y - x
pdiff <- (y-x)/x

hu_ind <- chl_oci$hu_ind
blend_ind <- chl_oci$blend_ind
ocx_ind <- !hu_ind & !blend_ind

tmp_df <- data.frame(nasa=x, rsg=y,
                     diff=diff, pdiff=pdiff,
                     algorithm = ifelse(hu_ind, "hu", ifelse(blend_ind, "blend", "ocx")),
                     stringsAsFactors = FALSE)

# p1 <- ggplot(tmp_df, aes(x=nasa, y=rsg, color=algorithm)) +
#     geom_point() +
#     scale_x_log10(limits=c(0.05,4)) +
#     scale_y_log10(limits=c(0.05,4))
# p2 <- ggplot(tmp_df) + geom_density(aes(x=diff, fill=algorithm), alpha=0.6)
p3 <- ggplot(tmp_df, aes(x=nasa, y=pdiff, color=algorithm)) +
    geom_point() +
    scale_x_log10(limits=c(0.05,4))

# print(p1 / p2)
print(p3)

sum(diff,na.rm=TRUE)/sum(x, na.rm=TRUE) * 100

sum(x != y)
sum(diff > 0.01)

tmp_df[diff > 0.01,]


# there are still some hu values different from nasa
# is it the color index? or the coefficients? or the conv_rrs_to_555 function? or the input itself?
tmp_df[tmp_df$algorithm=="hu" & tmp_df$diff < -0.005,]


stop()


CI_bound1 <- 0.15
CI_bound2 <- 0.2
chl_oci <- oci(rrs, wvs$blues, wvs$green, coefs, use_443nm=TRUE,
               sensor = ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)),
               CI_coef_version = 1, CI_bound1=CI_bound1, CI_bound2=CI_bound2)
chl_valid <- chl_oci$oci_chl
hu_ind <- chl_oci$hu_ind
blend_ind <- chl_oci$blend_ind
ocx_ind <- !hu_ind & !blend_ind
chl_ocx <- ocx(rrs, wvs$blues, wvs$green, coefs, use_443nm=TRUE)
chl_hu <- hu(rrs, get_ci_bands("modis"), get_ci_coefs(1))
#***
blend_fn <- function(chl_hu, chl_ocx, CI_bound1=0.15, CI_bound2=0.2) {
    a <- (chl_hu - CI_bound1)/(CI_bound2 - CI_bound1)
    b <- (CI_bound2 - chl_hu)/(CI_bound2 - CI_bound1)
    blend_chl <- a * chl_ocx + b * chl_hu
    return(list(blend_chl, a, b))
}

tmp_df <- data.frame(nasa=x, oci=chl_valid, ocx=chl_ocx, hu=chl_hu, blend=blend_fn(chl_hu, chl_ocx, CI_bound1, CI_bound2)[[1]],
                     algorithm = ifelse(hu_ind, "hu", ifelse(blend_ind, "blend", "ocx")),
                     stringsAsFactors = FALSE)
xlims <- c(0.1, 0.3)
ylims <- c(-0.3, 0.3)
p1 <- ggplot(tmp_df, aes(x=nasa, y=(oci-nasa)/nasa, color=algorithm)) +
    geom_point() +
    scale_x_log10(limits=xlims) +
    scale_y_continuous(limits=ylims)
p2 <- ggplot(tmp_df, aes(x=nasa, y=(ocx-nasa)/nasa, color=algorithm)) +
    geom_point() +
    scale_x_log10(limits=xlims) +
    scale_y_continuous(limits=ylims)
p3 <- ggplot(tmp_df, aes(x=nasa, y=(hu-nasa)/nasa, color=algorithm)) +
    geom_point() +
    scale_x_log10(limits=xlims) +
    scale_y_continuous(limits=ylims)
p4 <- ggplot(tmp_df, aes(x=nasa, y=(blend-nasa)/nasa, color=algorithm)) +
    geom_point() +
    scale_x_log10(limits=xlims) +
    scale_y_continuous(limits=ylims)
p1 / p2 / p3 / p4

tmp_df[tmp_df$algorithm=="blend" & tmp_df$nasa < 0.2,]





# assume your ocx component is accurate
# work backward - if all the difference is coming from hu, how much?

tmp_df$a = blend_fn(chl_hu, chl_ocx, CI_bound1, CI_bound2)[[2]]
tmp_df$b = blend_fn(chl_hu, chl_ocx, CI_bound1, CI_bound2)[[3]]

# a * chl_ocx + b * chl_hu
tmp_df$aocx = tmp_df$a * chl_ocx
tmp_df$bhu = tmp_df$b * chl_hu

tmp_df[tmp_df$algorithm=="blend" & tmp_df$nasa < 0.2,]


