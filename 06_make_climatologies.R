# Stephanie.Clay@dfo-mpo.gc.ca
# 7 Oct 2020

# Create 8day, monthly, or annual climatologies based on L3b files.
# 
# Input: files created by 05_make_composites.R
#
# Currently only compatible with 4km-resolution files.
#
# For 8day, monthly, or annual climatologies:
# Average over the appropriate days for each year, then average over years
# (this gives each year equal weight).

library(ncdf4)
library(oceancolouR)
library(stringr)
library(dplyr)
library(Rfast) # for colMedians and std
library(fst)

sensors <- c("MODIS", "SeaWiFS", "VIIRS-SNPP")
variables <- c("CHL_OCX", "CHL_POLY4", "CHL_GSM_GS", "CHL_EOF", "PAR", "SST")#, "RRS")

composite <- "monthly" # 8day, monthly, or annual

base_input_path <- "/mnt/data3/claysa"
output_path <- file.path(base_input_path, "climatologies")

# for monthly or annual composites, which type of mean should be used to first summarize the daily values across a selected month, or year?
# (after these are calculated, another set of stats are calculated across the full series of years, and those are all saved to .fst)
# options: geometric_mean or arithmetic_mean
which_mean <- "geometric_mean"


#*******************************************************************************

all_regions <- list("CHL_OCX"=c("PANCAN"),
                    "CHL_POLY4"=c("NWA", "NEP"),
                    "CHL_GSM_GS"=c("NWA", "NEP"),
                    "CHL_EOF"=c("GoSL"),
                    "CHL1"=c("PANCAN"),
                    "CHL2"=c("PANCAN"),
                    "PAR"=c("PANCAN"),
                    #"RRS"=c("PANCAN"),
                    "SST"=c("PANCAN"))

all_years <- list("MODIS"=2003:2020,
                  "SeaWiFS"=1997:2010,
                  "VIIRS-SNPP"=2012:2020,
                  "OLCI-A"=2016:2020,
                  "OLCI-B"=2018:2020)

if (composite=="8day") {
    num_loops <- 46
} else if (composite=="monthly") {
    num_loops <- 12
} else if (composite=="annual") {
    num_loops <- 1
}


#*******************************************************************************

for (variable in variables) {
    
    print(variable)
    
    if (startsWith(variable, "CHL")) {
        input_variable <- "chlor_a"
        var_units <- "mg m^-3"
    } else if (variable=="SST") {
        input_variable <- "sst"
        var_units <- "degrees Celsius"
    } else if (variable=="PAR") {
        input_variable <- "par"
        var_units <- "Einstein m^-2 d^-1"
    }# else if (variable=="RRS") {
    #     input_variable <- paste0("Rrs_", waves)
    #     var_units <- "sr^-1"
    # }
    
    regions <- all_regions[[variable]]
    
    
    for (region in regions) {
        
        print(region)
        
        tmp_num_pix <- num_pix[[region]][["4km"]]
        bins <- (function(v) get(data(list = v, package = "oceancolouR", envir = new.env())))(paste0(tolower(region), "_bins_4km"))
        
        # # for netcdf output, create bin dimension
        # dim_bindata <- ncdim_def(name="binDataDim",
        #                          units = "",
        #                          vals = 1:length(bins),
        #                          unlim = FALSE,
        #                          create_dimvar = FALSE)
        
        for (sensor in sensors) {
            
            print(sensor)
            
            years <- all_years[[sensor]]
            
            full_input_path <- file.path(base_input_path, sensor, variable, region, "composites")
            
            for (i in 1:num_loops) {
                
                # gather files across the time series for this week/month
                file_pattern <- paste0(sensor, "_", variable, "_", region)
                all_files <- sort(list.files(path=full_input_path, pattern=file_pattern, recursive=TRUE))
                all_files <- all_files[grepl(paste0(composite, "_composite", ifelse(composite=="annual", "", paste0("_",str_pad(i,width=2,side="left",pad="0")))), all_files)]
                
                # read data from input files
                df_mean <- matrix(NA, nrow=tmp_num_pix, ncol=length(all_files))
                df_num_obs <- matrix(NA, nrow=tmp_num_pix, ncol=length(all_files))
                for (file_ind in 1:length(all_files)) {
                    df_mean[,file_ind] <- read_fst(file.path(full_input_path,all_files[file_ind]))[,which_mean]
                    df_num_obs[,file_ind] <- read_fst(file.path(full_input_path,all_files[file_ind]))[,"num_obs"]
                }
                
                # calculate stats per pixel for this composite across all years
                df_mean[!is.finite(df_mean)] <- NA
                df_num_obs[!is.finite(df_num_obs)] <- NA
                
                num_obs <- rowSums(is.finite(df_num_obs))
                geometric_mean <- rep(NA, tmp_num_pix)
                geometric_mean[num_obs>0] <- apply(df_mean[num_obs>0,], 1, geoMean, na.rm = TRUE)
                arithmetic_mean <- rep(NA, tmp_num_pix)
                arithmetic_mean[num_obs>0] <- rowMeans(df_mean[num_obs>0,], na.rm = TRUE)
                sat_median <- rep(NA, tmp_num_pix)
                sat_median[num_obs>0] <- rowMedians(df_mean[num_obs>0,], na.rm = TRUE)
                sat_sd <- rep(NA, tmp_num_pix)
                sat_sd[num_obs>0] <- rowVars(df_mean[num_obs>0,], std = TRUE, na.rm = TRUE)
                sat_sd[!is.finite(sat_sd)] <- NA
                
                output_df <- data.frame(bin = bins,
                                        geometric_mean = geometric_mean,
                                        arithmetic_mean = arithmetic_mean, 
                                        median = sat_median,
                                        sd = sat_sd,
                                        num_years = num_obs,
                                        num_obs_total = rowSums(df_num_obs, na.rm=TRUE),
                                        stringsAsFactors = FALSE)
                
                
                #***************************************************************
                # SAVE AS .FST (dataframe)
                
                output_fname <- paste0(sensor, "_", variable, "_", region, "_",
                                       ifelse(which_mean=="geometric_mean", "geomean_", "mean_"),
                                       composite, "_climatology_", paste0(range(years), collapse="-"),
                                       ifelse(composite=="annual", "", paste0("_",str_pad(i,width=2,side="left",pad="0"))), ".fst")
                
                write_fst(output_df, path=file.path(output_path, output_fname), compress=100)
                
                
                # #***************************************************************
                # # SAVE AS NETCDF
                # 
                # output_fname <- gsub(".fst", ".nc", output_fname)
                # output_var <- list()
                # for (i in 2:ncol(all_data)) {
                #     vname <- colnames(all_data)[i]
                #     output_var[[i-1]] <- ncvar_def(name=vname,
                #                                    units=var_units,
                #                                    dim=list(dim_bindata),
                #                                    missval=NA,
                #                                    longname=vname)
                # }
                # # create new output netcdf
                # ncout <- nc_create(filename=output_fname, vars=output_var, force_v4=TRUE)
                # # put variables in file
                # for (i in 1:length(output_var)) {
                #     ncvar_put(ncout, output_var[[i]], vals=all_data[,i+1])
                # }
                # nc_close(ncout)
                # 
                # #***************************************************************
                
            }
            
        } # sensor loop
        
    } # region loop
    
} # variable loop
