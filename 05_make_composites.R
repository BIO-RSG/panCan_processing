# Stephanie.Clay@dfo-mpo.gc.ca
# 7 Oct 2020

# Create 8day, monthly, or annual composites based on L3b files.
# 
# Input: daily L3b netCDF files in panCanadian dataset.
#
# Currently only compatible with 4km-resolution files.

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
            
            full_input_path <- file.path(base_input_path, sensor, variable, region)
            output_path <- file.path(full_input_path, "composites")
            dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
            
            for (y in years) {
                
                print(y)
                
                for (i in 1:num_loops) {
                    
                    if (composite=="annual") {
                        doy_vec <- str_pad(1:366, width=3, side="left", pad="0")
                    } else if (composite=="8day") {
                        doy_vec <- str_pad(days_vector(year=y, week=i), width=3, side="left", pad="0")
                    } else if (composite=="monthly") {
                        doy_vec <- str_pad(days_vector(year=y, month=i), width=3, side="left", pad="0")
                    }
                    
                    files <- sort(unlist(sapply(paste0(y, doy_vec, ".L3b"), list.files, path=full_input_path, recursive=TRUE)))
                    
                    # get a dataframe of stats for each pixel across a number of days in a week/month/year
                    
                    if (length(files)==0) {
                        
                        output_df <- data.frame(bin = bins,
                                                geometric_mean=double(tmp_num_pix),
                                                arithmetic_mean=double(tmp_num_pix), 
                                                median=double(tmp_num_pix),
                                                sd=double(tmp_num_pix),
                                                num_obs=rep(0,tmp_num_pix),
                                                stringsAsFactors = FALSE)
                        
                    } else {
                        
                        # get the day of each existing netcdf filename
                        if (sensor %in% c("OLCI-A", "OLCI-B")) {
                            file_days <- as.numeric(format(as.Date(sapply(files, substr, start=10, stop=17), format="%Y%m%d"), "%j"))
                        } else {
                            file_days <- as.numeric(sapply(files, substr, start=11, stop=13))
                        }
                        
                        # read daily data into dataframe
                        data_composite <- matrix(NA, nrow=tmp_num_pix, ncol=length(doy_vec))
                        dv_num <- as.numeric(doy_vec)
                        for (k in 1:length(dv_num)) {
                            if (dv_num[k] %in% file_days) {
                                d <- nc_open(filename = file.path(full_input_path, files[which(file_days==dv_num[k])]))
                                data_composite[,k] <- ncvar_get(nc = d, varid = input_variable)
                                nc_close(d)
                            } else {
                                data_composite[,k] <- rep(NA, tmp_num_pix)
                            }
                        }
                        
                        # calculate stats per pixel for this composite
                        data_composite[!is.finite(data_composite)] <- NA
                        num_obs <- rowSums(is.finite(data_composite))
                        geometric_mean <- rep(NA, tmp_num_pix)
                        geometric_mean[num_obs>0] <- apply(data_composite[num_obs>0,], 1, geoMean, na.rm = TRUE)
                        arithmetic_mean <- rep(NA, tmp_num_pix)
                        arithmetic_mean[num_obs>0] <- rowMeans(data_composite[num_obs>0,], na.rm = TRUE)
                        sat_median <- rep(NA, tmp_num_pix)
                        sat_median[num_obs>0] <- rowMedians(data_composite[num_obs>0,], na.rm = TRUE)
                        sat_sd <- rep(NA, tmp_num_pix)
                        sat_sd[num_obs>0] <- rowVars(data_composite[num_obs>0,], std = TRUE, na.rm = TRUE)
                        sat_sd[!is.finite(sat_sd)] <- NA
                        
                        output_df <- data.frame(bin = bins,
                                                geometric_mean = geometric_mean,
                                                arithmetic_mean = arithmetic_mean, 
                                                median = sat_median,
                                                sd = sat_sd,
                                                num_obs = num_obs,
                                                stringsAsFactors = FALSE)
                        
                    }
                    
                    
                    #***********************************************************
                    # SAVE AS .FST (dataframe)
                    
                    output_fname <- paste0(sensor, "_", variable, "_", region, "_", y, "_", composite, "_composite",
                                           ifelse(composite=="annual", "", paste0("_",str_pad(i,width=2,side="left",pad="0"))), ".fst")
                    
                    write_fst(output_df, path=file.path(output_path, output_fname), compress=100)
                    
                    # #***********************************************************
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
                    # #***********************************************************
                    
                } # composite loop
                
            } # year loop
            
        } # sensor loop
        
    } # region loop
    
} # variable loop
