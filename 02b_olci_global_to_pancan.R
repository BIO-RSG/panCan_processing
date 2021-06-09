# Stephanie.Clay@dfo-mpo.gc.ca
# 11 Sep 2020

# Subset GlobColour OLCI files to the panCanadian grid.

library(pbapply)
library(parallel)
mult_num_cl <- 6 # or whatever works best, leaving at least one processing core free
library(lubridate)
library(ncdf4)
# library(sinkr) # for isin.convert
library(data.table)
library(dplyr)
library(oceancolouR)

waves <- c(754,779,865,885,1020)
nwaves <- c(400,412,443,490,510,555,560,620,670,674,681,709)

sensor <- "OLCI-B"

variable <- "CHL1" # CHL1, CHL2, CHL-OC5, or RRS

years <- 2020

days <- 281:324

input_path <- "/mnt/data2/claysa/"
output_path <- "/mnt/data3/claysa/"


#*******************************************************************************
# FUNCTIONS

subset_to_pancan <- function(i, input_file, input_var_name, pancan_colsrows) {
    
    input_file <- input_file[i]
    input_var_name <- input_var_name[i]
    
    nc <- nc_open(input_file)
    var_values <- ncvar_get(nc, input_var_name)
    # Longitudinal index of the bins stored in the product, zero based and beginning at west
    olci_cols <- ncvar_get(nc, "col")
    olci_rows <- nc$dim$row$vals
    nc_atts <- ncatt_get(nc, varid=0)
    nc_close(nc)
    
    start_time <- nc_atts$start_time
    olci_colsrows <- paste0(olci_cols, "x", olci_rows)
    
    # from global 4km-resolution raster, match olci to panCan grid
    matched_df <- left_join(data.frame(colsrows=pancan_colsrows, stringsAsFactors = FALSE),
                            data.frame(colsrows=olci_colsrows,
                                       var=var_values,
                                       stringsAsFactors = FALSE),
                            by="colsrows") %>%
        dplyr::select(., -colsrows)
    
    # # QUICK VIEW OF RESULTING RASTER
    # library(latticeExtra)
    # library(maptools)
    # data("wrld_simpl")
    # olci_rast <- var_to_rast(data.frame(bin=bins_out, var=matched_df$var))
    # tmp_extent <- extent(c(xmn=-140, xmx=-40, ymn=38, ymx=70))
    # spplot(crop(olci_rast, tmp_extent)) + layer(sp.polygons(wrld_simpl))
    
    # get variable for time dimension
    start_time <- as_datetime(start_time, format="%Y%m%dT%H%M%SZ")
    start_time <- as.numeric(as.POSIXct(start_time))
    
    return(list(matched_df=matched_df,
                start_time=start_time))
}



# input file is multiple files
# output file is only one file, with just "RRS", no wavebands
convert_multiple_files <- function(j, files, file_dates, full_input_path, full_output_path, input_var_name, pancan_colsrows,
                                   output_short_var_name, output_long_var_name, dim_bindata, var_units) {
    
    input_file <- sort(files[grep(paste0("L3b_", file_dates[j], "_"), files)])
    
    output_file <- paste0(gsub("_NRRS[[:digit:]]+_", "_RRS_", input_file[1]))
    output_file <- gsub(".nc", "_panCan.nc", output_file)
    
    # add paths
    input_file <- paste0(full_input_path, input_file)
    output_file <- paste0(full_output_path, output_file)
    
    # sort input_file by input_var_name
    input_file <- input_file[as.numeric(sapply(input_var_name, grep, input_file))]
    
    
    # Adjust variable names to get flags as well
    input_var_name <- paste0(rep(input_var_name, each=2), rep(c("_mean", "_flags"), length(input_var_name)))
    output_short_var_name <- paste0(rep(output_short_var_name, each=2), rep(c("", "_flags"), length(output_short_var_name)))
    output_long_var_name <- paste0(rep(c("", "Flags for "), length(output_long_var_name)), rep(output_long_var_name, each=2))
    input_file <- rep(input_file, each=2)
    
    sub_olci <- lapply(1:length(input_var_name),
                       FUN=subset_to_pancan,
                       input_file=input_file,
                       input_var_name=input_var_name,
                       pancan_colsrows=pancan_colsrows)
    
    start_time <- sub_olci[[1]]$start_time
    if (is.character(start_time)) {
        start_time <- as_datetime(start_time, format="%Y%m%dT%H%M%SZ")
        start_time <- as.numeric(as.POSIXct(start_time))
    }
    
    # create dimensions
    dim_time <- ncdim_def(name = "time",
                          units = "seconds since 1970-01-01 00:00:00.000 +0000",
                          vals = start_time,
                          unlim = TRUE,
                          calendar = "standard")
    
    output_var <- list()
    
    for (n in 1:length(output_short_var_name)) {
        
        # create base output chlorophyll variable
        output_var[[n]] <- ncvar_def(name=output_short_var_name[n],
                                     units=var_units,
                                     dim=list(dim_time, dim_bindata),
                                     missval=NA,
                                     longname=output_long_var_name[n])
        
    }
    
    # create new output netcdf
    ncout <- nc_create(filename=output_file,
                       vars=output_var,
                       force_v4=TRUE)
    
    for (n in 1:length(input_var_name)) {
        
        # put variables in file
        ncvar_put(ncout, output_var[[n]], vals=sub_olci[[n]]$matched_df$var)
        
    }
    
    # close file
    nc_close(ncout)
    
}



convert_single_file <- function(j, files, full_input_path, full_output_path, input_var_name, pancan_colsrows,
                                output_short_var_name, output_long_var_name, dim_bindata, var_units) {
    
    
    if (length(input_var_name)==1) {
        
        input_file <- paste0(full_input_path, files[j])
        output_file <- paste0(full_output_path, paste0(substr(files[j], 1, nchar(files[j])-3), "_panCan.nc"))
        
        sub_olci <- subset_to_pancan(1, input_file=input_file,
                                     input_var_name=input_var_name,
                                     pancan_colsrows=pancan_colsrows)
        
        # create dimensions
        dim_time <- ncdim_def(name = "time",
                              units = "seconds since 1970-01-01 00:00:00.000 +0000",
                              vals = sub_olci$start_time,
                              unlim = TRUE,
                              calendar = "standard")
        
        # create base output chlorophyll variable
        output_var <- ncvar_def(name=output_short_var_name,
                                units=var_units,
                                dim=list(dim_time, dim_bindata),
                                missval=NA,
                                longname=output_long_var_name)
        
        # create new output netcdf
        ncout <- nc_create(output_file, output_var, force_v4=TRUE)
        
        # put variables in file
        ncvar_put(ncout, output_var, vals=sub_olci$matched_df$var)
        
        # close file
        nc_close(ncout)
        
    } else {
        
        input_file <- paste0(full_input_path, files[j])
        output_file <- paste0(full_output_path, paste0(substr(files[j], 1, nchar(files[j])-3), "_panCan.nc"))
        
        input_file <- rep(input_file, length(input_var_name))
        
        sub_olci <- lapply(1:length(input_var_name),
                           FUN=subset_to_pancan,
                           input_file=input_file,
                           input_var_name=input_var_name,
                           pancan_colsrows=pancan_colsrows)
        
        start_time <- sub_olci[[1]]$start_time
        if (is.character(start_time)) {
            start_time <- as_datetime(start_time, format="%Y%m%dT%H%M%SZ")
            start_time <- as.numeric(as.POSIXct(start_time))
        }
        
        # create dimensions
        dim_time <- ncdim_def(name = "time",
                              units = "seconds since 1970-01-01 00:00:00.000 +0000",
                              vals = start_time,
                              unlim = TRUE,
                              calendar = "standard")
        
        output_var <- list()
        
        for (n in 1:length(output_short_var_name)) {
            
            # create base output chlorophyll variable
            output_var[[n]] <- ncvar_def(name=output_short_var_name[n],
                                         units=var_units,
                                         dim=list(dim_time, dim_bindata),
                                         missval=NA,
                                         longname=output_long_var_name[n])
            
        }
        
        # create new output netcdf
        ncout <- nc_create(filename=output_file,
                           vars=output_var,
                           force_v4=TRUE)
        
        for (n in 1:length(input_var_name)) {
            
            # put variables in file
            ncvar_put(ncout, output_var[[n]], vals=sub_olci[[n]]$matched_df$var)
            
        }
        
        # close file
        nc_close(ncout)
        
    }
    
}


#*******************************************************************************
# SUBSET L3B OLCI IN DATA2/ TO L3B PANCAN OLCI IN DATA3/

# depending on variable name, set up the variable name from the input file, and the names and units for the output file
if (variable=="CHL1") {
    input_var_name <- c("CHL1_mean", "CHL1_flags")
    output_short_var_name <- c("chlor_a", "flags")
    output_long_var_name <- c("Chlorophyll-a Concentration", "Flags")
    var_units <- c("mg m^-3", NULL)
} else if (variable=="CHL2") {
    input_var_name <- c("CHL2_mean", "CHL2_flags")
    output_short_var_name <- c("chlor_a", "flags")
    output_long_var_name <- c("Chlorophyll-a Concentration, Neural Network algorithm", "Flags")
    var_units <- c("mg m^-3", NULL)
} else if (variable=="CHL-OC5") {
    input_var_name <- c("CHL-OC5_mean", "CHL-OC5_flags")
    output_short_var_name <- c("chlor_a", "flags")
    output_long_var_name <- c("Chlorophyll-a Concentration, OC5 algorithm", "Flags")
    var_units <- c("mg m^-3", NULL)
} else if (variable=="RRS") {
    input_var_name <- c(paste0("NRRS", nwaves),
                        paste0("RRS", waves))
    output_short_var_name <- paste0("Rrs_", c(nwaves, waves))
    output_long_var_name <- c(paste0("Normalized Remote Sensing Reflectance at ", nwaves, "nm"),
                              paste0("Remote Sensing Reflectance at ", waves, "nm"))
    var_units <- rep("sr^-1", length(input_var_name))
}

# # CREATE .CSV CONTAINING COL/ROW VALUES BASED ON PANCAN BINS - need this to match to OLCI col/row
# # takes ~7-8 minutes to run
# data("pancan_bins_4km", package="oceancolouR")
# pancan_isin <- lapply(pancan_bins_4km, isin.convert)
# pancan_cols <- as.numeric(sapply(pancan_isin, "[[", 2))
# pancan_rows <- as.numeric(sapply(pancan_isin, "[[", 3))
# fwrite(data.frame(col=pancan_cols,
#                   row=pancan_rows,
#                   stringsAsFactors = FALSE),
#        file="data/panCan_ColsRows_from_ISIN.csv",
#        quote=FALSE)


# Note that you must subtract 1 because globcolour OLCI is zero-indexed
pancan_colsrows <- fread("data/panCan_ColsRows_from_ISIN.csv")
pancan_cols <- pancan_colsrows$col - 1
pancan_rows <- pancan_colsrows$row - 1
pancan_colsrows <- paste0(pancan_cols, "x", pancan_rows)

# create bin dimension, same for each pancan file
dim_bindata <- ncdim_def(name="binDataDim",
                         units = "",
                         vals = 1:529797,
                         unlim = FALSE,
                         create_dimvar = FALSE)


for (i in 1:length(years)) {
    
    full_input_path <- paste0(input_path, sensor, "/", variable, "/", years[i], "/")
    full_output_path <- paste0(output_path, sensor, "/", variable, "/PANCAN/", years[i], "/")
    get_dir(full_output_path)
    
    files <- list.files(path=full_input_path, pattern=".nc")
    files <- files[grep(variable, files)]
    files <- sort(files)
    
    if (length(files)==0) {next}
    
    # restrict files to selected days
    file_dates <- sapply(strsplit(files, "_"), "[[", 2)
    file_doys <- as.numeric(sapply(1:length(file_dates), function(x) format(as.Date(file_dates[x], format="%Y%m%d"), "%j")))
    doy_ind <- file_doys %in% days
    files <- files[doy_ind]
    file_dates <- unique(file_dates[doy_ind])
    
    if (length(files)==0) {
        next
    } else if (length(files)==1) {
        num_cl <- 1
    } else {
        num_cl <- mult_num_cl
    }
    
    
    # FOR ONE/MULTIPLE VARIABLES IN A SINGLE FILE, WRITTEN TO A SINGLE FILE
    # EXAMPLE: GLOBCOLOUR OLCI CHL1
    
    if (grepl("CHL", variable)) {
        
        # Create clusters for parallel processing.
        cl <- makeCluster(num_cl)
        # Load necessary variables or custom functions into cluster.
        clusterExport(cl,
                      varlist = c("files", "full_input_path", "full_output_path", "input_var_name", "pancan_colsrows",
                                  "output_short_var_name", "output_long_var_name",
                                  "dim_bindata", "var_units", "convert_single_file", "subset_to_pancan"),
                      envir = environment())
        # Load necessary libraries into cluster.
        clusterEvalQ(cl, c(library(ncdf4), library(dplyr), library(lubridate)))
        # Like lapply, but with the clusters variable cl as the first argument.
        pblapply(1:length(files), convert_single_file,
                 files=files, full_input_path=full_input_path, full_output_path=full_output_path, input_var_name=input_var_name,
                 pancan_colsrows=pancan_colsrows, output_short_var_name=output_short_var_name,
                 output_long_var_name=output_long_var_name, dim_bindata=dim_bindata, var_units=var_units,
                 cl=cl)
        # Stop parallel processing and return processing power to other operations.
        stopCluster(cl)
        
    } else if (grepl("RRS", variable)) {
        
        # FOR MULTIPLE VARIABLES IN MULTIPLE FILES, WRITTEN TO A SINGLE FILE
        # Example: GLOBCOLOUR OLCI RRS
        
        # Create clusters for parallel processing.
        cl <- makeCluster(num_cl)
        # Load necessary variables or custom functions into cluster.
        clusterExport(cl,
                      varlist = c("files", "file_dates", "full_input_path", "full_output_path", "input_var_name", "pancan_colsrows",
                                  "output_short_var_name", "output_long_var_name", 
                                  "dim_bindata", "var_units", "convert_multiple_files", "subset_to_pancan"),
                      envir = environment())
        # Load necessary libraries into cluster.
        clusterEvalQ(cl, c(library(ncdf4), library(dplyr), library(lubridate)))
        # Like lapply, but with the clusters variable cl as the first argument.
        pblapply(1:length(file_dates), convert_multiple_files,
                 files=files, file_dates=file_dates, full_input_path=full_input_path, full_output_path=full_output_path, input_var_name=input_var_name,
                 pancan_colsrows=pancan_colsrows, output_short_var_name=output_short_var_name,
                 output_long_var_name=output_long_var_name, dim_bindata=dim_bindata, var_units=var_units,
                 cl=cl)
        # Stop parallel processing and return processing power to other operations.
        stopCluster(cl)
        
    }
    
}


# # VIEW TEST IMAGE
# path = "/mnt/data3/claysa/OLCI-B/CHL1/PANCAN/2019/"
# fname = "L3b_20190724__296881628_4_AV-OLB_CHL1_DAY_00_panCan.nc"
# test = nc_open(paste0(path, fname))
# chlor_a = ncvar_get(test, "chlor_a")
# nc_close(test)
# plot_rast_from_bin(vec=log10(as.numeric(chlor_a)))

