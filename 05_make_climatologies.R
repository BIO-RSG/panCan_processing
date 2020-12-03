# Stephanie.Clay@dfo-mpo.gc.ca
# 7 Oct 2020

# Create daily, monthly, or annual climatologies based on L3b files.
#
# For monthly or daily climatologies:
# Average over the appropriate days for each year, then average over years
# (this gives each year equal weight).

library(ncdf4)
library(oceancolouR)
library(stringr)
library(dplyr)
# library(Rfast) # for colMedians and std
library(fst)

sensor <- "MODIS" # MODIS, SeaWiFS, or VIIRS-SNPP
variable <- "CHL_OCX" # CHL_OCX, CHL_POLY4, CHL_GSM_GS, PAR, or SST
region <- "PANCAN" # PANCAN, NWA, or NEP

composite <- "monthly" # daily, monthly, or annual

base_input_path <- "/mnt/data3/claysa"
output_path <- file.path(base_input_path, "climatologies")

# for monthly or annual composites, which type of mean should be used to first summarize the daily values across a selected month, or year?
# (after these are calculated, another set of stats are calculated across the full series of years, and those are all saved to .fst)
# options: geometric_mean or arithmetic_mean
which_mean <- "geometric_mean"


#*******************************************************************************

all_years <- list("MODIS"=2003:2020,
                  "SeaWiFS"=1997:2010,
                  "VIIRS-SNPP"=2012:2020)

tmp_num_pix <- num_pix[[region]]

if (region=="PANCAN") {
    data("pancan_bins_4km")
    bins <- pancan_bins_4km
} else if (region=="NWA") {
    data("nwa_bins_4km")
    bins <- nwa_bins_4km
} else if (region=="NEP") {
    data("nep_bins_4km")
    bins <- nep_bins_4km
}

if (variable=="CHL_GSM_GS") {
    input_variable <- "chl_GSM_GS"
} else if (startsWith(variable, "CHL")) {
    input_variable <- "chlor_a"
} else if (variable=="SST") {
    input_variable <- "sst"
} else if (variable=="PAR") {
    input_variable <- "par"
}

if (composite=="daily") {
    num_loops <- 366
} else if (composite=="monthly") {
    num_loops <- 12
} else if (composite=="annual") {
    num_loops <- 1
}

#*******************************************************************************

for (i in 1:num_loops) {
    
    if (composite=="daily") {
        doy_vec <- str_pad(i, width=3, side="left", pad="0")
    } else if (composite=="monthly") {
        doy_vec <- str_pad(days_vector(year=2004, month=i), width=3, side="left", pad="0")
    } else if (composite=="annual") {
        doy_vec <- str_pad(1:366, width=3, side="left", pad="0")
    }
    
    full_input_path <- file.path(base_input_path, sensor, variable, region)
    
    all_files <- sort(unlist(sapply(paste0(doy_vec, ".L3b"), list.files, path=full_input_path, recursive=TRUE)))
    
    years <- all_years[[sensor]]
    
    # Collect data
    data_yrs <- list()
    
    for (j in 1:length(years)) {
        
        year <- years[j]
        
        print(year)
        
        files <- all_files[grepl(paste0(year, "/"), all_files)]
        
        # either get a vector of pixel data for this day, or a dataframe of stats
        # for each pixel across a number of days (for a month, or year)
        
        if (composite=="daily") {
            
            if (length(files)==0) {
                data_yrs[[j]] <- data.frame(bin = bins,
                                            variable = rep(NaN, tmp_num_pix),
                                            stringsAsFactors = FALSE)
            } else {
                d <- nc_open(filename = paste0(full_input_path, files))
                data_yrs[[j]] <- data.frame(bin = bins,
                                            variable = ncvar_get(nc = d, varid = input_variable),
                                            stringsAsFactors = FALSE)
                nc_close(d)
            }
            
        } else {
            
            if (length(files)==0) {
                
                data_yrs[[j]] <- data.frame(bin = bins,
                                            geometric_mean=double(tmp_num_pix),
                                            arithmetic_mean=double(tmp_num_pix), 
                                            # median=double(tmp_num_pix), 
                                            # sd=double(tmp_num_pix),
                                            num_obs=rep(0,tmp_num_pix),
                                            stringsAsFactors = FALSE)
                
            } else {
                
                # get the day of each existing netcdf filename
                if (sensor %in% c("OLCI-A", "OLCI-B")) {
                    file_days <- as.numeric(format(as.Date(sapply(files, substr, start=10, stop=17), format="%Y%m%d"), "%j"))
                } else {
                    file_days <- as.numeric(sapply(files, substr, start=11, stop=13))
                }
                
                data_composite <- matrix(nrow=tmp_num_pix, ncol=length(doy_vec))
                dv_num <- as.numeric(doy_vec)
                for (k in 1:length(dv_num)) {
                    if (dv_num[k] %in% file_days) {
                        d <- nc_open(filename = file.path(full_input_path, files[which(file_days==dv_num[k])]))
                        data_composite[,k] <- ncvar_get(nc = d, varid = input_variable)
                        nc_close(d)
                    } else {
                        data_composite[,k] <- rep(NaN, tmp_num_pix)
                    }
                }
                
                num_obs <- rowSums(is.finite(data_composite))
                geometric_mean <- rep(NaN, tmp_num_pix)
                geometric_mean[num_obs>0] <- apply(data_composite[num_obs>0,], 1, geoMean, na.rm = TRUE)
                arithmetic_mean <- rep(NaN, tmp_num_pix)
                arithmetic_mean[num_obs>0] <- rowMeans(data_composite[num_obs>0,], na.rm = TRUE)
                # median <- apply(data_composite, 1, median, na.rm = TRUE)
                # sd <- apply(data_composite, 1, sd, na.rm = TRUE)
                
                data_yrs[[j]] <- data.frame(bin = bins,
                                            geometric_mean = geometric_mean,
                                            arithmetic_mean = arithmetic_mean, 
                                            # median = median, 
                                            # sd = sd,
                                            num_obs = num_obs,
                                            stringsAsFactors = FALSE)
                
            }
            
        }
        
    } # finish looping through years
    
    
    # summarize data across available years
    if (composite=="daily") {
        
        all_data <- do.call(rbind, data_yrs) %>%
            dplyr::group_by(bin) %>% 
            dplyr::summarise(geometric_mean = geoMean(variable, na.rm = TRUE),
                             arithmetic_mean = mean(variable, na.rm = TRUE),
                             median = median(variable, na.rm = TRUE),
                             sd = sd(variable, na.rm = TRUE),
                             num_obs_day = sum(is.finite(variable), na.rm = TRUE)) %>%
            dplyr::ungroup()
        
    } else {
        
        all_data <- do.call(rbind, data_yrs) %>%
            dplyr::rename(variable = sym(which_mean)) %>%
            dplyr::group_by(bin) %>%
            dplyr::summarise(geometric_mean = geoMean(variable, na.rm = TRUE),
                             arithmetic_mean = mean(variable, na.rm = TRUE),
                             median = median(variable, na.rm = TRUE),
                             sd = sd(variable, na.rm = TRUE),
                             num_obs_total = sum(num_obs, na.rm = TRUE),
                             num_obs = sum(is.finite(variable), na.rm = TRUE)) %>%
            dplyr::ungroup()
        
        colnames(all_data)[7] <- paste0("num_obs_", ifelse(composite=="monthly", "month", "year"))
        
    }
    
    
    output_fname <- paste0(sensor, "_", variable, "_", region, "_climatology_", paste0(range(years), collapse="-"), "_", composite, ".fst")
    if (composite=="daily") {
        output_fname <- gsub(".fst", paste0("_", str_pad(i, width=3, side="left", pad="0"), ".fst"), output_fname)
    } else if (composite=="monthly") {
        output_fname <- gsub(".fst", paste0("_", str_pad(i, width=2, side="left", pad="0"), ".fst"), output_fname)
    }
    
    write_fst(all_data, path=file.path(output_path, output_fname), compress=100)
    
}
