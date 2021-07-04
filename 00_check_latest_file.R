# Stephanie.Clay@dfo-mpo.gc.ca
# 2020-10-25

# List the latest files that have been downloaded and processed to the panCanadian
# grid so you know where to pick up the processing.
# (Not the same as the script used to check which files are missing from the entire
# dataset, 03_check_for_missing_panCan_files.R)

# MODIS, SeaWiFS, VIIRS-SNPP, OLCI-A, OLCI-B
sensors <- c("MODIS", "SeaWiFS", "VIIRS-SNPP", "OLCI-A", "OLCI-B")

# last year downloaded for each sensor above (vector, same length as "sensors")
years <- c(2021, 2010, 2021, 2020, 2020)

# SENSOR VARIABLES:
#   MODIS/SeaWiFS/VIIRS-SNPP: CHL, PAR, RRS, CHL_POLY4, CHL_GSM_GS
#   OLCI-A: CHL1, CHL2, CHL-OC5, RRS
#   OLCI-B: CHL1
# REGIONS: PANCAN for all but CHL_POLY4 and CHL_GSM_GS, which are NWA and NEP

# variables/regions for each sensor
vr <- list("MODIS"=c("CHL_OCX/PANCAN", "PAR/PANCAN", "RRS/PANCAN", "SST/PANCAN",
                     "CHL_POLY4/NWA", "CHL_POLY4/NEP", "CHL_GSM_GS/NWA", "CHL_GSM_GS/NEP"),
           "SeaWiFS"=c("CHL_OCX/PANCAN", "PAR/PANCAN", "RRS/PANCAN", "SST/PANCAN",
                       "CHL_POLY4/NWA", "CHL_POLY4/NEP", "CHL_GSM_GS/NWA", "CHL_GSM_GS/NEP"),
           "VIIRS-SNPP"=c("CHL_OCX/PANCAN", "PAR/PANCAN", "RRS/PANCAN", "SST/PANCAN",
                          "CHL_POLY4/NWA", "CHL_POLY4/NEP", "CHL_GSM_GS/NWA", "CHL_GSM_GS/NEP"),
           "OLCI-A"=c("CHL1/PANCAN", "CHL2/PANCAN", "CHL-OC5/PANCAN", "RRS/PANCAN"),
           "OLCI-B"=c("CHL1/PANCAN"))


download_path <- "/mnt/data2/claysa"
process_path <- "/mnt/data3/claysa"


#*******************************************************************************

last_day <- data.frame(sensor = rep(sensors, sapply(vr, length)),
                       year = rep(years, sapply(vr, length)),
                       variable = as.character(sapply(strsplit(unlist(vr[names(vr) %in% sensors]), "/"), "[[", 1)),
                       region = as.character(sapply(strsplit(unlist(vr[names(vr) %in% sensors]), "/"), "[[", 2)),
                       downloaded = NA,
                       processed = NA,
                       stringsAsFactors = FALSE)

for (i in 1:nrow(last_day)) {
    
    sensor <- last_day[i,"sensor"]
    year <- last_day[i,"year"]
    variable <- last_day[i,"variable"]
    region <- last_day[i,"region"]
    
    # CHL_GSM_GS and CHL_POLY4 are only in the processed folders (created from downloaded RRS)
    if (variable %in% c("CHL_GSM_GS", "CHL_POLY4")) {
        dday <- NA
    } else {
        dpath <- file.path(download_path, sensor, ifelse(variable=="CHL_OCX", "CHL", variable), year)
        dfiles <- list.files(dpath)
        if (length(dfiles)==0) {
            dday <- NA
        } else {
            if (grepl("OLCI", sensor)) {
                ddays <- sort(as.numeric(format(as.Date(sapply(dfiles, substr, start=5, stop=12), format="%Y%m%d"), "%j")))
            } else if (grepl("SST", variable)) {
                ddays <- sapply(1:length(dfiles), function(x) format(as.Date(strsplit(dfiles, "[.]")[[x]][2], format="%Y%m%d"), "%j"))
            } else {
                ddays <- sort(as.numeric(sapply(dfiles, substr, start=6, stop=8)))
            }
            dday <- ddays[length(ddays)]
        }
    }
    
    ppath <- file.path(process_path, sensor, variable, region, year)
    pfiles <- list.files(ppath)
    
    if (length(pfiles)==0) {
        pday <- NA
    } else {
        if (grepl("OLCI", sensor)) {
            pdays <- sort(as.numeric(format(as.Date(sapply(pfiles, substr, start=5, stop=12), format="%Y%m%d"), "%j")))
        } else {
            pdays <- sort(as.numeric(sapply(pfiles, substr, start=6, stop=8)))
        }
        pday <- pdays[length(pdays)]
    }
    
    last_day[i,] <- c(sensor, year, variable, region, dday, pday)
    
}

cat("Most recent downloaded/processed day of year for each sensor/year/variable/region combination:\n\n")
print(last_day)

cat(paste0("\nCurrent year: ", format(Sys.Date(), "%Y")))
cat(paste0("\nCurrent day of year: ", format(Sys.Date(), "%j")))
