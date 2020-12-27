# Stephanie.Clay@dfo-mpo.gc.ca
# 21 Dec 2020

# Extract in situ chla and satellite Rrs matchups from the netCDF file created
# by the Python L2 matchup script, and use them to create a training set for the
# EOF chlorophyll-a algorithm.

# Reduce matchups by Julien Laliberte's criteria:
#   non-flagged pixels, 3x3 box, no matches with any negative Rrs, at least 6 pixels valid, use the median
#   +- 3hr window

# Note: this works for matchups with Rrs and chla in a 5x5 box.
# (If they're 3x3, adjust the get_var function below)


#*******************************************************************************
# ADJUST THESE VARIABLES

sensor <- "MODIS"       # MODIS, SeaWiFS, or VIIRS-SNPP
region <- "GoSL"        # NWA, NEP, or GoSL
years <- 1997:2019      # years of in situ data, used only in the output csv filename
max_dist <- 10000       # in metres
max_depth <- 10         # in metres
max_timediff <- 3       # in hours
matchup_file <- "../satellite_validation/02_matchups/chl_in_situ_gosl_1997-2019.txt_extNA_modis_matches.nc"


#*******************************************************************************
# LOAD LIBRARIES AND GET THE DATA

if (!file.exists(matchup_file)) {
    stop("Missing input matchup file.")
}

library(ncdf4)
library(oceancolouR)
library(geodist)
library(dplyr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(patchwork)

# get Rrs or chla, using Julien's criteria:
# non-flagged pixels, 3x3 box, no matches with any negative Rrs, at least 6 pixels valid, use the median
get_var <- function(var) {
    var <- matrix(var, nrow=5)
    var <- as.numeric(var[2:4,2:4])
    var <- var[is.finite(var)]
    var <- var[var >= 0]
    if (length(var) < 6) {return(NA)}
    return(median(var, na.rm=TRUE))
}


lambda <- all_lambda[[sensor]]

nc <- nc_open(matchup_file)
rrs <- list()
for (i in 1:length(lambda)) {
    rrs[[i]] <- ncvar_get(nc, paste0("processing_data/Rrs_", lambda[i]))
}
islat <- ncvar_get(nc, "in_situ_data/in_situ_latitude")
islon <- ncvar_get(nc, "in_situ_data/in_situ_longitude")
chla <- ncvar_get(nc, "in_situ_data/in_situ_chl")
id <- ncvar_get(nc, "in_situ_data/in_situ_id")
depth <- ncvar_get(nc, "in_situ_data/in_situ_depth")
datetime <- ncvar_get(nc, "in_situ_data/in_situ_datetime")
nctime <- ncvar_get(nc, "processing_data/nctime")
sat_chla <- ncvar_get(nc, "processing_data/chlor_a")
nclat <- ncvar_get(nc, "processing_data/latitude")
nclon <- ncvar_get(nc, "processing_data/longitude")
nc_close(nc)


# REDUCE VARIABLES IN 5X5 BOX TO SINGLE POINT
window_size <- dim(nclat)[1]
# (need in situ and satellite lats/lons, and geodist to get accurate distance)
# NOTE: USING THE CENTRE PIXEL LAT/LON, even though the median of the 3x3 matrix is used for Rrs
nclat <- apply(nclat, MARGIN=2, FUN="[[", ...=ceiling(window_size/2))
nclon <- apply(nclon, MARGIN=2, FUN="[[", ...=ceiling(window_size/2))
final_rrs <- lapply(1:length(rrs), FUN=function(i) apply(rrs[[i]], MARGIN=2, FUN=get_var))
final_rrs <- do.call(cbind, final_rrs)
sat_chla <- apply(sat_chla, MARGIN=2, FUN=get_var)


# MERGE ALL VARIABLES INTO A SINGLE DATAFRAME
df <- data.frame(cbind(final_rrs, chla, id, depth, datetime, islon, islat, nclon, nclat, nctime, sat_chla), stringsAsFactors = FALSE)
colnames(df) <- c(paste0("Rrs_", lambda), "chla", "id", "depth", "datetime", "islon", "islat", "nclon", "nclat", "nctime", "sat_chla")


#*******************************************************************************
# ADD MORE COLUMNS, FORMAT, SPLIT INTO TRAINING AND TEST SETS

# CALCULATE DISTANCE, AND USE ONLY CLOSEST MATCHUP
df <- df %>%
    dplyr::mutate(dist = as.numeric(geodist(x=df %>% dplyr::select(nclon, nclat),
                                            y=df %>% dplyr::select(islon, islat),
                                            paired = TRUE,
                                            measure="geodesic"))) %>%
    dplyr::arrange(., dist) %>%
    dplyr::group_by(., id) %>%
    dplyr::summarize_all(., .funs = dplyr::first)


# CALCULATE TIME DIFFERENCE BETWEEN IN SITU SAMPLE AND SATELLITE PASS
df <- df %>%
    dplyr::mutate(t1 = as_datetime(as.numeric(nctime)),
                  t2 = as_datetime(datetime, tz="UTC")) %>%
    dplyr::mutate(timediff = abs(sapply(1:nrow(df), function(i) time_length(interval(t1[i], t2[i]), unit="hour"))))


# MAKE SURE THEY'RE IN THE RIGHT FORMAT
cols <- grepl("Rrs_[0-9]{3}$", colnames(df)) | colnames(df) %in% c("chla", "depth", "dist", "timediff", "sat_chla")
df[,cols] <- sapply(df[,cols], as.numeric)
df$datetime <- as.character(df$datetime)
df$id <- as.integer(df$id)
df <- df %>% as.data.frame()


# REMOVE DATA WITH INVALID RRS OR CHLA
rrs_chl_inds <- grepl("Rrs_[0-9]{3}$", colnames(df)) | colnames(df)=="chla"
all_finite <- apply(is.finite(as.matrix(df[,rrs_chl_inds])), MARGIN=1, sum, na.rm=TRUE)==sum(rrs_chl_inds)
df <- df[all_finite,]


# CREATE TRAINING INDEX BY DEPTH, DISTANCE, AND TIME DIFFERENCE
training_ind <- df$depth <= max_depth & df$dist <= max_dist & df$timediff <= max_timediff


# GET SATELLITE CHL FOR OCX COMPARISON
ocx_chl <- df$sat_chla[!training_ind]


# REMOVE UNNECESSARY COLUMNS (keep separate metadata df for output csv)
metadata_df <- df %>% dplyr::select(id, depth, datetime)
df <- df %>% dplyr::select(-id, -depth, -islon, -islat, -nclon, -nclat, -dist, -datetime, -nctime, -t1, -t2, -timediff, -sat_chla)


# APPLY INDEX for use in the EOF training set (the rest should be used in the test set)
training_df <- df[training_ind,]
test_df <- df[!training_ind,]



#*******************************************************************************
# WRITE TO OUTPUT

write.csv(dplyr::bind_cols(metadata_df[training_ind,], training_df),
          file=paste0("03c_EOF_training_set_", region, "_", sensor, "_", paste0(range(years), collapse="-"), ".csv"),
          row.names=FALSE,
          quote=FALSE)



stop()
#*******************************************************************************
# COMPARE ALGORITHMS BELOW



# EOF CHL
test_eof_chl <- eof_chl(rrs=test_df, training_set=training_df)
logy <- log10(test_eof_chl)
logy[!is.finite(logy)] <- NA
logx <- log10(test_df$chla)
logx[!is.finite(logx)] <- NA
eof_lm <- lm(logy ~ logx)
eof_stats <- round(get_lm_stats(eof_lm), 3)
plot_eof <- ggplot(data.frame(in_situ_chl = test_df$chla, satellite_chl = test_eof_chl, stringsAsFactors = FALSE)) +
    geom_point(aes(x=in_situ_chl, y=satellite_chl), alpha=0.5) +
    geom_abline(intercept=0, slope=1) +
    geom_abline(intercept=eof_stats$Intercept, slope=eof_stats$Slope, color="red", linetype="dotted") +
    theme_bw() +
    ggtitle("EOF") +
    scale_x_log10(limits=c(0.1,30)) +
    scale_y_log10(limits=c(0.1,30)) +
    annotation_custom(tableGrob(t(eof_stats), cols = NULL, theme=ttheme_minimal(base_size=10, padding=unit(c(1,1), "mm"))),
                      xmin=-Inf, xmax=0.01, ymin=0.8, ymax=Inf)



# POLY4 CHL
lambdas <- get_lambda(ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)), use_443nm = FALSE)
best_alg_coefs <- get_coefs(ifelse(sensor=="VIIRS-SNPP", "viirs", tolower(sensor)), "nwa", "poly4")
poly4_chl <- ocx(rrs=as.matrix(test_df), blues=lambdas$blues, green=lambdas$green, coefs=best_alg_coefs)
logy <- log10(poly4_chl)
logy[!is.finite(logy)] <- NA
logx <- log10(test_df$chla)
logx[!is.finite(logx)] <- NA
poly4_lm <- lm(logy ~ logx)
poly4_stats <- round(get_lm_stats(poly4_lm), 3)
plot_poly4 <- ggplot(data.frame(in_situ_chl = test_df$chla, satellite_chl = poly4_chl, stringsAsFactors = FALSE)) +
    geom_point(aes(x=in_situ_chl, y=satellite_chl), alpha=0.5) +
    geom_abline(intercept=0, slope=1) +
    geom_abline(intercept=poly4_stats$Intercept, slope=poly4_stats$Slope, color="red", linetype="dotted") +
    theme_bw() +
    ggtitle("POLY4") +
    scale_x_log10(limits=c(0.1,30)) +
    scale_y_log10(limits=c(0.1,30)) +
    annotation_custom(tableGrob(t(poly4_stats), cols = NULL, theme=ttheme_minimal(base_size=10, padding=unit(c(1,1), "mm"))),
                      xmin=-Inf, xmax=0.01, ymin=0.8, ymax=Inf)




# OCX (from satellite)
logy <- log10(ocx_chl)
logy[!is.finite(logy)] <- NA
logx <- log10(test_df$chla)
logx[!is.finite(logx)] <- NA
ocx_lm <- lm(logy ~ logx)
ocx_stats <- round(get_lm_stats(ocx_lm), 3)
plot_ocx <- ggplot(data.frame(in_situ_chl = test_df$chla, satellite_chl = ocx_chl, stringsAsFactors = FALSE)) +
    geom_point(aes(x=in_situ_chl, y=satellite_chl), alpha=0.5) +
    geom_abline(intercept=0, slope=1) +
    geom_abline(intercept=ocx_stats$Intercept, slope=ocx_stats$Slope, color="red", linetype="dotted") +
    theme_bw() +
    ggtitle("OCx") +
    scale_x_log10(limits=c(0.1,30)) +
    scale_y_log10(limits=c(0.1,30)) +
    annotation_custom(tableGrob(t(ocx_stats), cols = NULL, theme=ttheme_minimal(base_size=10, padding=unit(c(1,1), "mm"))),
                      xmin=-Inf, xmax=0.01, ymin=0.8, ymax=Inf)



# GSM_GS
library(parallel)
library(pbapply)
library(compiler)
rrs <- as.matrix(test_df)[,grepl("Rrs_[0-9]{3}$", colnames(df))]
exp_file <- "03b_gsm_exponents_2019.csv"
gtype <- "gs"
num_cl <- 10
gsm_cmp <- cmpfun(gsm)
exps <- read.csv(exp_file, stringsAsFactors = FALSE)
exps <- as.numeric(unlist(exps[exps$region=="NWA" & exps$sensor==sensor & exps$gtype==gtype,4:6]))
chl_exp <- exps[1]
adg_exp <- exps[2]
bbp_exp <- exps[3]
wvs <- all_lambda[[sensor]]
rrs <- rrs/(0.52 + 1.7*rrs)
num_cl <- min(num_cl, detectCores()-1)

cl <- makeCluster(num_cl)
clusterExport(cl, c('gsm_cmp','gsm_model','rrs','wvs','chl_exp','adg_exp','bbp_exp','gtype'))
iops <- pbapply(X=rrs,
                MARGIN=1,
                FUN=gsm_cmp,
                lambda=wvs,
                adg_exp=adg_exp,
                bbp_exp=bbp_exp,
                chl_exp=chl_exp,
                gtype=gtype,
                cl=num_cl)
stopCluster(cl)
iops <- t(iops)
gsm_gs_chl <- as.numeric(iops[,1])

logy <- log10(gsm_gs_chl)
logy[!is.finite(logy)] <- NA
logx <- log10(test_df$chla)
logx[!is.finite(logx)] <- NA
gsm_gs_lm <- lm(logy ~ logx)
gsm_gs_stats <- round(get_lm_stats(gsm_gs_lm), 3)
plot_gsm_gs <- ggplot(data.frame(in_situ_chl = test_df$chla, satellite_chl = gsm_gs_chl, stringsAsFactors = FALSE)) +
    geom_point(aes(x=in_situ_chl, y=satellite_chl), alpha=0.5) +
    geom_abline(intercept=0, slope=1) +
    geom_abline(intercept=gsm_gs_stats$Intercept, slope=gsm_gs_stats$Slope, color="red", linetype="dotted") +
    theme_bw() +
    ggtitle("GSM_GS") +
    scale_x_log10(limits=c(0.1,30)) +
    scale_y_log10(limits=c(0.1,30)) +
    annotation_custom(tableGrob(t(gsm_gs_stats), cols = NULL, theme=ttheme_minimal(base_size=10, padding=unit(c(1,1), "mm"))),
                      xmin=-Inf, xmax=0.01, ymin=0.8, ymax=Inf)




# png(filename=paste0(region, "_", sensor, "_", paste0(range(years), collapse="-"), "_matchups.png"), width=1000, height=1000)
plot_ocx + plot_poly4 + plot_gsm_gs + plot_eof
# dev.off()

