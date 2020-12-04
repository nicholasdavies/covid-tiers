library(data.table)
library(ggplot2)
library(lubridate)
library(here)
library(cowplot)
library(readxl)
library(sn)
library(qs)

# Sys.setenv(MP_DUPLICATE_LIB_OK=TRUE) # RCB needs to override this environment variable to run

#
# SETUP
#

# set up covidm
cm_path = "~/Documents/covidm_MTPs/covidm_for_fitting/";
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 2;
source(paste0(cm_path, "/R/covidm.R")) # RCB needs to move (backed up) Makevars file into .R repo for this to work with -fopenmp
source("~/Documents/covidm_MTPs/locate.R");
cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)
source("~/Documents/covidm_MTPs/distribution_fit.R");
source("~/Documents/covidm_MTPs/spim_output.R");



#
# DATA
#

nhs_regions = popUK[, unique(name)]
all_data = qread("~/Dropbox/uk_covid_data/processed-data-2020-10-23.qs")
# for more recent data:
# all_data = qread("~/Dropbox/uk_covid_data/processed-data-2020-11-27.qs")
ld = all_data[[1]]
sitreps = all_data[[2]]
virus = all_data[[3]]
sero = all_data[[4]]