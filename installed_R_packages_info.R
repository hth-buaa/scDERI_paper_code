##########################################
#                                        #
#   R Package List for Reproducibility   #
#                                        #
##########################################


# List my R-Pack information for other researchers to reproduce my work
setwd("/mnt/i/scDERI_Paper/")
installed_packages_info = as.data.frame(installed.packages(), stringsAsFactors = FALSE)
# install.packages("writexl")
library(writexl)
write_xlsx(installed_packages_info, "installed_packages_info.xlsx")
