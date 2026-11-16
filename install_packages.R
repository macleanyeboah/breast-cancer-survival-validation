# Install the packages this project needs (versions used: see README).
pkgs <- c("survival", "rstpm2", "randomForestSRC", "ggplot2", "scales", "jsonlite")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing)
