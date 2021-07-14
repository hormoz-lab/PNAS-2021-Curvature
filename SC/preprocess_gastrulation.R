preprocess_gastrulation <- function(data_path) {
  
  dat <- readRDS(paste(data_path, 'corrected_pcas.rds', sep='/'))
  write.csv(dat$all[,order(apply(dat$all, 2, var), decreasing=TRUE)], paste(data_path, 'PCA.csv', sep='/'), row.names=T)
  
}
