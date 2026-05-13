# Utility function to build the calculator input data from standard genomics data
build_inputs <- function(x) {
  res <- matrix(NA, nrow = nrow(x), ncol = 7)
  rownames(res) <- rownames(x)
  colnames(res) <- c("ipss_m","age","asxl1_kras","srsf2_nras","nras_runx1","atrx","jak2")
  for(i in 1:nrow(res)) {
    res[i,"ipss_m"] <- as.numeric(x[i,"ipss_m"])
    res[i,"age"] <- as.numeric(x[i,"age"])
    if(!is.na(x[i,"asxl1"])&&!is.na(x[i,"kras"])) {
      if(x[i,"asxl1"]==1&&x[i,"kras"]==1) {
        res[i,"asxl1_kras"] <- 1
      } else {
        res[i,"asxl1_kras"] <- 0
      }
    }
    if(!is.na(x[i,"srsf2"])&&!is.na(x[i,"nras"])) {
      if(x[i,"srsf2"]==1&&x[i,"nras"]==1) {
        res[i,"srsf2_nras"] <- 1
      } else {
        res[i,"srsf2_nras"] <- 0
      }
    }
    if(!is.na(x[i,"nras"])&&!is.na(x[i,"runx1"])) {
      if(x[i,"nras"]==1&&x[i,"runx1"]==1) {
        res[i,"nras_runx1"] <- 1
      } else {
        res[i,"nras_runx1"] <- 0
      }
    }
    res[i,"atrx"] <- as.numeric(x[i,"atrx"])
    res[i,"jak2"] <- as.numeric(x[i,"jak2"])
  }
  res <- as.data.frame(res)
  res
}
