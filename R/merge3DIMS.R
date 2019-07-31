#' @export
#' @title Merge a list of 2D MSImagingExperiments to 3D
#' @description Merges MSImagingExperiments to 3D using the 'run' variable of each
#' @param ims3d list of processed 3D IMS dataset
#' @return 3D MSImagingExperiment with z-dimension added

merge3DIMS <- function(ims3d){
  ims3d <- combine(ims3d[1:length(ims3d)])
  zdim <- run(ims3d)
  levels(zdim) <- seq(0,nlevels(zdim), 1) #add z-dimension to dataset based on sequentially ordered runs
  coord(ims3d)$z <- as.numeric(as.character(zdim))
  return(ims3d)
}
