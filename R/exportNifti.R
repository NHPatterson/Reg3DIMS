#' @export
#' @title Export data in Nifti format
#' @description Exports reconstructed 3D arrays which can be opened and visualized in FIJI ImageJ, Slicer or other scientific imaging software. Function writes *.nii files to the current working directory.
#' @param ims3d combined, processed 3D IMS dataset
#' @param xyz_um spatial resolution information in microns for all three dimensions
#' @param name_prefix character string to be prepended to exported nifty image filenames

exportNifti <- function(ims3d, xyz_um = c(20,20,20), name_prefix = 'my3d_data'){
  for(mz in Cardinal::mz(ims3d)){
    outarray <- slice(ims3d, mz=mz)
    pixdim(outarray) <- xyz_um
    pixunits(outarray) <- 'um'
    writeNifti(outarray, file=paste0(name_prefix,'_mz_',round(mz,3),'.nii'))
  }
}
