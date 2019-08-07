require(Cardinal)
require(RNiftyReg)
require(EBImage)
require(Reg3DIMS)


save(MuBr_seq2DIMS, file='./data/MuBr_seq2DIMS.RData')

##This procedure requires that 'run' desgination within the MSImageExperiment
#be ordered alphanumerically in sequence along the z-axis
#in the case of the test data, the 'run' indicates the slide# and sec# on the slide in sequence
#care should be taken here to use leading zeroes to avoid sorting issues when there are more than 10
#i.e. sorting will be 1,10,2 rather than 1,2,3... if there aren't leader zeroes (01,02,03)

#do PCA to get template imge
MuBr_seq2DIMS.PCA <- PCA(MuBr_seq2DIMS, ncomp=4)

#generate a list of PCA images from Principal Component 1 for registration
#reg_images <- get3DTemplateImgs(MuBr_seq2DIMS.PCA, canvas_xy = c(700,700), column = 'PC1')

#generate a list of m/z images for registration
reg_images <- get3DTemplateImgs(MuBr_seq2DIMS, canvas_xy = c(700,700), column = 807)

##develop registrations
#This procedure allows maximum flexibility for the registration scheme by allowing multiple registration degrees of freedom (rigid, affine, nonlinear) to be run in sequence.
#the reg_seq argument takes registration scope names: rigid (rotation & translation), affine(scaling,rotation,translation, shearing), non-linear (localized deformation).
#we can chain these schemes together and this is especially relevant for nonlinear transformations where the images should be initially well aligned before localized warping.
#after selecting a registration chain, scope-specific arguments are passed as a list with the scope name in the reg_seq_args argument. see ?niftyreg.linear and ?niftyreg.nonlinear for these arguments. We provide good default arguments in the examples below but these may need to be tuned for the inidividual application.

tforms_rig_aff <- reg3DIMS(reg_images,reg_seq = c('rigid','affine'), 
                          mid_point = 13, 
                          flip_positions = c(16), flip_type = c(1),
                          reg_seq_args = list(
                            rigid = list(nLevels=5, maxIterations = 10, interpolation = 0),
                            affine = list(nLevels=5, maxIterations = 10, interpolation = 0))) 

plotSeq3Dims(tforms_rig_aff)

tforms_rig_nl <- reg3DIMS(reg_images,reg_seq = c('rigid','nonlinear'), 
                   mid_point = 13, 
                   flip_positions = c(16), flip_type = c(1),
                   reg_seq_args = list(
                     rigid = list(nLevels=5, maxIterations = 10, interpolation = 0),
                     nonlinear = list(nBins = 128, finalSpacing = c(60,60,60), spacingUnit=c('voxel'),interpolation = 0))) 

plotSeq3Dims(tforms_rig_nl)

#apply transformations across the whole data set
MuBr_seq2DIMS.rigAff <- transform3DIMS(MuBr_seq2DIMS, tforms_rig_aff)
MuBr_seq2DIMS.rigNl  <- transform3DIMS(MuBr_seq2DIMS, tforms_rig_nl)

#merge 3D data, adding z dimension to dataset
MuBr_3DIMS.rigAff <- merge3DIMS(MuBr_seq2DIMS.rigAff)
MuBr_3DIMS.rigNl  <- merge3DIMS(MuBr_seq2DIMS.rigNl)

#export data to nifti for visualization outside of R
#visualization can be done inside of R but 
#other software offers better options for 3D scientific images
setwd('D:/temp/3d mouse brain')
exportNifti(MuBr_3DIMS.rigAff, xyz_um = c(20,20,20), name_prefix = 'mubr_3d_test_rig_aff')
exportNifti(MuBr_3DIMS.rigNl, xyz_um = c(20,20,20), name_prefix = 'mubr_3d_test_rig_nl')
