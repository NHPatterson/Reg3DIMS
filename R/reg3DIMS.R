#' @export
#' @title Register serial 3D IMS data
#' @description This function registers a sequence of 2D images saved as a list of 2D matrices into a 3D reconstruction, returning a list of transformations
#' @param reg_images list of 2D matrices representing template images
#' @param mid_point depth from which registration should proceed
#' @param flip_position optional, depth at which the image is mirror, length must match length of flip_type
#' @param flip_type optional, 1 (horizontral) or 2(vertical), length must match length of flip_position
#' @param ... arguments passed to niftyreg.linear
#' @return named list of NiftyReg transformations per 2D section

reg3DIMS <- function(reg_images, reg_seq = c('rigid','nonlinear'), mid_point = 5,flip_positions = c(), flip_type = c(), reg_seq_args = list(
  rigid = list(nLevels=5, maxIterations = 10, interpolation = 0),
  nonlinear = list(nBins = 128,finalSpacing = c(20,20,20), spacingUnit=c('voxel'),interpolation = 0)),
  plot = T){
  
  if(length(flip_positions) != length(flip_type)){
    stop('flip_positions is not equal to the number of flip_types \n
         Each position that is to be flipped must have a flip type(1 (horizontal) or 2 (vertical))')
  }
  
  if(any(flip_type < 1) & any(flip_type > 2)){
    stop('flip_type must be specified as 1 or 2 \n
         1  = (horizontal) or 2 =  (vertical))')
  }
  
  
  im_seq <- names(reg_images)
  
  backward <- seq(mid_point - 1, 1, by= -1)
  
  forward  <-seq(mid_point + 1,length(reg_images),by=1)
  
  dfs <- list()
  for(idx in backward){
    if(idx == max(backward)){
      dfs[[idx]] <- data.frame(idx = idx, source_im_tformed = F, source_im = im_seq[idx], 
                               target_im_tformed = F, target_im = im_seq[idx+1], stringsAsFactors = F)
      
    }else{
      dfs[[idx]] <- data.frame(idx = idx, source_im_tformed = F, source_im = im_seq[idx], 
                               target_im_tformed = T, target_im = im_seq[idx+1], stringsAsFactors = F)
      
    }
  }
  
  tform_df_backward <- do.call(rbind, dfs)
  tform_df_backward <- tform_df_backward[order(-tform_df_backward$idx),]
  
  dfs <- list()
  for(idx in forward){
    if(idx == min(forward)){
      dfs[[idx]] <- data.frame(idx = idx, source_im_tformed = F, source_im = im_seq[idx], 
                               target_im_tformed = F, target_im = im_seq[idx-1], stringsAsFactors = F)
      
    }else{
      dfs[[idx]] <- data.frame(idx = idx, source_im_tformed = F, source_im = im_seq[idx], 
                               target_im_tformed = T, target_im = im_seq[idx-1], stringsAsFactors = F)
      
    }
  }
  
  tform_df_forward <- do.call(rbind, dfs)
  
  tform_df_full <- rbind(tform_df_backward, tform_df_forward)
  
  tform_df_full$flip_source <- F
  tform_df_full$flip_type <- NA
  
  flip_idx <- which(tform_df_full$source_im %in% im_seq[flip_positions])
  
  tform_df_full$flip_source[flip_idx] <- T
  
  tform_df_full$flip_type[flip_idx] <- flip_type
  
  tforms <- list()
  for(row in seq_len(nrow(tform_df_full))){
    source_im <- tform_df_full[row,'source_im']
    target_im <- tform_df_full[row,'target_im']
    
    if(tform_df_full[row,'flip_source']){
      
      if(tform_df_full[row,'flip_type'] == 1){
        aff <- buildAffine(translation = c(0, 0, 0), scales = c(1, -1, 1),
                           skews = c(0, 0,0), angles = 0,
                           source = reg_images[[source_im]],
                           target = reg_images[[target_im]],
                           anchor =  "centre")
        f_type <- 'horizontally'
      }else{
        aff <- buildAffine(translation = c(0, 0, 0), scales = c(-1, 1, 1),
                           skews = c(0, 0,0), angles = 0,
                           source = reg_images[[source_im]],
                           target = reg_images[[target_im]],
                           anchor =  "centre")
        f_type <- 'vertically'
        
      }
      
      reg_images[[source_im]] <-applyTransform(aff,reg_images[[source_im]] ,0)
      tforms[[source_im]]$init <- aff
      
      
      cat(source_im,' was flipped ', f_type,'\n')
      
    }
    cat('Registration of ', source_im, ' to ', target_im, '\n')
    #add control over registration sequencing
    for(scope in reg_seq){
      
      reg_seq_args[[scope]]$source<- reg_images[[source_im]]
      reg_seq_args[[scope]]$target<- reg_images[[target_im]]
      reg_seq_args[[scope]]$scope <- scope
      
      #initialize with previous transformations
      #weird issue here, but ignore for now
      if(length(reg_seq) > 1 & which(reg_seq == scope) > 1){
        if(scope %in% c('rigid','affine')){
          reg_seq_args[[scope]]$init <- forward(tform)
          
        }else if(scope == 'nonlinear'){
          reg_seq_args[[scope]]$init <- forward(tform)
          
        }
      }
      tform <- do.call(niftyreg, reg_seq_args[[scope]])
    }
    
    tforms[[source_im]]$nifty <- tform
    
    tform_im_mat <- tforms[[source_im]]$nifty$image[1:nrow(tforms[[source_im]]$nifty$image),
                                                    1:ncol(tforms[[source_im]]$nifty$image)]
    
    reg_images[[source_im]] <- tform_im_mat

    bg <- EBImage::normalize(reg_images[[target_im]])
    
    overlay <- rgbImage(red = bg, green = EBImage::normalize(tform_im_mat))
    
    if(plot == T){
      display(overlay, method='raster')
      
      text(x = 20, y = 20, label = 
             paste0(source_im, ' to ', target_im), adj = c(0,1), col = "white", cex = 1)
    }

    
    tforms[[source_im]]$overlay <- overlay
  }
  
  tforms[[im_seq[mid_point]]]$init <- 'mid_point'
  tforms <- tforms[order(names(tforms))]
  return(tforms) 
  }

