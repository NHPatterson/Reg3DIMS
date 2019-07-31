#' @export
#' @title Register serial 3D IMS data
#' @description This function registers a sequence of 2D images saved as a list of 2D matrices into a 3D reconstruction, returning a list of transformations
#' @param reg_images list of 2D matrices representing template images
#' @param mid_point depth from which registration should proceed
#' @param flip_position optional, depth at which the image is mirror, length must match length of flip_type
#' @param flip_type optional, 1 (horizontral) or 2(vertical), length must match length of flip_position
#' @return named list of NiftyReg transformations per 2D section

reg3DIMS <- function(reg_images, mid_point = 5,flip_positions = c(), flip_type = c()){
  
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
    
    
    rigid <- niftyreg(reg_images[[source_im]],
                      reg_images[[target_im]],
                      scope=c("rigid"),
                      nLevels = 5L,maxIterations = 10,
                      verbose=F, interpolation=0)
    
    tforms[[source_im]]$nifty <- niftyreg(reg_images[[source_im]],
                                  reg_images[[target_im]],
                                  init=forward(rigid),
                                  scope=c("affine"),nLevels = 3,
                                  maxIterations = 10,verbose=F,
                                  interpolation=0)
    
    tform_im_mat <- tforms[[source_im]]$nifty$image[1:nrow(tforms[[source_im]]$nifty$image),
                                         1:ncol(tforms[[source_im]]$nifty$image)]
    
    reg_images[[source_im]] <- tform_im_mat
    
    bg <- EBImage::normalize(reg_images[[target_im]])
    
    overlay <- paintObjects(bg, 
                            toRGB(EBImage::normalize(tform_im_mat)), 
                            col=c("green"), opac=c(0.3), thick=F)
    
    display(overlay,method='raster')
    
  }
  
  tforms[[im_seq[mid_point]]]$init <- 'mid_point'
  
  return(tforms) 
}
  
  # 
  # ri_backward <- list()
  # ri_forward <- list()
  # 
  # # if(mid_point %in% flip_positions){
  # #   ft <- flip_df$flip_type[which(flip_df$flip_pos == mid_point)]
  # #   
  # #   if(ft == 1){
  # #     aff <- buildAffine(translation = c(0, 0, 0), scales = c(1, -1, 1),
  # #                        skews = c(0, 0,0), angles = 0,
  # #                        source = reg_images[[position-1]],
  # #                        target = reg_images[[position]],
  # #                        anchor =  "centre")
  # #   }else{
  # #     aff <- buildAffine(translation = c(0, 0, 0), scales = c(-1, 1, 1),
  # #                        skews = c(0, 0,0), angles = 0,
  # #                        source = reg_images[[position-1]],
  # #                        target = reg_images[[position]],
  # #                        anchor =  "centre")
  # #   }
  # #   reg_images[[mid_point]] <-applyTransform(aff,reg_images[[mid_point]] ,0)
  # #   
  # #   ri_backward[[mid_point]] <- T
  # #   ri_backward[[mid_point]]$init_tform <- aff
  # #   
  # # }
  # 
  # for(position in backward){
  # 
  #   if(position == max(backward)){
  # 
  #     cat('target: ', position, '\n')
  #     cat('source: ', position - 1, '\n')
  #     
  #     ri_backward[[position]] <- reg_images[[mid_point]]
  #     rigid <- niftyreg(reg_images[[position-1]],
  #                       reg_images[[position]],
  #                       scope=c("rigid"),
  #                       nLevels = 5L,maxIterations = 100,
  #                       verbose=F, interpolation=0)
  #     
  #     ri_backward[[position - 1]] <- niftyreg(reg_images[[position-1]],
  #                                             reg_images[[position]],
  #                                             init=forward(rigid),
  #                                             scope=c("affine"),nLevels = 3,
  #                                             maxIterations = 20,verbose=F,
  #                                             interpolation=0)
  #     
  #   } else {
  #     if(position %in% flip_positions){
  #       
  #       ft <- flip_df$flip_type[which(flip_df$flip_pos == position)]
  #       
  #       if(ft == 1){
  #         aff <- buildAffine(translation = c(0, 0, 0), scales = c(1, -1, 1),
  #                            skews = c(0, 0,0), angles = 0,
  #                            source = reg_images[[position+1]],
  #                            target = reg_images[[position]],
  #                            anchor =  "centre")
  #       }else{
  #         aff <- buildAffine(translation = c(0, 0, 0), scales = c(-1, 1, 1),
  #                            skews = c(0, 0,0), angles = 0,
  #                            source = reg_images[[position+1]],
  #                            target = reg_images[[position]],
  #                            anchor =  "centre")
  #       }
  #       
  #       
  #       
  #       out <-applyTransform(aff,reg_images[[position+1]] ,0)
  #       
  #       rigid <- niftyreg(out, ri_backward[[position]]$image, scope=c("rigid"),nLevels = 5L,maxIterations = 100,verbose=F,interpolation=0)
  #       
  #       ri_backward[[position + 1]] <- niftyreg(out, ri_backward[[position]]$image, scope=c("affine"),init=forward(rigid),nLevels = 3L,maxIterations = 20,verbose=F,interpolation=0)
  #       
  #       ri_backward[[position + 1]]$init_tform <- aff
  #       
  #       
  #     }else{
  #       cat('target: ', position, '\n')
  #       cat('source: ', position - 1, '\n')
  #       
  #       rigid <- niftyreg(reg_images[[position-1]],
  #                         ri_backward[[position]]$image,
  #                         scope=c("rigid"),nLevels = 5L,
  #                         maxIterations = 100,verbose=F,
  #                         interpolation=0)
  #       
  #       ri_backward[[position - 1]] <- niftyreg(reg_images[[position-1]],
  #                                               ri_backward[[position]]$image,
  #                                               init=forward(rigid),
  #                                               scope=c("affine"),
  #                                               nLevels = 3L, maxIterations = 20,
  #                                               verbose=F, interpolation=0)
  #     }
  #     
  #   }
  #   
  #   
  # }
  # 
  # for(position in forward){
  #   
  #   if(position == min(forward)){
  #     cat('target: ', position, '\n')
  #     cat('source: ', position + 1, '\n')
  #     
  #     
  #     
  #     ri_forward[[position]] <- reg_images[[mid_point]]
  #     
  #     rigid <- niftyreg(reg_images[[position+1]], reg_images[[position]] ,scope=c("rigid"),nLevels = 5L,maxIterations = 100,verbose=F,interpolation=0)
  #     
  #     ri_forward[[position + 1]] <- niftyreg(reg_images[[position+1]], reg_images[[position]] ,scope=c("affine"),init=forward(rigid),nLevels = 3L,maxIterations = 20,verbose=F,interpolation=0)
  #   } else {
  #     
  #     if(position %in% flip_positions){
  #       
  #       ft <- flip_df$flip_type[which(flip_df$flip_pos == position)]
  #       
  #       if(ft == 1){
  #         
  #         aff <- buildAffine(translation = c(0, 0, 0), scales = c(1, -1, 1), 
  #                            skews = c(0, 0,0), angles = 0, 
  #                            source = reg_images[[position+1]],
  #                            target = reg_images[[position]], 
  #                            anchor =  "centre")
  #       }else{
  #         
  #         aff <- buildAffine(translation = c(0, 0, 0), scales = c(-1, 1, 1), 
  #                            skews = c(0, 0,0), angles = 0, 
  #                            source = reg_images[[position+1]],
  #                            target = reg_images[[position]], 
  #                            anchor =  "centre")
  #       }
  #       
  #       
  #       out <-applyTransform(aff,reg_images[[position+1]] ,0)
  #       rigid <- niftyreg(out, ri_forward[[position]]$image, scope=c("rigid"),nLevels = 5L,maxIterations = 100,verbose=F,interpolation=0)
  #       
  #       ri_forward[[position + 1]] <- niftyreg(out, ri_forward[[position]]$image, scope=c("affine"),init=forward(rigid),nLevels = 3L,maxIterations = 20,verbose=F,interpolation=0)
  #       
  #       ri_forward[[position + 1]]$init_tform <- aff
  #       
  #       
  #       
  #     }else{
  #       cat('target: ', position, '\n')
  #       cat('source: ', position + 1, '\n')
  #       
  #       rigid <- niftyreg(reg_images[[position+1]], ri_forward[[position]]$image, scope=c("rigid"),nLevels = 5L,maxIterations = 100,verbose=F,interpolation=0)
  #       
  #       ri_forward[[position + 1]] <- niftyreg(reg_images[[position+1]], ri_forward[[position]]$image, scope=c("affine"),init=forward(rigid),nLevels = 3L,maxIterations = 20,verbose=F,interpolation=0)
  #     }
  #     
  #     
  #   }
  #   
  #   
  # }
  # tform_list <- c(ri_backward[1:mid_point], ri_forward[mid_point+1:length(ri_forward)])
  # 
  # tform_list[sapply(tform_list, is.null)] <- NULL
  # 
  # tform_list[mid_point] <- 'Untransformed'
  # 
  # names(tform_list) <- levels(pd_all$sample)
  # 
  # return(tform_list)
  #}

