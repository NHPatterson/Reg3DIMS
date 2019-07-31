#' @export
#' @title Get 3D template images for registration
#' @description Generate a named list  of template images based on Cardinal IMS or analysis data (m/z images, PCA, Spatial Shrunken Centroids clustering)
#' @param ims Cardinal data of class: MSImagingExperiment, PCA2, SpatialShrunkenCentroids2
#' @param canvas_xy standardized canvas size for all images, must accomodate all IMS pixels
#' @param column name of column to extract if character, or m/z value if numeric 
#' @return named list of images

get3DTemplateImgs <- function(ims, canvas_xy= c(700,700), column='PC1'){
  
  
  if(class(ims)[1] == "SpatialShrunkenCentroids2"){
    if(column != 'all'){
      if(!column %in% colnames(ims@resultData[[1]]$probabilities)){
        stop(paste0('The specified column (',column,') for data selection is not in the dataset (',paste(colnames(ims@resultData[[1]]$probabilities),collapse=' '),')'))
      }
    }
    
    resultdf <- data.frame(sample = run(ims),  
                           x=  coord(ims)$x, 
                           y = coord(ims)$y, 
                           resultData(ims)[[1]]$probabilities)
  }
  
  
  if(class(ims)[1] == "PCA2"){
    if(!column %in% colnames(ims@resultData[[1]]$scores) ){
      stop(paste0('The specified column (',column,') for data selection is not in the dataset (',paste(colnames(ims@resultData[[1]]$scores),collapse=' '),')'))
    }
    
    resultdf <- data.frame(sample = run(ims), 
                           x=  coord(ims)$x, 
                           y = coord(ims)$y, 
                           resultData(ims)[[1]]$scores[,column])
    
    colnames(resultdf)[4] <- column
  }
  
  if(class(ims)[1] == "MSImagingExperiment"){
    if(is.numeric(column) == F){
      stop(paste0('The specified column (',column,') for data selection is not numeric and not an m/z'))
    }
    resultdf <- data.frame(sample = run(ims), 
                           x = coord(ims)$x, 
                           y = coord(ims)$y,
                           mz = iData(ims)[features(ims, mz=column),])
  }
  
  
  dfs <- split(resultdf[c('x','y')], resultdf$sample)
  
  dfs <- lapply(dfs,function(x){
    x$x <- x$x - (min(x$x) - 1)
    x$y <- x$y - (min(x$y) - 1)
    return(x)
  })
  
  dfs <- do.call(rbind, dfs)
  
  sample_xy_min <- aggregate(dfs[c('x','y')], by=list(resultdf$sample), FUN=min)
  sample_xy_max <- aggregate(dfs[c('x','y')], by=list(resultdf$sample), FUN=max)
  
  xyminmax <- data.frame(x_min = sample_xy_min$x,
             y_min = sample_xy_min$y,
             x_max = sample_xy_max$x,
             y_max = sample_xy_max$y)
  
  if(max(sample_xy_max$x) > canvas_xy[1]){
    stop(paste0('specified canvas x size (',canvas_xy[1],') is too small to fit all IMS data extents (',max(sample_xy_max$x),')'))
  }
  
  if(max(sample_xy_max$y) > canvas_xy[2]){
    stop(paste0('specified canvas y size (',canvas_xy[2],') is too small to fit all IMS data extents (',max(sample_xy_max$y),')'))
  }
  
  if(class(ims)[1] == "MSImagingExperiment"){
    if(all(dim(imageData(ims)[1,1,,]) ==
         c(max(xyminmax$x_max),max(xyminmax$y_max)))){
    
    image_matrix <- imageData(ims)[features(ims, mz=column),,,]
    # image_matrix[!is.na(image_matrix)] <- 1
    image_matrix[is.na(image_matrix)] <- 0

    reg_images <- lapply(seq(dim(image_matrix)[1]), function(x) {
      
      im_xmax <- which(rowSums(image_matrix[1,,]) == 0) - 1
      im_ymax <- which(colSums(image_matrix[1,,]) == 0) - 1
      
      immat <- matrix(0, nrow=canvas_xy[1] , ncol=canvas_xy[2])

      x_padding = floor((canvas_xy[1] - xyminmax$x_max[x]) / 2) + 1
      y_padding = floor((canvas_xy[2] - xyminmax$y_max[x]) / 2) + 1
      
      x_width <- x_padding + xyminmax$x_max[x] - 1
      y_width <- y_padding + xyminmax$y_max[x] - 1

      immat[x_padding:x_width,
            y_padding:y_width] <- image_matrix[x,
                                               xyminmax$x_min[x]:xyminmax$x_max[x],
                                               xyminmax$y_min[x]:xyminmax$y_max[x]]

      return(immat)
    })
    return(reg_images)
    }
  }
  
  resultdf$x <- dfs$x
  resultdf$y <- dfs$y
  reg_images <- list()
  counter <- 0
  for(sample in levels(factor(resultdf$sample))){
    
    cat('Collecting resultSet image from : ', sample,
        '   Progress : ', round((counter / nlevels(factor(resultdf$sample)))*100,2),'%','\r')
    
    
    image_matrix <- matrix(0, nrow=canvas_xy[1] , ncol=canvas_xy[2])
    
    sampledf <- resultdf[resultdf$sample %in% sample,]
    x_padding = floor((canvas_xy[1] - max(sampledf$x)) / 2)
    y_padding = floor((canvas_xy[2] - max(sampledf$y)) / 2)
    sampledf$x <- sampledf$x + x_padding
    sampledf$y <- sampledf$y + y_padding
    
    for(i in 1:nrow(sampledf)){
      if(is.numeric(column)){
        
        image_matrix[sampledf[i,'x'],sampledf[i,'y']] <- sampledf[i,'mz']
        
      }else{
        
        image_matrix[sampledf[i,'x'],sampledf[i,'y']] <- sampledf[i,column]
        
      }
    }
    
    reg_images[[sample]] <- image_matrix
    
    counter <- counter + 1
    
  }
  
  return(reg_images)
}

