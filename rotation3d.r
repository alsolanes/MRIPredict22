rotation_matrix_3d <- function(degrees)
{
  a = degrees*pi/180
  R.x <- matrix( c(     1,       0,       0,  0,
                        0,  cos(a), -sin(a),  0,
                        0, sin(a),   cos(a),  0,
                        0,       0,       0,  1), 4)
  R.y <- matrix( c(cos(a),       0,  -sin(a), 0, 
                        0,       1,       0,  0,  
                   sin(a),       0,   cos(a), 0,
                        0,       0,       0,  1), 4)
  R.z <- matrix( c(cos(a), -sin(a),       0,  0,
                   sin(a),  cos(a),       0,  0,
                        0,       0,       1 , 0,
                        0,       0,       0,  1), 4)
  list(R.x=R.x, R.y=R.y, R.z=R.z)
}
translate_matrix_3d <- function(x,y,z)
{
  matrix(c(1, 0, 0, x,
           0, 1, 0, y,
           0, 0, 1, z,
           0, 0, 0, 1), 4)
}
rotate_coordinates_old <- function(degrees, axe=1, coordinates,translate_xyz)
{
  # 1st translate to origin (center=c(0,0,0)), 2nd rotate along an axe, 3rd translate back to original position
  # out = translate_matrix_3d(-translate_xyz[1], -translate_xyz[2], -translate_xyz[3]) %*% matrix(c(coordinates,1),ncol=1)
  # out = rotation_matrix_3d(degrees)[[axe]] %*% out 
  # translate_matrix_3d(translate_xyz[1], translate_xyz[2], translate_xyz[3]) %*% out
    translate_matrix_3d(translate_xyz[1], translate_xyz[2], translate_xyz[3]) %*%
          rotation_matrix_3d(degrees)[[axe]] %*%
          translate_matrix_3d(-translate_xyz[1], -translate_xyz[2], -translate_xyz[3]) %*%
          matrix(c(coordinates,1),ncol=1)
}
plot_brain = function(brain) {#, id_subject){
  .require('reshape2')
  .require('rgl')
  M = melt(brain>0.001)
  M = M[M$value==TRUE,]
  plot3d(M$Var1,M$Var2,M$Var3,size=10,transparency=0.5, xlab = 'X', ylab = 'Y', zlab = 'Z')
}
rotate_brain = function(brain, degrees, axe=3) {
  size_brain=dim(brain)
  rotated_brain = array(rep(0, prod(size_brain*3)), dim = size_brain * 3)
  for(i in 1:size_brain[1]) {
    for (j in 1:size_brain[2]) {
      for (k in 1:size_brain[3]) {
        if (brain[i,j,k] == TRUE) {
          rotated = rotate_coordinates(degrees, axe = axe, coordinates = c(i,j,k), translate_xyz = round(c(dim(brain)/2))) + size_brain*1.5
          rotated_brain[rotated[1],rotated[2],rotated[3]] = 1
        }
      }
    }
  }
  rotated_brain
}
rotate_brain_mc = function(brain, degrees, axe=3, n_cores = 2) {
  size_brain=dim(brain)
  rotated_brain = array(rep(0, prod(size_brain*3)), dim=size_brain*3)
  translate_xyz = round(c(dim(brain)/2))
  brain2=pbmclapply(expand.grid(i=1:size_brain[1],j=1:size_brain[2],k=1:size_brain[3]), function(brain) {
    coords = rotate_coordinates(degrees, axe = axe, coordinates = c(i,j,k), translate_xyz = translate_XYZ) + 1.5*size_brain
    rotated_brain[rotated[1],rotated[2],rotated[3]] = 1
    brain
  },brain,mc.cores = n_cores, max.vector.size = 22020096000
  )
}
split_brain = function(brain, degrees, axe=1) {
  # image corresponds to:
  # 1 from x axis rotation
  # 2 from y axis rotation
  # 3 from z axis rotation
  # will receive the mask
  # find the median zone (accounting for how many voxels are included in the corresponding plane)
  # will split the brain in two

  rotated_brain = rotate_brain(brain, degrees, axe)
  # define a translation to center the brain (to 0,0,0) and then rotate it

  k_cut = index_median(rotated_brain, axe = axe)
  switch(as.character(axe),
         "1"={
           rotated_brain[c(1:k_cut),,]=0
         },
         "2"={
           rotated_brain[,c(1:k_cut),]=0
         },
         "3"={
           rotated_brain[,,c(1:k_cut)]=0
         })
  list(rotated_brain = rotated_brain, z_cut = k_cut)
}
get_rotated_coordinates = function(degrees, axe, coordinates = c(i, j, k), brain) {
  size_brain = dim(brain)
  rotate_coordinates(degrees, axe = axe, coordinates = coordinates, translate_xyz = round(c(dim(brain)/2))) + 1.5*size_brain
}
index_median = function(brain, axe = 3) {
  z_median=unlist(lapply(1:dim(brain)[axe], function(i){ switch(as.character(axe),
                                                                "1"=sum(brain[i,,]),
                                                                "2"=sum(brain[,i,]),
                                                                "3"=sum(brain[,,i]))}
  ))
  median(which(z_median>0))
}
cut_brain = function(brain, degrees=0, axe=1){
  size_brain = dim(brain)
  z_cut = split_brain(brain, degrees, axe)$z_cut
  translate_XYZ = round(c(size_brain/2))
  for(i in 1:size_brain[1]) {
    for (j in 1:size_brain[2]) {
      for (k in 1:size_brain[3]) {
        coords = rotate_coordinates(degrees, axe = axe, coordinates = c(i,j,k), translate_xyz = translate_XYZ) + 1.5*size_brain
        if(coords[axe] >= z_cut) {
          brain[i,j,k] = FALSE
        }
      }
    }
  }
  list(brain = brain, z_cut = z_cut)
}
cut_brain_mc = function(brain, degrees=0, axe=1, n_cores=2) {
  size_brain = dim(brain)
  z_cut = split_brain(brain, degrees, axe)$z_cut
  translate_XYZ = round(c(size_brain/2))
  brain2=pbmclapply(expand.grid(i=1:size_brain[1],j=1:size_brain[2],k=1:size_brain[3]), function(brain) {
    coords = rotate_coordinates(degrees, axe = axe, coordinates = c(i,j,k), translate_xyz = translate_XYZ) + 1.5*size_brain
    if(coords[axe] >= z_cut) {
      brain[i,j,k] = FALSE
    }
    brain
  },brain,mc.cores = n_cores, max.vector.size = 22020096000
  )
}
plot_Brain2D = function(X, y_cut=round(dim(X)[2]/2)) {

  if (y_cut=="")
    image.default(X)
  else
    image.default(X[,y_cut,])
  # Close the file
  dev.off()
}
rotate_coordinates = function(list_coordinates, iter=1, return_3d=FALSE, img_3d=c()) {
  vectors = list(c(0,1,0),   # 0º eix X
                 c(0,1,1),   # 45º
                 c(0,0,1),   # 90º
                 c(0,-1,1),  #
                 c(0,-1,0),  # 180
                 c(0,-1,-1), #
                 c(0,0,-1),  # 270#
                 c(0,1,-1),   # 315
                 c(1,0,0),   # 0º eix Y
                 c(1,0,1),   # 45º
                 #c(0,0,1),   # 90º
                 c(-1,0,1),  #
                 c(-1,0,0),  # 180
                 c(-1,0,-1), #
                 #c(0,0,-1),  # 270#
                 c(1,0,-1),   # 315
                 #c(1,0,0),   # 0º eix Z
                 c(1,1,0),   # 45º
                 #c(0,1,0),   # 90º
                 c(-1,1,0),  #
                 #c(-1,0,0),  # 180
                 c(-1,-1,0), #
                 #c(0,-1,0),  # 270#
                 c(1,-1,0)   # 315
  )
  
  vect = vectors[[iter]]
  v1 = vect / sqrt(2)

  U = list_coordinates
  C = apply(U, 1, function(u) {sum(u * v1)})
  if (runif(1) > 0.5) {
    M = C > median(C)
  } else {
    M = C >= median(C)
  }
  if (return_3d){
    mask_ijk = U
    total = nrow(mask_ijk)
    img_3d = img_3d*0
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    for(i in seq_len(nrow(mask_ijk))) {
      if(M[i])
        img_3d[mask_ijk[i,1], mask_ijk[i,2], mask_ijk[i,3]] = TRUE
      else
        img_3d[mask_ijk[i,1], mask_ijk[i,2], mask_ijk[i,3]] = FALSE
      setTxtProgressBar(pb, i)
    } # masked mri
    close(pb)
    img_3d
  }
  else{
    cbind(U,M)
  }
}