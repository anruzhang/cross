
library(rTensor)
library(matrixcalc)
library(MASS)

mynorm = function(x, type){if(type=="2"){spectral.norm(x)}else{norm(x,type)}}
tcross <- function(B, Arm, Indices, c_t){
  # Cross: efficient tensor completion
  
  # B: body measurements; tensor class input; required
  # arm: list of d matrices representing arm measurements
  # Indices: list of d sets of indices for locations of observed core tensor, defaut value: the top {1, ..., r} 
  
  try(if(missing("B")) stop("missing argument: B is required as tensor type."))
  try(if(missing("Arm")) stop("invalid input: Arm is required as list of matrices."))
  try(if(missing("c_t")) c_t = 3)
  try(if(class(B) != "Tensor") stop("invalid input: S should be of tensor type."))
  try(if(! is.list(Arm)) stop("invalid input: Arm should be of a list of matrices."))
  m = dim(B)
  d = length(m)
  try(if(d != length(Arm)) stop("invalid input: order of S and arm are incompatible."))
  
  
  g = rep(0, d); p = rep(0, d) #Second step.
  for (t in 1:d){
    try(if(! is.matrix(Arm[[t]])) stop("invalid input: Arm should be of a list of matrices."))
    p[t] = dim(Arm[[t]])[1]; g[t] = dim(Arm[[t]])[2]
    try(if(p[t] < m[t]) stop("invalid input: core tensor and arm matrix's dimentions are incompatible."))
  }
  if(missing("Indices")) {
    Indices = list();
    for (t in 1:d){
      Indices = c(Indices, list(1:m[t]))
    }
  }
  
  hatR = list()  # Initialization for expanding matrix hatR[]
  hatr=rep(0, d)
  for (t in 1:d)
  {
    arm_matrix = Arm[[t]]
    body_matrix = k_unfold(B, t)@data
    U_body = svd(body_matrix)$u
    joint_matrix = arm_matrix[Indices[[t]],]
    V_arm = svd(arm_matrix)$v
    joint_matrix = t(U_body) %*% joint_matrix %*% V_arm
    arm_matrix = arm_matrix %*% V_arm
    
    for (s in min(g[t], m[t]):1){
      if(!is.singular.matrix(as.matrix(joint_matrix[1:s,1:s]))){
        if (mynorm(arm_matrix[, 1:s] %*% ginv(joint_matrix[1:s, 1:s]), "2") <= c_t*sqrt(p[t])/(sqrt(m[t]))){
          hatr[t] = s; break;
        }
      }
    }
    if (hatr[t] >0 ){
      hatRt = arm_matrix[, 1:hatr[t]] %*% ginv(joint_matrix[1:hatr[t], 1:hatr[t]]) %*% t(U_body[,1:hatr[t]])
    } else {
      hatRt = matrix(0, p[t], m[t]);
    }
    hatR = c(hatR, list(hatRt))
  }
  
  hatA = B;
  for (t in 1:d){
    hatA = ttm(hatA, hatR[[t]], t)
  }
  
  return(list(hatA, hatr, hatR))
}

tcross_cv <- function(B, Arm, Indices){ # cross validation without specifying c_t
  # Cross-validation c_t range: 1.5 - 3
  try(if(missing("B")) stop("missing argument: B is required as tensor type."))
  try(if(missing("Arm")) stop("invalid input: Arm is required as list of matrices."))
  try(if(class(B) != "Tensor") stop("invalid input: S should be of tensor type."))
  try(if(! is.list(Arm)) stop("invalid input: Arm should be of a list of matrices."))
  m = dim(B)
  d = length(m)
  try(if(d != length(Arm)) stop("invalid input: order of S and arm are incompatible."))
  if(missing("Indices")) {
    Indices = list();
    for (t in 1:d){
      Indices = c(Indices, list(1:m[t]))
    }
  }
  
  c_t_value = seq(2.2, 3.2, 0.2); # Candidate parameter sets
  cross_error = rep(0, length(c_t_value)) 
  m_train = ceiling(m*4/5) # 5-fold CV
  rmax = 5 # Repeat times for each tuning parameter value
  
  for (r in 1:rmax){
    Indices_train=list() # Indices set of training set
    for(i in 1:d){
      Indices_train = c(Indices_train, list(sort(sample(Indices[[i]], m_train[i]))))
    }
    for (l in 1:length(c_t_value)){
      c_t = c_t_value[l]
      B_train = B[Indices_train[[1]], Indices_train[[2]], Indices_train[[3]]]
      X_train = tcross(B_train, Arm, Indices_train, c_t); X_train = X_train[[1]]
      cross_error[l] = cross_error[l] + fnorm(B - X_train[Indices[[1]], Indices[[2]], Indices[[3]]])^2 -
        fnorm(B[Indices_train[[1]], Indices_train[[2]], Indices_train[[3]]] - X_train[Indices_train[[1]], Indices_train[[2]], Indices_train[[3]]])^2
    }
  }
  c_t_0 = c_t_value[which.min(cross_error)]
  return(c(tcross(B, Arm, c_t = c_t_0), c_t_0))
}  
### End of the Algorithm ###
