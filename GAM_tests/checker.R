

element_name_R <- "R0.txt"
element_name_CPP <- "GLM_mixedferegression_R0_.txt"

element_name_R <- "R1.txt"
element_name_CPP <- "GLM_mixedferegression_R1_.txt"

element_name_R <- "solution_lapl_f.txt"
element_name_CPP <- "apply_return_laplacian_coefs_.txt"

element_name_R <- "solution_f.txt"
element_name_CPP <- "apply_return_solution_coefs_.txt"

element_name_R <- "AS.txt"
element_name_CPP <- "mixedferegression_AS_.txt"

element_name_R <- "Q.txt"
element_name_CPP <- "mixedferegression_Q_.txt"

element_name_R <- "WeightMatrix.txt"
element_name_CPP <- "pseudoData_weightsmatrix_0.txt"

element_name_R <- "b.txt"
element_name_CPP <- "GLM_mixedfereg_righthandside_.txt"

element_name_R <- "Z.txt"
element_name_CPP <- "z_6.txt"

element_name_R <- "mu.txt"
element_name_CPP <- "mu_0.txt"

element_name_R <- "beta_hat.txt"
element_name_CPP <- "betahat_0.txt"

element_name_CPP <- "regressionestimates_desmatprod.txt"
element_name_CPP_2 <- "regressionestimates_obs_meno_fnhat.txt"



setwd("/home/alb/Scrivania/PACS/Git_Folder/debugging_output/R/")
element_R <- read.table(element_name_R, header = FALSE, sep = ",", dec = ".", row.names = NULL)

setwd("/home/alb/Scrivania/PACS/Git_Folder/debugging_output/CPP/")
element_CPP <- read.table(element_name_CPP, header = FALSE, sep = ",", dec = ".", row.names = NULL)
#element_CPP_2 <- read.table(element_name_CPP_2, header = FALSE, sep = ",", dec = ".", row.names = NULL)


#matrice <- as.matrix(element_CPP)

#vector <- as.matrix(element_CPP_2)
#dim(matrice)
#dim(vector)
#result <- matrice %*% vector


dim(element_R)
dim(element_CPP)
#element_CPP <- element_CPP[,1:100]

#element_R <- data.frame(element_R[1:100,])

length(element_CPP[,1]) == length(element_R[,1])
row_dim = length(element_CPP[,1])

length(element_CPP[1,]) == length(element_R[1,])
col_dim = length(element_CPP[1,]) 

result = matrix(nrow = row_dim, ncol = col_dim)

at_least_one = 0

for(i in 1:row_dim){
  for(j in 1:col_dim){
    
    result[i,j] = (abs(element_CPP[i,j] - element_R[i,j]) > 10^-15)
    
  }
  at_least_one = at_least_one + sum(result[i,])
}

at_least_one

#element_CPP[1,1]
#element_R[1,1]


#element_CPP[1:10,1]
#element_R[1:10,1]
#result[1:10,1]

format(matrice, digits=16)
format(element_CPP[result],digits=16)
format(element_R[result],digits=16)

format((element_CPP[result]-element_R[result]),digits=16)





#### PLOT 3D ######

library(rgl)

element_name_R <- "beta_hat.txt"
element_name_desmat <- "covariates.txt"
element_name_response <- "link_response.txt"

element_name_CPP <- "betahat_0.txt"
element_name_CPP <- "compute_mu_betahat_0.txt"

element_name_CPP_3 <- "regressionestimates_obs_meno_fnhat.txt"


setwd("/home/alb/Scrivania/PACS/Git_Folder/debugging_output/R/")
  beta_hat_R <- read.table(element_name_R, header = FALSE, sep = ",", dec = ".", row.names = NULL)
covariates <- read.table(element_name_desmat, header = FALSE, sep = ",", dec = ".", row.names = NULL)
y <- read.table(element_name_response, header = FALSE, sep = ",", dec = ".", row.names = NULL)


setwd("/home/alb/Scrivania/PACS/Git_Folder/debugging_output/CPP/")
beta_hat_CPP <- read.table(element_name_CPP, header = FALSE, sep = ",", dec = ".", row.names = NULL)
obs_centered_CPP <- read.table(element_name_CPP_3, header = FALSE, sep = ",", dec = ".", row.names = NULL)



tmp <- inv.link(obs_centered_CPP)
datapoints <- data.frame(covariates[,1],covariates[,2],tmp, y)
names(datapoints) <- c("X1","X2","pseudo_obs","y")


datapoints <- data.frame(covariates[,1],covariates[,2],obs_centered_CPP)
names(datapoints) <- c("X1","X2","y")

#plot3d(datapoints$X1, datapoints$X2, datapoints$y, type="s", size=0.75, lit=FALSE)

idxs <- rep(TRUE, length(datapoints[,1]))
idxs[30] <- FALSE
datapoints2 <- datapoints[idxs,]

plot(datapoints2$y)
#identify(datapoints2$y)
idxs[21] = FALSE


idxs <- (abs(datapoints$y) < 80)

datapoints3 <- datapoints[idxs,]
plot(datapoints3$y)


#plot3d(datapoints3$X1, datapoints3$X2, datapoints3$y, type="s", size=0.75, lit=FALSE)

xrange <- range(datapoints3$X1)  
yrange <- range(datapoints3$X2)
res = 16
newdata <- expand.grid(x = seq(xrange[1], xrange[2], length.out = res),
                       y = seq(yrange[1], yrange[2], length.out = res)) 

names(newdata) <- c("X1", "X2") 
#newdata[["y"]] <- beta_hat_CPP[1,1]*newdata$X1 + beta_hat_CPP[2,1]*newdata$X2
newdata[["y"]] <- beta_hat_R[1,1]*newdata$X1 + beta_hat_R[2,1]*newdata$X2

#newdata


df2mat <- function(p, xvar = NULL, yvar = NULL, zvar = NULL) {
  if (is.null(xvar)) 
    xvar <- names(p)[1]
  if (is.null(yvar)) 
    yvar <- names(p)[2]
  if (is.null(zvar)) 
    zvar <- names(p)[3]  
  
  x <- unique(p[[xvar]])  
  y <- unique(p[[yvar]])  
  z <- matrix(p[[zvar]], nrow = length(y), ncol = length(x))  
  m <- list(x, y, z)
  names(m) <- c(xvar, yvar, zvar)  
  m
}

interleave <- function(v1, v2)  as.vector(rbind(v1,v2))

mpgrid_list <- df2mat(newdata)
# Make the plot with the data points
plot3d(datapoints3$X1, datapoints3$X2, datapoints3$y, type="s", size=0.5, lit=FALSE)
# Add the corresponding predicted points (smaller)
#spheres3d(datapoints3$X1, datapoints3$X2, newdata$y, alpha=0.4, type="s", size=0.5, lit=FALSE)
# Add line segments showing the error
#segments3d(interleave(datapoints3$X1, datapoints3$X2), interleave(datapoints3$X2, datapoints3$X2), interleave(datapoints3$y, newdata$y),alpha=0.4, col="red")
# Add the mesh of predicted values
surface3d(mpgrid_list$X1, mpgrid_list$X2, mpgrid_list$y, alpha=0.4, front="lines", back="lines")
