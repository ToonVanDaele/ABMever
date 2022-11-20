# Functions for population matrices

#--------------------------------------------------------
# Create a projection matrix based on demographic parameters
#
# The matrix model is female only with a post-breeding sensus and 3 age classes.
#
# S = Ssurvival for each ageclass (vector of length 3)
# F = fertility for each ageclass (vector of length 3)
# H = hunting for each ageclass (vector of length 3)
#
# return list with two 3x3 projection matrices, with and without hunting
#
set_projmat_post <- function(S, F, H){

  mat <- matrix(c(F[1]*S[1]*0.5, F[2]*S[2]*0.5,	F[3]*S[3]*0.5,
                  S[1],          0,  	          0,
                  0,             S[2],	        S[3]),
                ncol = 3, byrow = TRUE)
  mat_h <- t(t(mat) * (1 - H))
  return(list(mat = mat, mat_h = mat_h))
}


#--------------------------------------------------------
# Execute projection
#
# wrapper for pop.projection
#
# Output is similar to pop.projection but in a data frame with long format
# column names: time, agecl, n
#
matrix_proj <- function(A, n, iterations){

  mm <- pop.projection(A = A, n = n, iterations = iterations)
  rownames(mm$stage.vectors) <- ageclasses

  # Process output matrix model
  mms <- mm$stage.vectors %>%
    t() %>%
    as.data.frame() %>%
    mutate(time = row_number(),
           sim = "mm") %>%
    pivot_longer(cols = all_of(ageclasses), names_to = "agecl", values_to = "n" )
  return(mms)
}
