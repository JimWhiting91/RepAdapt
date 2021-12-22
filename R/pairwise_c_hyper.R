#' Calculate C-hyper
#'
#' @description \code{pairwise_c_hyper} returns the pairwise C-score using
#' the hypergeometric approach. Setting na.rm = T assumes that missing values
#' are true negatives. Input array should contain binary data where 1's represent
#' genes that are adapted, while 0's represent genes that are not adapted.
#'
#' @param input numeric array with genes in rows and lineages in columns
#' @param na.rm boolean
#'
#' @return length-one numeric
#'
#' @examples
#' array1 <- array (0,c(5000,20))
#' array1[cbind(1:20,1:20)] <- 1
#' array1[1:5,1:3] <- 1
#' array1[6:10,4:6] <- 1
#' pairwise_c_hyper (array1)
#'
#' @export


pairwise_c_hyper <- function (input, na.rm = F){

  numcol <- ncol (input)

  results_c_hyper <- array (NA, c (numcol,numcol))

  for(loop1 in 1:(numcol - 1)){
    for(loop2 in loop1:numcol){
      if (loop1 != loop2){
        ax <- sum (input[,loop1], na.rm = na.rm)
        ay <- sum (input[,loop2], na.rm = na.rm)
        g0 <- nrow (input)

        sd_hyp <- sqrt((ax*ay)*(g0-ax)*(g0-ay)/(g0^2*(g0-1)))

        exp_hyp <- ax * ay / g0

        obs_hyp <- sum (input[,loop1] == 1 & input[,loop2] == 1, na.rm = na.rm)

        if (sd_hyp != 0){
          results_c_hyper[loop1,loop2] <- (obs_hyp - exp_hyp) / sd_hyp
        } else {
          results_c_hyper[loop1,loop2] <- 0
          warning ('Some pairwise contrasts have no shared adapted loci')
        }

      }
    }
  }

  mean (results_c_hyper,na.rm = T)

}
