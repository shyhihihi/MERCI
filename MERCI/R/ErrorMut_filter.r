#' Detect the potential sequencing error
#'
#' This function allows stratify each variant as error or non-error(detected) based on FPR
#' @param mut a SNV
#' 		  e the random sequencing error, default=0.001/3
#' 		  FPR False positive rate, default=0.0005
#' @return a value of 'detected' or 'error'
#' @export
#' @examples ErrorMut_filter(mut, FPR=0.0005)

ErrorMut_filter <- function(mut, e=0.001/3, FPR=0.0005)
{
	n = as.numeric(mut[5]) ;
	s = as.numeric(mut[4]) ;
	P = 1 - sum(dbinom(0:(s-1), n, e));
	if(P < FPR)
	{ return("detected") } else{
	  return("error")	
	}
}
