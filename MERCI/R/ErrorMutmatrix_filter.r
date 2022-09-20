#' estimate the potential sequencing error for a variant matrix
#'
#' This function allows stratify each variant in the input matrix as error or non-error(detected) based on FPR
#' @param MT_matrix a variant matrix (data.frame)
#' 		  FPR False positive rate, default=0.0005
#' @return a data.frame with the final column showing the estimated resutlts
#' @export
#' @examples ErrorMutmatrix_filter(MT_matrix, FPR=0.0005)

ErrorMutmatrix_filter <- function(MT_matrix, FPR=0.0005)
{
	t_matrix <- MT_matrix ;
	t_matrix$detection <- apply(t_matrix, 1, function(x) ErrorMut_filter(x, FPR=FPR) ) 
	return(t_matrix) ;
}
