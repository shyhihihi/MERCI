#' The weights (obs(mi)/exp(mi) ratios) MT mutations
#'
#' This function calculates the weights of MT mutations in a given cell.
#' @param muts the mitochondrial mutation list.
#' @param t_cell a cell name or barcode.
#' @param muts_inf the annotation of muts, also as the output of mut_inf_extr function.
#' @param mutMatrix_processed the processed vaf matrix of mtSNV generated by MTmutMatrix_refined function.
#' @param e sequencing error, default=0.001.
#' @return a a list weights with same length of muts.
#' @export
#' @examples w_cal(muts, t_cell, muts_inf, mutMatrix_processed)

w_cal <- function(muts, t_cell, muts_inf, mutMatrix_processed, MT_coverage, e=0.001)
{
	W <- sapply(muts, function(t_mut) {
		t_pos <- strsplit(t_mut, '_')[[1]][2] ;
		
		s <- muts_inf[t_mut, 'vaf_po_receiver'] ;
		Pi0 <- (1-s)*(e/3)+s*(1-e) ;
		
		mut_vaf <- mutMatrix_processed[t_mut, t_cell] ;
		ni <- MT_coverage[t_pos, t_cell] ;
		expected_mi <- 0 ;

		if(is.na(mut_vaf))
		{
			mi <- NA ;
		} else {
			mi <- mut_vaf * ni ;
			expected_mi <- ni*Pi0 ;
		}

		wi <- mi/expected_mi ;
	})
	
	return(W) ;
}

