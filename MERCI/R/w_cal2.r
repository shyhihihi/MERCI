#' The weights (Likelihood ratios) MT mutations
#'
#' This function calculates the weights of MT mutations in a given cell.
#' @param muts the mitochondrial mutation list.
#' @param t_cell a cell name or barcode.
#' @param muts_inf the annotation of muts, also as the output of mut_inf_extr function.
#' @param mutMatrix_processed the processed vaf matrix of mtSNV generated by MTmutMatrix_refined function.
#' @param e sequencing error, default=0.001.
#' @return a a list weights with same length of muts.
#' @export
#' @examples w_cal2(muts, t_cell, muts_inf, mutMatrix_processed)

w_cal2 <- function(muts, t_cell, muts_inf, mutMatrix_processed, MT_coverage, e=0.001)
{
	library(zipfR) ;
	e <- 0.001 ;
	W <- sapply(muts, function(t_mut) {
		t_pos <- strsplit(t_mut, '_')[[1]][2] ;
		mut_vaf <- mutMatrix_processed[t_mut, t_cell] ;
		ni <- MT_coverage[t_pos, t_cell] ;
		
		if(is.na(mut_vaf))
		{
			mi <- NA ;
			wi <- NA ;
		}else {
			mi <- mut_vaf * ni ;
			if(mi==0)
			{wi <- 0 ; } else {
				PiC_e <- muts_inf[t_mut, 'vaf_po_receiver'] ;
				pi0_e <- (1-PiC_e)*(e/3)+PiC_e*(1-e) ;

				Mi <- muts_inf[t_mut, 's_reads_receiver'] ;
				Ni <- muts_inf[t_mut, 't_reads_receiver'] ;

				PiC_d <- sqrt(Mi*(Ni-Mi)/Ni^3) ;

				Pi0_d <- (1-4/3*e) * PiC_d  ;

				s1 <- min(pi0_e+1.96*Pi0_d, 1) ;
				s2 <- max(pi0_e-1.96*Pi0_d, 0) ;
				#s3 <- min(pi1_e+1.96*Pi1_d, 1) ;
		
				#numerator <- Rbeta(s3, mi+1, ni-mi+1) - Rbeta(s2, mi+1, ni-mi+1) ;
				numerator <- 1 - Rbeta(s2, mi+1, ni-mi+1) ;
				denominator <- Rbeta(s1, mi+1, ni-mi+1) ;

				wi <- numerator/denominator ;
				}
			}
		return(wi)
	})
	return(W) ;
}

