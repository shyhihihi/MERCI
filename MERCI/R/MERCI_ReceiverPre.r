#' Predict receivers
#'
#' This function predicts the receivers based on adding a 'AND' gate to DNA and RNA rank.
#' @param MTvar_stat the output of MERCI_MTvar_cal function, it contains the DNA rank score for input mixed cells.
#' @param MT_constitute_stat the output of MERCI_LOO_MT_est or MERCI_MT_est function. it contains the RNA rank scores for input mixed cells.
#' @param The rank cutoff (top rank %), the cells with both DNA and RNA rank > this cutoff will be considered as receivers, default=50 
#' @return A data.frame with input mixed cells as rows, cell_names and prediction labels as columns.
#' @export
#' @examples MERCI_ReceiverPre(MTvar_stat, MT_constitute_stat, top_rank=50) ;

MERCI_ReceiverPre <- function(MTvar_stat, MT_constitute_stat, top_rank=50)
{
	cutoff <- top_rank *0.01 ;

	rank1 <- MTvar_stat$MTvar_rank ;
	names(rank1) <- rownames(MTvar_stat) ;

	rank2 <- MT_constitute_stat$W_Donor_rank ;
	names(rank2) <- rownames(MT_constitute_stat) ;
 
	if(length(rank1)!=length(rank2))
	{
		print('Error: MTvar_stat and MT_constitute_stat have different row numbers!')
		return(NULL) ;
	}

	break_point <- ceiling(quantile(1:length(rank1), cutoff)) ;
	Donor_Vrank <- sort(rank1, decreasing=TRUE)  ;
	Donor_MTrank <- sort(rank2, decreasing=TRUE) ;

	pa.list1 <- names(Donor_Vrank)[1:break_point] ;
	pa.list2 <- names(Donor_MTrank)[1:break_point] ;
	

	all.pa <- names(rank1) ;
	pre.positive <- Reduce(intersect, list(pa.list1, pa.list2)) ;

	
	label <- ifelse(all.pa%in%pre.positive, 'Receiver', 'non-Receiver') ;
	t.results <- cbind.data.frame(cell=all.pa, prediction=label, stringsAsFactors=FALSE) ;
	rownames(t.results) <- all.pa ; 
	return(t.results) ;
}

