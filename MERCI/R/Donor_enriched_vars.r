#' Identification of donor cell enriched mtSNVs
#'
#' This function identifys the donor cell enriched mtSNVs
#' @param mutMatrix_processed a preprocessed variant VAF matrix with rows are mtSNVs, columns are cells.
#' @param donor_cells the list of cell names for mitochondrial donor cells.
#' @param reciever_cells the list of cell names for supposed receiver cells, they can be pure reference non-receivers (regular MERCI) or mixed cells with unknown ratio of receirves vs non-receivers (MERCI_LOO).
#' @param Nmut_min the minimum requirement for mutated cell counts, only the mtSNVs with mutated cell counts > Nmut_min will be further calculated to determine if it is a donor cell enriched mtSNV.
#' @param pvalue Only return the mtSNVs with statistical p <= pvalue, defualt=0.05.
#' @param qvalue Only return the mtSNVs with ajusted p value from BH method <= qvalue, default=0.1.
#' @param OR Only return the mtSNVs with odds ratio (T cells vs reciever_cells) > OR, default=1.
#' @return a list mtSNVs that is considered as donor cell enriched
#' @export
#' @examples Donor_enriched_vars(mutMatrix_processed, donor_cells, reciever_cells, pvalue=0.05, qvalue=0.1, OR=1)

Donor_enriched_vars <- function(mutMatrix_processed, donor_cells, reciever_cells, Nmut_min=5, pvalue=0.05, qvalue=0.1, OR=1)
{
	t_varMatrix <- mutMatrix_processed ;
	mutNum <- apply(t_varMatrix, 1, function(x) sum(x>0, na.rm=TRUE)) ;
	t_index <- which(mutNum >= Nmut_min) ;
	t_varMatrix  <- t_varMatrix[t_index, ] ; 
	
	print( paste('There are', nrow(t_varMatrix), 'MTvariats that have mutated Cell number > ', Nmut_min) ) ;
	
	Cluster_cells <- intersect(donor_cells, colnames(t_varMatrix) ) ;	
	print( paste(length(Cluster_cells), "ref donor cells were used to calculate donor cell enriched mtSNVs!") ) ;
	Cluster_mutMatrix <- t_varMatrix[, Cluster_cells] ;
	Other_cells <- intersect(reciever_cells, colnames(t_varMatrix) ) ;
	print( paste(length(Other_cells), "ref non-receiver cells were used to calculate donor cell enriched mtSNVs!") ) ;
	Other_mutMatrix <- t_varMatrix[, Other_cells] ;
	
	AllVariants <- rownames(t_varMatrix) ;
	t_stat <- sapply(AllVariants, function(tmp_variant) {
		CT1_vaf <- Cluster_mutMatrix[tmp_variant, ] ;
		CT2_vaf <- Other_mutMatrix[tmp_variant, ] ;
		
		mut_n1 <- sum(CT1_vaf>0, na.rm=TRUE) ;
		total_n1 <- sum(CT1_vaf>=0, na.rm=TRUE) ;
		
		mut_n2 <- sum(CT2_vaf>0, na.rm=TRUE) ;
		total_n2 <- sum(CT2_vaf>=0, na.rm=TRUE) ;
		
		t_matrix <- matrix( c(mut_n1, total_n1-mut_n1, mut_n2, total_n2-mut_n2), 2, 2) ;
		tmp1 <- chisq.test(t_matrix) ;
		tmp2 <- fisher.test(t_matrix) ;
		p_value <- tmp1$p.value ;
		OR <- tmp2$estimate ;
	
		mutMedianAF1 <- median(CT1_vaf[CT1_vaf>0], na.rm=TRUE) ;
		mutMedianAF2 <- median(CT2_vaf[CT2_vaf>0], na.rm=TRUE) ;
		mutMedianAF1[is.na(mutMedianAF1)] <- 0 ;
		mutMedianAF2[is.na(mutMedianAF2)] <- 0 ;

		t_results <- c(mutNumber_inScluster=mut_n1, mutNumber_inOther=mut_n2, freq_inScluster=mut_n1/total_n1 , freq_inOther=mut_n2/total_n2, p_value=p_value, OR, MedianAF_inScluster=mutMedianAF1, MedianAF_inOther=mutMedianAF2) ;
		return(t_results) ;	
	})
	
	t_stat <- t(t_stat) ;
	q_values <- p.adjust(t_stat[, 'p_value'], method='BH') ;

	t_stat <- cbind.data.frame(t_stat, q_values, stringsAsFactors = FALSE) ;
	q_orders <- order(t_stat[, 'q_values']) ;
	t_stat <- t_stat[q_orders, ] ;
	
	t_stat$catogary <- "Non_sig" ;
	t_index1 <- which(t_stat[, 'p_value'] <= pvalue & t_stat[, 'q_values'] <= qvalue & t_stat[, 'odds ratio'] > OR ) ;
	t_stat$catogary[t_index1] <- "donor_enrich" ;
	t_index2 <- which(t_stat[, 'p_value'] <= pvalue & t_stat[, 'q_values'] <= qvalue & t_stat[, 'odds ratio'] < OR ) ;
	t_stat$catogary[t_index2] <- "donor_depleted" ;

	Donor_enriched_vars <- rownames(t_stat)[t_stat$catogary=="donor_enrich"] ;
	#Donor_enriched_vars <- rownames(t_stat)[t_stat$catogary=="Non_sig"] ;
	return(Donor_enriched_vars) ;
}

