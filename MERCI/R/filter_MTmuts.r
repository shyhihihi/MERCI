#' Preprocess of mtSNV matrix
#'
#' This function remove the SNVs with observaRate less than a pre-defined cutoff
#' @param varMatrix a transformed variant VAF matrix with rows are mtSNVs, columns are cells.
#' @param MT_coverage The coverage information of mitochondrial genome for cells.
#' @param donor_cells the list of cell names for mitochondrial donor cells.
#' @param reciever_cells the list of cell names for candidate receipt cells of MT.
#' @param min_d the minimum read coverage required for each mtSNV in different cells, the mutation record will be considered as 'unobservable' and set as NA if the the coverage at the mutation locus < min_d, default=5.
#' @param min_observeRate the minimum requirement for observatiion Rate (see parameter min_d), eahc mtSNV will be removed if its observeRate < min_observeRate, default=10%.
#' @return a processed VAF matrix
#' @export
#' @examples filter_MTmuts(varMatrix, MT_coverage, donor_cells, reciever_cells)

filter_MTmuts <- function(varMatrix=mutMatrix_Hcover, MT_coverage, donor_cells, reciever_cells, min_d=5, min_observeRate= 0.1)
{
	t_cells <- colnames(varMatrix) ;
	if(!all(t_cells %in% colnames(MT_coverage)))
	return("error, some cells have no coverage information!")
	
	t_coverage <- MT_coverage[, t_cells] ;
	
	t_pos <- sapply(strsplit(rownames(varMatrix), "_"), function(x) x[2]) ;
	t_coverMatrix <- t_coverage[t_pos, ] ; 
	indexMatrix <- t_coverMatrix < min_d ;
	filter_matrix <- varMatrix ;
	filter_matrix[indexMatrix] <- NA ;
	filter_matrix2 <- filter_matrix[, c(donor_cells, reciever_cells)] ;

	All_variants <- rownames(filter_matrix2) ;
	cellNum <- ncol(filter_matrix2) ;
	total_unRate <- apply(filter_matrix2, 1, function(x) sum(is.na(x))/cellNum ) ;
	
	doN <- length(donor_cells) ;
	reN <- length(reciever_cells) ;
	t.cell_inf <- cbind.data.frame(Cell_id=c(donor_cells, reciever_cells), label=c(rep('donor_cells', doN), rep('reciever_cells', reN)), stringsAsFactors=FALSE) ;
	cell_clusterLabel <- t.cell_inf$label ;
	
	celltype_unRate <- sapply(All_variants, function(tmp_variant) {
		tapply(filter_matrix2[tmp_variant, ], cell_clusterLabel, function(x) { sum(is.na(x))/length(x)} ) ;
	}) ;
	
	celltype_unRate <- t(celltype_unRate) ;
	options(stringsAsFactors=FALSE) ;
	unobserveRate_stat <- cbind.data.frame(total_unRate, celltype_unRate) ;
	observeRate_stat <- 1 - unobserveRate_stat ;
	colnames(observeRate_stat)[1] <- 'total' ;
	
	t_index <- which(apply(observeRate_stat, 1, function(x) {
		all( x >= min_observeRate )
	}) )
	selected_variants <- All_variants[t_index] ;
	New_mutMatrix <- filter_matrix[selected_variants, ] ;
	return(New_mutMatrix) ;
}

