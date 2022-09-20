#' read in mtSNV of 10x scRNA-seq data
#'
#' This function read in MT variants table (*.MT_variants.txt) generated by MERCI-mtSNP, only works for 10x_scRNA-seq data type.
#' @param varFile The path to .MT_variants.txt file generated by MERCI-mtSNP.
#' @param minReads Only mtSNVs in Cells with MT read counts > minReads will be retained.
#' @return the mtSNV annotation table for cells with MT read counts more than a user-defined number
#' @export
#' @examples readMTvar_10x(varFile='./X.MT_variants.txt', minReads=500)

####读入突变数据，并进行基本过滤
readMTvar_10x <- function(varFile, minReads=500)
{
	MT_variants <- read.csv(varFile, header=TRUE, sep="\t", fill=TRUE, quote="", check.names=FALSE, row.names = NULL) ;
	
	totalCells <- length(unique(MT_variants$Cell)) ;
	print(paste('a total of', totalCells, 'Cells have at least one MT variant!') ) ;
	
	MT_variants <- MT_variants[MT_variants$Cell_reads > minReads, ] ; 
	r.Cells <- length(unique(MT_variants$Cell)) ;
	print(paste('a total of', r.Cells, 'Cells have >', minReads, 'mapped MT reads !') ) ;
	return(MT_variants) ;
}


