#' Plot coverage distribution
#'
#' This function generate a visualization of coverage distrubution for selected cells
#' @param cells a list of cell names
#' @param coverage the mitochondrial read coverage information of cells
#' @return a boxplot showing the coverage distribution of selected cells
#' @export
#' @examples coverage_plot(cells, coverage)

coverage_plot <- function(cells, coverage)
{
	t_coverage <- coverage[, cells] ;
	depth_list <- c(1, 5, 10, 20, 30, 40, 50, 100) ;
	len <- nrow(t_coverage) ;
	
	t_results <- sapply(depth_list, function(x) {
		tmp <- apply(t_coverage, 2, function(y) sum(y>= x)/len )
		tmp
	})
	
	colnames(t_results) <- paste(depth_list, "x", sep="") ;
	boxplot(t_results, ylab='Proportion of MT genome', xlab='Read coverage', main='The coverage distribution of cells') ;
}
