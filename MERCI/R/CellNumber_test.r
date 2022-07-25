#' Significance estimation for MERCI positive calls
#'
#' This function tests if enough true-receivers are included in the input mixed cells based on the statistic Rcm.
#' @param MTvar_stat the output of MERCI_MTvar_cal function, it contains the DNA rank score for input mixed cells.
#' @param MT_constitute_stat, the output of MERCI_LOO_MT_est or MERCI_MT_est function. it contains the RNA rank scores for input mixed cells.
#' @param top_ranks a list of rank cutoffs (top rank %), default is a list contains 10 cutoffs from top 100% to top 10%. 
#' @param Number_R The times of randomly permutation, default=1000. 

#' @return a barlpot showing the Rcm values at different rank cutoffs, and a data.frame with rank cutoffs as rows, associated statistics as columns (p_value, fdr, preN, random_medianN, random_maxN, random_minN, Rcm)
#' @export
#' @examples CellNumber_test(MTvar_stat, MT_constitute_stat, Number_R=1000) ;

CellNumber_test <- function(MTvar_stat, MT_constitute_stat, top_ranks=seq(100, 10, by=-10), Number_R=1000)
{
	perf_list_test <- c();
	for(i in 1:length(top_ranks) )
	{
		top_rank <- top_ranks[i] ;
		t.MTreceiver_pre <- MERCI_ReceiverPre(MTvar_stat, MT_constitute_stat, top_rank) ;
		#new version of Random_perfs_test
		MTreceiver_evl <- Random_perfs_test(t.MTreceiver_pre, top_rank=top_rank, Number_R=Number_R) ;
		perf_list_test <- cbind(perf_list_test, MTreceiver_evl) ;
		colnames(perf_list_test)[i] <- paste('rank ', top_rank, '%', sep='') ;
		print(top_rank) ;
	}
	
	p_value <- apply(perf_list_test, 2, function(x) {
		Ob.n <- x[1] ; 
		R.ns <- x[-1] ;
		sum(R.ns >= Ob.n)/length(R.ns) ;
	})
	medianN <-  apply(perf_list_test[-1, ], 2, median) ;
	preN <- perf_list_test[1, ] ;
	maxN <- apply(perf_list_test[-1, ], 2, max) ;
	minN <- apply(perf_list_test[-1, ], 2, min) ;
	fdr <- p.adjust(p_value, method='BH') ;
	
	Ratio <- preN/maxN ;
	t.results <- cbind.data.frame(p_value, fdr, preN, random_medianN=medianN, random_maxN=maxN, random_minN=minN, Rcm=Ratio) ;
	
	barplot(Ratio, las=2, names.arg=rownames(t.results), ylim=c(0,2), ylab='Rcm value', col='peachpuff')
	abline(h=1, lty='dotted', lwd=3, col = 'red') ;

	return(t.results) ;
}

