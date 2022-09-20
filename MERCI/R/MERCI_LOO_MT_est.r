#' Decovolution to estimate MT constitute with Leave-one-out (LOO) strategy
#'
#' This function performas SVR to estiamte the relative abundance for donor-transfered and intrinsic mitochondria.
#' @param cell_exp the gene expression profiles with rows are gene sumbols and columns are cells, the MT genes must be included.
#' @param donor_cells the list of cell names for mitochondrial donor cells.
#' @param reciever_cells the list of cell names for supposed MT reciever cells, ideally this list should contain a unknown ratio of true receivers and non-receivers.
#' @param organism the species of the analyzed sample data,  only supports 'Human' and 'mouse', default='Human'. If s.markers are provided, this parameter can be set as NULL.
#' @param s.markers the list of MT genes used for deconvolution, if not provided, then all the MT genes included in the cell_exp will be used, default=NULL.
#' @param epsilon epsilon in the insensitive-loss function of ε-SVR, default=0.1.
#' @return a a data.frame with the estimated relative abundance of mitochondria from donor and receiver itself, and the RNA ranks.
#' @export
#' @examples MERCI_LOO_MT_est(cell_exp, reciever_cells, donor_cells, organism='mouse') ;

MERCI_LOO_MT_est <- function(cell_exp, reciever_cells, donor_cells, organism='Human', s.markers=NULL, epsilon=0.1)
{
	print(paste('organism is', organism)) ;
	if(organism=='Human')
	{
		MT_genes <- rownames(cell_exp)[grep('^MT-', rownames(cell_exp))] ;
	}

	if(organism=='mouse')
	{
		MT_genes <- rownames(cell_exp)[grep('^mt-', rownames(cell_exp))] ;
	}
	
	if(is.null(s.markers))
	{
		s.markers <- MT_genes ;
		print(paste('A total of', length(s.markers), 'MT genes are used as features!'))
	}

	#using mean expression as reference expressioon
	reciever_cells <- intersect(reciever_cells, colnames(cell_exp)) ;
	donor_cells <- intersect(donor_cells, colnames(cell_exp)) ;
	print(paste('There are', length(reciever_cells), 'supposed reciever cells with RNA expression information!') ) ;
	print(paste('There are', length(donor_cells), 'supposed donor cells with RNA expression information!'))
	
	D_exp <- cell_exp[, donor_cells] ;
	marker_D <- D_exp[s.markers, ] ;
	ref_D <- apply(marker_D, 1, mean, na.rm=TRUE) ;
	
	##leave one out when calculate the ref_R
	R_exp <- cell_exp[, reciever_cells] ;
	marker_R <- R_exp[s.markers, ] ;
	#svr
	library(e1071) ;
	MTmarker_exp <- marker_R ;
	cell_w <- sapply(1:ncol(MTmarker_exp), function(x) {
		t_exp <- MTmarker_exp[, x] ;
		t.marker_R <- MTmarker_exp[, -x] ;
		ref_R <- apply(t.marker_R, 1, mean, na.rm=TRUE) ;
		data <- cbind.data.frame(X1=ref_D, X2=ref_R, Y=t_exp, stringsAsFactors=FALSE) ;
		modelsvm = svm(Y~X1+X2, data, kernel='linear', epsilon=epsilon) ;
		#Calculate parameters of the SVR model
		W = t(modelsvm$coefs) %*% modelsvm$SV
		return(W)
	})

	cell_w <- t(cell_w) ;
	dimnames(cell_w) <- list(colnames(MTmarker_exp), c('W_Donor', 'W_Receiver')) ;

	t.cell_w <- cell_w ;
	t.cell_w <- apply(t.cell_w, 1, function(x) {
		if(min(x)>= 0 )
		{
			return(x/sum(x)) 
		}else{
			return( (x-min(x))/(max(x)-min(x)) )
		}
	} )

	t.cell_w <- t(t.cell_w) ;

	par(mfrow = c(1, 2))
	boxplot(cell_w[,'W_Donor'], cell_w[,'W_Receiver'], main='MT weight', names=c('MT_Donor_index','MT_Receiver_index')) ;
	boxplot(t.cell_w[,'W_Donor'], t.cell_w[,'W_Receiver'], main='MT weight', names=c('MT_Donor_frac','MT_Receiver_frac')) ;

	Cell_weights <- cbind.data.frame(cell_w, t.cell_w) ;
	colnames(Cell_weights) <- c('Donor_MT_ind', 'Receiver_MT_ind', 'Donor_MT_frac', 'Receiver_MT_frac') ;

	W_Donor_rank <- cell_w[,'W_Donor'] ;
	W_Donor_rank[order(W_Donor_rank)] <- 1:length(W_Donor_rank) ;

	W_Receiver_rank <- cell_w[,'W_Receiver'] ;
	W_Receiver_rank[order(W_Receiver_rank, decreasing=TRUE)] <- 1:length(W_Receiver_rank) ;

	t.results <- cbind.data.frame(Cell_weights, W_Donor_rank, W_Receiver_rank, stringsAsFactors=FALSE) ;
	return(t.results)
}

