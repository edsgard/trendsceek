#' Set the mark distribution of a point pattern
#'
#' \code{set_marks} sets the mark distribution of a point pattern
#'
#' @param pp A point pattern
#' @param gene.marks A matrix (genes x cells) with as many columns as points (cells) in the point-pattern.
#' @param log.fcn A log-function with the gene.marks is to be transformed
#' @param pseudo.count A numeric specifying a pseudo-count to be added if a log-fcn is supplied.
#' 
#' @return A point pattern with marks having been set
#'
#' @examples
#' a_mat = matrix(rnorm(100, 0, 1), ncol = 10)
#' marx = matrix(rnorm(20, 0, 1), ncol = 10)
#'
#' pp = pos2pp(a_mat)
#' pp = set_marks(pp, marx)
#' 
#' @export
set_marks <- function(pp, gene.marks, log.fcn = NA, pseudo.count = 1){

    ##log
    if(!is.logical(log.fcn)){
        gene.marks = log.fcn(gene.marks + pseudo.count)
    }

    ##set marks    
    pp[['marks']] = as.data.frame(t(gene.marks))
    
    return(pp)
}

#' Subset the mark distributions of a point pattern
#'
#' \code{pp_select} subsets the mark distributions of a point pattern
#'
#' @param pp A point pattern
#' @param sel.genes A character vector of mark distributions (genes) to be kept
#' 
#' @return A point pattern retaining only the marks of the input sel.genes
#'
#' @examples
#' a_mat = matrix(rnorm(100, 0, 1), ncol = 10)
#' marx = matrix(rnorm(20, 0, 1), ncol = 10)
#'
#' pp = pos2pp(a_mat)
#' pp = set_marks(pp, marx)
#' pp = pp_select(pp, 1)
#' 
#' @export
pp_select <-  function(pp, sel.genes){

    marx_orig = pp[['marks']]
    marx = marx_orig[, sel.genes]
    
    ##handle drop of rownames if a single gene selected (TBD: possibly use drop = FALSE in the subsetting instead)
    if(length(sel.genes) == 1 & !is.numeric(marx_orig)){
        names(marx) = rownames(marx_orig)
    }

    ##set marks
    pp[['marks']] = marx
    
    return(pp)
}

#' Convert positions to point-pattern
#'
#' \code{pos2pp} converts 2-dimensional positions to a point-pattern object
#'
#' @param pos_mat A numeric matrix where each row corresponds to a point with x positions in the first column and y-positions in the second column.
#' 
#' @return A point pattern
#'
#' @examples
#' a_mat = matrix(rnorm(100, 0, 1), ncol = 2)
#' 
#' pp = pos2pp(a_mat)
#' 
#' @export
pos2pp <- function(pos_mat){
    x.range = range(pos_mat[, 1])
    y.range = range(pos_mat[, 2])
    pp = spatstat::ppp(x = pos_mat[, 1], y = pos_mat[, 2], x.range, y.range)

    return(pp)
}

#' Filter genes on being expressed
#'
#' \code{genefilter_exprmat} filter genes on being expressed in a number of cells.
#'
#' @param exprmat A numeric matrix (genes x cells)
#' @param min.expr A numeric specifying the minimum expression required in a cell
#' @param min.ncells.expr An integer specifying the minimum number of cells a gene needs to be expressed in with an expression level above 'min.expr'
#' 
#' @return A matrix where genes not passing the expression criteria have been removed
#'
#' @examples
#' a_mat = matrix(rnorm(100, 0, 1), nrow = 10, dimnames = list(1:10, 1:10))
#' 
#' filt_mat = genefilter_exprmat(a_mat, min.expr = 0, min.ncells.expr = 3)
#' print(nrow(filt_mat))
#' 
#' @export
genefilter_exprmat <- function(exprmat, min.expr, min.ncells.expr){
    gene2nsamples.expr = apply(exprmat, 1, function(jgene.expr, expr.cutoff){length(which(jgene.expr > expr.cutoff));}, expr.cutoff = min.expr)
    rm.genes = names(gene2nsamples.expr)[which(gene2nsamples.expr < min.ncells.expr)]
    exprmat_filt = exprmat[setdiff(rownames(exprmat), rm.genes), ]
    return(exprmat_filt)
}

#' Perform dimensionality reduction using t-SNE
#'
#' \code{trend.tsne} runs t-SNE, including a PCA pre-processing step
#'
#' @param counts A numeric matrix (genes x cells).
#' @param tsne.k An integer with the dimension of the resulting embedding.
#' @param init.dims An integer with the number of dimensions to be reduced to by the initial PCA
#' before t-SNE is run.
#' @param perp.frac A numeric specifying the fraction of samples to be used as neighbors by t-SNE
#' (perplexity parameter).
#' @param max.iter An integer with the maximum number of iterations to perform.
#' @param epoch An integer with the number of iterations between printing update messages.
#' @param pseudo.count A numeric with the pseudo.count to be used before logging (log10) the read
#' counts.
#' 
#' @return A matrix with the positions of the cells in the t-SNE embedding
#'
#' @examples
#' data('scialdone')
#' counts = scialdone[['counts']]
#' counts_norm = deseq_norm(counts, min.count = 1)
#' tsne_res = trend.tsne(counts_norm[1:500, ], tsne.k = 2, init.dims = 100, perp.frac = 0.2,
#' max.iter = 400, epoch = 50)
#' 
#' @export
trend.tsne <- function(counts, tsne.k = 2, init.dims = 100, perp.frac = 0.2, max.iter = 500, epoch = 50, pseudo.count = 1){

    
    ##*#####
    ##PCA
    ##*#####
    
    ##filter constant rows
    counts = rm.const.vec(counts, col.rm = FALSE)

    ##log
    if(min(counts, na.rm = TRUE) < 0){ ##if COMBAT normalized then there are negative values
        pseudo.count = abs(min(counts, na.rm = TRUE)) + 1
    }
    log.counts = log10(counts + pseudo.count)
    
    #pca
    pca = stats::prcomp(t(log.counts), scale=T, center=T)
    pca.basis = pca[['x']]
    
    
    ##*#####
    ##tSNE
    ##*####
        
    ##get dimred data
    X = pca.basis
    
    ##set params
    n.samples = nrow(X)
    perp = round(perp.frac * n.samples)        
    j.init.dims = min(init.dims, n.samples)
    cat(sprintf('nsamples: %i\n', n.samples))
    cat(sprintf('perp: %i\n', perp))
    cat(sprintf('n.dims input: %i\n', j.init.dims))
    
    ##tSNE
    tsne.res = tsne::tsne(X, perplexity = perp, k = tsne.k, max_iter = max.iter, epoch = epoch, initial_dims = j.init.dims, whiten = FALSE)
    colnames(tsne.res) = 1:tsne.k
    rownames(tsne.res) = rownames(X)
    colnames(tsne.res) = paste('d', colnames(tsne.res), sep = '')

    return(tsne.res)
}

rm.const.vec <- function(data.mat, col.rm = TRUE, row.rm = TRUE){

    #rm cols with constant variance
    if(col.rm){
        const.cols = which(apply(data.mat, 2, function(x){stats::var(x) == 0}))
        if(length(const.cols) != 0){
            warning('There were constant columns. These were removed.')
            data.mat = data.mat[, setdiff(1:ncol(data.mat), const.cols)]
        }
    }

    #rm rows (genes) with constant variance
    if(row.rm){
    
        const.rows = which(apply(data.mat, 1, function(x){stats::var(x) == 0}))
        if(length(const.rows) != 0){
            warning('There were constant rows. These were removed.')
            data.mat = data.mat[setdiff(1:nrow(data.mat), const.rows), ]
        }
    }

    return(data.mat)
}

winsorize <- function(x, n.win = 1){
###Winsorize extreme values
###Ex: Two most extreme values (from each side):
###t(apply(ed, 1, winsorize, 2))

    n.vals = length(x)
    fraction = n.win / n.vals
    if(length(fraction) != 1 || fraction < 0 || fraction > 0.5){
        stop("bad value for 'fraction'")
    }

    win.sorted.ind = order(x)
    x[win.sorted.ind[1:n.win]] = x[win.sorted.ind[n.win + 1]]
    x[win.sorted.ind[(n.vals - n.win + 1):n.vals]] = x[win.sorted.ind[n.vals - n.win]]

    ##fraction-based
    if(0){
        lim = stats::quantile(x, probs = c(fraction, 1 - fraction))

        ##set extreme values
        x[ x < lim[1] ] = lim[1]
        x[ x > lim[2] ] = lim[2]
    }
    
    return(x)
}


#' Normalize read counts using DESeq2
#'
#' \code{deseq_norm} normalizes read counts using size-factors estimated by DESeq2
#'
#' @param counts A numeric matrix with read count expression values (genes x cells)
#' @param min.count A numeric specifying the minimum total read count across all cells that a transcript need to have. Transcripts with lower expression level are removed from size-factor estimation.
#' 
#' @return A numeric matrix with normalized expression values (genes x cells)
#'
#' @examples
#' data('scialdone')
#' counts = scialdone[['counts']]
#' counts_norm = deseq_norm(counts, min.count = 1)
#' 
#' @export
deseq_norm <- function(counts, min.count){

    ##create deseq dataset
    colData = as.data.frame(colnames(counts), stringsAsFactors = FALSE)
    colnames(colData) = 'cellName'
    dds = DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ 1)

    ##filter genes on being expressed
    dds = dds[ rowSums(DESeq2::counts(dds)) > min.count, ]

    ##estimate size factors
    dds = DESeq2::estimateSizeFactors(dds)
    ##size.factors = DESeq2::sizeFactors(dds)

    ##divide with size factor
    counts_norm = DESeq2::counts(dds, normalized = TRUE)

    return(counts_norm)
}
