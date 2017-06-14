
#' Calculate gene variability stats
#'
#' \code{calc_varstats} calculates the variability of each gene, conditioned on their expression level
#'
#' @param counts A numeric matrix with read count expression values (genes x cells)
#' @param ercc.raw A numeric matrix with read count expression values of ERCC spike-in transcripts. If not available then input 'counts' also for this argument.
#' @param n.ercc.min An integer specifying the minimum number of ERCC transcripts with any read mapped to it. Samples with less ERCC transcripts are removed from the fitting against the ERCC data.
#' @param min.count An integer specifying the minimum read count for which a transcript is considered to be expressed
#' @param min.prop A numeric of the proportion of ERCC transcripts that needs to be expressed (above 'min.count'). Samples with lower proportion of expressed ERCC transcripts are removed from the ERCC fitting.
#' @param min.cv2 A numeric specifying the minimum squared coefficient of variation (CV^2) of the ERCC transcript expression. The quantile of the expression distribution of transcripts with higher CV^2 is used to set a minumum expression threshold for genes to be used for ERCC fitting (see argument 'quant.cutoff' below).
#' @param quant.cutoff A numeric specifying the quantile of the expression distribution of ERCC transcripts. Transcripts with mean expression below the expression at this quantile are excluded from ERCC fitting.
#' @param min.biol.cv2 A numeric specifying the minimum squared coefficient of variation of real biological data.
#' @param n.win An integer specifying the number of the extremest expression values in the expression distribution of a sample to be set to the 'n.win' highest value (Winsorizing).
#' @param method A character specifying if generalized linear model ('glm') or robust linear regression ('rlm') should be used.
#' @param min.mean.limit An integer, transcripts with mean expression below this limit are excluded from ERCC fitting.
#' @param counts.pseudo.count An integer with the pseudo-count to be added before taking the logarithm
#' @param ercc.pseudo.count An integer with the pseudo-count to be added to the ERCC data before taking the logarithm
#' 
#' @return A list with variability test statistics for all genes and results from fitting a glm or rlm regression model to the ERCC expression data
#'
#' @examples
#' data('scialdone')
#' counts = scialdone[['counts']]
#' counts_filt = genefilter_exprmat(counts, min.expr = 5, min.ncells.expr = 3)
#' 
#' quantile.cutoff = 0.9 ##filter out the most lowly expressed genes from the fitting
#' vargenes_stats = calc_varstats(counts_filt, counts_filt, quant.cutoff = quantile.cutoff)
#' 
#' @export
calc_varstats <- function(counts, ercc.raw, n.ercc.min = 2, min.count = 2, min.prop = 0.5, min.cv2 = 0.3, quant.cutoff = 0.9, min.biol.cv2 = 0.5^2, n.win = 2, method = 'glm', min.mean.limit = 0, counts.pseudo.count = 1, ercc.pseudo.count = 0){            
    
    ###
    #Filter samples such that all agree with samples in counts
    ###
    pass.samples = colnames(counts)
    counts = counts[, pass.samples]
    ercc.raw = ercc.raw[, intersect(colnames(ercc.raw), pass.samples)]
    
    
    ###
    #Filter ERCC samples on presence of spike-in
    ###

    #Params: n.ercc.min, min.count, min.prop
    ercc.filt = ercc.filter(ercc.raw, n.ercc.min, min.count, min.prop)
    dim(ercc.filt)

    
    ###
    #Normalize counts
    ###
    ##Real data
    if(min(counts, na.rm = TRUE) < 0){ ##if COMBAT normalized then there are negative values
        counts.pseudo.count = abs(min(counts, na.rm = TRUE)) + 1
    }
    real.size.factors = DESeq2::estimateSizeFactorsForMatrix(counts + counts.pseudo.count)
    expr = t( t(counts) / real.size.factors )

    ##ERCC
    if(min(ercc.filt, na.rm = TRUE) < 0){
        ercc.pseudo.count = abs(min(ercc.filt, na.rm = TRUE)) + 1
    }
    ercc.size.factors = DESeq2::estimateSizeFactorsForMatrix(ercc.filt + ercc.pseudo.count)
    ercc = t( t(ercc.filt) / ercc.size.factors )

    
    #######
    #Estimate tech noise
    #######

    #Params: min.cv2, quant.cutoff
    
    #Get stats
    ercc.stats = get.ercc.stats(ercc)    

    ##Filter NAs
    ercc.stats = ercc.stats[!is.na(ercc.stats[, 'mean']), ]
    
    #Filter genes on min mean expression
    min.mean = unname( stats::quantile( ercc.stats[which(ercc.stats[, 'cv2'] > min.cv2 ), 'mean'], quant.cutoff ) )
    if(min.mean < min.mean.limit){
        min.mean = min.mean.limit
    }
    pass.genes = rownames(ercc.stats)[which(ercc.stats[, 'mean'] >= min.mean)]
    print(length(pass.genes))

    #Fit to ERCC
    ercc.fit.res = fit.tech(ercc.stats, ercc.size.factors, real.size.factors, pass.genes, method)
    print(ercc.fit.res['expl.var'])
    
    
    ######
    #Test for var conditioned on tech var
    ######

    #Params: min.biol.cv2, n.win
    
    ##Stats for all cells
    win.expr = t(apply(expr, 1, winsorize, n.win))
    real.stats = get.real.stats(win.expr, ercc.fit.res, min.biol.cv2, method)

    ##Filter NAs
    real.stats = real.stats[which(!is.na(real.stats[, 'cv2'])), ]

    ##Store all results in list
    res.list = list(real.stats = real.stats, ercc.stats = ercc.stats, ercc.fit.res = ercc.fit.res)
    
    return(res.list)
}

ercc.filter <- function(ercc.raw, n.ercc.min, min.count, min.prop){
    
    #Rm samples with no ERCC in the rpkm annot (they didn't have ERCC spiked in to begin with either)
    ercc.ind = which(apply(ercc.raw, 2, function(j.ercc){length(which(!is.na(j.ercc))) > 0}))
    pass.samples = names(ercc.ind) 
    ercc.raw = ercc.raw[, pass.samples]
    ncol(ercc.raw)

    #Rm samples with reads mapped to less than n ERCC transcripts
    ercc.ind = which(apply(ercc.raw, 2, function(j.ercc){length(which(j.ercc != 0)) > n.ercc.min}))
    pass.samples = names(ercc.ind)
    length(pass.samples)
    ercc.raw = ercc.raw[, pass.samples]
        
    #For each cell, proportion of expressed transcripts > min.rpkm should be greater than min.prop
    ercc.ind = which(apply(ercc.raw, 2, function(j.ercc, min.count){j.ercc = j.ercc[which(j.ercc != 0)]; prop.expr = length(which(j.ercc > min.count)) / length(j.ercc); return(prop.expr >= min.prop)}, min.count = min.count))
    pass.samples = names(ercc.ind)    
    length(pass.samples)

    #apply filter
    ercc.raw = ercc.raw[, pass.samples]    

    ###
    #Filter ercc transcripts that are all zero
    ###
    ercc.ind = which(apply(ercc.raw, 1, function(j.ercc){length(which(j.ercc > min.count)) >= n.ercc.min}))
    pass.genes = names(ercc.ind)
    ercc.raw = ercc.raw[pass.genes, ]

    return(ercc.raw)
}

get.ercc.stats <- function(ercc){

    mean = rowMeans(ercc)
    var = genefilter::rowVars(ercc)
    cv2 = var / mean^2

    #bind
    ercc.stats = cbind(mean, var, cv2)

    return(ercc.stats)
}


fit.tech <- function(ercc.stats, ercc.size.factors, real.size.factors, pass.genes, method = 'glm'){

    if(method == 'glm'){
        fit.res = fit.tech.glm(ercc.stats, ercc.size.factors, real.size.factors, pass.genes)
    }
    if(method == 'rlm'){
        fit.res = fit.tech.rlm(ercc.stats, ercc.size.factors, real.size.factors, pass.genes)
    }

    return(fit.res)
}

fit.tech.glm <- function(ercc.stats, ercc.size.factors, real.size.factors, pass.genes){

    
    #Fit
    fit = statmod::glmgam.fit(cbind(a0 = 1, a1.tilde = 1 / ercc.stats[pass.genes, 'mean']), ercc.stats[pass.genes, 'cv2'])

    #get coefs
    xi = mean(1 / ercc.size.factors)
    a0 = unname(fit$coefficients["a0"])
    a1.tilde = fit$coefficients["a1.tilde"]
    a1 = unname(a1.tilde - xi)

    #psi + a1*theta
    pass.samples = names(ercc.size.factors)
    psia1theta = mean(1 / real.size.factors[pass.samples]) + a1 * mean(ercc.size.factors / real.size.factors[pass.samples])

    #Explained variance
    residual.var = stats::var(log(stats::fitted.values(fit)) - log(ercc.stats[pass.genes, 'cv2']) )
    total.var = stats::var(log(ercc.stats[pass.genes, 'cv2']))
    expl.var = 1 - (residual.var / total.var)

    #bind
    fit.res = c(xi, a0, a1.tilde, a1, psia1theta, residual.var, total.var, expl.var)
    names(fit.res) = as.character(expression(xi, a0, a1.tilde, a1, psia1theta, residual.var, total.var, expl.var))

    return(fit.res)
}

fit.tech.rlm <- function(ercc.stats, ercc.size.factors, real.size.factors, pass.genes){
    
    ercc.stats = as.data.frame(ercc.stats[pass.genes, ])
    log10.mean = log10(ercc.stats[, 'mean'])
    log10.cv2 = log10(ercc.stats[, 'cv2'])
    fit = MASS::rlm(log10.cv2 ~ log10.mean, data = ercc.stats)
    a0 = fit[['coefficients']][1]
    a1 = fit[['coefficients']][2]

    ##Explained variance    
    pseudo.count = abs(min(stats::fitted.values(fit))) + 0.1 ##predicted values can become negative
    residual.var = stats::var(log(stats::fitted.values(fit) + pseudo.count) - log(ercc.stats[pass.genes, 'cv2'] + pseudo.count) )
    total.var = stats::var(log(ercc.stats[pass.genes, 'cv2'] + pseudo.count))
    expl.var = 1 - (residual.var / total.var)

    ##store
    fit.res = c(a0, a1, residual.var, total.var, expl.var)
    names(fit.res) = as.character(expression(a0, a1, residual.var, total.var, expl.var))

    return(fit.res)    
}

get.real.stats <- function(expr, tech.fit.res, min.biol.cv2, method = 'glm'){
    
    #basic summary stats
    gene.means = rowMeans(expr)
    gene.vars = matrixStats::rowVars(expr)
    gene.cv2 = gene.vars / gene.means^2

    if(method == 'glm'){
        
        ##tech fit res
        a0 = tech.fit.res['a0']
        a1.tilde = tech.fit.res['a1.tilde']
        psia1theta = tech.fit.res['psia1theta']
        
        ##other test-statistica
        m = ncol(expr)
        cv2 = a0 + min.biol.cv2 + a0 * min.biol.cv2
        testDenom = (gene.means * psia1theta + gene.means^2 * cv2) / ( 1 + cv2/m )
        gene.var.quants = gene.vars * (m-1) / testDenom
        
        ##get p-value from a chi-sq dist
        p = 1 - stats::pchisq(gene.var.quants, m-1)

        ##mult-test correction
        p.adj = stats::p.adjust(p, "BH")

        ##bind res
        var.teststat = gene.var.quants
        test.res = cbind(gene.means, gene.vars, gene.cv2, min.biol.cv2, var.teststat, p, p.adj)
    }
    if(method == 'rlm'){
        a0 = tech.fit.res['a0']
        a1 = tech.fit.res['a1']
        pred.y = a0 + a1 * log10(gene.means)
        obs.y = log10(gene.cv2)
        var.ratio = 10^obs.y / 10^pred.y
        
        ##store
        var.teststat = var.ratio
        test.res = cbind(gene.means, gene.vars, gene.cv2, min.biol.cv2, var.teststat)
    }
    colnames(test.res)[1:3] = c('mean', 'var', 'cv2')
    
    #order by test var
    test.res = test.res[order(test.res[, 'var.teststat'], decreasing = TRUE), ]

    ##add rank
    rank = 1:nrow(test.res)
    test.res = cbind(test.res, rank)
    
    return(test.res)
}


#' Plot squared coefficient of variation against mean expression
#'
#' \code{plot_cv2vsmean} plots the squared coefficient of variation of all genes against the mean expression using the output from \code{calc_varstats} as input
#'
#' @param vargenes_stats A list with variability test statistics output by \code{calc_varstats}
#' @param sel.genes A character vector with names of genes to be highlighed in the plot as points
#' @param sel.highlight A logical specifying whether 'sel.genes' should be highlighed
#' @param plot.ercc.points A logical specifying whether the transcripts used for fitting should be shown as points
#' @param method A character specifying if a generalized linear model ('glm') or robust linear regression ('rlm') was used.
#'
#' @examples
#' data('scialdone')
#' counts = scialdone[['counts']]
#' counts_filt = genefilter_exprmat(counts, min.expr = 5, min.ncells.expr = 3)
#'
#' ##Calculate gene variability stats
#' quantile.cutoff = 0.9 ##filter out the most lowly expressed genes from the fitting
#' vargenes_stats = calc_varstats(counts_filt, counts_filt, quant.cutoff = quantile.cutoff)
#'
#' ##Select subset of top most variable genes
#' topvar.genes = rownames(vargenes_stats[['real.stats']])[1:500]
#'
#' ##Plot
#' plot_cv2vsmean(vargenes_stats, topvar.genes, plot.ercc.points = FALSE)

#' @export
plot_cv2vsmean <- function(vargenes_stats, sel.genes, sel.highlight = TRUE, plot.ercc.points = TRUE, method = 'glm'){

    real.stats = vargenes_stats[['real.stats']]
    ercc.stats = vargenes_stats[['ercc.stats']]
    tech.fit.res = vargenes_stats[['ercc.fit.res']]

    
    #get gene summary stats
    gene.means = real.stats[, 'mean']
    gene.cv2 = real.stats[, 'cv2']
    min.biol.cv2 = unique(real.stats[, 'min.biol.cv2'])
    
    ercc.gene.means = ercc.stats[, 'mean']
    ercc.gene.cv2 = ercc.stats[, 'cv2']
    
    ##Limits
    if(min(gene.means) < 0){
        gene.means = gene.means + abs(min(gene.means)) + 1
    }
    
    xlim = range(log10(gene.means))
    xlim[1] = floor(xlim[1])
    xlim[2] = ceiling(xlim[2])
    
    ylim = range(log10(gene.cv2))
    ylim[1] = floor(ylim[1])
    ylim[2] = ceiling(ylim[2])
    
    #labels
    xlab = "average normalized read count"
    ylab = "squared coefficient of variation (CV^2)"

    #color map
    n.genes = length(gene.means)
    gene2color = rep('blue', n.genes)
    names(gene2color) = names(gene.means)
    if(sel.highlight){
        gene2color[sel.genes] = '#C0007090'
    }

    
    ##*###
    ##Plot the real data
    ##*###
    if(sel.highlight){
        ##plot(log10(gene.means), log10(gene.cv2), pch=20, cex=.2, col = gene2color, xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n')
        graphics::smoothScatter(log10(gene.means), log10(gene.cv2), pch=20, cex=.2, col = gene2color, xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n', nrpoints = 0)
        graphics::points(log10(gene.means[sel.genes]), log10(gene.cv2[sel.genes]), pch = 20, cex = .5, col = 'red', xaxt = 'n', yaxt = 'n')
    }else{
        graphics::smoothScatter(log10(gene.means), log10(gene.cv2), pch=20, cex=.2, col = gene2color, xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n')
    }    
    x = xlim[1]:xlim[2]
    y = ylim[1]:ylim[2]
    graphics::axis(1, x, as.character(10^x))
    graphics::axis(2, y, as.character(10^y), las = 2)

    
    ##*###
    ##Plot technical noise fit
    ##*###
    if(method == 'glm'){
        
        ##get tech fit res
        a0 = tech.fit.res['a0']
        a1.tilde = tech.fit.res['a1.tilde']
        psia1theta = tech.fit.res['psia1theta']

        xg = 10^seq(xlim[1], xlim[2], length.out = 1000)
        cv2.pred = (a1.tilde / xg) + a0
        graphics::lines(log10(xg), log10(cv2.pred), col = "#FF000080", lwd = 3)

        ##Add a curve showing the expectation for the chosen biological CV^2 threshold
        graphics::lines(log10(xg), log10(psia1theta/xg + a0 + min.biol.cv2), lty="dashed", col="#C0007090", lwd=3 )
    }
    if(method == 'rlm'){
        a0 = tech.fit.res['a0']
        a1 = tech.fit.res['a1']
        
        x.mean.expr = seq(xlim[1], xlim[2], length.out=1000)
        log10.cv2.pred = a0 + a1 * x.mean.expr #the fit was done in log-space
        graphics::lines(x.mean.expr, log10.cv2.pred, col = "#FF000080", lwd = 3)
    }

    
    ##*###
    ##Add the normalised ERCC transcripts as points
    ##*###
    if(plot.ercc.points){
        graphics::points(log10(ercc.gene.means), log10(ercc.gene.cv2), pch=20, cex=1, col="black")
    }

}    
