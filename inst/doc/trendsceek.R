## ------------------------------------------------------------------------
##Load library
library('trendsceek')

## ------------------------------------------------------------------------
##Set parameters
##approximate number of cells
nsim_cells = 300

##low and high expression levels
low_expr = c(10, 10)
high_expr = c(20, 50)

##Randomly position the cells in 2d space
pp = sim_pois(nsim_cells)

##Add gene expression levels as mark distributions
pp = add_markdist_step(pp, low_expr, high_expr)
pp = add_markdist_hotspot(pp, low_expr, high_expr)
pp = add_markdist_streak(pp, low_expr, high_expr)

## ----sim-scatter-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_pp_scatter(pp, log_marks = FALSE, scale_marks = FALSE, pal.direction = -1)


## ----sim-density-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_pp_density(pp, log_marks = FALSE)

## ----sim-trendsceek, results='hide'--------------------------------------
##set parameters
nrand = 100
ncores = 1

##run
trendstat_list = trendsceek_test(pp, nrand, ncores)

##extract significant genes
alpha = 0.1 ##Benjamini-Hochberg
sig_list = extract_sig_genes(trendstat_list, alpha)

## ----sim-trendstats-plot, fig.width = 7, fig.height = 7, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
sig_genes = sig_list[['markcorr']][, 'gene']
plot_trendstats(trendstat_list, sig_genes)

## ----real-data-----------------------------------------------------------

## ------------------------------------------------------------------------
data('scialdone')
counts = scialdone[['counts']]

## ------------------------------------------------------------------------
##Pass genes having at least 3 cells with minimum 5 in read count
min.ncells.expr = 3
min.expr = 5
counts_filt = genefilter_exprmat(counts, min.expr, min.ncells.expr)
dim(counts_filt)

## ------------------------------------------------------------------------
quantile.cutoff = 0.9 ##filter out the most lowly expressed genes from the fitting
method = 'glm' ##For (robust) linear regression set to 'rlm'
vargenes_stats = calc_varstats(counts_filt, counts_filt, quant.cutoff = quantile.cutoff, method = method)

## ------------------------------------------------------------------------
n.topvar = 500
topvar.genes = rownames(vargenes_stats[['real.stats']])[1:n.topvar]

## ----real-varstats-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot.ercc.points = FALSE
plot_cv2vsmean(vargenes_stats, topvar.genes, plot.ercc.points = plot.ercc.points)

## ------------------------------------------------------------------------
min.count = 1
counts_norm = deseq_norm(counts, min.count)

## ------------------------------------------------------------------------
counts_sub = counts_norm[topvar.genes, ]
dim(counts_sub)

## ------------------------------------------------------------------------
##params
tsne.k = 2 #Number of tSNE dimensions
init.dims = 100 #Number of PCA-dimensions to reduce to before tSNE
perp.frac = 0.2 #Fraction of cells 
max.iter = 300 #Maximum number of iterations
epoch = 50 #print status every epoch iteration

##run tsne including initial pca
tsne_res = trend.tsne(counts_sub, tsne.k, init.dims, perp.frac, max.iter, epoch)

## ----create-point-pattern------------------------------------------------

## ------------------------------------------------------------------------
##Convert tSNE cell positions to point pattern
pp = pos2pp(tsne_res)

##Set marks as the logged normalized gene expression 
log.fcn = log10
pp = set_marks(pp, counts_sub, log.fcn = log.fcn)

## ------------------------------------------------------------------------
##Subset on top variable genes
n.top2plot = 10
topvar.genes = rownames(vargenes_stats[['real.stats']])[1:n.top2plot]
pp2plot = pp_select(pp, topvar.genes)


## ----real-tsne-scatter-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_pp_scatter(pp2plot, log_marks = FALSE, scale_marks = FALSE, pal.direction = -1)

## ----real-run-trendsceek-------------------------------------------------

## ----real-trendsceek, results='hide'-------------------------------------
##set parameters
nrand = 100
ncores = 1

##run
trendstat_list = trendsceek_test(pp2plot, nrand, ncores)

## ------------------------------------------------------------------------
head(trendstat_list[['supstats_wide']])

## ------------------------------------------------------------------------
alpha = 0.05 ##Benjamini-Hochberg
sig_list = extract_sig_genes(trendstat_list, alpha)

## ------------------------------------------------------------------------
lapply(sig_list, nrow)

## ------------------------------------------------------------------------
##get significant genes
sig_genes = sig_list[['markcorr']][, 'gene']

## ----real-markcorr-trendstats-plot, fig.width = 7, fig.height = 7, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_trendstats(trendstat_list, sig_genes)

## ------------------------------------------------------------------------
pp_sig = pp_select(pp, sig_genes)

## ----real-markcorr-scatter-notscaled-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_pp_scatter(pp_sig, log_marks = FALSE, scale_marks = FALSE, pal.direction = -1)

## ----real-markcorr-scatter-scaled-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_pp_scatter(pp_sig, log_marks = FALSE, scale_marks = TRUE, pal.direction = -1)


## ----real-markcorr-density-scaled-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_pp_density(pp_sig, log_marks = FALSE)

## ----real-celltest, results='hide', message=FALSE, warning=FALSE---------
nrand = 100
cellpeaks_siggenes = cellsceek_test(pp_sig, nrand = nrand)
sig_cells = get_sigcells(cellpeaks_siggenes)


## ----real-celltest-sigcells----------------------------------------------
head(sig_cells, 10)

## ----real-sigcells-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_pp_density(pp_sig, log_marks = FALSE, cells2highlight = sig_cells)

## ------------------------------------------------------------------------
##get significant genes
sig_genes = sig_list[['markvario']][, 'gene']

## ----real-markvario-trendstats-plot, fig.width = 7, fig.height = 7, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_trendstats(trendstat_list, sig_genes)

## ------------------------------------------------------------------------
pp_sig = pp_select(pp, sig_genes)

## ----real-markvario-scatter-notscaled-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_pp_scatter(pp_sig, log_marks = FALSE, scale_marks = FALSE, pal.direction = -1)

## ----real-markvario-scatter-scaled-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_pp_scatter(pp_sig, log_marks = FALSE, scale_marks = TRUE, pal.direction = -1)


## ----real-markvario-density-scaled-plot, fig.width = 5, fig.height = 5, include = TRUE, echo = TRUE, fig.align = 'left', fig.cap = ''----
plot_pp_density(pp_sig, log_marks = FALSE)

