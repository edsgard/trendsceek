# trendsceek
Identify genes with spatial expression trends in single-cell gene-expression data 

## Installation
First, install the package dependencies which are available on bioconductor but
not on CRAN:
```R
source("http://www.bioconductor.org/biocLite.R")
deps = c('BiocParallel')
new_deps = deps[!(deps %in% installed.packages()[,"Package"])]
if(length(new_deps) != 0){biocLite(new_deps)}
```

Installation can then be done via the devtools package:
```R
library('devtools')
devtools::install_github('edsgard/trendsceek')
```

## Tutorial
Once you've installed trendsceek you'll be able to follow the
vignette-tutorial. You can open it by:
```R
vignette('trendsceek')
```

## Minimal example
```R
     library('trendsceek')

     ##create synthetic dataset
     pp = sim_pois(300)
     low_expr = c(10, 10)
     high_expr = c(20, 50)
     pp = add_markdist_hotspot(pp, low_expr, high_expr)

     ##run trendsceek
     trendstat_list = trendsceek_test(pp, nrand = 100, ncores = 1)
     head(trendstat_list[['supstats_wide']])

     ##show significant genes
     sig_list = extract_sig_genes(trendstat_list, alpha = 0.1)
     sig_genes = sig_list[['markcorr']][, 'gene']
     print(sig_genes)
     plot_trendstats(trendstat_list, sig_genes)
     pp_sig = pp_select(pp, sig_genes)
     plot_pp_scatter(pp_sig, log_marks = FALSE, scale_marks = TRUE, pal.direction = -1)

     ##cells located in high-expressing regions of the significant genes
     cellpeaks_siggenes = cellsceek_test(pp_sig)
     sig_cells = get_sigcells(cellpeaks_siggenes)
     plot_pp_density(pp_sig, log_marks = FALSE, cells2highlight = sig_cells)	 
```

## Function reference manual
To get help for specific functions you can use ?fcn, for example:
```R
library('trendsceek')
?trendsceek_test
```

The complete function reference manual for all functions can be found
at "doc/refman.pdf" within the installed library directory (to find
your R library directories you can call
.libPaths() from within R). You can
also view the latest version by:
```R
browseURL('https://github.com/edsgard/trendsceek/tree/master/inst/doc/refman.pdf')
```

## Citation
If you use trendsceek, please cite it as follows:

Edsgärd D. and Sandberg R., Identification of spatial expression
trends in single-cell gene expression data, <em>Nature Methods</em>, 2017<br>
