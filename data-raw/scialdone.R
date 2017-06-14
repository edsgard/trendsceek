

cloud.dir = '/Volumes/Data/cloud/btsync/work/rspd'
data.dir = file.path(cloud.dir, 'data/external/scialdone_2016')

counts.rds = file.path(data.dir, 'counts.rds')
meta.rds = file.path(data.dir, 'meta.rds')

main <- function(){

    meta.mat = readRDS(meta.rds)
    counts = readRDS(counts.rds)
    
    ##get 481 E6.5 cells used in Figure 2A
    cluster3.cells = meta.mat[which(meta.mat[, 'embryoStage'] == 'E6.5' & meta.mat[, 'cluster'] == 'turquoise'), 'cellName']
    length(cluster3.cells) ##481    

    ##Subset cells
    counts_filt = counts[, cluster3.cells]
    rownames(meta.mat) = meta.mat[, 'cellName']
    meta.mat = meta.mat[cluster3.cells, ]
    dim(counts_filt) ##481 x 41388
    dim(meta.mat) #481

    ##*###
    ##ensg2mgi mapping
    ##*###
    ensg2mgi.rds = file.path(cloud.dir, 'data/annot/ensembl', 'ensg2mgi.rds')
    ensg2mgi = readRDS(ensg2mgi.rds)
    feats = base::merge(ensg2mgi, as.matrix(rownames(counts_filt)), by.x = 'ensg', by.y = 1, all.y = TRUE)

    ##if ensg2mgi is one-to-many then collapse into one string with multiple mgi
    feats = feats %>% group_by(ensg) %>% summarise(base::paste(gene, collapse = ';')) %>% as.matrix
    colnames(feats) = c('ensg', 'mgi')

    ##set to ensg in cases where mgi == NA
    na_ind = which(feats[, 'mgi'] == 'NA')
    genes = feats[, 'mgi']
    genes[na_ind] = feats[na_ind, 'ensg']
    length(na_ind) #817
    length(grep('^ENSMUSG*', genes)) #817

    feats = cbind(feats, genes)
    colnames(feats) = c('ensg', 'mgi', 'gene')

    ##for non-unique mgis append ensg to the mgi
    feats = as.data.frame(feats, stringsAsFactors = FALSE) %>% group_by(gene) %>% mutate(n_ensg = n()) %>% mutate(new_gene = ifelse(n_ensg >= 2, base::paste(gene, ensg, sep = ';'), gene)) %>% ungroup %>% dplyr::select(ensg, mgi, new_gene) %>% dplyr::rename(gene = new_gene) %>% as.data.frame(stringsAsFactors = FALSE)
    
    ##ensure that ordering is the same
    rownames(feats) = feats[, 'ensg']
    feats = feats[rownames(counts_filt), ]
    rownames(counts_filt) = feats[, 'gene']
    
    ##Dump
    scialdone = list(counts = counts_filt, meta = meta.mat)
    devtools::use_data(scialdone)
    
}
