#'
#' This function performs pseudobulk analysis on scRNAseq data by random sampling
#'
#' @param sce: SCE object, takes raw counts
#' @param cellIDs: IDs of cells of interest belonging to cluster(s)
#' @param obs: obs/metadata variable to group by - must column name in be in sce colData
#' @param replicates: number of artificial replicates to split each condition into, recommend at least 500 cells per replicate after splitting
#' @param contrasts: DESeq contrast for results - take a vector of 2, eg. c('treatment','control')
#' @param alpha: DESeq alpha parameter (significance level for DE testing)
#'
#' @return: list of Pseudobulk raw counts, DESeq2 object, VST-transformed matrix, DE results, colData design table
#' @export
#' 
runPseudoBulk <- function(sce, 
                          cellIDs = NULL, 
                          obs,
                          contrasts,
                          replicates = NULL,
                          alpha = 0.05
                          ){
    # subset sce object to cluster(s) of interest                          
    if(is.null(cellIDs) == TRUE){
        sce <- sce
    } else if(!is.null(cellIDs) == TRUE){
        sce <- sce[,colnames(sce) %in% cellIDs]
    }

    # create vector of obs to contrast
    groups <- as.character(unique(sce@colData[[obs]]))

    # set replicate number, automatically sets 500 per group if none specified
    if(is.null(replicates) == TRUE){
        reps <- ncol(sce) / (500 * length(groups))
    } else if(!is.null(replicates) == TRUE){
        reps <- replicates
    }

    # indices of cells split by obs to contrast
    cond_index <- vector(mode = "list", length = length(groups))
    for(i in seq_along(groups)){
        cond_index[[i]] <- rownames(sce@colData[sce@colData[[obs]] == groups[i],])
    }

    # subset sce object
    sce_list <- vector(mode = "list", length = length(groups))
    for(j in seq_along(groups)){
        sce_list[[j]] <- sce[,colnames(sce) %in% cond_index[[j]]]
    }

    # create counts matrices and wrangle to create x replicates per obs condition
    cond_list <- vector(mode = "list", length = length(groups))
    for(k in seq_along(groups)){
        cond_list[[k]] <- data.frame(as.matrix(sce_list[[k]]@assays@data$counts),
                                     row.names = rownames(sce_list[[k]]))

        # randomly assign cells to number of artificial replicates specified                           
        colnames(cond_list[[k]]) <- sample(1:reps, ncol(cond_list[[k]]), replace = T)

        # collapse columns assigned to same artifical replicate
        cond_list[[k]] <- sapply(1:reps, 
                                 function(x)  rowSums(cond_list[[k]][names(cond_list[[k]]) %in% x]))
    }

    # stitch counts together
    bulk_counts <- do.call(cbind, cond_list)

    # assign column names
    colnames(bulk_counts) <- 1:ncol(bulk_counts)

    # create colData annotation with pseudobulk samples
    design <- data.frame(row.names = colnames(bulk_counts),
                         contrast = rep(groups, each = reps)) # note this only takes one condition for now..

    # create DESeq object
    dds <- DESeqDataSetFromMatrix(countData = bulk_counts, 
                                  colData = design,
                                  design = ~ contrast)

    # filter genes with low counts
    dds <- dds[rowSums(counts(dds)) > reps*2,] 

    # DEseq normalisation and vst transformation
    dds <- DESeq(dds, betaPrior = TRUE)
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
    countmatrix <- assay(vsd)

    # DE analysis
    DEresult <- results(dds, alpha = alpha, contrast = c('contrast', contrasts[1], contrasts[2]))
    DEresult <- as.data.frame(DEresult)
    DEresult <- DEresult[order(DEresult$padj),]
    DEresult$significance <- as.factor(DEresult$padj < alpha)

    ## RETURN key outputs
    output <- list('Pseudobulk_raw_counts' = bulk_counts,
                   'DESeq_object' = dds, 
                   'Transformed_matrix' = countmatrix,
                   'DE_result' = DEresult,
                   'Design_table' = design)
    return(output)
}
#'
#' 
#' This function iterates 'pseudobulk' DE analysis using DESeq2 pipeline \
#' on scRNAseq data for genes of interest.
#' It is essentially bootstrapping for scRNAseq pseudobulk analysis.
#'
#' @param sce: SCE object, takes raw counts
#' @param cellIDs: IDs of cells of interest belonging to cluster(s)
#' @param obs: obs/metadata variable to group by - must column name in be in sce colData
#' @param replicates: number of artificial replicates to split each condition into, recommend at least 500 cells per replicate after splitting
#' @param contrast: DESeq contrast for results - take a vector of 2, eg. c('treatment','control')
#' @param alpha: significance level 
#' @param genes: vector of genes of interest, if NULL, iterate and test for all
#' @param iterations: number of times to permute random sampling and results function
#' 
#' @return A list of LogFoldChange, adjusted P value, and p value statistics for genes of interest
#' @export
#' 
iteratePseudoBulk <- function(sce, 
                              cellIDs = NULL, 
                              obs,
                              contrasts,
                              replicates = NULL,
                              alpha = 0.05,
                              genes = NULL, 
                              iterations = 10
                              ){
    LFC_container <- vector(mode = 'list', length = iterations)
    Padj_container <- vector(mode = 'list', length = iterations)
    Pval_container <- vector(mode = 'list', length = iterations)

    # run resultsPseudoBulk function for each iteration
    for(x in 1:iterations){
    tmp <- resultsPseudoBulk(sce_int = sce, cellIDs_int = cellIDs, obs_int = obs, 
                             contrasts_int = contrasts, replicates_int = replicates,
                             alpha_int = alpha, genes_int = genes)

    # populate container                         
    LFC_container[[x]] <- as.matrix(tmp$`LFC`)
    Padj_container[[x]] <- as.matrix(tmp$`Padj`)
    Pval_container[[x]] <- as.matrix(tmp$`Pval`)

    # print progress
    counter <- paste('DESeq modelling round', x, 'complete')
    print(counter)
    }

    ## collapse iterated values into matrix for LFC and Padj
    LFC_mat <- do.call(cbind, LFC_container)
    Padj_mat <- do.call(cbind, Padj_container)
    Pval_mat <- do.call(cbind, Pval_container)

    # calculate statistics for LFC (basic penalty implemented for multiple re-sampling)
    LFC_rowMean <- rowMeans(LFC_mat)
    LFC_pvalue <- apply(LFC_mat, 1, function(x){t.test(x, mu = 0, alternative = 'two.sided')$p.value * iterations})
    LFC_ci95 <- t(apply(LFC_mat, 1, function(x){t.test(x, mu = 0, alternative = 'two.sided')$conf.int}))
    colnames(LFC_ci95) <- c('conf.int.lo','conf.int.up')
    LFC_stat <- cbind(LFC_rowMean, LFC_pvalue, LFC_ci95)

    # calculate statistics for Padj
    Padj_rowMean <- rowMeans(Padj_mat)
    FDR <- apply(Padj_mat, 1, function(x){t.test(x, mu = alpha, alternative = 'less')$p.value})
    Padj_ci95 <- t(apply(Padj_mat, 1, function(x){t.test(x, mu = alpha, alternative = 'less')$conf.int}))
    colnames(Padj_ci95) <- c('conf.int.lo','conf.int.up')
    Padj_stat <- cbind(Padj_rowMean, FDR, Padj_ci95)

    # calculate statistics for Pval
    Pval_rowMean <- rowMeans(Pval_mat)
    FDR_pval <- apply(Pval_mat, 1, function(x){t.test(x, mu = alpha, alternative = 'less')$p.value})
    Pval_ci95 <- t(apply(Pval_mat, 1, function(x){t.test(x, mu = alpha, alternative = 'less')$conf.int}))
    colnames(Pval_ci95) <- c('conf.int.lo','conf.int.up')
    Pval_stat <- cbind(Pval_rowMean, FDR_pval, Pval_ci95)

    # RETURN output
    stats <- list('LFC_stat' = LFC_stat,
                  'Padj_stat' = Padj_stat,
                  'Pval_stat' = Pval_stat)
    return(stats)
}
#'
#' This is an internal function used by iteratePseudoBulk
#' @return A list of LogFoldChanges, adjusted P values, and p values for genes of interest
#' @export
#' 
resultsPseudoBulk <- function(sce_int, 
                              cellIDs_int = NULL, 
                              obs_int,
                              contrasts_int,
                              replicates_int = NULL,
                              alpha_int = 0.05,
                              genes_int = NULL
                              ){
    # subset sce object to cluster(s) of interest                          
    if(is.null(cellIDs_int) == TRUE){
        sce <- sce_int
    } else if(!is.null(cellIDs_int) == TRUE){
        sce <- sce_int[,colnames(sce_int) %in% cellIDs_int]
    }

    # create vector of obs to contrast
    groups <- as.character(unique(sce@colData[[obs_int]]))

    # set replicate number, automatically sets 500 per group if none specified
    if(is.null(replicates_int) == TRUE){
        reps <- ncol(sce) / (500 * length(groups))
    } else if(!is.null(replicates_int) == TRUE){
        reps <- replicates_int
    }

    # indices of cells split by obs to contrast
    cond_index <- vector(mode = "list", length = length(groups))
    for(i in seq_along(groups)){
        cond_index[[i]] <- rownames(sce@colData[sce@colData[[obs_int]] == groups[i],])
    }

    # subset sce object
    sce_list <- vector(mode = "list", length = length(groups))
    for(j in seq_along(groups)){
        sce_list[[j]] <- sce[,colnames(sce) %in% cond_index[[j]]]
    }

    # create counts matrices and wrangle to create x replicates per obs condition
    cond_list <- vector(mode = "list", length = length(groups))
    for(k in seq_along(groups)){
        cond_list[[k]] <- data.frame(as.matrix(sce_list[[k]]@assays@data$counts),
                                     row.names = rownames(sce_list[[k]]))

        # randomly assign cells to number of artificial replicates specified                           
        colnames(cond_list[[k]]) <- sample(1:reps, ncol(cond_list[[k]]), replace = T)

        # collapse columns assigned to same artifical replicate
        cond_list[[k]] <- sapply(1:reps, 
                                 function(x)  rowSums(cond_list[[k]][names(cond_list[[k]]) %in% x]))
    }

    # stitch counts together
    bulk_counts <- do.call(cbind, cond_list)

    # assign column names
    colnames(bulk_counts) <- 1:ncol(bulk_counts)

    # create colData annotation with pseudobulk samples
    design <- data.frame(row.names = colnames(bulk_counts),
                         contrast = rep(groups, each = reps)) # note this only takes one condition for now..

    # create DESeq object
    dds <- DESeqDataSetFromMatrix(countData = bulk_counts, 
                                  colData = design,
                                  design = ~ contrast)

    # filter genes with low counts
    dds <- dds[rowSums(counts(dds)) > reps*2,] 

    # DEseq normalisation and vst transformation
    dds <- DESeq(dds, betaPrior = TRUE, quiet = TRUE)

    # DE analysis
    DEresult <- results(dds, alpha = alpha_int, contrast = c('contrast', contrasts_int[1], contrasts_int[2]))
    DEresult <- as.data.frame(DEresult)
    DEresult <- DEresult[order(DEresult$padj),]

    # account for no genes specified
    if(is.null(genes_int) == TRUE){
        genelist <- rownames(DEresult)
    } else if(!is.null(genes_int) == TRUE){
        genelist <- genes_int
    }                             
                                 
    # retrieve L2FC and Padj
    LFC <- DEresult[genelist,c('log2FoldChange')]
    names(LFC) <- genelist                          
    padj <- DEresult[genelist,c('padj')] 
    names(padj) <- genelist
    pval <- DEresult[genelist,c('pvalue')] 
    names(pval) <- genelist
                                 
    ## RETURN output
    res_output <- list('LFC' = LFC,
                       'Padj' = padj,
                       'Pval' = pval)
    return(res_output)
}