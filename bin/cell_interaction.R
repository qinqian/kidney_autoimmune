library(furrr)
library(purrr)
library(future)
library(Seurat)
library(dplyr)
library(arrow)
library(SeuratObject)
library(presto)
library(ggplot2)
library(patchwork)


coloc_one_type = function(index_type, adj, y, nperm = 100, max_dist=30, compartments=NULL, verbose=TRUE) {
    if (verbose) message(index_type)
    types = unique(y)
    i_index = which(y == index_type)
    i_shuffle = setdiff(seq_len(length(y)), i_index)

    X = adj[i_index, ] %*% Matrix::sparse.model.matrix(~0+y) %>% as.matrix()
    colnames(X) = gsub('^y', '', colnames(X))
    freq = (colSums(X) / nrow(X))[types]

    freq_perm = map(seq_len(nperm), function(i) {
        set.seed(i)
        yperm = y
        if (is.null(compartments)) {
            yperm[i_shuffle] = sample(y[i_shuffle])
        } else {
            ## shuffle inside compartments, to preserve total composition within compartment
            .x = split(i_shuffle, compartments[i_shuffle]) %>%
                map(function(.i) {
                    ## CAUTION: if .i is a single number, sample will interpret it as 1:.i
                    if (length(.i) == 1) {
                        res = .i
                    } else {
                        res = sample(.i) ## shuffle non-index cells inside hub            
                    }
                   
                    names(res) = .i
                    return(res)
                }) %>%
                reduce(c)
            yperm[as.integer(names(.x))] <- y[.x]
        }

        X = adj[i_index, ] %*% Matrix::sparse.model.matrix(~0+yperm) %>% as.matrix() #%>% prop.table(1)
        colnames(X) = gsub('^yperm', '', colnames(X))
        (colSums(X) / nrow(X))[types]    
    }) %>%
        purrr::reduce(rbind2)

    stats = tibble(
        type = types,
        freq,
        zscore = (freq - apply(freq_perm, 2, mean)) / apply(freq_perm, 2, sd),
        pval = exp(pnorm(-zscore, log.p = TRUE, lower.tail = TRUE)), ## one-tailed
        fdr = p.adjust(pval)
    ) %>%
        cbind(dplyr::rename(data.frame(t(apply(freq_perm, 2, quantile, c(.025, .975)))), q025 = `X2.5.`, q975 = `X97.5.`)) %>% ## 95% CI
        subset(type != index_type) %>%
        dplyr::mutate(index_type = index_type) %>%
        dplyr::select(index_type, type, everything()) %>%
        arrange(fdr)

    return(stats)    
}


coloc_all_types = function(index_types, coords, y, nperm = 100, nsteps=1, max_dist=30, compartments=NULL, parallel=TRUE, verbose=TRUE) {
    if (parallel & length(index_types) > 1) {
        plan(multicore)
    } else {
        plan(sequential)
    }

    ## Define neighbors
    ## NOTE: max_dist only refers to directly adjacent neighbors
    adj = spatula::getSpatialNeighbors(coords, return_weights = TRUE)
    adj@x[adj@x > max_dist] = 0
    adj = Matrix::drop0(adj)
    adj@x = rep(1, length(adj@x))
   
    ## If nsteps>1, consider not only your adjacent neighbors
    ##   but also your neighbor's neighbors etc.
    if (nsteps > 1) {
        adj = adj + Matrix::Diagonal(n = nrow(adj)) ## add self
        for (iter in seq_len(nsteps - 1)) {
            adj = adj %*% adj
        }
        ## Ignore weights. Only care if cell is a neighbor or not
        adj@x = rep(1, length(adj@x))
       
        ## Remove self as neighbor
        adj = adj - Matrix::Diagonal(n = nrow(adj))
        adj = Matrix::drop0(adj)
    }
       
    index_types %>%
        future_map(coloc_one_type, adj, y, nperm, max_dist, compartments, verbose, .options = furrr::furrr_options(seed = 1)) %>%
        rbindlist()  
}

