




createSeurateObj10x <- function(dataPath = NULL, genome="mm9",sampleName="",minGene=1500, mincells =3){
  if(is.null(dataPath)){
    stop("Please provide the dataPath")
  }
  
  
  fname <- file.path(dataPath,"outs/filtered_gene_bc_matrices/",genome)
  
  if(!dir.exists(fname)){
    stop("Please check the provided path is a valid 10X folder")
  }
  
  
  sobj <- Read10X(fname)
  
  cnames <- colnames(sobj)
  rnames <- rownames(sobj)
  sobj <- as(sobj, "dgCMatrix")
  colnames(sobj) <- cnames
  rownames(sobj) <- rnames
  
  sobj <- CreateSeuratObject(raw.data = sobj,
                             scale.factor = 1e6,
                             min.cells = mincells,
                             min.genes = minGene, 
                             project = sampleName, 
                             do.scale = F, 
                             do.center = F, 
                             names.field = 2,
                             names.delim = "\\-")
  
  
  sobj <- SetIdent(sobj,ident.use = sampleName)
  mito.genes <- grep("^mt-", rownames(sobj@data), value = T)
  
  if(length(mito.genes) > 0){
    percent.mito <- colSums(as.matrix(sobj@raw.data[mito.genes, colnames(sobj@data)])) / colSums(as.matrix(sobj@raw.data[, colnames(sobj@data)]))
    sobj <- AddMetaData(sobj, percent.mito, "percent.mito")
  } 
  
  sobj <- SubsetData(sobj,subset.name = "percent.mito",accept.high = 0.1)
  sobj
}



removeMitoandRiboGenes <- function(sobj){
  
  pos_mito <- grep("^MT-",rownames(sobj@data),value = T,ignore.case = T)
  pos_ribo <- grep("^RPS",rownames(sobj@data),value = T,ignore.case = T)
  pos_ribo2 <- grep("^RPL",rownames(sobj@data),value = T,ignore.case = T)
  pos_hb <- grep("^HB.-.*",rownames(sobj@data),value = T,ignore.case = T)
  
  sobj@data <- sobj@data[!rownames(sobj@data) %in% c(pos_mito,pos_ribo,pos_ribo2,pos_hb),]
  
  #sobj@raw.data <- sobj@raw.data[, colnames(sobj@data)]
  
  # remove genes not expressed
  nbExpr <- rowSums(as.matrix(sobj@data))
  
  sobj@data <- sobj@data[nbExpr >0,]
  
  #sobj@data <- LogNormalize(sobj@raw.data,scale.factor = 1e6)
  
  sobj
}


removeMitoandRiboGenes_raw <- function(sobj){
  
  pos_mito <- grep("^MT-",rownames(sobj@raw.data),value = T,ignore.case = T)
  pos_ribo <- grep("^RPS",rownames(sobj@raw.data),value = T,ignore.case = T)
  pos_ribo2 <- grep("^RPL",rownames(sobj@raw.data),value = T,ignore.case = T)
  pos_hb <- grep("^HB.-.*",rownames(sobj@raw.data),value = T,ignore.case = T)
  
  sobj@raw.data <- sobj@raw.data[!rownames(sobj@raw.data) %in% c(pos_mito,pos_ribo,pos_ribo2,pos_hb),]
  
  #sobj@raw.data <- sobj@raw.data[, colnames(sobj@data)]
  
  # remove genes not expressed
  nbExpr <- rowSums(as.matrix(sobj@raw.data))
  
  sobj@raw.data <- sobj@raw.data[nbExpr >0,]
  
  #sobj@data <- LogNormalize(sobj@raw.data,scale.factor = 1e6)
  
  sobj
}


SeuratToSummarizedExperiment <- function(sobj){
  require(SingleCellExperiment)
  
  counts <- as.matrix(sobj@raw.data[rownames(sobj@data), colnames(sobj@data)])
  
  rData <- data.frame(gene_names=rownames(counts))
  cData <- sobj@meta.data
  
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    rowData = rData,
    colData = cData
  )
  
  return(sce)
}






#' Title
#'
#' @param genesSymbols
#' @param species
#'
#' @return
#' @export
#'
#' @examples
GeneNameToEntrez <- function(genesSymbols, species = "mouse") {
  if (species  == "human") {
    require(org.Hs.eg.db)
    xx <- as.list(org.Hs.egALIAS2EG)
    xx <- xx[!is.na(xx)]
  }
  else{
    if (species  == "mouse") {
      require(org.Mm.eg.db)
      xx <- as.list(org.Mm.egALIAS2EG)
      xx <- xx[!is.na(xx)]
    }
  }
  
  pos <- match(genesSymbols, names(xx))
  genesSymbols <- genesSymbols[!is.na(pos)]
  df <-
    data.frame(
      ensemblID = genesSymbols,
      entrezID = "",
      stringsAsFactors = FALSE
    )
  
  for (i in 1:nrow(df)) {
    df$entrezID[i] <- xx[[genesSymbols[i]]][1]
  }
  return(df)
}


#' Title
#'
#' @param deResult
#'
#' @return
#' @export
#'
#' @examples
doGOAnalysis <- function(genesID,keyType = "ENTREZID") {
  require(clusterProfiler)
  go_res <- list()
  for (trm in c("BP")){#, "CC", "MF")) {
    go_res[[trm]] <-
      enrichGO(genesID,keyType = keyType,
               'org.Mm.eg.db',
               ont = trm,minGSSize = 3,
               pvalueCutoff = 1)
  }
  return(go_res)
}


#' Title
#'
#' @param DEgenes
#'
#' @return
#' @export
#'
#' @examples
doGOAnalysis.vector <- function(DEgenes) {
  treatment_DE <- unique(unlist(DEgenes))
  ## Convert gene names to EntrezID
  treatment_DE.entrez <- GeneNameToEntrez(treatment_DE)
  treatment_DE.entrez <-
    subset(treatment_DE.entrez, entrezID != "NA")
  
  
  go_res <- list()
  for (trm in c("BP", "CC", "MF")) {
    go_res[[trm]] <-
      enrichGO(
        treatment_DE.entrez$entrezID,
        'org.Mm.eg.db',
        ont = trm,
        pvalueCutoff = 0.05,
        readable = TRUE
      )
  }
  return(go_res)
}

#' Title
#'
#' @param go_res
#' @param top
#' @param col
#' @param title
#'
#' @return
#' @export
#'
#' @examples
plotGoBars <-
  function(go_res,
           top = 20,
           col = "red",
           title = "top 20 enriched GO") {
    tmp <- as.data.frame(go_res)
    
    if (nrow(tmp) > 0) {
      tmp$p.adjust <- as.double(as.character(tmp$pvalue))
      
      ord <- order(tmp$pvalue)
      tmp <- head(tmp[ord,], top)
      
      ord <- order(tmp$pvalue, decreasing = T)
      tmp <- tmp[ord,]
      tmp$Description <- factor(tmp$Description, levels = tmp$Description)
      
      
      ggplot(tmp, aes(x = Description, y = -log10(pvalue))) + 
        geom_bar(stat = "identity", fill = col) +
        coord_flip() +
        theme_bw() + ggtitle(title)
    }
  }



#' Title
#'
#' @param sobj
#' @param ngenes
#' @param ttl
#' @param rowName
#' @param pc.use
#' @param disp.min
#' @param disp.max
#'
#' @return
#' @export
#'
#' @examples
plotPCtopGenes <-
  function(sobj,
           ngenes,
           ttl,
           rowName = FALSE,
           pc.use = 1,
           disp.min = -3,
           disp.max = 3) {
    #data.use=PCTopCells(sobj, pc.use = pc.use, cells.use = NULL, do.balanced = TRUE, label.columns = F, num.genes = ngenes, do.return = TRUE)
    cells.use = sobj@cell.names
    genes.use = rev(PCTopGenes(sobj, pc.use, ngenes, do.balanced = TRUE))
    
    if("pca.rot" %in% slotNames(sobj)){
      cells.ordered = cells.use[order(sobj@pca.rot[cells.use, pc.use])]  
    }else{
      cells.ordered = cells.use[order(sobj@dr$pca@cell.embeddings[cells.use, pc.use])]  
    }
    
    data.use = sobj@scale.data[genes.use, cells.ordered]
    data.use = MinMax(data.use, min = disp.min, max = disp.max)
    
    ## Prepare for plotting
    
    col_anno <- data.frame(CellType = sobj@ident)
    
    
    if (length(levels(sobj@ident)) > 20) {
      TypesCols = ggsci::pal_d3(palette = "category20c")(20)
      
      n <- length(levels(sobj@ident)) - 20
      TypesCols <- c(TypesCols, ggsci::pal_rickandmorty()(n))
    } else{
      TypesCols = ggsci::pal_d3(palette = "category20c")(length(levels(sobj@ident)))
    }
    
    
    names(TypesCols) <- levels(sobj@ident)
    ann_colors = list(CellType = TypesCols)
    
    pheatmap(
      data.use,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = rowName,
      show_colnames = FALSE,
      annotation_colors = ann_colors,
      annotation_col = col_anno,
      fontsize = 6,
      main = ttl
    )
  }



plotPCSigGenes <-
  function(sobj,
           ngenes,
           ttl,
           rowName = FALSE,
           pc.use = 1,
           disp.min = -3,
           disp.max = 3,
           pval = 0.05) {
    #data.use=PCTopCells(sobj, pc.use = pc.use, cells.use = NULL, do.balanced = TRUE, label.columns = F, num.genes = ngenes, do.return = TRUE)
    cells.use = sobj@cell.names
    genes.use = PCASigGenes(
      object = sobj,
      pcs.use = pc.use,
      pval.cut = pval,
      use.full = F,
      max.per.pc = ngenes
    )
    
    if (length(genes.use) == 0)
      return(NULL)
    cells.ordered = cells.use[order(sobj@dr$pca@cell.embeddings[cells.use, pc.use])]
    data.use = sobj@scale.data[genes.use, cells.ordered]
    data.use = MinMax(data.use, min = disp.min, max = disp.max)
    
    ## Prepare for plotting
    
    col_anno <- data.frame(CellType = sobj@ident)
    
    TypesCols <-
      brewer.pal(8, "Dark2")[1:length(levels(sobj@ident))]
    names(TypesCols) <- levels(sobj@ident)
    ann_colors = list(CellType = TypesCols)
    
    pheatmap(
      data.use,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = rowName,
      show_colnames = FALSE,
      annotation_colors = ann_colors,
      annotation_col = col_anno,
      fontsize = 6,
      main = ttl
    )
  }



#' applies one of the Seurat implemented algorithms on the data
#'
#' @param sobj Seurat object
#' @param algo The algorithm to use 1,2 or 3
#' @param res The clustering resolution (<1 less clusters, >1 more clusters)
#' @param pcs The PCs to use
#'
#' @return A Seurat object in which the column "algo_#_res_#" is added to the data.info data frame.
#' @export
#'
#' @examples
doSeuratClustering <- function(sobj,
                               algo = 2,
                               res = 1,
                               pcs = 1:15,
                               k = 10) {
  sobj <-
    FindClusters(
      sobj,
      dims.use  = pcs,
      resolution = res,
      k.param = k,
      print.output = 0,
      algorithm = algo,
      save.SNN = FALSE
    )
  
  sobj <-
    BuildClusterTree(
      sobj,
      reorder.numeric = TRUE,
      do.reorder = TRUE,
      do.plot = F,
      pcs.use = pcs
    )
  col <- paste0("res.", res)
  sobj@data.info[, col] <- NULL
  col <- paste0("algo_", algo, "_res_", res)
  cnames <- c(colnames(sobj@data.info), col)
  sobj@data.info <- cbind(sobj@data.info, sobj@data.info$tree.ident)
  colnames(sobj@data.info) <- cnames
  sobj@data.info$tree.ident <- NULL
  return(sobj)
}


doKmeansClustering <- function(sobj,
                               pcs = 1:10,
                               max.nc = 20) {
  nb <-
    NbClust(
      sobj@pca.rot[, pcs],
      diss = NULL,
      distance = "euclidean",
      min.nc = 2,
      max.nc = 15,
      method = "kmeans",
      index = "all",
      alphaBeale = 0.1
    )
  
  nbclus <-  median(nb$Best.nc[1,])
  
  km <- kmeans(sobj@pca.rot[, pcs], centers = nbclus, nstart = 20)
  
  sobj2 <- SetIdent(sobj, ident.use = km$cluster)
  message("Ordering")
  sobj2 <- OrderClusters(sobj2)
  
  cname <- paste0("algo_kmeans_", space, "_", criteria)
  cname <- c(colnames(sobj@data.info), cname)
  sobj@data.info <- cbind(sobj@data.info, sobj2@ident)
  rm(sobj2)
  colnames(sobj@data.info) <- cname
  
  return(sobj)
}



#' Clustering using the Mclust method
#'
#' @param sobj Seurat object
#' @param space which projection space to use. Supports "pca" or "tSNE"
#' @param nClust the number of clustets to test for. Default c(1:40)
#' @param method by default is NULL. The method is detected autmatically using the BIC method.
#' The methods supported are the one supported in the mclust package.
#' @param plot.BIC Plot the BIC estimation plot or not
#' @param pcs.touse is "pca" is selected in the method parameter, what are the PCs to use.
#'
#' @return A Seurat object with the `data.info` table update with the clustering results.
#' @export
#'
#' @examples
doMclust <- function(sobj,
                     space = "pca",
                     nClust = 2:40,
                     method = NULL,
                     plot.criteria = FALSE,
                     pcs.touse = 1:3,
                     criteria = c("BIC", "ICL")) {
  require(mclust)
  require(stringr)
  if (space == "pca") {
    coords <- sobj@pca.rot[, pcs.touse]
  }
  
  if (space == "tSNE") {
    coords <- sobj@tsne.rot
  }
  
  if (criteria == "BIC") {
    if (is.null(method)) {
      message("Calculating BIC")
      mclust_clustersBIC <- mclustBIC(coords, G = nClust)
      print(summary(mclust_clustersBIC))
      if (plot.criteria) {
        plot(mclust_clustersBIC,
             legendArgs = list(
               x = "topright",
               ncol = 2,
               cex = 0.5,
               inset = c(0.1, 0)
             ))
      }
    }
    
    message("Clustering")
    mclust_clusters <- Mclust(coords, x = mclust_clustersBIC)
    
  } else {
    if (criteria == "ICL") {
      message("Calculating ICL")
      mclust_clustersICL <- mclustICL(coords, G = nClust)
      print(summary(mclust_clustersICL))
      if (plot.criteria) {
        plot(mclust_clustersICL,
             legendArgs = list(
               x = "topright",
               ncol = 2,
               cex = 0.5,
               inset = c(0.1, 0)
             ))
      }
    }
    top3Models <- pickBIC(mclust_clustersICL)
    top3Models <-
      str_match(string = names(top3Models), pattern = "(.+),(.+)")
    minpos <- which(top3Models[, 3] == min(top3Models[, 3]))
    message("Clustering")
    mclust_clusters <-
      Mclust(coords, G = top3Models[minpos, 3], modelNames = top3Models[minpos, 2])
    
  }
  
  #plot(mclust_clusters, what="classification", col=distinctColorPalette(mclust_clusters$G))
  #sm <- summary(mclust_clusters, parameters = TRUE)
  
  sobj2 <-
    SetIdent(sobj, ident.use = mclust_clusters$classification)
  message("Ordering")
  sobj2 <- OrderClusters(sobj2)
  
  cname <- paste0("algo_mclust_", space, "_", criteria)
  cname <- c(colnames(sobj@data.info), cname)
  sobj@data.info <- cbind(sobj@data.info, sobj2@ident)
  rm(sobj2)
  colnames(sobj@data.info) <- cname
  
  return(sobj)
}



calcl_kPCA <- function(mat){
  
  require(distances)
  require(RSpectra)
  
  
  K <- distances(mat)
  K <- distance_matrix(K)^2
  gamma <- 1/max(K)
  K <- exp(-gamma * K)
  K <- as.matrix(K)
  
  ## Normalize K
  
  nr <- dim(K)[1]
  A <- matrix(1/nr, nr, nr) 
  Ktild <- K - A %*% K- K %*% A + A %*% K %*% A          
  
  rm(A)
  
  pca <- eigs_sym(A = Ktild,30,which = "LM")
  
  rm(Ktild)
  transform <- K %*% pca$vectors
  res <- list(KPCs=transform)
  return(res)
}



doRaceID <- function(sobj, k = NULL) {
  source("D:/tools/RaceID-master/RaceID_class.R")
  
  sc <- SCseq(as.data.frame(sobj@scale.data))
  
  
  sc <-
    clustexp(
      sc,
      metric = "pearson",
      cln = 0,
      do.gap = TRUE,
      clustnr = 20,
      B.gap = 50,
      SE.method = "Tibs2001SEmax",
      SE.factor = .25,
      bootnr = 50,
      rseed = 17000
    )
  
  cname <- paste0("algo_RaceID")
  cname <- c(colnames(sobj@data.info), cname)
  sobj@data.info <- cbind(sobj@data.info, sc@kmeans$kpart)
  colnames(sobj@data.info) <- cname
  
  return(sobj)
}


doDBScan <- function(sobj, pcs.touse = 1:15) {
  dbscan_clusters <-
    dbscan(sobj@pca.rot[, pcs.touse], minPts = 10, eps = 9)
}



calcCSPA <- function(sobj) {
  algos <- grep("algo", colnames(sobj@data.info))
  if (length(algos) == 0) {
    stop("Please run some clustering algorithms than call this method")
  }
  
  CSPA <-
    matrix(0,
           nrow = nrow(sobj@data.info),
           ncol = nrow(sobj@data.info))
  
  for (a in algos) {
    clusRes <- sobj@data.info[, a]
    m <-
      apply(matrix(clusRes), 1, function(i) {
        as.numeric(i == clusRes)
      })
    CSPA <- CSPA + m
  }
  
  colnames(CSPA) <- rownames(sobj@data.info)
  rownames(CSPA) <- rownames(sobj@data.info)
  CSPA <- CSPA / length(algos)
  return(CSPA)
}






#' Title
#'
#' @param sobj
#' @param ngenes
#' @param ttl
#' @param rowName
#' @param pc.use
#'
#' @return
#' @export
#'
#' @examples
GetRiboinPCtopGenes <-
  function(sobj,
           ngenes,
           ttl,
           rowName = FALSE,
           pc.use = 1) {
    genes.use = rev(PCTopGenes(sobj, pc.use, ngenes))
    
    ribo <- grep("^Rp", genes.use, ignore.case = TRUE)
    return(length(ribo))
  }



#' Title
#'
#' @param sobj
#' @param ident1
#' @param perc.expr
#' @param ident2
#'
#' @return
#' @export
#'
#' @examples
getMASTmarkers <-
  function(sobj,
           ident1,
           perc.expr = 0.05,
           ident2 = NULL) {
    require(MAST)
    ## Change all the other labels to two
    clusters <- as.character(sobj@ident)
    pos1 <- which(clusters == ident1)
    if (is.null(ident2)) {
      pos2 <- which(clusters != ident1)
    }
    else{
      pos2 <- which(clusters == ident2)
    }
    
    clusters[pos1] = "Treatment"
    clusters[pos2] = "Control"
    clusters <- factor(clusters)
    #thr1 <- max(ceiling(length(pos1) * perc.expr),3)
    #thr2 <- max(ceiling(length(pos2) * perc.expr),3)
    
    
    
    
    #nbExpr_Treatment <- apply(sobj@scale.data[,pos1],1,function(x) sum(x>0))
    #nbExpr_control <- apply(sobj@scale.data[,pos2],1,function(x) sum(x>0))
    
    #pos <- which( nbExpr_Treatment >= thr1 | nbExpr_control >= thr2 )
    #sobj@scale.data <- sobj@scale.data[pos,]
    
    geneToUse <- rownames(sobj@data)
    cellToUse <- colnames(sobj@data)
    
    
    if ("nGene" %in% colnames(sobj@data.info)) {
      nGeneOn <- sobj@data.info$nGene
    } else{
      mat <- sobj@raw.data[geneToUse, cellToUse]
      mat <- c@raw.data[geneToUse, cellToUse]
    }
    
    
    colMetaData <-
      data.frame(Condition = clusters,
                 ncells = 1,
                 nGeneOn = nGeneOn)
    colMetaData$cngeneson <-
      scale(colMetaData$nGeneOn, center = TRUE, scale = F)
    
    # if(min(as.matrix(sobj@scale.data)) < 0){
    #   mat <- as.matrix(sobj@scale.data)
    # }else{
    #   mat <-  as.matrix(log2(expm1(sobj@data)+1))  #as(mat, "dgCMatrix")
    # }
    
    mat <-  as.matrix(log2(expm1(sobj@data) + 1))
    #mat <- Log2Norm(data = mat,scale_factor = 1e6)
    #colnames(mat) <- cellToUse
    #rownames(mat) <- geneToUse
    rm(sobj)
    gc()
    
    mastDat <- FromMatrix(mat, cData = colMetaData)
    
    #mastDat <- mastDat[which(freq(mastDat)> perc.expr),]
    
    zlmCond <- zlm.SingleCellAssay(~ Condition, mastDat)
    coefAndCI <- summary(zlmCond, logFC = TRUE)
    coefAndCI <- coefAndCI$datatable
    coefAndCI <-
      coefAndCI[contrast != '(Intercept)' & component == "logFC",]
    coefAndCI[, contrast := abbreviate(contrast)]
    
    zlm.lr <- lrTest(zlmCond, 'Condition')
    zlm.lr_pvalue <- melt(zlm.lr[, , 'Pr(>Chisq)'])
    zlm.lr_pvalue <-
      zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]
    setorder(zlm.lr_pvalue, value)
    
    pos <- match(zlm.lr_pvalue$primerid, coefAndCI$primerid)
    
    res <-
      cbind(as.character(zlm.lr_pvalue$primerid),
            zlm.lr_pvalue$value,
            coefAndCI$ci.hi[pos])
    colnames(res) <- c("Gene", "pvalue", "log2FC")
    res <- as.data.frame(res)
    return(res)
  }




my.MASTDETest <- function(
  object,
  cells.1,
  cells.2,
  min.cells = 3,
  genes.use = NULL,
  latent.vars = NULL,
  assay.type = "RNA",
  min.freq = 0.2,
  ...
) {
  
  require(data.table)
  
  # Check for MAST
  if (!'MAST' %in% rownames(x = installed.packages())) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }
  data.test <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  
  
  tmp_sobj <- SubsetData(object, cells.use = c(cells.1, cells.2))
  newIdent <- c(rep("Group1",length(cells.1)), rep("Group2",length(cells.2)))
  tmp_sobj <- SetIdent(tmp_sobj, ident.use = newIdent)
  
  ae <- AverageDetectionRate(tmp_sobj)
  ae <- subset(ae, Group1 > min.freq | Group2 > min.freq)
  
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = ae))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(ae)]
  
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  
  # if (length(x = latent.vars) > 0) {
  #   my.latent <- scale(x = my.latent)
  # }
  
  coldata <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  coldata[cells.1, "group"] <- "Group1"
  coldata[cells.2, "group"] <- "Group2"
  coldata$group <- factor(x = coldata$group)
  coldata$wellKey <- rownames(x = coldata)
  latent.vars <- c("condition", latent.vars)
  countdata.test <- log2(expm1(data.test[genes.use, rownames(x = coldata)])+1)
  fdat <- data.frame(rownames(x = countdata.test))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(x = countdata.test),
    cData = coldata,
    fData = fdat
  )
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond
  fmla <- as.formula(
    object = paste0(" ~ ", paste(latent.vars, collapse = "+"))
  )
  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  summaryCond <- summary(object = zlmCond, doLRT = 'conditionGroup2')
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(
    summaryDt[contrast=='conditionGroup2' & component=='H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
    summaryDt[contrast=='conditionGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid'
  ) #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'BH')]
  
  colnames(fcHurdle) <- c("GeneName","pvalue","logFC","logFC.highEstimate","logFC.lowEstimate","FDR")
  
  fcHurdle <- fcHurdle[,c("GeneName","logFC","logFC.highEstimate","logFC.lowEstimate","pvalue","FDR")]
  
  
  fcHurdle <- setorder(fcHurdle,"FDR")
  # p_val <- summaryDt[summaryDt[, "component"] == "H", 4]
  # genes.return <- summaryDt[summaryDt[, "component"] == "H", 1]
  # p_val <- subset(summaryDt, component == "H")[, 4]
  # genes.return <- subset(summaryDt, component == "H")[, 1]
  #to.return <- data.frame(p_val, row.names = genes.return)
  
  
  nbCells_cond1 <- apply(countdata.test[fcHurdle$GeneName,cells.1],1,function(x) {sum(x>0)})
  nbCells_cond2 <- apply(countdata.test[fcHurdle$GeneName,cells.2],1,function(x) {sum(x>0)})
  
  
  fcHurdle$Cond1_totalCells <- length(cells.1)
  fcHurdle$Cond2_totalCells <- length(cells.2)
  
  fcHurdle$Cond1_ExpressedCells <- nbCells_cond1
  fcHurdle$Cond2_ExpressedCells <- nbCells_cond2
  
  data.1 = apply(as.matrix(countdata.test[fcHurdle$GeneName, cells.1, drop = F]), 1, ExpMean)
  data.2 = apply(as.matrix(countdata.test[fcHurdle$GeneName, cells.2, drop = F]), 1, ExpMean)
  
  log2FC = (data.2 -data.1)
  
  
  fcHurdle$empirical_logFC <-log2FC
  
  return(fcHurdle)
}


doGOForGroup <- function(sobj, clus, period){
  

  message(paste("cluster",clus,"period",period))
  smp_saline <- rownames(subset(sobj@meta.data, Period == period & treatment == "Saline" & tree.ident== clus))
  smp_cocaine <- rownames(subset(sobj@meta.data, Period == period & treatment == "Cocaine" & tree.ident== clus))
  
  
  if(length(smp_saline)>1 & length(smp_cocaine)>1){
    res = gseaEnrichment.MAST(object =  sobj,
                                                    cells.1 = smp_saline,
                                                    cells.2 =  smp_cocaine,
                                                    # we already regressed out the UMI bias
                                                    latent.vars = NULL 
                              )
  }else{
    return(NULL)
  }
  
  return(res)
  
}


gseaEnrichment.MAST <- function(
  object,
  cells.1,
  cells.2,
  min.cells = 3,
  genes.use = NULL,
  latent.vars = NULL,
  assay.type = "RNA",
  ...
){
  require(data.table)
  
  # Check for MAST
  if (!'MAST' %in% rownames(x = installed.packages())) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }
  data.test <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.test))
  # check that the gene made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(x = data.test)]
  my.latent <- FetchData(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.raw = TRUE
  )
  
  # if (length(x = latent.vars) > 0) {
  #   my.latent <- scale(x = my.latent)
  # }
  
  coldata <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  coldata[cells.1, "group"] <- "Group1"
  coldata[cells.2, "group"] <- "Group2"
  coldata$group <- factor(x = coldata$group)
  coldata$wellKey <- rownames(x = coldata)
  latent.vars <- c("condition", latent.vars)
  countdata.test <- data.test[genes.use, rownames(x = coldata)]
  fdat <- data.frame(rownames(x = countdata.test))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(x = countdata.test),
    cData = coldata,
    fData = fdat
  )
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond
  fmla <- as.formula(
    object = paste0(" ~ ", paste(latent.vars, collapse = "+"))
  )
  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  
  
  boots <- bootVcov1(zlmCond, R=100)
}

##

getMASTmarkers.group <-
  function(sobj, pos1, pos2, perc.expr = 0.3) {
    require(MAST)
    require(data.table)
    require(reshape2)
    
    clusters <- rep("", length(pos1) + length(pos2))
    
    clusters[1:length(pos1)] = "Treatment"
    clusters[(length(pos1) + 1):length(clusters)] = "Control"
    clusters <- factor(clusters)
    
    #thr1 <- ceiling(length(pos1) * perc.expr)
    #thr2 <- ceiling(length(pos2) * perc.expr)
    
    
    #pos <- which( nbExpr_Treatment >= thr1 | nbExpr_control >= thr2 )
    #sobj@scale.data <- sobj@scale.data[pos,]
    #nGeneOn <- apply(sobj@data[,c(pos1,pos2)], 2, function(x) length(x[x>0]))
    
    colMetaData <- data.frame(
      Condition = clusters,
      nGene = sobj@meta.data[c(pos1, pos2),]$nGene #,
      #sample=sobj@data.info[c(pos1,pos2),]$Sample
      )
      
      colMetaData$nGene <-
        scale(colMetaData$nGene, center = TRUE, scale = F)
      
      
      expressedGenes <-
        which(apply(sobj@data[, c(pos1, pos2)], 1, function(x)
          sum(x > 0)) > 0)
      
      mat <-
        sobj@data[expressedGenes, c(pos1, pos2)]
      
      data.1 = apply(as.matrix(mat[, pos1, drop = F]), 1, expMean)
      data.2 = apply(as.matrix(mat[, pos2, drop = F]), 1, expMean)
      
      log2FC = (data.1 -data.2)
      
      #nbExpr_Treatment <- apply(mat[,pos1],1,function(x) sum(x>0)/length(pos1))
      #nbExpr_control <- apply(mat[,pos2],1,function(x) sum(x>0)/length(pos2))
      
      
      mat <-  as.matrix(log2(expm1(mat) + 1))
      
      mastDat <- FromMatrix(mat, cData = colMetaData)
      
      zlmCond <- zlm.SingleCellAssay(~ Condition + nGene, mastDat)
      coefAndCI <- summary(zlmCond, logFC = TRUE)
      coefAndCI <- coefAndCI$datatable
      coefAndCI <- coefAndCI[contrast != '(Intercept)' & component == "logFC",]
      coefAndCI[, contrast := abbreviate(contrast)]
      
      zlm.lr <- lrTest(zlmCond, 'Condition')
      zlm.lr_pvalue <- melt(zlm.lr[, , 'Pr(>Chisq)'])
      zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]
      setorder(zlm.lr_pvalue, value)
      
      pos <- match(zlm.lr_pvalue$primerid, coefAndCI$primerid)
      
      res <- cbind(
          as.character(zlm.lr_pvalue$primerid),
          zlm.lr_pvalue$value,
          coefAndCI$ci.hi[pos]
        )
      
      colnames(res) <- c("Gene", "pvalue", "log2FC")
      res <- as.data.frame(res)
      res$Gene <- as.character(res$Gene)
      res$pvalue <- as.double(as.character(res$pvalue))
      res$log2FC <- as.double(as.character(res$log2FC))
      res$emp_log2FC <- log2FC[as.character(res$Gene)]
      return(res)
  }


getDEgenes <- function(sobj, pos1, pos2, perc.expr = 0.3) {
  require(MAST)
  require(data.table)
  require(reshape2)
  
  
  
  subSobj <- SubsetData(sobj, cells.use = c(pos1, pos2))
  
  newIdent <- rep("", length(pos1) + length(pos2))
  
  
  p1 <- match(pos1, colnames(subSobj@data))
  p2 <- match(pos2, colnames(subSobj@data))
  
  newIdent[p1] <- "treatment"
  newIdent[p2] <- "control"
  
  
  subSobj <- SetIdent(subSobj, ident.use = newIdent)
  
  
  mrk <- FindMarkers(subSobj, ident.1 = "treatment")
  
  return(mrk)
  
}


plot_ggjoy <-
  function(markers,
           clusters,
           expr,
           clustersColors,
           title,
           plotType = "boxplot") {
    ## order cells by cluster
    
    orderedData <- matrix(0, nrow = length(markers), ncol = 0)
    
    clusters <- c()
    
    for (cls in levels(sobj@ident)) {
      pos <- which(sobj@ident == cls)
      clusCells <- sobj@data[markers, pos]
      
      ph <-
        pheatmap(
          clusCells,
          cluster_rows = F,
          cluster_cols = T,
          silent = T
        )
      
      clusCells <- clusCells[, ph$tree_col$order]
      
      orderedData <- cbind(orderedData, clusCells)
      clusters <- c(clusters, rep(cls, ncol(clusCells)))
    }
    
    ph <-
      pheatmap(
        clusCells,
        cluster_rows = T,
        cluster_cols = F,
        silent = T
      )
    
    clusCells <- clusCells[ph$tree_row$order,]
    
    colnames(clusCells) <- 1:ncol(clusCells)
    
    clusCells.mlt <- melt(clusCells)
    
    clusCells.mlt$cluster <- rep(clusters, each = length(markers))
    
    
    
  }


## Plot violin plots for marker genes per cluster
plotViolinPerGene <-
  function(markers,
           clusters,
           expr,
           clustersColors,
           title,
           plotType = "boxplot") {
    toUse <- expr[markers,]
    
    df <-
      data.frame(
        gene = factor(rep(markers, ncol(expr)), levels = markers) ,
        expr = unlist(as.list(expr[markers,])),
        cluster = rep(clusters, each = length(markers))
      )
    
    
    if (plotType == "boxplot") {
      g <-
        ggplot(df, aes(x = cluster, y = expr)) + geom_boxplot(aes(fill = cluster)) + theme_bw() +
        facet_wrap( ~ gene , nrow = length(markers), scales = "free") +
        scale_fill_manual(values = clustersColors) + ggtitle(title) +
        theme(
          axis.text.y = element_text(face = "bold", size = 10),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 12),
          strip.text.x = element_text(size = 8, face = "bold")
        )
    } else{
      g <-
        ggplot(df, aes(x = cluster, y = expr)) + geom_violin(aes(fill = cluster), scale = "width") + theme_bw() +
        facet_wrap( ~ gene, nrow = length(markers), scales = "free") +
        scale_fill_manual(values = clustersColors) + ggtitle(title) +
        theme(
          axis.text.y = element_text(face = "bold", size = 10),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.title = element_text(face = "bold", size = 12),
          strip.text.x = element_text(size = 8, face = "bold")
        )
    }
    
    g
  }



SelectCellsInCluster <- function(neuro_threeDTsne_pos, clus) {
  tmp <- subset(neuro_threeDTsne_pos, cluster == clus)
  ## Normally the cells in each cluster should be in the same location in brain
  location <- table(tmp$location) / nrow(tmp)
  
  ## if we have more than 2 locations, filter the outlier one.
  if (length(location) > 1) {
    m = min(location)
    if (m < 40) {
      toUse <- names(location)[which(location != m)]
      ## Get only the cells that represent tha majority of this location
      tmp <- subset(tmp, location == toUse)
    }
    else{
      ## Skip this cluster
      tmp <- NULL
    }
  }
  
  return(tmp)
}


getDescendants <- function(tree, node, curr = NULL) {
  if (is.null(curr))
    curr <- vector()
  daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
  curr <- c(curr, daughters)
  w <- which(daughters >= length(tree$tip))
  if (length(w) > 0)
    for (i in 1:length(w))
      curr <- getDescendants(tree, daughters[w[i]], curr)
  return(curr)
}

doReorder <-
  function(object,
           pytho,
           reorder.numeric = TRUE,
           genes.use = NULL,
           pcs.use = FALSE) {
    old.ident.order = sort(unique(object@ident))
    data.tree = object@cluster.tree[[1]]
    all.desc = getDescendants(data.tree, data.tree$Nnode + 2)
    all.desc = old.ident.order[all.desc[all.desc <= (data.tree$Nnode + 1)]]
    object@ident = factor(object@ident, levels = all.desc, ordered = TRUE)
    if (reorder.numeric) {
      object = SetIdent(object, object@cell.names, as.integer(object@ident))
      object@data.info[object@cell.names, "tree.ident"] = as.integer(object@ident)
    }
    object = BuildClusterTree(object,
                              genes.use,
                              pcs.use,
                              do.plot = FALSE,
                              do.reorder = FALSE)
  }



OrderClusters <- function(object) {
  library(ape)
  ## Rename clusters according to their name
  data.avg = AverageExpression(object)
  data.dist <- dist(t(data.avg))
  data.tree = as.phylo(hclust(data.dist))
  object@cluster.tree[[1]] = data.tree
  object <- doReorder(object, data.tree)
  object
}

GetAssayData <-
  function(object,
           assay.type = "RNA",
           slot = "data") {
    if (assay.type == "RNA") {
      if (slot == "raw.data") {
        to.return <- object@raw.data
      } else if (slot == "data") {
        to.return <- object@data
      } else if (slot == "scale.data") {
        if (length(x = object@scale.data) == 0) {
          stop("Object@scale.data has not been set. Run ScaleData() and then retry.")
        }
        to.return <- object@scale.data
      }
      #note that we check for this to avoid a long subset for large matrices if it can be avoided
      if (length(x = object@cell.names) == ncol(to.return)) {
        return(to.return)
      }
      return(to.return[, object@cell.names])
    }
    if (!(assay.type %in% names(object@assay))) {
      stop(paste(assay.type, "data has not been added"))
    }
    if (!(slot %in% slotNames(eval(expr = parse(
      text = paste0("object@assay$", assay.type)
    ))))) {
      stop(paste(slot, "slot doesn't exist"))
    }
    to.return <-
      (eval(expr = parse(
        text = paste0("object@assay$", assay.type, "@", slot)
      )))
    if (length(x = object@cell.names) == ncol(x = to.return)) {
      return(to.return)
    }
    return(to.return[, object@cell.names])
  }



WhichCells_my <- function(object,
                          ident = NULL,
                          ident.remove = NULL,
                          cells.use = NULL,
                          subset.name = NULL,
                          accept.low = -Inf,
                          accept.high = Inf,
                          accept.value = NULL,
                          max.cells.per.ident = Inf,
                          random.seed = 1) {
  set.seed(seed = random.seed)
  cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  ident <- SetIfNull(x = ident, default = unique(x = object@ident))
  ident <- setdiff(x = ident, y = ident.remove)
  if (!all(ident %in% unique(x = object@ident))) {
    bad.idents <- ident[!(ident %in% unique(x = object@ident))]
    stop(paste("Identity :", bad.idents, "not found.   "))
  }
  cells.to.use <- character()
  for (id in ident) {
    cells.in.ident <- object@ident[cells.use]
    cells.in.ident <-
      names(x = cells.in.ident[cells.in.ident == id])
    cells.in.ident <- cells.in.ident[!is.na(x = cells.in.ident)]
    if (length(x = cells.in.ident) > max.cells.per.ident) {
      cells.in.ident <-
        sample(x = cells.in.ident, size = max.cells.per.ident)
    }
    cells.to.use <- c(cells.to.use, cells.in.ident)
  }
  cells.use <- cells.to.use
  if (!is.null(x = subset.name)) {
    subset.name <- as.character(subset.name)
    data.use <- FetchData(object = object,
                          vars.all = subset.name,
                          cells.use = cells.use)
    if (length(x = data.use) == 0) {
      stop(paste("Error : ", id, " not found"))
    }
    subset.data <- data.use[, subset.name, drop = F]
    if (!is.null(x = accept.value)) {
      pass.inds <- which(x = subset.data == accept.value)
    } else {
      pass.inds <-
        which(x = (subset.data > accept.low) &
                (subset.data < accept.high))
    }
    cells.use <- rownames(x = data.use)[pass.inds]
  }
  return(cells.use)
}


AverageExpression_my <- function(object,
                                 genes.use = NULL,
                                 return.seurat = FALSE,
                                 add.ident = NULL,
                                 use.scale = FALSE,
                                 use.raw = FALSE,
                                 show.progress = TRUE,
                                 ...) {
  ident.orig <- object@ident
  orig.levels <- levels(x = object@ident)
  ident.new <- c()
  if (!is.null(x = add.ident)) {
    new.data <- FetchData(object = object, vars.all = add.ident)
    new.ident <- paste(object@ident[rownames(x = new.data)],
                       new.data[, 1],
                       sep = '_')
    object <- SetIdent(
      object = object,
      cells.use = rownames(x = new.data),
      ident.use = new.ident
    )
  }
  if (return.seurat) {
    assays.use <- c("RNA", names(x = object@assay))
  } else {
    assays.use <- "RNA"
  }
  slot.use <- "data"
  fxn.average <- function(x)
    mean(expm1(x))
  if (use.scale) {
    slot.use <- "scale.data"
    fxn.average <- mean
  }
  if (use.raw) {
    slot.use <- "raw.data"
    fxn.average <- mean
  }
  data.return <- list()
  for (i in 1:length(x = assays.use)) {
    data.use <- GetAssayData(object = object,
                             assay.type = assays.use[i],
                             slot = slot.use)
    
    
    genes.assay <- genes.use
    if (length(x = intersect(x = genes.use, y = rownames(x = data.use))) <
        1) {
      genes.assay <- rownames(x = data.use)
    }
    data.all <- data.frame(row.names = genes.assay)
    for (j in levels(x = object@ident)) {
      temp.cells <- WhichCells_my(object = object, ident = j)
      genes.assay <-
        unique(x = intersect(x = genes.assay, y = rownames(x = data.use)))
      if (length(x = temp.cells) == 1) {
        data.temp <- as.matrix(data.use[genes.assay, temp.cells])
      }
      if (length(x = temp.cells) > 1) {
        data.temp <- apply(X = as.matrix(data.use[genes.assay, temp.cells]),
                           MARGIN = 1,
                           FUN = fxn.average)
      }
      data.all <- cbind(data.all, data.temp)
      colnames(x = data.all)[ncol(x = data.all)] <- j
      if (show.progress) {
        print(paste0("Finished averaging ", assays.use[i], " for cluster ", j))
      }
      if (i == 1) {
        ident.new <-
          c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
      }
    }
    names(x = ident.new) <- levels(x = object@ident)
    data.return[[i]] <- data.all
    names(x = data.return)[i] <- assays.use[[i]]
  }
  if (return.seurat) {
    toRet <- CreateSeuratObject(
      raw.data = data.return[[1]],
      project = "Average",
      min.cells = 0,
      min.genes = 0,
      is.expr = 0,
      ...
    )
    #for multimodal data
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        toRet <- SetAssayData(
          object = toRet,
          assay.type = names(x = data.return)[i],
          slot = "raw.data",
          new.data = data.return[[i]]
        )
      }
    }
    toRet <- SetIdent(
      object = toRet,
      cells.use = toRet@cell.names,
      ident.use = ident.new[toRet@cell.names]
    )
    toRet@ident <- factor(
      x = toRet@ident,
      levels = as.character(x = orig.levels),
      ordered = TRUE
    )
    
    # finish setting up object if it is to be returned
    
    toRet <- NormalizeData(toRet, display.progress = show.progress)
    toRet <- ScaleData(toRet, display.progress = show.progress)
    
    return(toRet)
  } else {
    return(data.return[[1]])
  }
}



buildExprSimTree <- function(object) {
  require(ape)
  
  data.avg <- AverageExpression_my(object)
  ## remove un-expessed genes
  totalExpr <- apply(data.avg, 1, sum)
  
  data.avg <- data.avg[which(totalExpr > 0),]
  
  data.dist <- dist(t(data.avg))
  data.tree = as.phylo(hclust(data.dist))
  
  object@cluster.tree[[1]] = data.tree
  object
}

calculateInsulationScore <- function(coclust, w = 5, order = T) {
  require(fastcluster)
  
  hc <- fastcluster::hclust(as.dist(1 - coclust))
  
  coclust <- coclust[hc$order, hc$order]
  insulation_score <- c()
  for (i in w:(nrow(coclust) - w)) {
    region <- coclust[(i - w + 1):i, (i + 1):(i + w)]
    score <- mean(region)
    
    insulation_score <- c(insulation_score, score)
  }
  return(insulation_score)
}



localMinima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}




## Find cluster modifier from Seurat

FindClusters.my <-
  function(object,
           genes.use = NULL,
           pc.use = NULL,
           k.param = 30,
           k.scale = 25,
           plot.SNN = FALSE,
           prune.SNN = 1 / 15,
           save.SNN = FALSE,
           reuse.SNN = FALSE,
           use.tsne = NULL,
           do.sparse = FALSE,
           modularity.fxn = 1,
           resolution = 0.8,
           algorithm = 1,
           n.start = 100,
           n.iter = 10,
           random.seed = 0,
           print.output = TRUE,
           temp.file.location = NULL) {
    # for older objects without the snn.k slot
    if (typeof(validObject(object, test = T)) == "character") {
      object@snn.k <- numeric()
    }
    
    snn.built <- FALSE
    if (.hasSlot(object, "snn.dense")) {
      if (length(object@snn.dense) > 1) {
        snn.built <- TRUE
      }
    }
    if (.hasSlot(object, "snn.sparse")) {
      if (length(object@snn.sparse) > 1) {
        snn.built <- TRUE
      }
    }
    
    if ((
      missing(genes.use) &&
      missing(pc.use) && missing(k.param) && missing(k.scale) &&
      missing(prune.SNN) && snn.built
    ) || reuse.SNN) {
      save.SNN <- TRUE
      if (reuse.SNN && !snn.built) {
        stop("No SNN stored to reuse.")
      }
      if (reuse.SNN &&
          (
            !missing(genes.use) || !missing(pc.use) || !missing(k.param) ||
            !missing(k.scale) || !missing(prune.SNN)
          )) {
        warning(
          "SNN was not be rebuilt with new parameters. Continued with stored SNN. To suppress this
          warning, remove all SNN building parameters."
        )
      }
      }
    # if any SNN building parameters are provided or it hasn't been built, build a new SNN
    else{
      object <- BuildSNN.my(
        object,
        genes.use,
        pc.use,
        k.param,
        k.scale,
        plot.SNN,
        prune.SNN,
        use.tsne,
        do.sparse,
        print.output
      )
    }
    
    # deal with sparse SNNs
    if (length(object@snn.sparse) > 1) {
      SNN.use <- object@snn.sparse
    } else {
      SNN.use <- object@snn.dense
    }
    for (r in resolution) {
      object <-
        Seurat:::RunModularityClustering(
          object,
          SNN.use,
          modularity.fxn,
          r,
          algorithm,
          n.start,
          n.iter,
          random.seed,
          print.output
        )
      object <- Seurat:::GroupSingletons(object, SNN.use)
      name <- paste("res.", r, sep = "")
      object <- StashIdent(object, name)
    }
    
    if (!save.SNN) {
      object@snn.sparse <- sparseMatrix(1, 1, x = 1)
      object@snn.dense <- matrix()
      object@snn.k <- integer()
    }
    return(object)
    }






#' Build SNN network using tSNE coordinates. Modified from Seurate package
#'
#' @param object
#' @param genes.use
#' @param pc.use
#' @param k.param
#' @param k.scale
#' @param plot.SNN
#' @param prune.SNN
#' @param use.tsne : True or FALSE
#' @param do.sparse
#' @param print.output
#'
#' @return
#' @export
#'
#' @examples
BuildSNN.my <-
  function(object,
           genes.use = NULL,
           pc.use = NULL,
           k.param = 10,
           k.scale = 10,
           plot.SNN = FALSE,
           prune.SNN = 1 / 15,
           use.tsne = NULL,
           do.sparse = FALSE,
           print.output = TRUE) {
    if (is.null(genes.use) && is.null(pc.use)) {
      genes.use <- object@var.genes
      data.use <- t(as.matrix(object@data[genes.use, ]))
    } else if (!is.null(use.tsne)) {
      data.use <- as.matrix(object@tsne.rot)
    } else if (!is.null(pc.use)) {
      data.use <- as.matrix(object@pca.rot[, pc.use])
    } else if (!is.null(genes.use) && is.null(pc.use)) {
      data.use <- t(as.matrix(object@data[genes.use, ]))
    } else {
      stop("Data error!")
    }
    
    n.cells <- nrow(data.use)
    if (n.cells < k.param) {
      warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.")
      k.param <- n.cells - 1
    }
    
    #find the k-nearest neighbors for each single cell
    my.knn <-
      FNN::get.knn(as.matrix(data.use), k = min(k.scale * k.param, n.cells - 1))
    nn.ranked <-
      cbind(1:n.cells, my.knn$nn.index[, 1:(k.param - 1)])
    nn.large <- my.knn$nn.index
    if (do.sparse) {
      w <-
        Seurat:::CalcSNNSparse(object,
                               n.cells,
                               k.param,
                               nn.large,
                               nn.ranked,
                               prune.SNN,
                               print.output)
    } else {
      w <-
        Seurat:::CalcSNNDense(object,
                              n.cells,
                              nn.large,
                              nn.ranked,
                              prune.SNN,
                              print.output)
    }
    if (plot.SNN) {
      if (length(object@tsne.rot) < 1) {
        warning("Please compute a tSNE for SNN visualization. See RunTSNE().")
      }
      else{
        net <- graph.adjacency(w,
                               mode = "undirected",
                               weighted = TRUE,
                               diag = FALSE)
        plot.igraph(
          net,
          layout = as.matrix(object@tsne.rot),
          edge.width = E(net)$weight,
          vertex.label = NA,
          vertex.size = 0
        )
      }
    }
    
    #only allow one of the snn matrix slots to be filled
    object@snn.k <- k.param
    if (do.sparse == TRUE) {
      object@snn.sparse <- w
      object@snn.dense <- matrix()
    } else {
      object@snn.dense <- w
      object@snn.sparse <- sparseMatrix(1, 1, x = 1)
    }
    return(object)
  }



topMarkerViolionPlot <-
  function(sobj,
           Markers,
           sigPV,
           neuro_clustersColors,
           ttl) {
    Markers.sig <- subset(Markers, pvalue < sigPV)
    ord <- order(Markers.sig$log2FC, decreasing = T)
    Markers.sig <- head(Markers.sig[ord,], 15)
    plotViolinPerGene(Markers.sig$Gene,
                      sobj@ident,
                      sobj@data,
                      neuro_clustersColors,
                      ttl)
  }



getVarGenes <-
  function(sobj,
           cutoff = 2,
           figname = "var_genes.pdf",
           fname = "var_genes.csv") {
    pdf(figname, height = 8, width = 11)
    
    sobj <-
      MeanVarPlot(
        sobj ,
        fxn.x = expMean,
        fxn.y = logVarDivMean,
        x.low.cutoff = 0.5,
        x.high.cutoff = 16,
        y.cutoff = cutoff,
        do.contour = F,
        cex.use = 0.2,
        cex.text.use = 0.5
      )
    dev.off()
    
    write.table(sobj@var.genes, file = fname, row.names = FALSE)
    
    return(sobj)
  }



addConvexHull <- function(coord) {
  if (ncol(coord) < 2) {
  }
  
  
  ch <- chull(cls1_pos$tSNE_1, cls1_pos$tSNE_2)
}

#' Plot the different projections of tSNE-plot aloong the the 3D space axis
#'
#' @param sobj Seurat object
#'
#' @return
#' @export
#'
#' @examples
plottSNEClustering <-
  function(sobj,
           clustersColors = NULL,
           clusLevels = NULL,
           contours = F, 
           label.size=4) {
    if (is.null(clustersColors)) {
      clustersColors <- distinctColorPalette(length(levels(sobj@ident)))
      names(clustersColors) <- levels(sobj@ident)
    }
    
    if (!is.null(clusLevels)) {
      sobj@ident <-
        factor(as.character(sobj@ident),
               levels = clusLevels,
               ordered = T)
      #sobj <- SetIdent(sobj,ident.use = ident)
    }
    
    
    boldtxt <- theme(text = element_text(face = "bold", size = 12))
    plt_12_all <-
      TSNEPlot(
        sobj,
        colors.use = clustersColors,
        do.return = T,
        dim.1 = 1,
        dim.2 = 2,
        do.label = T, 
        label.size = label.size
      ) + boldtxt
    
    if (ncol(sobj@tsne.rot) >= 3) {
      plt_13_all <-
        TSNEPlot(
          sobj,
          colors.use = clustersColors,
          do.return = T,
          dim.1 = 1,
          dim.2 = 3,
          do.label = T,
          label.size = label.size
        ) + boldtxt
      plt_23_all <-
        TSNEPlot(
          sobj,
          colors.use = clustersColors,
          do.return = T,
          dim.1 = 2,
          dim.2 = 3,
          do.label = T,
          label.size = label.size
        ) + boldtxt
      plot_grid(
        plt_12_all + theme(legend.position = "none") ,
        plt_13_all + theme(legend.position = "none"),
        plt_23_all,
        ncol = 2
      )
    } else{
      plt_12_all
    }
  }



plotCoClustheatmap <-
  function(sobj,
           coclust,
           clustersColors = NULL,
           cells_order = NULL) {
    ord <- order(sobj@ident)
    
    cells_cluster_anno <- data.frame(Cluster = sobj@ident[ord])
    
    if (is.null(clustersColors)) {
      clustersColors <- distinctColorPalette(length(levels(sobj@ident)))
      names(clustersColors) <- levels(sobj@ident)
    }
    
    ann_colors <- list(Cluster = clustersColors)
    
    if (is.null(cells_order)) {
      ph <- pheatmap(
        coclust,
        cluster_cols = T,
        cluster_rows = T,
        annotation_col = cells_cluster_anno,
        show_colnames = F,
        show_rownames = F,
        annotation_colors = ann_colors
      )
    } else{
      ph <-
        pheatmap(
          coclust[cells_order, cells_order],
          cluster_cols = F,
          cluster_rows = F,
          annotation_col = cells_cluster_anno,
          show_colnames = F,
          show_rownames = F,
          annotation_colors = ann_colors
        )
    }
    
    
    #tmp <- coclust[ord,ord]
    #colnames(tmp) <- NULL
    #rownames(tmp) <- NULL
    #levelplot(t(tmp[rev(ph$tree_row$order) ,rev(ph$tree_col$order)]), par.settings = BuRdTheme)
  }


## CEF format used by backspin

ConvertToCef <- function(sobj,  genes=NULL, fout) {
  
  if(is.null(genes)){
    expr <- round(expm1(sobj@data))  
  }else{
    expr <- round(expm1(sobj@data[genes,]))
  }
  
  
  colnames(expr) <- gsub("\\.+$", "", colnames(expr))
  colnames(expr) <- gsub("\\.", "_", colnames(expr))
  
  write(file = fout,
        x = paste("CEF\t0\t1\t1", nrow(expr), ncol(expr), "0", sep = "\t"))
  write(file = fout,
        paste('\t', c("cell"), paste(colnames(expr), collapse = "\t"), sep = '\t'),
        append = T)
  x <- cbind(rep("", nrow(expr)), as.matrix(expr))
  write(file = fout, paste(c("gene"), sep = "\t"), append = T)
  write.table(
    file = fout,
    x,
    append = T,
    col.names = F,
    row.names = T,
    quote = F,
    sep = "\t"
  )
}


## MNF


myNMF <-
  function(data,
           prefix = "NMF",
           cluster = 3,
           top = 1500,
           nrun = 100,
           norm = F,
           ncores = 8,
           algorithm = "brunet",
           mode = "real",
           seed = 123211) {
    if (mode == "estim") {
      r <- 2:cluster
      estim.r <-
        nmf(
          data,
          r,
          algorithm,
          .opt = paste0("vp", ncores),
          nrun = nrun,
          seed = seed,
          maxIter = 5000
        )
      res <- estim.r
    } else if (mode == "compare") {
      res.multi.method <-
        nmf(
          data,
          cluster,
          list("brunet", "lee", "ns"),
          nrun = nrun,
          .opt = paste0("vtp", ncores),
          seed = seed
        )
      # compare(res.multi.method)
      # print(compare(res.multi.method))
      res <- res.multi.method
    } else if (mode == "real") {
      # option 't' will toggle error track function
      res <-
        nmf(
          data,
          cluster,
          algorithm,
          .opt = paste0("vtp", ncores),
          nrun = nrun,
          seed = seed,
          maxIter = 5000
        )
    }
    
    return(res)
  }


getClustersOrignalGroupComp <- function(sobj) {
  info <-
    data.frame(cluster = sobj@ident,
               treatment = sobj@data.info$orig.ident)
  
  p <-
    ggplot(info, aes(x = cluster, fill = treatment)) + geom_bar(position = position_dodge()) + facet_grid(treatment ~
                                                                                                            .)
  p
}



CompareClusters <- function(sobj, clus1, clus2) {
  comp_res <- list()
  if (clus1 == clus2) {
    stop("Clus1 should be different than clus2")
  }
  
  if (!clus1 %in% levels(sobj@ident)) {
    stop(paste0("Couldn't find cells in", clus1))
  }
  
  if (!clus2 %in% levels(sobj@ident)) {
    stop(paste0("Couldn't find cells in", clus2))
  }
  
  inClus1 <- names(sobj@ident[sobj@ident == clus1])
  #inClus1 <- grep("Saline",inClus1,value = T)
  inClus2 <- names(sobj@ident[sobj@ident == clus2])
  #inClus2 <- grep("Saline",inClus2,value = T)
  
  if (length(inClus1) > 1 & length(inClus2) > 1) {
    n <- paste0(clus1, "_vs_", clus2)
    comp_res[[n]] <-
      getMASTmarkers.group(sobj = sobj,
                           pos1 = inClus2,
                           pos2 = inClus1)
    D1_clusters_Markers.sig <- subset(comp_res[[n]], pvalue < 1e-3)
    
    ttl <- paste("DE genes between cluster", clus1, "and", clus2)
    plotViolinPerGene(
      head(D1_clusters_Markers.sig$Gene, 15),
      sobj@ident,
      sobj@data,
      NAc_D1_clustersColors_reclust,
      ttl,
      "violon"
    )
  }
  
  return(comp_res)
}

### GO analysis

plotGoBars <-
  function(go_res,
           top = 20,
           col = "red",
           title = "top 10 enriched GO") {
    tmp <- as.data.frame(go_res)
    ord <- order(tmp$qvalue)
    tmp <- head(tmp[ord,], top)
    
    ord <- order(tmp$qvalue, decreasing = T)
    tmp <- tmp[ord,]
    tmp$Description <-
      factor(tmp$Description, levels = tmp$Description)
    
    boldText = theme(
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 8)
    )
    
    g <-
      ggplot(tmp, aes(x = Description, y = -log10(qvalue))) + geom_bar(stat = "identity", fill = col) + coord_flip() +
      theme_bw() + ggtitle(title) + boldText
    g
  }



DOgsea <- function(deResults){
  
  
  require(annotables)
  require(clusterProfiler)
  
  data("grcm38")
  gogmt <- read.gmt("D:/genomes/mm9/GO_gmt/Mus_musculus_GSEA_GO_sets_bp_ids_highquality_April_2015.gmt")
  
  pos <- match(toupper(deResults$GeneName), toupper(grcm38$symbol) )
  
  entrezIDs <- grcm38$entrez[pos]
  
  fc <- deResults$logFC[!is.na(entrezIDs)]
  
  #fc[is.nan(fc)] <- deResults$empirical_logFC[is.nan(fc)]
  
  names(fc ) <- entrezIDs[!is.na(entrezIDs)]
  
  fc <- fc[!is.nan(fc)]
  
  fc <- sort(fc,decreasing = T)
  
  gsea <- GSEA(fc,TERM2GENE = gogmt,pvalueCutoff = 1)
  
  return(gsea)
}


doGOAnalysis <- function(deResult) {
  DEgenes <- deResult$geneID
  treatment_DE <- unique(unlist(DEgenes))
  ## Convert gene names to EntrezID
  treatment_DE.entrez <- EnsemblToEntrez(treatment_DE)
  treatment_DE.entrez <-
    subset(treatment_DE.entrez, entrezID != "NA")
  
  
  go_res <- list()
  for (trm in c("BP", "CC", "MF")) {
    go_res[[trm]] <-
      enrichGO(
        treatment_DE.entrez$entrezID,
        'org.Mm.eg.db',
        ont = trm,
        pvalueCutoff = 0.05
      )
  }
  return(go_res)
}



GeneNameToEnsemble <- function(genesSymbols, species = "mouse") {
  if (species  == "human") {
    require(org.Hs.eg.db)
    xx <- as.list(org.Hs.egSYMBOL)
    xx <- xx[!is.na(xx)]
  }
  else{
    if (species  == "mouse") {
      require(org.Mm.eg.db)
      xx <- as.list(org.Mm.egSYMBOL2EG)
      xx <- xx[!is.na(xx)]
    }
  }
  
  pos <- match(toupper(genesSymbols), toupper(names(xx)))
  genesSymbols <- genesSymbols[!is.na(pos)]
  df <-
    data.frame(
      entrezID = genesSymbols,
      geneSymbol = "",
      stringsAsFactors = FALSE
    )
  
  for (i in 1:nrow(df)) {
    df$geneSymbol[i] <- xx[[genesSymbols[i]]][1]
  }
  return(df)
}


scale_rows = function(x) {
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale) {
  if (!(scale %in% c("none", "row", "column"))) {
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(
    scale,
    none = mat,
    row = scale_rows(mat),
    column = t(scale_rows(t(mat)))
  )
  return(mat)
}




FeaturePlot3D <- function(sobj, marker, dolog = FALSE) {
  sobj_3D <- as.data.frame(sobj@tsne.rot)
  sobj_3D$cluster <- factor(sobj@ident)
  
  if (!marker %in% rownames(sobj@data)) {
    stop(paste(marker, "not expressed"))
  }
  
  if (dolog) {
    sobj_3D$expr <- log(sobj@data[marker,] + 1)
  } else{
    sobj_3D$expr <- sobj@data[marker,]
  }
  
  
  plot_ly(
    sobj_3D,
    x = ~ tSNE_1,
    y =  ~ tSNE_2,
    z = ~ tSNE_3,
    marker = list(
      size = 2,
      color = ~ expr,
      colorscale = c('#FFFFCC', '#800026'),
      showscale = TRUE
    )
  ) %>%
    add_markers() %>%
    layout(title = paste(marker, "expression"))
}

plottSNEClustering3D <- function(sobj, cluster_colors) {
  sobj_3D <- as.data.frame(sobj@tsne.rot)
  sobj_3D$cluster <- factor(sobj@ident)
  
  plot_ly(
    sobj_3D,
    x = ~ tSNE_1,
    y =  ~ tSNE_2,
    z = ~ tSNE_3,
    color = ~ cluster,
    colors = cluster_colors,
    marker = list(size = 2)
  )
}

doConcensusClustering <- function(sobj, pcs = 1:15) {
  for (algo in c(1:3)) {
    for (res in c(0.6, 1, 1.5, 2)) {
      message(paste("Trying algorithm", algo, "Resoltuion:", res))
      sobj <-
        doSeuratClustering(sobj,
                           algo = algo,
                           res = res,
                           pcs = pcs)
    }
  }
  
  for (space in c("pca")) {
    for (criteria in c("BIC", "ICL")) {
      message(paste("Trying mclust in ", space, "Using:", criteria))
      sobj <-
        doMclust(
          sobj,
          space =  space,
          pcs.touse = pcs,
          criteria =  criteria,
          nClust = 2:20
        )
    }
  }
  
  #message("Kmean clustering")
  #sobj <- doKmeansClustering(sobj,pcs = pcs, max.nc = 20)
  
  coclust <-  calcCSPA(sobj)
  return(coclust)
  
}



selectBestClustering <- function(sobj, pcs = 1:15) {
  for (algo in c(1:3)) {
    for (res in c(0.6, 1, 1.5, 2, 3)) {
      message(paste("Trying algorithm", algo, "Resoltuion:", res))
      sobj <-
        doSeuratClustering(sobj,
                           algo = algo,
                           res = res,
                           pcs = pcs)
    }
  }
  
  for (space in c("pca")) {
    for (criteria in c("BIC", "ICL")) {
      message(paste("Trying mclust in ", space, "Using:", criteria))
      sobj <-
        doMclust(
          sobj,
          space =  space,
          pcs.touse = pcs,
          criteria =  criteria,
          nClust = 2:20
        )
    }
  }
  
  #message("Kmean clustering")
  #sobj <- doKmeansClustering(sobj,pcs = pcs, max.nc = 20)
  
  algos <- grep("algo", colnames(sobj@data.info))
  
  clustRes <- sobj@data.info[, algos]
  
  dunnIndex <- calcDunnIdent(sobj, pcs)
  
  res <- list(clustRes = clustRes, separation = dunnIndex)
  return(res)
}



calcDunnIdent <- function(sobj, pcs = 1:15, pattern="algo") {
  require(distances)
  require(clValid)
  
  D <- distances(sobj@pca.rot[, pcs])
  D <- as.matrix(D)
  
  
  algos <- grep(pattern, colnames(sobj@data.info))
  if (length(algos) == 0) {
    stop("Please run some clustering algorithms than call this method")
  }
  
  
  separation <- data.frame()
  
  
  for (a  in algos) {
    dn <- dunn(D, clusters = as.numeric(sobj@data.info[, a]))
    df <- data.frame(method = a, dunn = dn)
    
    separation <- rbind(separation, df)
  }
  
  return(separation)
}


calcDunnIdent_varGene <- function(sobj) {
  require(distances)
  require(clValid)
  
  D <- distances(t(sobj@scale.data[sobj@var.genes,]))
  D <- as.matrix(D)
  
  
  algos <- grep("algo", colnames(sobj@data.info))
  if (length(algos) == 0) {
    stop("Please run some clustering algorithms than call this method")
  }
  
  
  separation <- data.frame()
  
  
  for (a  in algos) {
    dn <- dunn(D, clusters = sobj@data.info[, a])
    df <- data.frame(method = a, dunn = dn)
    
    separation <- rbind(separation, df)
  }
  
  return(separation)
}


plotMarkers.pheatmap <-
  function(sobj,
           markers,
           clustersColors = NULL,
           markerClusters = NULL,
           clusterCols = F) {
    require(ComplexHeatmap)
    require(circlize)
    ord <- order(as.numeric(as.character(sobj@ident)))
    
    cells_cluster_anno <- data.frame(Cluster = sobj@ident[ord])
    rownames(cells_cluster_anno) <- colnames(sobj@data)
    s
    if (is.null(clustersColors)) {
      clustersColors <- distinctColorPalette(length(levels(sobj@ident)))
      names(clustersColors) <- levels(sobj@ident)
    }
    
    ann_colors <- list(Cluster = clustersColors)
    
    
    if (!all(markers %in%  rownames(sobj@data))) {
      stop("Check that all of the provided genes are in sobj")
    }
    
    
    
    
    tmp <- pheatmap:::scale_rows(sobj@data[markers, ord])
    
    pheatmap(
      log2(sobj@data[markers, ord] + 1),
      cluster_cols = F,
      cluster_rows = T,
      annotation_col = cells_cluster_anno,
      show_colnames = F,
      show_rownames = T,
      annotation_colors = ann_colors
    )
  }





plotMarkersHeatmap <-
  function(sobj,
           markers,
           scale_rows=TRUE,
           scale_min=-2, scale_max=2,
           clustersColors = NULL,
           markerClusters = NULL,
           show_gene_names = T,
           clusterCols = F, cluster_rows=F,
           identlevels = NULL, title=NULL) {
    require(ComplexHeatmap)
    require(circlize)
    #ord <- order(as.numeric(as.character(sobj@ident)))
    
    
    if (!is.null(identlevels)) {
      id <- as.character(sobj@ident)
      id <- factor(id, levels = identlevels, ordered = T)
      ord <- order(id)
    } else{
      ord <- order(sobj@ident)
      identlevels <- levels(sobj@ident)
    }
    
    
    cells_cluster_anno <- data.frame(Cluster = factor(
      as.character(sobj@ident[ord]),
      levels = identlevels,
      ordered = T
    ))
    rownames(cells_cluster_anno) <- colnames(sobj@data)
    
    if (is.null(clustersColors)) {
      clustersColors <- distinctColorPalette(length(levels(sobj@ident)))
      names(clustersColors) <- levels(sobj@ident)
    }
    
    ann_colors <- list(Cluster = clustersColors[identlevels])
    
    if (!is.null(markerClusters)) {
      markerClusters <-
        factor(as.character(markerClusters), levels = unique(c(identlevels,markerClusters)))
    }
    
    
    # cellOrder <- ord
    #
    # minOrd <- 0
    # for(cls in unique(sobj@ident[ord])){
    #   pos <- which(sobj@ident == cls)
    #
    #   mat <- sobj@data[markers,pos]
    #
    #   d <- dist(t(mat))
    #   sub_ord <- hclust(d)$order
    #
    #   sub_ord <- minOrd + sub_ord
    #   minOrd <- sub_ord
    #
    #   cellOrder[pos] <- sub_ord
    # }
    
    
    
    if (!all(markers %in%  rownames(sobj@data))) {
      stop("Check that all of the provided genes are in sobj")
    }
    
    
    if(scale_rows){
      tmp <- pheatmap:::scale_rows(sobj@data[markers, ord]) 
      
      ha_column = HeatmapAnnotation(df = cells_cluster_anno, col = ann_colors)
      hcol = colorRamp2(seq(scale_min, scale_max, length.out = 100), colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100))
      
    }else{
      tmp <- sobj@data[markers, ord]  
      
      ha_column = HeatmapAnnotation(df = cells_cluster_anno, col = ann_colors)
      q = quantile(tmp,0.99)
      hcol = colorRamp2(seq(0, q, length.out = 100), colorRampPalette(rev(brewer.pal(
        n = 11, name = "RdYlBu")))(100))
    }
    
    tmp <- as.matrix(tmp)
    #tmp[tmp < -2] <- -2
    #tmp[tmp > 2] <- 2
    #tmp <- sobj@data[markers,ord]
    #pheatmap(log2(sobj@data[markers,ord]+1),cluster_cols = F,cluster_rows = T,
    #               annotation_col = cells_cluster_anno,show_colnames = F,show_rownames = T,
    #               annotation_colors = ann_colors)
    
    
    
    if (is.null(markerClusters)) {
      H <- Heatmap(
        tmp,
        col = hcol,
        top_annotation = ha_column,
        cluster_columns = clusterCols,
        show_column_names = F,
        show_row_names = show_gene_names,
        show_row_dend = F,
        cluster_rows = cluster_rows,
        column_title = title,
        name = "log(Expr)"
      )
    } else{
      H <- Heatmap(
        tmp,
        col = hcol,
        split = markerClusters,
        top_annotation = ha_column,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = show_gene_names,
        show_row_dend = F,
        cluster_rows = cluster_rows,
        column_title = title,
        name = "log(Expr)"
      )
    }
    
    draw(H)
    ## Add the clusters line separators
    nbCells <- table(sobj@ident)[unique(sobj@ident[ord])]
    
    nbSplits <-
      ifelse(is.null(markerClusters), 1, length(unique(markerClusters)))
    
    for (i in 1:length(nbCells)) {
      sm <- sum(nbCells[1:i])
      
      for (slice in 1:nbSplits) {
        decorate_heatmap_body("log(Expr)", {
          x = sm / ncol(sobj@data)
          grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 2))
          grid.rect(gp = gpar(col = "black", fill = NA))
        }, slice = slice)
      }
    }
}


plotMarkersHeatmap_v3 <-
  function(sobj,
           markers,
           scale_rows=TRUE,
           scale_min=-2, scale_max=2,
           clustersColors = NULL,
           markerClusters = NULL,
           show_gene_names = T,
           clusterCols = F, 
           cluster_rows=F,
           slotToUse="data",
           identlevels = NULL, title=NULL) {
    
    require(ComplexHeatmap)
    require(circlize)
    require(RColorBrewer)
    
    
    if (!is.null(identlevels)) {
      id <- as.character(Idents(sobj))
      id <- factor(id, levels = identlevels, ordered = T)
      ord <- order(id)
    } else{
      ord <- order(Idents(sobj))
      identlevels <- levels(Idents(sobj))
    }
    
    
    cat = gsub("_.*$","",as.character(Idents(sobj)[ord]))
    cells_cluster_anno <- data.frame(Cluster = factor(as.character(Idents(sobj)[ord]),levels = identlevels,ordered = T),
                                     Category = factor(cat, levels = unique(gsub("_.*$","",identlevels)),ordered = T))
    rownames(cells_cluster_anno) <- colnames(sobj)
    
    if (is.null(clustersColors)) {
      clustersColors <- distinctColorPalette(length(levels(Idents(sobj))))
      names(clustersColors) <- levels(Idents(sobj))
    }
    
    caterogiesColors = ggsci::pal_d3("category20c")(length(unique(cat)))
    names(caterogiesColors) = unique(gsub("_.*$","",identlevels))
    
    ann_colors <- list(Cluster = clustersColors[identlevels], Category =caterogiesColors)
    
    if (!is.null(markerClusters)) {
      markerClusters <-
        factor(as.character(markerClusters), levels = unique(c(identlevels,markerClusters)))
    }
    
    
    # cellOrder <- ord
    #
    # minOrd <- 0
    # for(cls in unique(sobj@ident[ord])){
    #   pos <- which(sobj@ident == cls)
    #
    #   mat <- sobj@data[markers,pos]
    #
    #   d <- dist(t(mat))
    #   sub_ord <- hclust(d)$order
    #
    #   sub_ord <- minOrd + sub_ord
    #   minOrd <- sub_ord
    #
    #   cellOrder[pos] <- sub_ord
    # }
    
    
    
    if (!all(markers %in%  rownames(sobj))) {
      stop("Check that all of the provided genes are in sobj")
    }
    
    
    if(scale_rows){
      tmp <- pheatmap:::scale_rows(sobj[['RNA']]@data[markers, ord]) 
      
      ha_column = HeatmapAnnotation(df = cells_cluster_anno, col = ann_colors,simple_anno_size_adjust = TRUE)
      hcol = colorRamp2(seq(scale_min, scale_max, length.out = 100), colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100))
      
    }else{
      tmp <- sobj[['RNA']]@data[markers, ord]
      
      ha_column = HeatmapAnnotation(df = cells_cluster_anno, col = ann_colors)
      q = quantile(tmp,0.99)
      hcol = colorRamp2(seq(0, q, length.out = 100), colorRampPalette(rev(brewer.pal(
        n = 11, name = "RdYlBu")))(100))
    }
    
    tmp <- as.matrix(tmp)
    #tmp[tmp < -2] <- -2
    #tmp[tmp > 2] <- 2
    #tmp <- sobj@data[markers,ord]
    #pheatmap(log2(sobj@data[markers,ord]+1),cluster_cols = F,cluster_rows = T,
    #               annotation_col = cells_cluster_anno,show_colnames = F,show_rownames = T,
    #               annotation_colors = ann_colors)
    
    
    
    if (is.null(markerClusters)) {
      H <- Heatmap(
        tmp,
        col = hcol,
        top_annotation = ha_column,
        cluster_columns = clusterCols,
        show_column_names = F,
        show_row_names = show_gene_names,
        show_row_dend = F,
        cluster_rows = cluster_rows,
        column_title = title,
        name = "log(Expr)"
      )
    } else{
      H <- Heatmap(
        tmp,
        col = hcol,
        split = markerClusters,
        top_annotation = ha_column,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = show_gene_names,
        show_row_dend = F,
        cluster_rows = cluster_rows,
        column_title = title,
        name = "log(Expr)"
      )
    }
    
    draw(H) #heatmap_legend_side = "right", annotation_legend_side = "bottom")
    ## Add the clusters line separators
    nbCells <- table(Idents(sobj))[unique(Idents(sobj)[ord])]
    
    nbSplits <-
      ifelse(is.null(markerClusters), 1, length(unique(markerClusters)))
    
    for (i in 1:length(nbCells)) {
      sm <- sum(nbCells[1:i])
      
      for (slice in 1:nbSplits) {
        decorate_heatmap_body("log(Expr)", {
          x = sm / ncol(sobj)
          grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 2))
          grid.rect(gp = gpar(col = "black", fill = NA))
        }, slice = slice)
      }
    }
  }



plotCombinedMarkers <- function(sobj_list,
                                sobj_mrks_list,
                                power.thr = 0.7,
                                pct.1.thr = 0.7,
                                pct.2.thr = 0.2,
                                is.pval = FALSE,
                                pval.thr = 1e-3,
                                nbGenes=5,
                                scale_min=-2,
                                scale_max=2,
                                logFC.thr=log(1.5),
                                clustersColors = NULL,
                                cols=NULL,scale_rows=TRUE,
                                additional_markers= NULL){
  
  require(dplyr)
  require(ComplexHeatmap)
  require(RColorBrewer)
  require(circlize)
  
  
  ## Get the list of genes
  markers <- c()
  
  for(i in 1:length(sobj_mrks_list)){
    
    if (is.pval) {
      mrks.sig <- sobj_mrks_list[[i]] %>% 
        group_by(cluster) %>%
        filter(pct.1 > pct.1.thr & pct.2 < pct.2.thr & p_val < pval.thr & avg_logFC > logFC.thr)
        
    } else{
      mrks.sig <- sobj_mrks_list[[i]] %>% 
        group_by(cluster) %>%
        filter(pct.1 > pct.1.thr & pct.2 < pct.2.thr & power > power.thr & avg_logFC > logFC.thr)
    }
    
    markers <- c(markers, mrks.sig$gene)
  }
  
  markers <- unique(markers)
  
  ## Plot the heatmaps
  
  if(is.null(names(sobj_list))){
    names(sobj_list) <- 1:length(sobj_list)
  }
  
  Heatmaps <- NULL

  merged_expr <- matrix(0,nrow = length(markers))
  rownames(merged_expr) <- markers
  
  
  legends <- list()
  
  for(i in 1:length(sobj_list)){
    sobj <- sobj_list[[i]]
    
    
    ord <- order(sobj@ident)
    identlevels <- levels(sobj@ident)
    
    mrks_in <- intersect(markers, rownames(sobj@data))
    mrks_notIn <- setdiff(markers, rownames(sobj@data))
    
    expr <- as.matrix(sobj@data[mrks_in,ord])
    
    # if there are some genes thare are not expressed in the Seurat object, set them to 0
    if(length(mrks_notIn) >0){
      tmp.expr <- matrix(0,nrow=length(mrks_notIn), ncol=length(ord))
      rownames(tmp.expr) <- mrks_notIn
      colnames(tmp.expr) <- colnames(sobj@data)[ord]
      
      expr <- rbind(expr, tmp.expr)
      rm(tmp.expr)
    }
    rm(mrks_in, mrks_notIn)
    
    ## Prepare annotation
    cells_cluster_anno <- data.frame(Cluster = factor(
      as.character(sobj@ident[ord]),
      levels = identlevels,
      ordered = T
    ))
    rownames(cells_cluster_anno) <- colnames(sobj@data)
    
    if (is.null(clustersColors) | class(clustersColors)!="list") {
      clus_cols <- generateColors(length(levels(sobj@ident)))
      names(clus_cols) <- levels(sobj@ident)
    }else{
      clus_cols <- clustersColors[[i]] 
      names(clus_cols) <- levels(sobj@ident)
    }
    
    ann_colors <- list(Cluster = clus_cols[identlevels])
    
    ## Scale the matrix if needed
    if(scale_rows){
      
      tmp <- pheatmap:::scale_rows(expr[markers,]) 
      
      ha_column = HeatmapAnnotation(df = cells_cluster_anno, col = ann_colors)
      hcol = colorRamp2(seq(scale_min, scale_max, length.out = 100), colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100))
      
    }else{
      tmp <- sobj@data[markers, ord]  
      
      ha_column = HeatmapAnnotation(df = cells_cluster_anno, col = ann_colors)
      q = quantile(tmp,0.99)
      hcol = colorRamp2(seq(0, q, length.out = 100), colorRampPalette(rev(brewer.pal(
        n = 11, name = "RdYlBu")))(100))
    }
    
    
    legends[[i]]  <- Legend(at = names(ann_colors$Cluster), 
                            title = paste("Clusters of",names(sobj_list)[i]),
                            background = ann_colors$Cluster,
                            type = "grid",
                            legend_gp = gpar(col=ann_colors$Cluster)
                            )
    
    # Create Heatmap
    H <- Heatmap(
      tmp,
      col = hcol,
      top_annotation = ha_column,
      cluster_columns = F,
      show_column_names = F,
      row_names_side="left",
      show_heatmap_legend  = FALSE,
      show_row_names = ifelse(i==1,TRUE,FALSE),
      show_row_dend = F,
      cluster_rows = FALSE,
      name = names(sobj_list)[i]
    )
    
    
    merged_expr <- cbind(merged_expr, tmp[markers,])
    
    if(is.null(Heatmaps)){
      Heatmaps <- H
    }else{
      Heatmaps <- Heatmaps + H
    }
  }
  
  tmp[is.nan(tmp)] <- scale_min
  tmp[is.na(tmp)] <- scale_min
  hc = hclust(dist(tmp))
  
  #draw(Heatmaps,row_order = hc$order)
  draw(Heatmaps, annotation_legend_list = legends, show_heatmap_legend =FALSE)
  
  ## Add Clusters borders
  for(n in 1:length(sobj_list)){
    sobj <- sobj_list[[n]]
    ord <- order(sobj@ident)
    
    nbCells <- table(sobj@ident)[unique(sobj@ident[ord])]
    
    nbSplits <- 1 #length(unique(sobj@ident[ord]))
      
    
    for (i in 1:length(nbCells)) {
      sm <- sum(nbCells[1:i])
      
      for (slice in 1:nbSplits) {
        decorate_heatmap_body(names(sobj_list)[n], {
          x = sm / ncol(sobj@data)
          grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 2))
          grid.rect(gp = gpar(col = "black", fill = NA))
        }, slice = slice)
      }
    }
  }
}




WhiteNoise_ll <- function(x, y) {
  options = gpOptions()
  #options$learnScales = FALSE
  options$kern$comp = list("white")
  kerntype = list(type = "cmpnd", comp = list("white"))
  trueKern = kernCreate(x, kerntype)
  K = kernCompute(trueKern, x)
  model = gpCreate(dim(x)[2], dim(y)[2], x, y, options)
  
  inithypers = log(model$kern$comp[[1]]$variance)
  model = gpExpandParam(model, inithypers) ## This forces kernel computation.
  ll_init = gpLogLikelihood(model) ## GP log-marginal likelihood for this model.
  
  model = gpOptimise(model, display = TRUE, iters = 500)
  opthypers = gpExtractParam(model, only.values = FALSE)
  
  opthypers = exp(opthypers)
  
  ll_opt = gpLogLikelihood(model)
  
  return(ll_opt)
}

rbfKernel_ll <- function(x, y) {
  options = gpOptions()
  #options$learnScales = T
  options$kern$comp = list("rbf", "white")
  
  model = gpCreate(dim(x)[2], dim(y)[2], x, y, options)
  
  inithypers = log(c(1 / 0.05, 1, model$kern$comp[[2]]$variance))
  model = gpExpandParam(model, inithypers) ## This forces kernel computation.
  ll_init = gpLogLikelihood(model) ## GP log-marginal likelihood for this model.
  
  model = gpOptimise(model, display = TRUE, iters = 500)
  #opthypers = gpExtractParam(model, only.values = FALSE)
  #opthypers = exp(opthypers);
  ll_opt = gpLogLikelihood(model)
  
  # Return the fitted function to use as the expression of the pathay
  xtest = matrix(1:nrow(x))
  meanVar = gpPosteriorMeanVar(model, xtest, varsigma.return = TRUE)
  res <- list(likelihood = ll_opt,
              posterioMean =  t(meanVar$mu))
  return(res)
}





Plot_PseudoTime_TSNE <- function(monoObj, sobj) {
  pd <- pData(monoObj)
  
  tSNE_pos <- sobj@tsne.rot
  tSNE_pos <-
    cbind(tSNE_pos, pd[colnames(sobj@scale.data), c("Pseudotime", "State")])
  
  
  ggplot(tSNE_pos, aes(x = tSNE_1, y = tSNE_2, color = Pseudotime)) + geom_point() + theme_bw() +
    scale_color_gradient(high = "red", low = "yellow")
}

plot_PseudoTime_heatmap <-
  function(monoObj,
           sobj,
           markers,
           show_rownames = FALSE,
           cluster_row = FALSE,
           max_scale = 3,
           min_scale = -3) {
    pd <- pData(monoObj)
    ord <- order(pd$Pseudotime)
    annotation_col <- data.frame(CellType = factor(sobj@ident))
    tmp <- scale(t(sobj@scale.data[markers, rownames(pd)[ord]]))
    tmp <- apply(tmp, 2, EnrichedHeatmap::default_smooth_fun)
    tmp <- t(tmp)
    colnames(tmp) <- rownames(pd)[ord]
    tmp[tmp > max_scale] <- max_scale
    tmp[tmp < min_scale] <- min_scale
    ph <- pheatmap(
      tmp,
      show_rownames = FALSE,
      scale = "none",
      show_colnames = FALSE,
      annotation_col = annotation_col,
      cluster_cols = F,
      cluster_rows = F
    )
    ph
  }


smooth_fun = function(x) {
  l = !is.na(x)
  if (sum(l) >= 2) {
    oe1 = try(x <-
                suppressWarnings(predict(locfit(x[l] ~ lp(
                  seq_along(x)[l], nn = 0.1, h = 0.8
                )), seq_along(x))), silent = TRUE)
    if (inherits(oe1, "try-error")) {
      oe2 = try(x <-
                  suppressWarnings(predict(
                    loess(x[l] ~ seq_along(x)[l], control = loess.control(surface = "direct")),
                    seq_along(x)
                  )))
      
      if (inherits(oe2, "try-error")) {
        stop("error when doing locfit or loess smoothing")
      } else {
        return(x)
      }
    } else {
      return(x)
    }
  } else {
    stop("Too few data points.")
  }
  return(x)
}




SNN <- function(data, outfile, k, distance) {
  if (missing(data)) {
    stop(paste("Input data missing.", help, sep = "\n"))
  }
  if (missing(outfile)) {
    stop(paste("Output file name missing.", help, sep = "\n"))
  }
  if (missing(k)) {
    k = 3
  }
  
  if (missing(distance)) {
    distance <-
      "euclidean"  # other distance options refer to dist() in R
  }
  #m<-as.data.frame(data)
  
  numSpl <- dim(data)[1]
  m <- dist(data, distance, diag = TRUE, upper = TRUE)
  x <- as.matrix(m)
  IDX <- t(apply(x, 1, order)[1:k,]) # knn list
  
  edges <- list()              # SNN graph
  for (i in 1:numSpl) {
    j <- i
    while (j < numSpl) {
      j <- j + 1
      shared <- intersect(IDX[i,], IDX[j,])
      if (length(shared) > 0) {
        s <- k - 0.5 * (match(shared, IDX[i,]) + match(shared, IDX[j,]))
        strength <- max(s)
        if (strength > 0)
          edges <- rbind(edges, c(i, j, strength))
      }
    }
  }
  return(edges)
}


getVariableGenesM3Drop <-
  function(sobj,
           min.genes = 1500,
           mt_method = "BH",
           threshold = 1e-30) {
    require(M3Drop)
    
    M3DNorm <- M3DropCleanData(
      as.matrix(sobj@raw.data),
      labels = colnames(sobj@raw.data),
      is.counts = TRUE,
      min_detected_genes = min.genes
    )
    
    M3Dnorma_fits <- M3DropDropoutModels(M3DNorm$data)
    
    var_genes <-
      M3DropDifferentialExpression(M3DNorm$data,
                                   mt_method = "fdr",
                                   mt_threshold = threshold)
    
    return(var_genes)
  }



plotBetweenTreatmentHeatmap <- function(sobj = NULL,
                                        cells1 = NULL,
                                        cells2 = NULL,
                                        types  = NULL,
                                        colors = NULL,
                                        genes.use = NULL,
                                        ttl = "",
                                        scale = "none", 
                                        useCheatmap=FALSE) {
  if (is.null(genes.use)) {
    genes.use <- rownames(sobj@data)
  }
  
  ## check if the cell names exist
  if (!all(cells1 %in% colnames(sobj@raw.data)) |
      !all(cells2 %in% colnames(sobj@data))) {
    stop("Some of the provided cell names do not exist in sobj@data")
  }
  
  
  if (length(types) != 2) {
    stop("Two types only can be specified in the 'types' variable")
  }
  
  
  
  expr <-
    matrix(sobj@data[genes.use, c(cells1, cells2)],
           nrow = length(genes.use),
           ncol = length(c(cells1, cells2)))
  colnames(expr) <- c(cells1, cells2)
  rownames(expr) <- genes.use
  
  
  if (nrow(expr) > 1) {
    #expr <- pheatmap:::scale_rows(expr)
    if (length(table(expr[, cells1])) > 1) {
      ph1 <- pheatmap(expr[, cells1], cluster_rows = FALSE, silent = TRUE)
    } else{
      ph1 <- list()
      ph1[["tree_col"]] <- list()
      ph1[["order"]] <- 1:length(cells1)
    }
    
    if (length(table(expr[, cells2])) > 1) {
      ph2 <- pheatmap(expr[, cells2], cluster_rows = FALSE, silent = TRUE)
    } else{
      ph2 <- list()
      ph2[["tree_col"]] <- list()
      ph2[["order"]] <- 1:length(cells2)
    }
    
    expr <-
      expr[, c(cells1[ph1$tree_col$order], cells2[ph2$tree_col$order])]
    
    #expr <- expr[rowSums(expr)> 0,]
  } else{
    rn <- rownames(expr)
    cn <- c(cells1[order(expr[, cells1])], cells2[order(expr[, cells2])])
    
    expr <- matrix(expr[, cn], nrow = 1)
    rownames(expr) <- rn
    colnames(expr) <- cn
  }
  
  
  anno <-
    data.frame(type = c(rep(types[1], length(cells1)), rep(types[2], length(cells2))))
  rownames(anno) <- c(cells1, cells2)
  
  
  if (is.null(colors)) {
    colors <- ggsci::pal_locuszoom()(2)
  }
  
  names(colors) <- types
  
  anno_col <- list(type = colors)
  
  
  if (nrow(expr) > 1) {
    
    
    if(useCheatmap){
      
      hname="expression"
      hcol = colorRamp2(seq(min(expr), quantile(expr,0.99), length.out = 100), 
                        colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100))
      
      if(scale=="row"){
        hname="z-score"
        rn <- rownames(expr)
        expr <- pheatmap:::scale_rows(expr)
        rownames(expr) <- rn
        
        hcol = colorRamp2(seq(-2, 2, length.out = 100), 
                          colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100))
      }
      
      ha_column = HeatmapAnnotation(df = anno, col = anno_col)
      
      H <- Heatmap(expr,
                   col = hcol,
                   top_annotation = ha_column,
                   cluster_columns = F,
                   cluster_rows = T,
                   show_column_names = F,
                   show_row_dend = F,
                   column_title = ttl,
                   name = hname)
      return(H)
      
    }else{
      col_breaks <- NA
      if (scale == "row") {
        col_breaks <- seq(-4, 4, length.out = 100)
      }
      ph <- pheatmap(
        expr,
        cluster_cols = FALSE,
        show_colnames = FALSE,
        annotation_col = anno,
        breaks = col_breaks,
        annotation_colors = anno_col,
        scale = scale,
        main = ttl
      )
      return(ph)  
    }
    
    
    
  } else{
    ha_column = HeatmapAnnotation(df = anno, col = anno_col)
    q <- quantile(expr, 0.99)
    hcol = colorRamp2(seq(0, q, length.out = 100), colorRampPalette(rev(brewer.pal(
      n = 9, name = "RdYlBu"
    )))(100))
    
    H <- Heatmap(
      expr,
      cluster_columns = F,
      cluster_rows = F,
      show_column_names = F,
      col = hcol,
      top_annotation = ha_column
    )
    draw(H)
    return(NULL)
  }
}

PredictCellsIdentity <-
  function(object,
           classifier,
           training.genes = NULL,
           training.classes = NULL,
           new.data = NULL,
           ...) {
    # build the classifier
    if (missing(classifier)) {
      classifier <-
        BuildRFClassifier(object,
                          training.genes = training.genes,
                          training.classes = training.classes,
                          ...)
    }
    # run the classifier on the new data
    features <- classifier$forest$independent.variable.names
    genes.to.add <- setdiff(features, rownames(new.data))
    data.to.add <-
      matrix(0,
             nrow = length(genes.to.add),
             ncol = ncol(new.data))
    rownames(data.to.add) <- genes.to.add
    new.data <- rbind(new.data, data.to.add)
    new.data <- new.data[features, ]
    new.data <- as.matrix(t(new.data))
    cat("Running Classifier ...", file = stderr())
    prediction <- predict(classifier, new.data)
    new.classes <- prediction$predictions
    return(new.classes)
  }




plotGlobalMarkersHeatmap <-
  function(sobj, markers, markerClusters, scale = TRUE) {
    require(ComplexHeatmap)
    require(circlize)
    require(stringr)
    require(RColorBrewer)
    
    if (!all(markers %in%  rownames(sobj@data))) {
      stop("Check that all of the provided genes are in sobj")
    }
    
    
    
    info <-
      str_match(NAc_cells_all@ident, pattern = "(\\w+)_(\\d+)")
    info <- as.data.frame(info, stringsAsFactors = FALSE)
    
    colnames(info) <- c("Ident", "Cell", "Cluster")
    info$Cluster <- as.numeric(info$Cluster)
    rownames(info) <- names(NAc_cells_all@ident)
    
    
    setorder(info, Cell, Cluster)
    ord <- order(as.character(sobj@ident))
    
    
    
    cells_cluster_anno <- data.frame(CellType = info[, 2],
                                     Cluster = info[, 3])
    
    rownames(cells_cluster_anno) <- rownames(info)
    
    
    cellType.col <- ggsci::pal_igv()(length(unique(info[, 2])))
    names(cellType.col) <- unique(info[, 2])
    
    Cluster.col <-
      ggsci::pal_d3(palette = "category20c")(length(unique(info[, 3])))
    names(Cluster.col) <- unique(info[, 3])
    
    ann_colors <- list(Cluster = Cluster.col,
                       CellType = cellType.col)
    
    
    if (scale) {
      tmp <- pheatmap:::scale_rows(sobj@data[markers, ord])
      tmp <- as.matrix(tmp)
      hcol = colorRamp2(seq(-4, 4, length.out = 100), colorRampPalette(rev(brewer.pal(
        n = 11, name = "RdYlBu"
      )))(100))
      hname <- "relative_expression"
    } else{
      tmp <- as.matrix(sobj@data[markers, ord])
      q <- quantile(tmp, 0.99)
      hcol = colorRamp2(seq(0, q, length.out = 100),
                        colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100))
      hname = "log(expresion+1)"
    }
    
    ha_column = HeatmapAnnotation(df = cells_cluster_anno, col = ann_colors)
    
    
    
    H <- Heatmap(
      tmp,
      col = hcol,
      split = markerClusters,
      top_annotation = ha_column,
      cluster_columns = F,
      show_column_names = F,
      show_row_names = T,
      cluster_rows = F,
      name = hname
    )
    
    
    draw(H)
    ## Add the clusters line separators
    nbCells <- table(sobj@ident)
    
    nbSplits <-
      ifelse(is.null(markerClusters), 1, length(unique(markerClusters)))
    
    for (i in 1:length(nbCells)) {
      sm <- sum(nbCells[1:i])
      
      for (slice in 1:nbSplits) {
        decorate_heatmap_body(hname, {
          x = sm / ncol(sobj@data)
          grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 1))
          grid.rect(gp = gpar(
            fill = "transparent",
            col = "black",
            lwd = 1
          ))
        }, slice = slice)
      }
    }
  }

PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

DotPlot_my <- function(object,
                       genes.plot,
                       cols.use = c("lightgrey", "blue"),
                       col.min = -2.5,
                       col.max = 2.5,
                       dot.min = 0,
                       dot.scale = 6,
                       group.by,
                       plot.legend = FALSE,
                       do.return = FALSE,
                       x.lab.rot = FALSE) {
  if (!missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <-
    data.frame(FetchData(object = object, vars.all = genes.plot))
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot %>% gather(key = genes.plot,
                          value = expression, -c(cell, id)) -> data.to.plot
  data.to.plot %>%
    group_by(id, genes.plot) %>%
    summarize(avg.exp = mean(expm1(x = expression)),
              pct.exp = PercentAbove(x = expression, threshold = 0)) -> data.to.plot
  data.to.plot %>%
    ungroup() %>%
    group_by(genes.plot) %>%
    mutate(avg.exp.scale = scale(x = avg.exp)) %>%
    mutate(avg.exp.scale = MinMax(data = avg.exp.scale,
                                  max = col.max,
                                  min = col.min)) ->  data.to.plot
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot,
                                    levels = rev(x = sub(
                                      pattern = "-",
                                      replacement = ".",
                                      x = genes.plot
                                    )))
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  p <-
    ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
    geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_gradient(low = cols.use[1], high = cols.use[2]) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}


## used for dubugging

plotSplits <- function(sobj, splits) {
  #splits <- splitClusters(NAc_Neuro_IN_sobj)
  
  NewIdent <- rep(0, length(sobj@ident))
  
  pos1 <- which(colnames(sobj@data) %in% splits$clus1Cells)
  pos2 <- which(colnames(sobj@data) %in% splits$clus2Cells)
  
  NewIdent[pos1] <- "1"
  NewIdent[pos2] <- "2"
  
  sobj <- SetIdent(sobj, ident.use = NewIdent)
  
  clusCols <- pal_aaas()(length(unique(NewIdent)))
  names(clusCols) <- unique(NewIdent)
  
  plotMarkersHeatmap(
    sobj,
    markers = unique(splits$markers$gene),
    clustersColors = clusCols,
    clusterCols = F
    #markerClusters = factor(NAc_Neuro_IN_markersClus, levels=c(levels(NAc_Neuro_IN_sobj@ident),"cmp_4_5_6", "global") )
  )
}

SetIfNull <- function(x, default) {
  if (is.null(x = x)) {
    return(default)
  } else {
    return(x)
  }
}


NoGrid <- function(...) {
  no.grid <- theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   ...)
  return(no.grid)
}

MarkerViolonPlot <- function(features,
                             sobj,
                             cell.ident,
                             ident.order=NULL,
                             do.sort = FALSE,
                             y.max = 8,
                             size.x.use = 16,
                             size.y.use = 16,
                             size.title.use = 20,
                             adjust.use = 1,
                             point.size.use = 1,
                             cols.use = NULL,
                             gene.names = NULL,
                             y.log = FALSE,
                             x.lab.rot = FALSE,
                             y.lab.rot = FALSE,
                             legend.position = "right",
                             remove.legend = FALSE,
                             inGrid = FALSE,
                             nrow = 3,
                             add.jitter = FALSE) {
  require(lemon)
  
  AllData <- data.frame()
  
  
  for (feature in features) {
    data <- data.frame(FetchData(object = sobj, vars = feature))
    colnames(data) <- c("Expression")
    data$gene <- feature
    
    set.seed(seed = 42)
    data$ident <- cell.ident
    if (do.sort) {
      data$ident <- factor(x = data$ident,
                           levels = names(x = rev(x = sort(
                             x = tapply(
                               X = data[, feature],
                               INDEX = data$ident,
                               FUN = mean
                             )
                           ))))
    }
    
    if(!is.null(ident.order)){
      data$ident <- factor(data$ident, levels = ident.order)
    }
    
    if (y.log) {
      noise <- rnorm(n = length(x = data[, feature])) / 200
      data[, "Expression"] <- data[, "Expression"] + 1
    } else {
      noise <- rnorm(n = length(x = data[, "Expression"])) / 100000
    }
    data[, "Expression"] <- data[, "Expression"] + noise
    y.max <- SetIfNull(x = y.max, default = max(data[, "Expression"]))
    
    AllData <- rbind(AllData, data)
  }
  
  AllData$gene <- factor(AllData$gene, levels = features)
  
  boldtxt <- theme(text = element_text(face = "bold", size = 12))
  
  plot <- ggplot(data = AllData,
                 mapping = aes(x = factor(x = ident),
                               y = Expression)) +
    geom_violin(
      scale = "width",
      adjust = adjust.use,
      alpha = 0.9,
      trim = TRUE,
      mapping = aes(fill = factor(x = ident))
    ) +
    theme_classic() +
    theme(
      panel.spacing = unit(-1, "mm"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.ticks = element_blank(),
      
      
      legend.position = legend.position,
      axis.title.x = element_text(
        face = "bold",
        colour = "black",
        size = size.x.use
      ),
      axis.title.y = element_text(
        face = "bold",
        colour = "black",
        size = size.y.use
      ),
      strip.background = element_blank(), 
      strip.placement = "outside",
      strip.text.y = element_text(colour = "Black",face = "bold", angle = 30, size = 8)
    ) +
    xlab("Cluster") +
    guides(fill = guide_legend(title = NULL)) +
    boldtxt
  
  if (add.jitter) {
    plot <-
      plot + ggbeeswarm::geom_quasirandom(alpha = 0.6, size = point.size.use)
    #geom_jitter(height = 0, size = point.size.use)
  }
  #geom_jitter(height = 0, size = point.size.use) +
  
  #NoGrid() +
  #ggtitle(feature) +
  #theme(plot.title = element_text(size = size.title.use, face = "bold", color="black"))+
  
  if (y.log) {
    plot <- plot + scale_y_log10()
  } else {
    #plot <- plot + ylim(min(data[, feature]), y.max)
  }
  if (feature %in% gene.names) {
    if (y.log) {
      plot <- plot + ylab(label = "Log Expression level")
    } else {
      plot <- plot + ylab(label = "Log Expression level")
    }
  } else {
    plot <- plot + ylab(label = "")
  }
  if (!is.null(x = cols.use)) {
    plot <- plot + scale_fill_manual(values = cols.use)
  }
  if (x.lab.rot) {
    plot <-
      plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  if (y.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90))
  }
  if (remove.legend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  
  if (inGrid) {
    plot <- plot + facet_rep_grid(gene ~.,scales = "free_y")
  } else{
    plot <- plot + facet_wrap(~ gene, nrow = nrow, scales = "free_y")
  }

  return(plot)
}


sparse_pca <-
  function(x,
           n_pcs,
           mu = NULL,
           s = NULL,
           center_scale = TRUE) {
    if (is.null(mu) && center_scale)
      mu <- colMeans(x)
    if (is.null(s) && center_scale)
      s <- apply(x, 2, sd, na.rm = TRUE)
    
    if (center_scale) {
      s[s == 0] <- min(s[s > 0])
      svd_res <- irlba::irlba(x, n_pcs, center = mu, scale = s)
    } else {
      svd_res <- irlba::irlba(x, n_pcs)
    }
    
    # compute explained variance
    n <- dim(x)[1]
    variance_sum <-
      sum(apply(x, 2, var, na.rm = TRUE) / (s ^ 2)) # sample variance sum
    var_pcs <- svd_res$d ^ 2 / (n - 1) / variance_sum
    
    return(
      list(
        x = svd_res$u %*% diag(svd_res$d),
        rotation = svd_res$v,
        sdev = svd_res$d / sqrt(n - 1),
        tot_var = variance_sum,
        var_pcs = var_pcs
      )
    )
  }



ValidateClusters.my <- function(object,
                                pc.use = NULL,
                                top.genes = 30,
                                min.connectivity = 0.01,
                                acc.cutoff = 0.9,
                                verbose = TRUE) {
  # probably should refactor to make cleaner
  if (length(x = object@snn) > 1) {
    SNN.use <- object@snn
  } else {
    stop("SNN matrix required. Please run BuildSNN() to save the SNN matrix in
         the object slot")
  }
  if (is.null(pc.use)) {
    stop("pc.use not set. Please choose PCs.")
  }
  num.clusters.orig <- length(x = unique(x = object@ident))
  still_merging <- TRUE
  if (verbose) {
    connectivity <-
      Seurat:::CalcConnectivity(object = object)
    end <- length(x = connectivity[connectivity > min.connectivity])
    progress <- end
    status <- 0
  }
  # find connectedness of every two clusters
  while (still_merging) {
    connectivity <- Seurat:::CalcConnectivity(object = object)
    merge.done <- FALSE
    while (!merge.done) {
      m <- max(connectivity, na.rm = TRUE)
      mi <- which(x = connectivity == m, arr.ind = TRUE)
      c1 <- rownames(x = connectivity)[mi[, 1]]
      c2 <- rownames(x = connectivity)[mi[, 2]]
      if (m > min.connectivity) {
        acc <- RunClassifier.my(
          object = object,
          group1 = c1,
          group2 = c2,
          pcs = pc.use,
          num.genes = top.genes
        )
        # if classifier can't classify them well enough, merge clusters
        if (acc < acc.cutoff) {
          object <- SetIdent(
            object = object,
            cells.use = WhichCells(object = object, ident = c1),
            ident.use = c2
          )
          if (verbose) {
            progress <-
              length(x = connectivity[connectivity > min.connectivity])
            print(
              paste0(
                sprintf("%3.0f", (1 - progress / end) * 100),
                "% complete --- merge clusters ",
                c1,
                " and ",
                c2,
                ", classification accuracy of ",
                sprintf("%1.4f", acc)
              )
            )
          }
          merge.done <- TRUE
        } else {
          if (verbose & status == 5) {
            print(
              paste0(
                sprintf("%3.0f", (1 - progress / end) * 100),
                "% complete --- Last 5 cluster comparisons failed to merge, ",
                "still checking possible merges ..."
              )
            )
            status <- 0
          }
          status <- status + 1
          connectivity[c1, c2] <- 0
          connectivity[c2, c1] <- 0
        }
      } else {
        still_merging <- FALSE
        break
      }
    }
  }
  if (verbose) {
    print(
      paste0(
        "100% complete --- started with ",
        num.clusters.orig,
        " clusters, ",
        length(x = unique(x = object@ident)),
        " clusters remaining"
      )
    )
  }
  return(object)
  }




ValidateClusters_CCA.my <- function(object,
                                pc.use = NULL,
                                top.genes = 30,
                                min.connectivity = 0.01,
                                acc.cutoff = 0.9,
                                verbose = TRUE) {
  # probably should refactor to make cleaner
  if (length(x = object@snn) > 1) {
    SNN.use <- object@snn
  } else {
    stop("SNN matrix required. Please run BuildSNN() to save the SNN matrix in
         the object slot")
  }
  if (is.null(pc.use)) {
    stop("pc.use not set. Please choose PCs.")
  }
  num.clusters.orig <- length(x = unique(x = object@ident))
  still_merging <- TRUE
  if (verbose) {
    connectivity <-
      Seurat:::CalcConnectivity(object = object)
    end <- length(x = connectivity[connectivity > min.connectivity])
    progress <- end
    status <- 0
  }
  # find connectedness of every two clusters
  while (still_merging) {
    connectivity <-
      Seurat:::CalcConnectivity(object = object)
    merge.done <- FALSE
    while (!merge.done) {
      m <- max(connectivity, na.rm = TRUE)
      mi <- which(x = connectivity == m, arr.ind = TRUE)
      c1 <- rownames(x = connectivity)[mi[, 1]]
      c2 <- rownames(x = connectivity)[mi[, 2]]
      if (m > min.connectivity) {
        acc <- RunClassifier_CCA.my(
          object = object,
          group1 = c1,
          group2 = c2,
          CCs =  pc.use,
          num.genes = top.genes
        )
        # if classifier can't classify them well enough, merge clusters
        if (acc < acc.cutoff) {
          object <- SetIdent(
            object = object,
            cells.use = WhichCells(object = object, ident = c1),
            ident.use = c2
          )
          if (verbose) {
            progress <-
              length(x = connectivity[connectivity > min.connectivity])
            print(
              paste0(
                sprintf("%3.0f", (1 - progress / end) * 100),
                "% complete --- merge clusters ",
                c1,
                " and ",
                c2,
                ", classification accuracy of ",
                sprintf("%1.4f", acc)
              )
            )
          }
          merge.done <- TRUE
        } else {
          if (verbose & status == 5) {
            print(
              paste0(
                sprintf("%3.0f", (1 - progress / end) * 100),
                "% complete --- Last 5 cluster comparisons failed to merge, ",
                "still checking possible merges ..."
              )
            )
            status <- 0
          }
          status <- status + 1
          connectivity[c1, c2] <- 0
          connectivity[c2, c1] <- 0
        }
      } else {
        still_merging <- FALSE
        break
      }
    }
  }
  if (verbose) {
    print(
      paste0(
        "100% complete --- started with ",
        num.clusters.orig,
        " clusters, ",
        length(x = unique(x = object@ident)),
        " clusters remaining"
      )
    )
  }
  return(object)
  }

ValidateSpecificClusters.my <- function(object,
                                        cluster1 = NULL,
                                        cluster2 = 1,
                                        pc.use = 2,
                                        top.genes = 30,
                                        acc.cutoff = 0.9) {
  acc <- RunClassifier.my(
    object = object,
    group1 = cluster1,
    group2 = cluster2,
    pcs = pc.use,
    num.genes = top.genes
  )
  print(paste0("Comparing cluster ",
               cluster1,
               " and ",
               cluster2,
               ": Acc = ",
               acc))
  if (acc < acc.cutoff) {
    object <- SetIdent(
      object = object,
      cells.use = WhichCells(object = object, ident = cluster1),
      ident.use = cluster2
    )
    print(paste("merge cluster", cluster1, "and", cluster2))
    merge.done <- TRUE
  }
  return(object)
}

# Train an SVM classifier and return the accuracy after 5 fold CV
#
# @param object     Seurat object
# @param group1     One identity to train classifier on
# @param group2     Second identity to train classifier on
# @param pcs        Vector of PCs on which to base genes to train classifier on.
#                   Pulls top num.genes genes associated with these PCs
# @param num.genes  Number of genes to pull for each PC
# @return           Returns the accuracy of the classifier after CV

RunClassifier.my <-
  function(object, group1, group2, pcs, num.genes) {
    d1 <- WhichCells(object = object, ident = group1)
    d2 <- WhichCells(object = object, ident = group2)
    y  <- as.numeric(x = object@ident[c(d1, d2)]) - 1
    x  <- data.frame(as.matrix(t(x = object@data[ICTopGenes(
      object = object,
      ic.use = pcs,
      num.genes =
        num.genes,
      do.balanced = T
    ), c(d1, d2)])))
    xv <- apply(X = x, MARGIN = 2, FUN = var)
    x  <- x[, names(x = which(xv > 0))]
    # run k-fold cross validation
    ctrl <- trainControl(method = "repeatedcv", repeats = 5)
    set.seed(seed = 1500)
    model <- train(
      x = x,
      y = as.factor(x = y),
      formula = as.factor(x = y) ~ .,
      method = "gbm",
      trControl = ctrl
    )
    acc <- model$results[, 2]
    return(acc)
  }



RunClassifier_CCA.my <-
  function(object, group1, group2, CCs, num.genes) {
    d1 <- WhichCells(object = object, ident = group1)
    d2 <- WhichCells(object = object, ident = group2)
    y  <- as.numeric(x = object@ident[c(d1, d2)]) - 1
    y= as.factor(y)
    x  <- data.frame(t(x = as.matrix(object@data[DimTopGenes(
      object = object,
      dim.use = CCs,
      reduction.type = "cca",
      num.genes =num.genes,
      do.balanced = T
    ), c(d1, d2)])))
    xv <- apply(X = x, MARGIN = 2, FUN = var)
    x  <- x[, names(x = which(xv > 0))]
    # run k-fold cross validation
    ctrl <- trainControl(method = "repeatedcv", repeats = 5)
    set.seed(seed = 1500)
    x$Class = y
    model <- caret::train(
      #data = x, 
      x = model.matrix(Class ~., data= x),
      y = factor(x$Class),
      #formula = factor(y) ~ .,
      method = "gbm",
      trControl = ctrl
    )
    #acc <- model$results[, 2]
    acc = max(model$results$Accuracy)
    return(acc)
  }


PlotMarketGenesHeatmap <-
  function(sobj,
           sobj_mrks,
           power.thr = 0.7,
           pct.1.thr = 0.7,
           pct.2.thr = 0.2,
           is.pval = FALSE,
           pval.thr = 1e-3,
           nbGenes=5,
           cols=NULL,
           additional_markers= NULL) {
    mrks.sig <- data.frame()
    
    for (cls in unique(sobj_mrks$cluster)) {
      if (is.pval) {
        tmp <-
          subset(sobj_mrks,
                 cluster == cls &
                   pct.1 > pct.1.thr &
                   pct.2 < pct.2.thr & p_val < pval.thr)
      } else{
        tmp <-
          subset(sobj_mrks,
                 cluster == cls &
                   pct.1 > pct.1.thr &
                   pct.2 < pct.2.thr & power > power.thr)
      }
      
      
      if (nrow(tmp) == 0) {
        tmp <- subset(sobj_mrks, cluster == cls)
      }
      tmp <- head(tmp, nbGenes)
      
      mrks.sig <- rbind(mrks.sig, tmp)
    }
    
    
    if(is.null(cols)){
      cols <- generateColors(length(levels(sobj@ident)))
        #pal_d3(palette = "category20c")(length(levels(sobj@ident)))  
    }
    
    if(is.null(names(cols))){
      names(cols) <- levels(sobj@ident)  
    }
    
    
    if(!is.null(additional_markers)){
      toUsegenes <- c(mrks.sig$gene, additional_markers)
      clus_names <- c(mrks.sig$cluster, rep("General",length(additional_markers)))
    }else{
      toUsegenes <- mrks.sig$gene
      clus_names <- mrks.sig$cluster
    }
    
    
    plotMarkersHeatmap(
      sobj,
      markers = toUsegenes,
      clustersColors = cols,
      markerClusters = clus_names,
      #identlevels = clus_lvls,
      clusterCols = FALSE
    )
    
  }


generateColors <- function(nbCols){
  require(ggsci)
  pallets <- list(pal_d3(),ggsci::pal_simpsons(),ggsci::pal_futurama(),ggsci::pal_ucscgb(), pal_jco())
  
  cols <- c()
  i=1
  while(nbCols > 0){
    if(nbCols > 10){
      cols <- c(cols,pallets[[i]](10))
      nbCols <- nbCols - 10
      i=i+1
    }else{
      cols <- c(cols,pallets[[i]](nbCols))
      nbCols <- 0
    }
  }
  
  #cols <- Shuffle(cols)
  cols
}

plotInteractive2DTsne <- function(sobj, clusCols = NULL) {
  df <- sobj@tsne.rot
  df$cluster <- as.character(sobj@ident)
  
  if (!is.null(clusCols)) {
    p <-
      plot_ly(
        data = df,
        x = ~ tSNE_1,
        y = ~ tSNE_2,
        color = ~ cluster,
        colors = clusCols
      )
  } else{
    p <- plot_ly(
      data = df,
      x = ~ tSNE_1,
      y = ~ tSNE_2,
      color = ~ cluster
    )
  }
  p
}



addSampleInfo <- function(sobj) {
  require(stringr)
  if("data.info" %in%  slotNames(sobj)){
    sobj@data.info$Sample <-
      str_match(rownames(sobj@data.info), ".+\\.(.+)")[, 2]
  }else{
    sobj@meta.data$Sample <-
      str_match(rownames(sobj@meta.data), ".+\\.(.+)")[, 2]
  }
  
  
  return(sobj)
}


SetSamplesLabels_metaData <- function(sobj) {
  sobj@meta.data$treatment <- ""
  sobj@meta.data$Period <- ""
  
  pos <-
    which(sobj@meta.data$Sample %in% c("NAc_Sample1", "NAc_Sample4"))
  sobj@meta.data$treatment[pos] <- "Cocaine"
  sobj@meta.data$Period[pos] <- "Maintenance"
  
  
  pos <-
    which(sobj@meta.data$Sample %in% c("NAc_Sample2", "NAc_Sample3"))
  sobj@meta.data$treatment[pos] <- "Saline"
  sobj@meta.data$Period[pos] <- "Maintenance"
  
  
  pos <-
    which(sobj@meta.data$Sample %in% c("NAc_Sample5", "NAc_Sample7"))
  sobj@meta.data$treatment[pos] <- "Saline"
  sobj@meta.data$Period[pos] <- "withdraw_48h"
  
  
  pos <-
    which(sobj@meta.data$Sample %in% c("NAc_Sample6", "NAc_Sample8"))
  sobj@meta.data$treatment[pos] <- "Cocaine"
  sobj@meta.data$Period[pos] <- "withdraw_48h"
  
  
  pos <-
    which(sobj@meta.data$Sample %in% c("NAc_Sample9", "NAc_Sample10"))
  sobj@meta.data$treatment[pos] <- "Cocaine"
  sobj@meta.data$Period[pos] <- "withdraw_15d"
  
  
  pos <-
    which(sobj@meta.data$Sample %in% c("NAc_Sample11", "NAc_Sample12"))
  sobj@meta.data$treatment[pos] <- "Saline"
  sobj@meta.data$Period[pos] <- "withdraw_15d"
  
  return(sobj)
}



#' Title
#'
#' @param sobj 
#'
#' @return
#' @export
#'
#' @examples
SetSamplesLabels <- function(sobj) {
  sobj@data.info$treatment <- ""
  sobj@data.info$Period <- ""
  
  pos <-
    which(sobj@data.info$Sample %in% c("NAc_Sample1", "NAc_Sample4"))
  sobj@data.info$treatment[pos] <- "Cocaine"
  sobj@data.info$Period[pos] <- "Maintenance"
  
  
  pos <-
    which(sobj@data.info$Sample %in% c("NAc_Sample2", "NAc_Sample3"))
  sobj@data.info$treatment[pos] <- "Saline"
  sobj@data.info$Period[pos] <- "Maintenance"
  
  
  pos <-
    which(sobj@data.info$Sample %in% c("NAc_Sample5", "NAc_Sample7"))
  sobj@data.info$treatment[pos] <- "Saline"
  sobj@data.info$Period[pos] <- "withdraw_48h"
  
  
  pos <-
    which(sobj@data.info$Sample %in% c("NAc_Sample6", "NAc_Sample8"))
  sobj@data.info$treatment[pos] <- "Cocaine"
  sobj@data.info$Period[pos] <- "withdraw_48h"
  
  
  pos <-
    which(sobj@data.info$Sample %in% c("NAc_Sample9", "NAc_Sample10"))
  sobj@data.info$treatment[pos] <- "Cocaine"
  sobj@data.info$Period[pos] <- "withdraw_15d"
  
  
  pos <-
    which(sobj@data.info$Sample %in% c("NAc_Sample11", "NAc_Sample12"))
  sobj@data.info$treatment[pos] <- "Saline"
  sobj@data.info$Period[pos] <- "withdraw_15d"
  
  return(sobj)
}



SetSamplesLabels_PFC <- function(sobj) {
  sobj@meta.data$treatment <- ""
  sobj@meta.data$Period <- ""
  
  pos <-
    which(sobj@meta.data$Sample %in% c("PFCSample1", "PFCSample4"))
  sobj@meta.data$treatment[pos] <- "Cocaine"
  sobj@meta.data$Period[pos] <- "Maintenance"
  
  
  pos <-
    which(sobj@meta.data$Sample %in% c("PFCSample2", "PFCSample3"))
  sobj@meta.data$treatment[pos] <- "Saline"
  sobj@meta.data$Period[pos] <- "Maintenance"
  
  
  pos <-
    which(sobj@meta.data$Sample %in% c("PFCSample5", "PFCSample7"))
  sobj@meta.data$treatment[pos] <- "Saline"
  sobj@meta.data$Period[pos] <- "withdraw_48h"
  
  
  pos <-
    which(sobj@meta.data$Sample %in% c("PFCSample6", "PFCSample8"))
  sobj@meta.data$treatment[pos] <- "Cocaine"
  sobj@meta.data$Period[pos] <- "withdraw_48h"
  
  
  pos <-
    which(sobj@meta.data$Sample %in% c("PFCSample9", "PFCSample10"))
  sobj@meta.data$treatment[pos] <- "Cocaine"
  sobj@meta.data$Period[pos] <- "withdraw_15d"
  
  
  pos <-
    which(sobj@meta.data$Sample %in% c("PFCSample11", "PFCSample12"))
  sobj@meta.data$treatment[pos] <- "Saline"
  sobj@meta.data$Period[pos] <- "withdraw_15d"
  
  return(sobj)
}




highlightClusters <- function(sobj,
                              clustersID = NA,
                              clus_cols = NULL) {
  if (is.null(sobj)) {
    stop("No Seurat object specificied")
  }
  
  if (is.na(clustersID) | is.null(clustersID)) {
    stop("No Clusters specified")
  }
  
  
  newIdent <-  as.character(sobj@ident)
  
  newIdent <-
    paste0("clust_", newIdent) #To avoid any problem if a clusters is named 2
  clustersID <- paste0("clust_", clustersID)
  
  newIdent[newIdent %in% clustersID] = "2"
  newIdent[newIdent != "2"] <- "1"
  
  sobj <- SetIdent(sobj, ident.use = newIdent)
  
  
  if (is.null(clus_cols)) {
    clus_cols <- c("gray", "red")
  }
  
  
  TSNEPlot(sobj, colors.use = clus_cols)
}





plotPseudoTimeHeatmap <-
  function(sobj,
           dpt,
           branches = 1:3,
           add.undeceided = FALSE,
           dir.branch= c(F,F,T),
           clusCols,
           heatcolor=viridis(10, direction = -1, option = "B"),
           order_by=1
           ) {
    NB_brances_genes <- data.frame()
    pt_orders <- list()
    
    sde_ord=NULL
    
    
    for (b in branches) {
      
      branch_pseudoTime <- destiny:::dpt_for_branch(dpt, branch_id = b)
      
      
      if(length(branches) >1){
        cells_in_brach <- which(dpt@branch[, 1L] == b)
      }else{
        cells_in_brach <- 1:nrow(dpt@branch)
      }
     
      # expr <- as.matrix(sobj@data[,cells_in_brach])
      expr <- as.matrix(sobj@data[sobj@var.genes, cells_in_brach])
      nbCells <- apply(expr, 1, function(x) {
        sum(x > 0)
      })
      
      expr <- expr[nbCells > length(cells_in_brach) / 10, ]
      
      sde <- switchde(expr,
                      branch_pseudoTime[cells_in_brach], verbose = TRUE)
      
      ord <- order(sde$qval, decreasing = F)
      head(sde[ord,])
      
      qval_cutoff <- 1e-3
      
      sde_sig <- head(subset(sde[ord,], qval < qval_cutoff), 200)
      NB_brances_genes <- rbind(NB_brances_genes, sde_sig)
      
      
      br_ord <- order(branch_pseudoTime[cells_in_brach], decreasing = dir.branch[b])
      
      
      if(b %in% names(pt_orders)){
        b <- paste0(b,"_1")
      }
      
      pt_orders[[as.character(b)]] <-
        colnames(sobj@data)[cells_in_brach][br_ord]
      
      if(b == order_by){
        sde_ord<- sde
      }
    }
    
    
    if (add.undeceided) {
      ## Get the cells in the decision points
      cells_in_brach <- which(is.na(dpt@branch[, 1L]))
      br_ord <- order(branch_pseudoTime[cells_in_brach])
      pt_orders[["undeceided"]] <-
        colnames(sobj@data)[cells_in_brach][br_ord]
    }
    
    
    ## Plot heatmap
    dynamicGenes <- unique(NB_brances_genes$gene)
    
    Pt_heatmaps_all <- NULL
    for (b in names(pt_orders)) {
      #cells_in_brach <- which(dpt@branch[, 1L] == b)
      
      pseudo_expr <- as.matrix(sobj@data[dynamicGenes, pt_orders[[b]]])
      
      rownames(pseudo_expr) <- NULL
      #pseudo_expr <- apply(pseudo_expr,2,EnrichedHeatmap::default_smooth_fun)
      
      
      
      cluster_anno <- data.frame(
        Cluster = factor(as.character(sobj@ident[pt_orders[[b]]]),
                         levels = levels(sobj@ident),ordered = T)
      )
      
      rownames(cluster_anno) <- pt_orders[[b]]
      
      
      ann_colors <-  list(
        Cluster = clusCols
      )
      
      
      ha_column = HeatmapAnnotation(df =cluster_anno,
                                         col= ann_colors)
      
      
      
      #pseudo_expr <- apply(pseudo_expr,2,EnrichedHeatmap::default_smooth_fun)
      
      
      if (b == names(pt_orders)[length(pt_orders)]) {
        rownames(pseudo_expr) <- dynamicGenes #c(sde_sig$gene,"Stmn2","Mfge8","Tk1")
      }
      
      #rownames(pseudo_expr) <- dynamicGenes #c(sde_sig$gene,"Stmn2","Mfge8","Tk1")
      H <- Heatmap(
        pseudo_expr,
        top_annotation = ha_column,
        show_row_dend = F,
        show_column_names = F,
        cluster_columns = F,
        col = heatcolor,
        name = paste("branch", b)
      )
      
      
      
      if (is.null(Pt_heatmaps_all)) {
        Pt_heatmaps_all <- H
      } else{
        Pt_heatmaps_all <- Pt_heatmaps_all + H
      }
    }
    
    draw(Pt_heatmaps_all)
    
    for (b in names(pt_orders)) {
      decorate_heatmap_body(paste("branch", b), {
        grid.rect(gp = gpar(col = "black", fill = NA))
      }, slice = 1)
    }
    
    
    if(length(Pt_heatmaps_all)==1){
      return(pseudo_expr)
    }
    #result <- c(pt_orders,Pt_heatmaps_all)
}


plotClustersMarkersSignal <- function(sobj, markers,cls_levels =NULL){
  
  require(reshape)
  require(ggjoy)
  
  if(is.null(sobj) | class(sobj) != "seurat"){
    stop("No Seurat object provided")
  }
  
  if(!all(markers %in% rownames(sobj@data))){
    stop("Check that all the gene names are in sobj@data")
  }
  
  ## Get the expression matrix
  mat <- as.matrix(sobj@data[markers,])
  
  clus_order <- order(sobj@ident)
  
  cell_ident <- as.character(sobj@ident)[clus_order]
  mat <- mat[,clus_order]
  mat <- as.data.frame(mat)
  
  mlt <- melt(mat)
  
  mlt$cluster <- rep(cell_ident,each=nrow(mat))
  mlt$gene <- rep(rownames(mat),ncol(mat))
  
  ## For each cluster get the cell order 
  
  cls_orders <- c()
  
  for(cls in unique(cell_ident)){
    
    tmp <- mat[,cell_ident == cls]
    
    ord <- hclust(dist(t(tmp)))$order
  
    cls_orders <- c(cls_orders,ord)
  }
  
  mlt$order <- rep(cls_orders,nrow(mat))
  
  mlt$gene <- factor(mlt$gene, levels=markers)
  
  if(!is.null(cls_levels)){
    mlt$cluster= factor(mlt$cluster,levels = cls_levels)
  }
  
  
  
  
  
  mlt$baseline = rep(0,nrow(mlt))
  ggplot(mlt, aes(x=order,y=baseline, height =value ,fill=cluster,colour=cluster)) + geom_joy(stat = "identity", scale = 1) + 
    #geom_bar(stat="identity",width = 1) +
    #coord_cartesian(ylim = c(0,6)) +
    facet_grid(gene ~cluster,scales = "free") + 
    theme_bw() + 
    theme(panel.spacing = unit(0, "lines"),axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),panel.background = element_rect(fill = "#ececec"), strip.text.y = element_text(angle = 0,face = "bold",colour = "black",hjust = 0),strip.background = element_blank()) + 
    scale_color_d3() +xlab("")  + 
    scale_fill_d3(palette = "category20c")
  
}



getDEGStats <- function(clust_DEG, ttl="", 
                        p.val=0.01,
                        FC=2, 
                        emp_FC=log2(1.5)){
  
  
  DE_genes_cpt <- data.frame()
  
  for(clus in names(clust_DEG)){
    for(period in names(clust_DEG[[clus]])){
      
      if(class(clust_DEG[[clus]][[period]]) == "data.frame"){
        pos <- which(is.na(clust_DEG[[clus]][[period]]$log2FC))
        clust_DEG[[clus]][[period]]$log2FC[pos] <- -clust_DEG[[clus]][[period]]$emp_log2FC[pos]
        
        deg_up <- subset(clust_DEG[[clus]][[period]], pvalue < p.val & log2FC > FC & emp_log2FC > emp_FC)
        deg_down <- subset(clust_DEG[[clus]][[period]], pvalue < p.val & log2FC < -FC &  emp_log2FC < -emp_FC)
        
        
        df <- data.frame(cellype=clus, period=period,
                         nbGenes=c(nrow(deg_up),-nrow(deg_down)),
                         direction=c("Up","Down"))
        
        DE_genes_cpt <- rbind(DE_genes_cpt,df)  
      }
    }
  }
  
  
  #DE_genes_cpt$cellype <- as.character(DE_genes_cpt$cellype)
  #DE_genes_cpt$cellype <- gsub("IN_12","Sst",DE_genes_cpt$cellype)
  #DE_genes_cpt$cellype <- gsub("IN_1","Chat",DE_genes_cpt$cellype)
  #DE_genes_cpt$cellype <- gsub("IN_9","Th",DE_genes_cpt$cellype)
  #DE_genes_cpt$cellype <- factor(DE_genes_cpt$cellype)
  
  p <- ggplot(DE_genes_cpt, aes(x=cellype,fill=direction))+ 
    geom_bar(data=subset(DE_genes_cpt,direction == "Up"),aes(y=nbGenes),width=0.3, stat = "identity", position = "dodge",color="black") +
    geom_bar(data=subset(DE_genes_cpt,direction == "Down"),aes(y=nbGenes),width=0.3,stat = "identity", position = "dodge",color="black") +
    geom_hline(yintercept = 0,colour = "black") + theme_bw() + scale_fill_lancet()+
    facet_grid(period~.) + ylab("") +
    theme(text = element_text(face = "bold",color = "black"),
          axis.text = element_text(face = "bold",color = "black"), 
          strip.background = element_blank(),
          strip.text = element_text(face = "bold",color = "black"),
          plot.title = element_text(hjust = 0.5)
          #panel.border = element_rect(colour = "black")
    ) +
    ggtitle(ttl)
  
  return(p)
}


M3DropDifferentialExpression.my <-  function (expr_mat, mt_method = "bon", mt_threshold = 0.05, suppress.plot = FALSE) {
  BasePlot <- M3Drop:::bg__dropout_plot_base(expr_mat, xlim = NA, suppress.plot = suppress.plot)
  MM <- M3Drop:::bg__fit_MM(BasePlot$p, BasePlot$s)
  if (!suppress.plot) {
    sizeloc <- M3Drop:::bg__add_model_to_plot(MM, BasePlot, lty = 1, 
                                     lwd = 2.5, col = "black", legend_loc = "topright")
  }
  DEoutput <- M3Drop:::bg__test_DE_K_equiv(expr_mat, fit = MM)
  sig <- which(p.adjust(DEoutput$pval, method = mt_method) < 
                 mt_threshold)
  DEgenes <- rownames(expr_mat)[sig]
  DEgenes <- DEgenes[!is.na(DEgenes)]
  if (!suppress.plot) {
    M3Drop:::bg__highlight_genes(BasePlot, DEgenes)
  }
  TABLE <- data.frame(Gene = DEgenes,
                      dropout = BasePlot$p[DEgenes],
                      log10_expr = BasePlot$s[DEgenes],
                      p.value = DEoutput$pval[sig], 
                      q.value = p.adjust(DEoutput$pval, method = mt_method)[sig])
  return(TABLE)
}



getMotifs <- function(peaks, only.groups = FALSE, bg = NULL){
  require(PWMEnrich)
  require(PWMEnrich.Mmusculus.background)
  require(BSgenome.Mmusculus.UCSC.mm9)
  data("PWMLogn.mm9.MotifDb.Mmus")
  
  peak_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9, 
                     seqnames(peaks),
                     start(peaks),
                     end(peaks))
  
  if(is.null(bg)){
    bg = PWMLogn.mm9.MotifDb.Mmus
  }
  
  res <- motifEnrichment(peak_seq, bg, group.only = only.groups)
  
  return(res)
}

bimodalirityTest <- function(bimodalData) {
  TEST <- dip.test(bimodalData)
  return(TEST$statistic[[1]])   # return(TEST$p.value[[1]])    to get the p value
}



getExpressedGenes <- function(sigGenes, expr){
  
  require(mclust)
  
  if(is.null(rownames(expr))){
    stop("rownames empty")
  }

  # This list will hold the names of cells in which each gene is expressed.
  CellsExpressedPerGene <- list()
  
  nGenes <- nrow(sigGenes)  
  for(gene in rownames(sigGenes)){
    # get gene expression
    x = as.numeric(expr[gene,])
    
    # cluster it into two
    res =Mclust(x,G = 2,verbose = F)
    
    # get the ID of the cluster with the expressed genes
    largestMean <- which.max(res$parameters$mean)
    clus = names(res$parameters$mean)[largestMean]
    
    ## get the cells in this cluster
    pos <- which(res$classification == clus)
  
    CellsExpressedPerGene[[gene]] <- colnames(expr)[pos]
    
    
  }
  
  return(CellsExpressedPerGene)
}



getCellsCorrelationPvalue <- function(exprs, Cells_cors, nbIter=1000){
  
  
  diag(Cells_cors) <- 0
  Cells_cor_smry <- summary(Matrix(Cells_cors,sparse = T))
  
  Cells_cor_smry$i <- rownames(Cells_cors)[Cells_cor_smry$i]
  Cells_cor_smry$j <- rownames(Cells_cors)[Cells_cor_smry$j]
  
  
  sampled_cor <- c()
  
  for(i in 1:nbIter){
    
    ## pick two genes 
    sampled_genes <- sample(colnames(exprs),2)
    tmp_cor <-  cor(exprs[,sampled_genes[1]], exprs[,sampled_genes[2]])
    sampled_cor <- c(sampled_cor, tmp_cor)
  }
  
  
  ## calculcate pvalues 
  
  Cells_cor_smry$pvalue <- sapply(Cells_cor_smry$x, function(x){
    if(x>=0){
      res <- length(sampled_cor[sampled_cor>x])/length(sampled_cor[sampled_cor>=0])
    }else{
      res <- length(sampled_cor[sampled_cor<x])/length(sampled_cor[sampled_cor<0])
    }
    res
  })
  
  Cells_cor_smry$qvalue <- p.adjust(Cells_cor_smry$pvalue,method = "BH")
  
  Cells_cor_smry
  
}

getCorrelationPvalue <- function(exprs, genes_cors, nbIter=1000){
  

  diag(genes_cors) <- 0
  varGenes_sig_cor_smry <- summary(Matrix(genes_cors,sparse = T))
  
  varGenes_sig_cor_smry$i <- rownames(genes_cors)[varGenes_sig_cor_smry$i]
  varGenes_sig_cor_smry$j <- rownames(genes_cors)[varGenes_sig_cor_smry$j]
  
  
  sampled_cor <- c()
  
  for(i in 1:nbIter){
    
    ## pick two genes 
    sampled_genes <- sample(rownames(exprs),2)
    tmp_cor <-  cor(exprs[sampled_genes[1],], exprs[sampled_genes[2],])
    sampled_cor <- c(sampled_cor, tmp_cor)
  }
  
  sampled_cor <- sampled_cor[!is.na(sampled_cor)]
  ## calculcate pvalues 
  
  varGenes_sig_cor_smry$pvalue <- sapply(varGenes_sig_cor_smry$x, function(x){
    if(x>=0){
      res <- length(sampled_cor[sampled_cor>x])/length(sampled_cor[sampled_cor>=0])
    }else{
      res <- length(sampled_cor[sampled_cor<x])/length(sampled_cor[sampled_cor<0])
    }
    res
  })
  
  varGenes_sig_cor_smry$qvalue <- p.adjust(varGenes_sig_cor_smry$pvalue,method = "BH")
  
  varGenes_sig_cor_smry
  
}


plotAnnotatedExpr <- function(Expr, ident){
  
  anno_col <- list(cluster = generateColors(20))
  names(anno_col$cluster) <-  unique(ident)
  
  ord <- order(ident)
  #ord <- names(ident)[ord]
  
  
  df <- data.frame(cluster = ident[ord])
  rownames(df) <- colnames(Expr)
  
  pheatmap(Expr[ord,ord],show_colnames = F,show_rownames = F,
           cluster_cols = F, cluster_rows =  F,
           annotation_col = df,annotation_colors = anno_col)
    
}


doWGCNA_clustering <- function(sobj, 
                               power=NULL,
                               deepSplit=4,
                               minClustSize=20
                               ){
  require(WGCNA)
  require(distances)
  suppressMessages(allowWGCNAThreads())
  
  
  if(is.null(power)){
    message("No power value was provided.")
    message("Automatic selection using soft-threshold")
    powers <- 1:30
    sft <- pickSoftThreshold(as.matrix(sobj@data[sobj@var.genes,]), 
                             powerVector = powers, verbose = 5)  
    power <- selectSoftThr(sft, 0.8)
  }
  
  
  A <- adjacency(sobj@data[sobj@var.genes,],power = power,distFnc = "distances")
  dissTOM = TOMdist(A)
  
  d2 = distances(dissTOM)
  d2  = as.dist(d2)
  #geneTree = flashClust(as.dist(dissTOM), method = "average");
  geneTree = hclust(d2,method = "ward.D2");
  # here we define the modules by cutting branches
  moduleLabelsManual1 = cutreeDynamic(dendro = geneTree, 
                                      distM = dissTOM, 
                                      method = "hybrid", verbose = 2,
                                      deepSplit = deepSplit,
                                      pamRespectsDendro = F, 
                                      minClusterSize = minClustSize)
  
  dynamicColors = labels2colors(moduleLabelsManual1)
  
  sobj <- SetIdent(sobj,ident.use = moduleLabelsManual1)
  return(sobj)
}


doWGCNA_Seurat_clustering <- function(sobj, 
                               power=NULL,
                               algo=1, res=1
){
  require(WGCNA)
  suppressMessages(allowWGCNAThreads())
  
  
  if(is.null(power)){
    message("No power value was provided.")
    message("Automatic selection using soft-threshold")
    powers <- 1:30
    sft <- pickSoftThreshold(as.matrix(sobj@data[sobj@var.genes,]), 
                             powerVector = powers, verbose = 5)  
    power <- selectSoftThr(sft, 0.8)
  }
  
  
  A <- adjacency(sobj@data[sobj@var.genes,],power = power)
  dissTOM = TOMdist(A)
  
  colnames(dissTOM) <- sobj@cell.names
  rownames(dissTOM) <- sobj@cell.names
  
  dissTOM = 1-dissTOM
  q = quantile(dissTOM,0.5)
  dissTOM[dissTOM>q] = 1
  #geneTree = flashClust(as.dist(dissTOM), method = "average");
  sobj <- FindClusters(sobj,resolution = res,distance.matrix = dissTOM,algorithm = algo)
  
  return(sobj)
}



getVariableGenes <- function(sobj, diptest.qval = 0.04,
                             expr.thr = 0.8, 
                             blacklist.genes = NULL,
                             cor.qval = 0.05){
  
  sobj <- FindVariableGenes(sobj,do.plot = F)
  inhib_varGenes <- sobj@hvg.info
  
  ######################################################
  ## Step1: Get the list of candidate variable genes  ##
  ######################################################
  
  inhib_varGenes$dipTest <- apply(sobj@data[rownames(inhib_varGenes),],1,bimodalirityTest)
  
  
  q_threshold <- quantile(inhib_varGenes$gene.dispersion.scaled, 0.95)
  
  varGenes <- subset(inhib_varGenes, gene.dispersion.scaled > q_threshold & dipTest > diptest.qval )
  
  varGenes_exrCells <- getExpressedGenes(varGenes, sobj@data)
  varGenes_exrCells_perc <- sapply(varGenes_exrCells, length)/ncol(sobj@data)
  
  varGenes_sig <- varGenes[varGenes_exrCells_perc < expr.thr,]
  
  varGenes_sig <- varGenes_sig[grep("mt-",rownames(varGenes_sig),invert = T),]
  
  
  if(!is.null(blacklist.genes)){
    varGenes_sig <- varGenes_sig[! rownames(varGenes_sig) %in% IEG$V1,]  
  }
  
  ################################################################################
  ## Step2: Get the genes that show a significant correlation/anti-correlation  ##
  ################################################################################
  
  varGenes_sig_cor <- cor(t(as.matrix(sobj@data[rownames(varGenes_sig),])))
  diag(varGenes_sig_cor) <- 0
  
  varGenes_sig_cor_smry <- getCorrelationPvalue(sobj@data, varGenes_sig_cor)
  
  
  varGenes_sig_cor_smry_sig <- subset(varGenes_sig_cor_smry, qvalue < cor.qval)
  
  varGenes_sig_cor_smry_sig <- unique(c(varGenes_sig_cor_smry_sig$i, varGenes_sig_cor_smry_sig$j ))
  
  return(varGenes_sig_cor_smry_sig)
}



clustSim_Hao <- function(sobj){
  
  library("reshape2")
  expr.matrix <- as.data.frame(as.matrix(sobj@data))
  expr.matrix$gene <- rownames(expr.matrix)
  expr.matrix <- melt(expr.matrix, id = "gene")
  tt <- sobj@meta.data
  tt$cell <- rownames(tt)
  expr.matrix2 <- expr.matrix %>% mutate(cell = variable) %>% left_join(tt[, c(5:7,12:13)], by = "cell")
  head(expr.matrix2)
  expr.matrix3 <- expr.matrix2 %>% group_by(gene, res.6) %>% dplyr::summarise(meanExp = mean(value)) %>% 
    ungroup
  expr.matrix3 <- as.data.frame(expr.matrix3)
  expr.matrix4 <- dcast(expr.matrix3, gene ~ res.6, value.var = "meanExp")
  head(expr.matrix4)
  rownames(expr.matrix4) <- expr.matrix4$gene
  expr.matrix4 <- as.matrix(expr.matrix4[, -1])
  expr.matrix4 <- t(scale(t(expr.matrix4), scale = T))  ## scale= T
  dim(expr.matrix4)
  # apply(expr.matrix4,1,sd) round(apply(expr.matrix4,1,sum),2)
  expr.matrix.mak <- expr.matrix4[unique(sobj.mak$gene), ]
  dim(expr.matrix.mak)  # 2,090 genes across 40 clusters
  
  
  d <- dist(as.matrix(t(expr.matrix.mak)))
  hcc <- hclust(d,method="ward.D")
}





Get_timecourse_DEG <- function(sobj, logFC_thr = log(1.5), pval_thr = 0.01, pct.1_thr = 0.4){
  
  
  ## check that the levels exist
  if(!all(levels(sobj@ident) %in% c("Maintenance","withdraw_48h","withdraw_15d"))){
    stop("Please make sure that the levels `Maintenance, withdraw_48h, withdraw_15d`  are defined")
  }
  
  
  idents_toUse <- names(which(table(sobj@ident)>=3))
  
  betweenPhases_DEG <- list()
  for(p1 in idents_toUse){
    for(p2 in idents_toUse){
      n <- paste0(p1,"_",p2)
      betweenPhases_DEG[[n]] <- FindMarkers(sobj,ident.1 = p1,ident.2 = p2,test.use = "bimod")
    }
  }
  
  betweenPhases_DEG_sig <- data.frame()
  
  for(p1 in idents_toUse){
    for(p2 in idents_toUse){
      n <- paste0(p1,"_",p2)
      #betweenPhases_DEG[[n]] <- FindMarkers(cluster35_saline,ident.1 = p1,ident.2 = p2,test.use = "bimod")
      
      if(!is.null(betweenPhases_DEG[[n]])){
        tmp <- subset(betweenPhases_DEG[[n]], avg_logFC > logFC_thr & p_val_adj < pval_thr & pct.1 > pct.1_thr )
        tmp$gene <- rownames(tmp)
        betweenPhases_DEG_sig <- rbind(betweenPhases_DEG_sig, tmp)  
      }   
    }
  }
  
  return(betweenPhases_DEG_sig)
  
}


RUVs_mat <-  function(x, cIdx, k, scIdx, round=TRUE, epsilon=1, tolerance=1e-8, isLog=FALSE) {
  
  require(foreach)
  require(bigmemory)
  require(irlba)
  require(Matrix)
  # if(!isLog && !all(.isWholeNumber(x))) {
  #   warning(paste0("The expression matrix does not contain counts.\n",
  #                  "Please, pass a matrix of counts (not logged)
  #                  or set isLog to TRUE to skip the log transformation"))
  # }
  
  if(isLog) {
    Y <- t(x)
  } else {
    message("log")
    Y <- t(log(x+epsilon))
  }
  
  scIdx <- scIdx[rowSums(scIdx > 0) >= 2, , drop = FALSE]
  Yctls <- matrix(0, prod(dim(scIdx)), ncol(Y))
  m <- nrow(Y)
  n <- ncol(Y)
  c <- 0
  
  
  Yctls = as.big.matrix(Yctls)
  Yctls_desc <- describe(Yctls)
  
  Y <- as.big.matrix(Y)
  Y_desc = describe(Y)
  message("doing parallel")
  foreach(ii = 1:nrow(scIdx),.combine = c) %dopar% {
    require(bigmemory)
    Yctls <- attach.big.matrix(Yctls_desc)
    Y <- attach.big.matrix(Y_desc)
    mu =  colMeans(Y[scIdx[ii, (scIdx[ii, ] > 0)], , drop = FALSE])
    
    for (jj in 1:(ncol(scIdx))) {
      if (scIdx[ii, jj] == -1)
        next
      c <- c + 1
      Yctls[c, ] <- Y[scIdx[ii, jj], , drop = FALSE] - mu
    }
    return(NULL)
  }
  
  # for (ii in 1:nrow(scIdx)) {
  #   mu =  colMeans(Y[scIdx[ii, (scIdx[ii, ] > 0)], , drop = FALSE])
  #
  #   for (jj in 1:(ncol(scIdx))) {
  #     if (scIdx[ii, jj] == -1)
  #       next
  #     c <- c + 1
  #     Yctls[c, ] <- Y[scIdx[ii, jj], , drop = FALSE] - mu
  #   }
  # }
  Yctls <- as.matrix(Yctls)
  Y <- as.matrix(Y)
  
  Yctls <- Yctls[rowSums(Yctls) != 0, ]
  Y <- rbind(Y, Yctls)
  
  sctl <- (m + 1):(m + nrow(Yctls))
  rm(Yctls)
  message("svd")
  svdRes <- irlba(Y[sctl, ], nu = 0, nv = k)
  k <- min(k, max(which(svdRes$d > tolerance)))
  a <- diag(as.vector(svdRes$d[1:k]), ncol=k, nrow=k) %*% t(as.matrix(svdRes$v[, 1:k]))
  colnames(a) <- colnames(Y)
  message('solve')
  W <- Y[, cIdx] %*% t(solve(a[, cIdx, drop = FALSE] %*% t(a[, cIdx, drop = FALSE]), a[, cIdx, drop = FALSE]))
  Wa <- W %*% a
  correctedY <- Y[1:m, ] - W[1:m, ] %*% a
  
  if(!isLog && all(.isWholeNumber(x))) {
    if(round) {
      correctedY <- round(exp(correctedY) - epsilon)
      correctedY[correctedY<0] <- 0
    } else {
      correctedY <- exp(correctedY) - epsilon
    }
  }
  
  W <- as.matrix(W[1:m,])
  colnames(W) <- paste("W", seq(1, ncol(W)), sep="_")
  return(list(W = W, normalizedCounts = t(correctedY)))
}  


RUVr <- function(x, cIdx, k, residuals, center=TRUE, round=TRUE, epsilon=1, tolerance=1e-8, isLog=FALSE) {
    
  
    require(Matrix)
    require(irlba)
  
    if(!isLog && !all(.isWholeNumber(x))) {
      warning(paste0("The expression matrix does not contain counts.\n",
                     "Please, pass a matrix of counts (not logged) or set isLog to TRUE to skip the log transformation"))
    }
    
    if(isLog) {
      Y <- t(x)
    } else {
      Y <- t(log(x+epsilon))
    }
    
    
    if(center) {
      E <- apply(residuals, 1, function(x) scale(x, center=TRUE, scale=FALSE))
    } else {
      E <- t(residuals)
    }
    m <- nrow(Y)
    n <- ncol(Y)
    svdWa <- irlba(E[, cIdx],nu = min(m, n), nv = min(m, n))
    
    k <- min(k, max(which(svdWa$d > tolerance)))
    W <- svdWa$u[, (1:k), drop = FALSE]
    alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
    correctedY <- Y - W %*% alpha
    if(!isLog && all(.isWholeNumber(x))) {
      if(round) {
        correctedY <- round(exp(correctedY) - epsilon)
        correctedY[correctedY<0] <- 0
      } else {
        correctedY <- exp(correctedY) - epsilon
      }
    }
    colnames(W) <- paste("W", seq(1, ncol(W)), sep="_")
    return(list(W = W, normalizedCounts = t(correctedY)))
}





plotDiffrentialExpression <- function(sobj, clus, period, genes=NULL){
  
  
  if(is.null(genes)){
    stop("Please provide a gene name")
  }
  
  if(! genes %in% rownames(sobj@data)){
    stop("Make sure that the provided gene name is exists in your data")
  }
  
  smp_saline <- rownames(subset(sobj@meta.data, Period == period & treatment == "Saline" & tree.ident== clus))
  smp_cocaine <- rownames(subset(sobj@meta.data, Period == period & treatment == "Cocaine" & tree.ident== clus))
  
  if(length(smp_saline) > 0 &length(smp_cocaine)>0){
    tmp_sobj <- SubsetData(sobj, cells.use=c(smp_saline, smp_cocaine))
    
    tmp_sobj <- SetIdent(tmp_sobj, ident.use  = tmp_sobj@meta.data$treatment)
    VlnPlot(tmp_sobj,genes)
  }
}



plotVolcanoPlot <- function(deg, pval, FC){
  
  
  deg$log2FC <- log2(exp(deg$avg_logFC))
  
  deg$isDE <- "NO"
  deg$gene = rownames(deg)
  
  pos <- which(abs(deg$avg_logFC) > log(FC) & deg$p_val_adj < pval)
  
  deg$isDE[pos] <- "Yes"
  
  deg$isDE <- factor(deg$isDE, levels = c("Yes","NO"))
  
  isDE <- subset(deg, isDE == "Yes")
  
  
  p1 = ggplot(deg, aes(x=log2FC, y= -log10(p_val_adj), color=isDE)) + geom_point() +
    theme_bw() + scale_color_manual(values = c("#925E9FFF","grey80")) + 
    geom_vline(xintercept = c(-1,1) * log2(FC), color="black",size=0.5)+
    geom_hline(yintercept = -log10(pval), color="black",size=0.5)
  
  if(nrow(isDE)>0){
    
    if(nrow(isDE)>10){
      ord <- order(abs(isDE$log2FC) * -log10(isDE$p_val_adj), decreasing = T)
      isDE <- head(isDE[ord,],10)
      
    }
    
   p1 = p1 + geom_text_repel(data=isDE, aes(x=log2FC, y= -log10(p_val_adj), label= gene))  
  }
  
  p1
  
}



plotAverageMclusResults <- function(mclus_res, clus_expr, period_lens){
  
  tmp <- split(rownames(mclus_res$data),mclus_res$classification)
  
  clus_means <- foreach( x =1:length(tmp),.combine = "rbind") %do% {
    
    if(length(tmp[[x]])> 1){
      df <- data.frame(group=x,
                       cell = 1:ncol(clus_expr),
                       mean = log2(colMeans(2^as.matrix(clus_expr[tmp[[x]],]))),
                       period = rep(c("Maintenance","withdraw_48h","withdraw_15d"), period_lens)
      )  
    }else{
      df <- data.frame(group=x,
                       cell = 1:ncol(clus_expr),
                       mean = as.numeric(clus_expr[tmp[[x]],]),
                       period = rep(c("Maintenance","withdraw_48h","withdraw_15d"),period_lens)
      ) 
    }
    
    return(df)
  }  
  
  clus_means$group <- factor(clus_means$group, levels = 1:length(tmp))
  clus_means$period <- factor(clus_means$period, levels = c("Maintenance","withdraw_48h","withdraw_15d"))
  
  
  #pdf(fname, width = 10,height = 20)
  p1 = ggplot(clus_means, aes(x=cell, y=mean, group=group, colour=group)) + geom_point() +
    theme_bw() +  facet_grid(group~., scales = "free_y") + 
    scale_color_npg() +
    #scale_color_manual(values = generateColors(20)) +
    geom_vline(xintercept = cumsum(period_lens)[-3],linetype="dashed") + ylab("log2(Cocain/saline)")  
  
  return(p1)
}


my_AddModuleScore <- function(
  object,
  genes.list = NULL,
  genes.pool = NULL,
  n.bin = 25,
  seed.use = 1,
  ctrl.size = 100,
  use.k = FALSE,
  enrich.name = "Cluster",
  random.seed = 1
) {
  set.seed(seed = random.seed) 
  genes.old <- genes.list
  if (use.k) {
    genes.list <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = genes.list)
  } else {
    if (is.null(x = genes.list)) {
      stop("Missing input gene list")
    }
    genes.list <- lapply(
      X = genes.list,
      FUN = function(x) {
        return(intersect(x = x, y = rownames(x = object@data)))
      }
    )
    cluster.length <- length(x = genes.list)
  }
  if (! all(Seurat:::LengthCheck(values = genes.list))) {
    warning(paste(
      'Could not find enough genes in the object from the following gene lists:',
      paste(names(x = which(x = ! Seurat:::LengthCheck(values = genes.list)))),
      'Attempting to match case...'
    ))
    genes.list <- lapply(
      X = genes.old,
      FUN = CaseMatch, match = rownames(x = object@data)
    )
  }
  if (! all(Seurat:::LengthCheck(values = genes.list))) {
    stop(paste(
      'The following gene lists do not have enough genes present in the object:',
      paste(names(x = which(x = ! Seurat:::LengthCheck(values = genes.list)))),
      'exiting...'
    ))
  }
  if (is.null(x = genes.pool)) {
    genes.pool = rownames(x = object@data)
  }
  data.avg <- apply(X = object@data[genes.pool,], MARGIN = 1, FUN = ExpMean)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = Hmisc::cut2(
    x = data.avg,
    m = round(x = length(x = data.avg) / n.bin)
  ))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector("list", cluster.length)
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(x = genes.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[genes.use[j]])],
          size = ctrl.size,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- c()
  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores <- rbind(
      ctrl.scores,
      apply(X = object@data[genes.use, ], MARGIN = 2, FUN = ExpMean)
    )
  }
  genes.scores <- c()
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    if (length(genes.use) == 1) {
      data.use <- t(as.matrix(object@data[genes.use, ]))
    } else {
      data.use <- object@data[genes.use, ]
    }
    genes.scores <- rbind(
      genes.scores,
      apply(X = data.use, MARGIN = 2, FUN = ExpMean)
    )
  }
  
  genes.scores.use <- genes.scores - ctrl.scores
  rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
  genes.scores.use <- t(x = as.data.frame(x = genes.scores.use))
  object <- AddMetaData(
    object = object,
    metadata = genes.scores.use,
    col.name = colnames(x = genes.scores.use)
  )
  return (object)
}



my_CellCycleScoring <- function (object, g2m.genes, s.genes, set.ident = FALSE) 
{
  enrich.name <- "Cell Cycle"
  genes.list <- list(S.Score = s.genes, G2M.Score = g2m.genes)
  object.cc <- my_AddModuleScore(object = object, genes.list = genes.list, n.bin = 50,
                              enrich.name = enrich.name, ctrl.size = min(vapply(X = genes.list, 
                                                                                FUN = length, FUN.VALUE = numeric(1))))
  cc.columns <- grep(pattern = enrich.name, x = colnames(x = object.cc@meta.data))
  cc.scores <- object.cc@meta.data[, cc.columns]
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores, 
                                                                 first = "S", second = "G2M", null = "G1") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      return(c(first, second)[which(x = scores == max(scores))])
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), 
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "S.Score", "G2M.Score", 
                               "Phase")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("S.Score", "G2M.Score", "Phase")]
  object <- AddMetaData(object = object, metadata = cc.scores)
  
  if (set.ident) {
    object <- StashIdent(object = object, save.name = "old.ident")
    object <- SetAllIdent(object = object, id = "Phase")
  }
  return(object)
}



removeCCGenes <- function(sobj, cc.genes = NULL){
  
  if(is.null(cc.genes)){
    stop("Please provide the list of Cell-cycle genes")
  }
  
  if(!is.character(cc.genes)){
    stop("cc.genes should be a charachter list")
  }
  
  sobj@raw.data <- sobj@raw.data[!rownames(sobj@raw.data) %in% cc.genes,]
  sobj@data <- sobj@data[!rownames(sobj@data) %in% cc.genes,]
  sobj
}


#From: https://github.com/ChristophH/in-lineage/blob/a00166e9798203e4bda03c4afcc03ad60c77c458/R/lib.R#L243

dim.red <- function(expr, max.dim, ev.red.th, plot.title=NA, do.scale.result=FALSE) {
  require(diffusionMap)
  cat('Dimensionality reduction via diffusion maps using', nrow(expr), 'genes and', ncol(expr), 'cells\n')
  if (sum(is.na(expr)) > 0) {
    dmat <- 1 - cor(expr, use = 'pairwise.complete.obs')
  } else {
    dmat <- 1 - cor(expr)
  }
  
  max.dim <- min(max.dim, nrow(dmat)/2)
  dmap <- diffuse(dmat, neigen=max.dim, maxdim=max.dim)
  ev <- dmap$eigenvals
  
  ev.red <- ev/sum(ev)
  evdim <- rev(which(ev.red > ev.red.th))[1]
  
  if (is.character(plot.title)) {
    plot(ev, ylim=c(0, max(ev)), main = plot.title)
    abline(v=evdim + 0.5, col='blue')
  }
  
  evdim <- max(2, evdim, na.rm=TRUE)
  cat('Using', evdim, 'significant DM coordinates\n')
  
  colnames(dmap$X) <- paste0('DMC', 1:ncol(dmap$X))
  res <- dmap$X[, 1:evdim]
  if (do.scale.result) {
    res <- scale(dmap$X[, 1:evdim])
  } 
  return(res)
}



plot_geneExpr_trajectory <- function(rmono, gene, do.log=T){
  
  
  if(!gene %in% rownames(exprs(rmono))){
    stop("Not all the genes are defined in rmono")
  }
  
  genes_expr <- as.matrix(exprs(rmono)[gene,])
  genes_expr <- as.data.frame(genes_expr)
  colnames(genes_expr) <- "gene_expr"
  
  if(do.log){
    genes_expr <- log1p(genes_expr)
  }
  
  pData(rmono) <- cbind(pData(rmono), genes_expr)
  
  p=plot_cell_trajectory(rmono, color_by = "gene_expr") + ggsci::scale_color_gsea() +
    ggtitle(gene) + theme(plot.title=element_text(hjust=0.5))
  
  p
}


getSmoothedPseudoTime_matrix <- function(cds_subset, 
                                         norm_method = c("log", "vstExprs"), 
                                         trend_formula = "~sm.ns(Pseudotime, df=3)", 
                                         scale_max=3, scale_min=-3,
                                         cores=1,pseudocount=1, scale_mat=TRUE){
 
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
                                         max(pData(cds_subset)$Pseudotime), 
                                         length.out = 100))
  
  m <- genSmoothCurves(cds_subset, cores = cores, 
                       trend_formula = trend_formula, 
                       relative_expr = T, 
                       new_data = newdata) 
  
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == 
      FALSE) {
    m = vstExprs(cds_subset, expr_matrix = m)
  }
  else if (norm_method == "log") {
    m = log2(m + pseudocount)
  }
  
  if(scale_mat){
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale_max] = scale_max
    m[m < scale_min] = scale_min  
  }
  m
}


# detectChange <- function(x,cpmType="Lepage"){
#   
#   cpm <- makeChangePointModel(cpmType=cpmType, ARL0=500)
#   
#   changepoints <- c()
#   
#   while (i < nrow(x)) {
#     i <- i + 1
#     #process each observation in turn
#     cpm <- processObservation(cpm,x[i,])
#     #if a change has been found, log it, and reset the CPM
#     if (changeDetected(cpm) == TRUE) {
#       print(sprintf("Change detected at observation %d", i))
#       #the change point estimate is the maximum D_kt statistic
#       Ds <- getStatistics(cpm)
#       tau <- which.max(Ds)
#       if (length(changepoints) > 0) {
#         tau <- tau + changepoints[length(changepoints)]
#       }
#       changepoints <- c(changepoints,tau)
#       #reset the CPM
#       cpm <- cpmReset(cpm)
#     }
#   }
#   
  
#  cpm <- processObservation(cpm,x[i])
#}



# From:
# https://github.com/mayer-lab/Mayer-et-al-2018_IntegratedAnalysis/blob/master/R/2_cell_assignment.R

# assign precursor cells to adult subtypes
MapCells=function(object, timepoint.dev="P10", timepoint.end="P56", num.k=10, thresh.require=0.9,map.col="AdultTypes1",new.col="Map1") {
  library(pbapply)
  tsne.dims=object@calc.params$RunTSNE$dims.use
  input.dims=GetCellEmbeddings(object,reduction.type = "cca.aligned",dims.use = tsne.dims)
  input.dist=as.matrix(dist(input.dims))
  cells.end=FastWhichCells(object,"DevStage",timepoint.end)
  cells.map=FastWhichCells(object,"DevStage",timepoint.dev)
  map.id=pbsapply(cells.map,function(x)
    names(which(sort(table(object@meta.data[names(sort(input.dist[x,cells.end])[1:num.k]),map.col]),decreasing = T)>=num.k*thresh.require))[1])
  nn.same=pbsapply(cells.map, function(x) length(intersect(names(sort(input.dist[x,])[1:num.k]),cells.map)))
  print(table(map.id))
  map.id[is.na(map.id)]="Unassigned"
  map.id[names(which(nn.same>=num.k*1))]="Unassigned"
  print(table(map.id))
  object@meta.data[cells.end,new.col]=object@meta.data[cells.end,map.col]
  object@meta.data[cells.map,new.col]=map.id
  return(object)
}


## Sub-sample sobj

SubSample_sobj <- function(sobj,nbCell=NULL,conserve_dist=TRUE){
  
  if(is.null(sobj) | class(sobj) !="seurat"){
    stop("No Seurat object is provided")
  }
  
  if(is.null(nbCell)){
    stop("You need specify the final subsample size")
  }
  
  if(length(sobj@cell.names) < nbCell){
    stop("nbCell is bigger than the number of available cells")
  }
  
  
  ## if we don't conserve the distribution just sample and return the object
  if(!conserve_dist){
    toUse <- sample(1:length(sobj@cell.names),size = nbCell)
    
    res_sobj <- SubsetData(sobj,cells.use = sobj@cell.names[toUse],subset.raw = T)
    
    return(res_sobj)
  }
  
  
  ## Otherwise get the distribution
  propotions <- table(sobj@ident)/length(sobj@ident)
  
  
  nbCells <- ceiling(propotions * nbCell)
  
  toUse <- c()

  
  for(clus in names(propotions)){
    
    cells_in_class <- which(sobj@ident == clus)
    
    toUse <- c(toUse, sample(cells_in_class,size = nbCells[clus]))
  }
  
  res_sobj <- SubsetData(sobj,cells.use = sobj@cell.names[toUse],subset.raw = T)
  
  return(res_sobj)
}



plotMultiClassROC <- function(sobj,RfClassifier){
  
  require(pROC)  
  
  ROCs_df <- data.frame()
  for(clus in levels(sobj@ident)){
    ids = as.character(sobj@ident)
    pos <- ids==clus
    ids[pos]=1
    ids[!pos]=0
    
    #predClass <- colnames(RfClassifier$predictions)[apply(RfClassifier$predictions,1,which.max)]
    # pos <- predClass==clus
    # predClass[pos]=1
    # predClass[!pos]=0
    
    
    predClass =RfClassifier$predictions[,clus]
    ROC1 <- smooth(roc(ids, as.numeric(predClass)))
    
    df = data.frame(specificity=ROC1$specificities, sensitivity=ROC1$sensitivities, cluster=clus)
    
    ROCs_df <- rbind(ROCs_df,df)
  }
  
  return(ROCs_df)
}


convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


AverageDetectionRateByCluster <- function (object, thresh.min = 0,clus=NULL) 
{
  if(is.null(clus)){
    stop("Please specify a cluster")
  }
  
  if(sum(object@ident == clus) ==0 ) {
    stop("Please check that specificied cluster name exist")
  }
  
  ident1 = clus
  ident2 = setdiff(unique(object@ident), ident1)
  
  tmp.cell1 = WhichCells(object = object, ident = ident1)
  tmp.cell2 = WhichCells(object = object, ident = ident2)
  
  data.temp1 <- apply(X = object@data[, tmp.cell1], MARGIN = 1, 
                     FUN = function(x) {
                       return(sum(x > thresh.min)/length(x = x))
                     })
  
  data.temp2 <- apply(X = object@data[, tmp.cell2], MARGIN = 1, 
                     FUN = function(x) {
                       return(sum(x > thresh.min)/length(x = x))
                     })
  
  
  data.all <- cbind(data.temp1, data.temp2)
  
  colnames(data.all) <- c("pct.1","pct.2")
  return(data.all)
}


getSpecificMarkers <- function(sobj,
                               groups,
                               gpval = 1e-3,
                               lpval = 1e-3){
  sobj_expr = expm1(sobj@data)
  
  AE_expression <- matrix(0,nrow = nrow(sobj_expr), ncol = length(levels(sobj@ident)))
  rownames(AE_expression)= rownames(sobj_expr)
  colnames(AE_expression) <- levels(sobj@ident)
  
  for(clus in levels(sobj@ident)){
    
    pos <- which(sobj@ident == clus)
    
    cls_mean  <- apply(sobj_expr[,pos],1,function(x){
      q99= as.numeric(quantile(x,probs = 0.99))
      x = x[x<q99]
      if(length(x)>0){
        mean(x)  
      }else{
        0
      }
    })
    
    AE_expression[,clus]= cls_mean
    print(paste0("Finished processing cluster: ",clus))
  }
  
  AE_detection <-  AverageDetectionRate(sobj)
  
  df = data.frame(cluster = levels(sobj@ident),
                  group = factor(groups))
  
  design = model.matrix( ~ -1 + group, data = df)
  rownames(design) = levels(sobj@ident)
  
  
  res =matrix(0,
              nrow = nrow(AE_detection),
              ncol = ncol(design))
  
  rownames(res) <- rownames(AE_detection)
  colnames(res) <- colnames(design)


  
  
  
  #exprs = m1^1.5 * m2
  exprs = AE_detection * AE_expression #log1p(m2)
  exprs = as.data.frame(exprs)#[toUse,]
  exprs$GENE = rownames(exprs)
  exprs = exprs[, c(ncol(exprs), 1:(ncol(exprs) - 1))]
  
  Results <- RN_calc(exprs, design)
  Results <- RN_select(Results, gpv_t = gpval, lpv_t = lpval)
  
  res = Results$selected
  res = res[, -c(1:3)]
  #tmp[tmp < 1] = 0
  #rs = rowSums(tmp)
  #tmp = tmp[rs > 0, ]
  res = data.matrix(res)
  res
}


getDiseaseAssociation_old <-
  function(disease_GWAS_genes,
           sobj,
           groups,
           gpval = 1e-3,
           lpval = 1e-3) {
    
  
    suppressMessages(require(RNentropy))
    #AE_expression <-  AverageExpression(sobj)
    
    sobj_expr = expm1(sobj@data)
    
    AE_expression <- matrix(0,nrow = nrow(sobj_expr), ncol = length(levels(sobj@ident)))
    rownames(AE_expression)= rownames(sobj_expr)
    colnames(AE_expression) <- levels(sobj@ident)
    
    for(clus in levels(sobj@ident)){
      
      pos <- which(sobj@ident == clus)
      
      cls_mean  <- apply(sobj_expr[,pos],1,function(x){
        q99= as.numeric(quantile(x,probs = 0.99))
        x = x[x<q99]
        if(length(x)>0){
          mean(x)  
        }else{
          0
        }
      })
      
      AE_expression[,clus]= cls_mean
      print(paste0("Finished processing cluster: ",clus))
    }
    
    AE_detection <-  AverageDetectionRate(sobj)
    
    df = data.frame(cluster = levels(sobj@ident),
                    group = factor(groups))
    
    design = model.matrix( ~ -1 + group, data = df)
    rownames(design) = levels(sobj@ident)
    disease_GWAS_genes_association = list()
    
    for (disease in names(disease_GWAS_genes)) {
      print(disease)
      disease_GWAS_genes_association[[disease]] <- matrix(0,
                                                          nrow = length(disease_GWAS_genes[[disease]]),
                                                          ncol = ncol(design))
      
      rownames(disease_GWAS_genes_association[[disease]]) <-
        disease_GWAS_genes[[disease]]
      colnames(disease_GWAS_genes_association[[disease]]) <-
        colnames(design)
      
      disease_GWAS_genes[[disease]] = intersect(disease_GWAS_genes[[disease]], rownames(AE_expression))
  
      
      m1 = AE_detection[disease_GWAS_genes[[disease]], ]
      #m1[m1<0.1]=0
      m2 = AE_expression[disease_GWAS_genes[[disease]], ]
      
      
    
      #exprs = m1^1.5 * m2
      exprs = m1 * m2 #log1p(m2)
      exprs = as.data.frame(exprs)#[toUse,]
      exprs$GENE = rownames(exprs)
      exprs = exprs[, c(ncol(exprs), 1:(ncol(exprs) - 1))]
      
      Results <- RN_calc(exprs, design)
      Results <- RN_select(Results, gpv_t = gpval, lpv_t = lpval)
      
      tmp = Results$selected
      tmp = tmp[, -c(1:3)]
      #tmp[tmp < 1] = 0
      #rs = rowSums(tmp)
      #tmp = tmp[rs > 0, ]
      tmp = data.matrix(tmp)
      
      disease_GWAS_genes_association[[disease]][rownames(tmp), colnames(tmp)] = disease_GWAS_genes_association[[disease]][rownames(tmp), colnames(tmp)] + tmp
    }
    
    return(disease_GWAS_genes_association)
  }



getDiseaseAssociation <-
  function(disease_GWAS_genes,
           sobj,
           groups,
           gpval = 1e-3,
           perpval = 1e-10,
           lpval = 1e-3,
           all_mrks= NULL, 
           AE_expression = NULL,
           AE_detection = NULL) {
    
    
    suppressMessages(require(RNentropy))
    #AE_expression <-  AverageExpression(sobj)
    
    sobj_expr = expm1(sobj[['RNA']]@data)
    
    if(is.null(AE_expression)){
      AE_expression <- matrix(0,nrow = nrow(sobj_expr), ncol = length(levels(sobj@active.ident)))
      rownames(AE_expression)= rownames(sobj_expr)
      colnames(AE_expression) <- levels(sobj@active.ident)
      
      for(clus in levels(sobj@active.ident)){
        
        pos <- which(sobj@active.ident == clus)
        
        cls_mean  <- apply(sobj_expr[,pos],1,function(x){
          q99= as.numeric(quantile(x,probs = 0.99))
          x = x[x<q99]
          if(length(x)>0){
            mean(x)  
          }else{
            0
          }
        })
        
        AE_expression[,clus]= cls_mean
        print(paste0("Finished processing cluster: ",clus))
      }  
    }
    
    
    if(is.null(AE_detection)){
      AE_detection <-  AverageDetectionRate(sobj)  
    }
    
  
    df = data.frame(cluster = levels(sobj@active.ident),
                    group = factor(groups))
    
    design = model.matrix( ~ -1 + group, data = df)
    rownames(design) = levels(sobj@active.ident)
    disease_GWAS_genes_association = list()
    
    
    #all_mrks$cluster  = factor(all_mrks$cluster, levels = colnames(m2))
  
    
    for (disease in names(disease_GWAS_genes)) {
      print(disease)
      
      disease_GWAS_genes[[disease]] = intersect(disease_GWAS_genes[[disease]], rownames(AE_expression))
      
      disease_GWAS_genes_association[[disease]] <- matrix(0,
                                                          nrow = length(disease_GWAS_genes[[disease]]),
                                                          ncol = ncol(design))
      
      rownames(disease_GWAS_genes_association[[disease]]) <-
        disease_GWAS_genes[[disease]]
      colnames(disease_GWAS_genes_association[[disease]]) <-
        colnames(AE_detection)
      
      
      
      
      m1 = AE_detection[disease_GWAS_genes[[disease]], ]
      #m1[m1<0.1]=0
      m2 = AE_expression[disease_GWAS_genes[[disease]], ]
      #m2[m2<1]=0
      
      
      
      #exprs = m1^1.5 * m2
      exprs = m2 #m1 * m2 
      exprs = as.data.frame(exprs)#[toUse,]
      exprs$GENE = rownames(exprs)
      exprs = exprs[, c(ncol(exprs), 1:(ncol(exprs) - 1))]
      
      Results <- RN_calc(exprs, design)
      
      
      percent_express = as.data.frame(100*m1)
      percent_express$GENE = rownames(m1)
      percent_express = percent_express[, c(ncol(percent_express), 1:(ncol(percent_express) - 1))]
      Results2 <- RN_calc(percent_express, design)
      #Results <- RN_select(Results, gpv_t = gpval, lpv_t = lpval)
      
      
      
      # mrks = subset(all_mrks, gene %in% disease_GWAS_genes[[disease]] & p_val_adj<lpval)
      # mrks$p_val_adj = 1
      # mrks$cluster  = factor(mrks$cluster, levels = colnames(m2))
      # 
      #loc_pval = reshape2::dcast(data = mrks, formula = gene~cluster,value.var = "p_val_adj",fill=0,drop = F)
      # 
      #rownames(loc_pval) = loc_pval$gene 
      #loc_pval$gene = NULL
      #loc_pval = data.matrix(loc_pval)
      # 
      #loc_pval[loc_pval<lpval]=1
      
      loc_pval = m1 * 0
      
      selected_genes = names(Results$gpv[Results$gpv > -log10(gpval)])
      selected_genes2 = names(Results2$gpv[Results2$gpv > -log10(perpval)])
      #tmp = Results$selected
      #tmp = tmp[, -c(1:3)]
      #tmp[tmp < 1] = 0
      #rs = rowSums(tmp)
      #tmp = tmp[rs > 0, ]
      
      selected_genes = intersect(selected_genes, rownames(loc_pval))
      selected_genes = intersect(selected_genes, selected_genes2)
      
      m1[m1>0.05]=1
      m1[m1<=0.05]=0
      
      
      # m2 = t(apply(m2,1,function(x){
      #   if(sum(x)>0){
      #     log2(x+1)/log2(max(x+1))  
      #   }else{
      #     x
      #   }
      # }))
      
      for(i in 1:nrow(m2)){
        gene = rownames(m2)[i]
        #print(paste0(gene,"\r"))
        
        if(max(as.numeric(m2[gene,])) ==0) { next}
        
        fc =log2(as.numeric(m2[gene,])+1)/log2(max(as.numeric(m2[gene,]))+1)
        
        if(min(fc)>0.5){
          m2[gene,]=1
        }else{
          detected <- SamSPECTRAL::kneepointDetection(vect=sort(as.numeric(m2[gene,])), PlotFlag=FALSE)
          thr = sort(as.numeric(m2[gene,]))[detected$MinIndex]
          pos <- which(as.numeric(m2[gene,]) < thr)
          m2[gene,pos] = 0 
          pos <- which(as.numeric(m2[gene,]) >= thr)
          m2[gene,pos] = 1
        }
      }
      
      #m2[m2 > 0.5]=1
      #m2[m2 <= 0.5]=0
      #tmp = m1[m1>0.1] * m2[m2>1]
      tmp = m1 * m2
      
      tmp = tmp[selected_genes,]
      
      rs = rowSums(tmp)
      tmp = tmp[rs>0 ,]
      # tmp = loc_pval[selected_genes,]
      # tmp = data.matrix(tmp)
      # tmp2 = m2[selected_genes, ]
      # tmp2[tmp2>1]=1
      # tmp2[tmp2<1]=0
      # tmp = tmp + tmp2
      
      # tmp = exprs[selected_genes,]
      # tmp = tmp[,-1]
      # tmp[tmp>0] = 1
      # 
      disease_GWAS_genes_association[[disease]] = tmp
    }
    
    return(disease_GWAS_genes_association)
}


KNNSmooth <- function(sobj, k=10, reduction.type="cca.aligned",dims.use=1:10){
  
  
  data.use <- GetCellEmbeddings(object = sobj,
                                reduction.type = reduction.type,
                                dims.use = dims.use)
  
  exprs  = as.matrix(expm1(sobj@data))
  smoothed = exprs * 0
  
  knn = RANN::nn2(data.use,k = k+1,searchtype = 'standard')
  
  
  smoothed = apply(knn$nn.idx,1,function(x){
    neighbors <- as.numeric(x)[-1]
    tmp =as.matrix(exprs[,neighbors])
    #rowMeans some how is very slow
    m = Rfast::rowmeans(tmp)
    m
  })
  
  colnames(smoothed) = colnames(exprs)
  rownames(smoothed) = rownames(exprs)
  
  smoothed = as(smoothed, "dgCMatrix")
  smoothed = log1p(smoothed)
  sobj = SetAssayData(sobj,  assay.type ="smoothed", "raw.data", smoothed)
  sobj
}



getExpressedClasses <- function(exprs){
  
  genes_class <- data.frame()
  for(gene in rownames(exprs)){
    
    # get the clusters in which the gene is expressed
    expressedin = colnames(exprs)[which(exprs[gene,]==1)]
    
    
    expressedin = unique(gsub("_.+","",expressedin))
    
    if(all( c("Exc","Inhib","Endo","Astro","Oligo","NF Oligo","OPC","Microglia") %in% expressedin)){
      df <- data.frame(gene= gene, class= "All")	
    }else{
      df <- data.frame(gene= gene, class= expressedin)
    } 	
    genes_class =rbind(genes_class, df)
  } 
  return(genes_class)
}


getExpressedClasses2 <- function(exprs){
  
  
  
  
  gaps = gsub("_\\d+","",colnames(exprs))
  gaps = factor(gaps)
  
  
  
  
  percents = tapply(colnames(exprs),gaps,
                    function(x){
                      rs = rowSums(exprs[,x])                          
                      rs/length(x)
                    })
  
  percents = do.call("rbind",percents)
  percents = t(percents)
  
  mx = apply(percents,1,max)
  cluster_specifc_genes = which(mx <0.5)
  
  mx = apply(percents,1,min)
  non_specifc_genes = which(mx >= 0.5)
  
  
  anno = percents
  anno[percents>=0.5]=1
  anno[percents<0.5]=0
  
  rs = rowSums(anno)
  
  # cell-type specific
  specific_genes = which(rs==1)
  
  
  #Consider the genes expressed in neuronal clusters also specific
  neuronal = which(rs==2 & percents[,"D1"]==1 & percents[,"D2"]==1)
  
  
  genes_class <- data.frame(gene =rownames(percents),
                            category=c("non-specific"),
                            stringsAsFactors = F
  )
  
  genes_class$category[cluster_specifc_genes] = "cluster-Specific"
  genes_class$category[non_specifc_genes] = "non-specific"
  genes_class$category[specific_genes] = "cellType-specific"
  genes_class$category[neuronal] = "cellType-specific"
  
  genes_class
}





my_RN_select = function (Results, gpv_t = 0.01, lpv_t = 0.01, method = "BH") 
{
  lpv_t <- -log10(lpv_t)
  gpv_t <- -log10(gpv_t)
  Results$gpv_bh <- -log10(p.adjust(10^-Results$gpv, method = method))
  true_rows <- (Results$gpv_bh >= gpv_t)
  design_b <- t(Results$design > 0)
  Results$lpv_sel <- data.frame(row.names = rownames(Results$lpv)[true_rows])
  boh <- list()
  for (d in seq_along(design_b[, 1])) {
    Results$lpv_sel <- cbind(Results$lpv_sel, apply(as.matrix(Results$lpv[true_rows, 
                                                                          design_b[d, , drop = F]]), 1,my_RN_select_lpv_row, 
                                                    lpv_t))
  }
  colnames(Results$lpv_sel) <- colnames(Results$design)
  lbl <- Results$res[, !sapply(Results$res, is.numeric), drop = FALSE]
  Results$selected <- cbind(lbl[true_rows, ], Results$gpv[true_rows], 
                            Results$gpv_bh[true_rows], Results$lpv_sel)
  colnames(Results$selected) <- c(names(which(!sapply(Results$res, 
                                                      is.numeric))), "GL_LPV", "Corr. GL_LPV", colnames(Results$lpv_sel))
  Results$selected <- Results$selected[order(Results$selected[, 
                                                              3], decreasing = TRUE), ]
  Results$lpv_sel <- NULL
  return(Results)
}


my_RN_select_lpv_row = function (x, lpv_t) 
{
  if (sum(abs(x) >= lpv_t)/length(x) >=0.6) {
    if (sum(x > 0)/length(x) >=0.6) {
      return(1)
    }
    else if (all(x < 0)) {
      return(-1)
    }
    else {
      return(NA)
    }
  }
  else {
    return(0)
  }
}



#Use a similar strategy as Hongkui zen's paper
FindSimilarClusters <- function(sobj1, sobj2, q.cutoff = 0.05){
  
  require(dplyr)
  
  # Train a random classifier to assign-cell type
  sobj1 = FindVariableGenes(sobj1)
  sobj2 = FindVariableGenes(sobj2)
  
  common_genes = intersect(sobj1@var.genes, sobj2@var.genes)
  
  sobj1_classifier  = BuildRFClassifier(sobj1,
                                        training.genes = common_genes,
                                        training.classes =sobj1@ident)
  
  preIdent <- predict(sobj1_classifier, as.matrix(t(sobj2@data[common_genes,])))
  sobj2@meta.data$Predicted_ident = preIdent
  
  
  # Get Average expression of the reference dataset (don't consider outliers)
  sobj1_avgExpr <- RobustAverageExpression(sobj1)	
  
  common_genes = intersect(rownames(sobj1@data), rownames(sobj2@data))
  
  # Calculate correlation between each cell and the best cluster
  sobjs_cells_cor <- matrix(0,nrow=ncol(sobj2@data), ncol=1)
  rownames(sobjs_cells_cor) = colnames(sobj2@data)
  colnames(sobjs_cells_cor) = "Correlations"
  
  sobjs_cells_cor = as.data.frame(sobjs_cells_cor)
  
  sobjs_cells_cor$Correlations <- sapply(colnames(sobj2@data), function(cellname){
    
    predId = obj2@meta.data[cellname,"Predicted_ident"]
    cr = cor(as.numeric(sobj2@data[common_genes,cellname]), as.numeric(sobj1_avgExpr[common_genes,predId]))
    cr		
  })
  
  # Calculate correlations z-score
  sobjs_cells_cor$zscore = scale(sobjs_cells_cor$Correlations, scale=T, center =T)
  
  sobjs_cells_cor$Cluster = sobj2@meta.data$Predicted_ident
  
  # Find clusters that don't show a high mean correlation
  
  lower_cutoff = qnorm(q.cutoff)
  
  
  # Calculate the mean z-score for each cluster
  
  cluster_mean = sobjs_cells_cor %>% 
    group_by(Cluster) %>%
    dplyr::summarize(Mean_zscore = mean(zscore, na.rm=TRUE))
  
  
  cluster_mean$Pass = "Yes"
  cluster_mean$Pass[cluster_mean$Mean_zscore < lower_cutoff] = "No"
  
  return(cluster_mean)
}



getDiseaseAssociationByGroup <-
  function(disease_GWAS_genes,
           sobj,
           groups,
           gpval = 1e-3,
           lpval = 1e-3) {
    
    
    
    AE_expression <-  AverageExpression(sobj)
    AE_detection <-  AverageDetectionRate(sobj)
    
    df = data.frame(cluster = levels(sobj@ident),
                    group = factor(groups))
    
    design = model.matrix( ~ -1 + group, data = df)
    rownames(design) = levels(sobj@ident)
    disease_GWAS_genes_association = list()
    
    for (disease in names(disease_GWAS_genes)) {
      print(disease)
      disease_GWAS_genes_association[[disease]] <- matrix(0,
                                                          nrow = length(disease_GWAS_genes[[disease]]),
                                                          ncol = ncol(design))
      
      rownames(disease_GWAS_genes_association[[disease]]) <-
        disease_GWAS_genes[[disease]]
      colnames(disease_GWAS_genes_association[[disease]]) <-
        colnames(design)
      
      disease_GWAS_genes[[disease]] = intersect(disease_GWAS_genes[[disease]], rownames(AE_expression))
      
      
      m1 = AE_detection[disease_GWAS_genes[[disease]], ]
      #m1[m1<0.1]=0
      m2 = AE_expression[disease_GWAS_genes[[disease]], ]
      
      quan = apply(m1, 1, median)
      
      #toUse = which(quan<0.5)
      exprs = m1 * m2
      exprs = as.data.frame(exprs)#[toUse,]
      exprs$GENE = rownames(exprs)
      exprs = exprs[, c(ncol(exprs), 1:(ncol(exprs) - 1))]
      
      Results <- RN_calc(exprs, design)
      Results <- my_RN_select(Results, gpv_t = gpval, lpv_t = lpval)
      
      tmp = Results$selected
      tmp = tmp[, -c(1:3)]
      #tmp[tmp < 1] = 0
      #rs = rowSums(tmp)
      #tmp = tmp[rs > 0, ]
      tmp = data.matrix(tmp)
      
      disease_GWAS_genes_association[[disease]][rownames(tmp), colnames(tmp)] = disease_GWAS_genes_association[[disease]][rownames(tmp), colnames(tmp)] + tmp
    }
    
    return(disease_GWAS_genes_association)
  }



RobustAverageExpression <- function(sobj, assay ="RNA",slot="data"){
  
  
  if(assay != "RNA"){
    sobj_expr <- GetAssayData(sobj, assay.type = assay, slot = slot)
    sobj_expr = expm1(sobj_expr)
  }else{
    sobj_expr = expm1(sobj[['RNA']]@data)    
  }
  
  AE_expression <- matrix(0,nrow = nrow(sobj_expr), ncol = length(levels(sobj@active.ident)))
  rownames(AE_expression)= rownames(sobj_expr)
  colnames(AE_expression) <- levels(sobj@active.ident)
  
  for(clus in levels(sobj@active.ident)){
    
    pos <- which(sobj@active.ident == clus)
    
    cls_mean  <- apply(sobj_expr[,pos],1,function(x){
      q99= as.numeric(quantile(x,probs = 0.99))
      x = x[x<=q99]
      if(length(x)>0){
        mean(x)  
      }else{
        0
      }
    })
    
    AE_expression[,clus]= cls_mean
    print(paste0("Finished processing cluster: ",clus))
  }
  
  #AE_expression = log2(AE_expression+1)
  AE_expression
}


RobustAverageDetection <- function(sobj, assay ="RNA",slot="data"){
  
  
  if(assay != "RNA"){
    sobj_expr <- GetAssayData(sobj, assay.type = assay, slot = slot)
    sobj_expr = expm1(sobj_expr)
  }else{
    sobj_expr = expm1(sobj[['RNA']]@data)  
  }
  
  AE_detection <- matrix(0,nrow = nrow(sobj_expr), ncol = length(levels(sobj@active.ident)))
  rownames(AE_detection)= rownames(sobj_expr)
  colnames(AE_detection) <- levels(sobj@active.ident)
  
  for(clus in levels(sobj@active.ident)){
    
    pos <- which(sobj@active.ident == clus)
    
    cls_percent  <- apply(sobj_expr[,pos],1,function(x){
      
      q99= as.numeric(quantile(x,probs = 0.99))
      x = x[x<q99]
      if(length(x)>0){
        sum(x>=1)/length(x)        
      }else{
        0
      }
      
    })
    
    AE_detection[,clus]= cls_percent
    print(paste0("Finished processing cluster: ",clus))
  }
  
  AE_detection
}



#https://www.r-bloggers.com/bigcor-large-correlation-matrices-in-r/
# Not very optimized but it should be much faster.

bigcor <- function(x, nblocks = 10, verbose = TRUE, method ="s", nb.cores=4)
{
  library(bigmemory, quietly = TRUE)
  library(bigstatsr, quietly = TRUE)
  
  
  NCOL <- ncol(x)
  
  ## test if ncol(x) %% nblocks gives remainder 0
  #if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
  
  
  blockIds  = cut(1:NCOL,nblocks)
  blockIds=as.numeric(blockIds)
  
  MAT = as_FBM(x)
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  corMAT <- FBM(NCOL, NCOL)
  
  ## split column numbers into 'nblocks' groups
  #SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(sort(unique(blockIds)), sort(unique(blockIds)))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  
  fct <- function(cMat,ind, mat, combs,splits, method){                  
    COMB <- combs[ind, ]        
    G1 <- which(splits == COMB[1])
    G2 <- which(splits == COMB[2])
    
    cMat[G1,G2] <- cor(mat[, G1], mat[, G2], method = method)    
    cMat[G2, G1] <- t(cMat[G1,G2])
    NULL    
  }
  
  big_apply(corMAT,a.FUN=fct,a.combine = 'c', 
            ind= rows_along(COMBS),
            mat = MAT,
            combs=COMBS, 
            splits=blockIds,        
            method=method,
            block.size=1, # We process a combination at a time (don't change this value)
            ncores=nb.cores)
  
  gc()
  
  res = corMAT[1:nrow(corMAT),1:ncol(corMAT)]
  colnames(res) = colnames(x)
  rownames(res) = colnames(x)
  return(res)
}



MetaNeighborUS_parallel <- function(var_genes, dat, i = 1, study_id, cell_type){
  library(bigmemory)
  library(bigstatsr)
  
  #dat    <- SummarizedExperiment::assay(dat, i = i)
  
  if(class(dat) !="matrix"){
    stop("dat should be a matrix")
  }
  
  samples <- colnames(dat)
  
  #check obj contains study_id
  if(length(study_id)!=length(samples)){
    stop('study_id length does not match number of samples')
  }
  
  #check obj contains cell_type
  if(length(cell_type)!=length(samples)){
    stop('cell_type length does not match number of samples')
  } 
  
  pheno <- as.data.frame(cbind(study_id,cell_type), stringsAsFactors = FALSE)
  pheno$StudyID_CT <- paste(pheno$study_id, pheno$cell_type, sep = "|")
  celltypes   <- unique(pheno$StudyID_CT)
  cell_labels <- matrix(0, ncol=length(celltypes), nrow=dim(pheno)[1])
  rownames(cell_labels) <-colnames(dat)
  colnames(cell_labels) <- celltypes
  
  for(i in seq_along(celltypes)){
    type <- celltypes[i]
    matching_celltype <- match(pheno$StudyID_CT, type)
    cell_labels[!is.na(matching_celltype),i]  <- 1
  }
  
  matching_vargenes <- match(rownames(dat), var_genes)
  matching_vargenes_count   <- sum(!is.na(matching_vargenes))
  
  if(matching_vargenes_count < 2){
    stop("matching_vargenes should have more than 1 matching genes!",
         call. = TRUE)
  } else if(matching_vargenes_count < 5) {
    warning("matching_vargenes should have more matching genes!", 
            immediate. = TRUE)
  }
  
  print("calculating correlation")
  cor_data    <- bigcor(dat[!is.na(matching_vargenes),],nblocks=30, method="s")
  rm(dat)
  gc()
  rank_data   <- cor_data*0
  rank_data[] <- rank(cor_data, ties.method = "average", na.last = "keep")
  rank_data[is.na(rank_data)] <- 0
  rank_data   <- rank_data/max(rank_data)
  sum_in      <- (rank_data) %*% cell_labels
  sum_all     <- matrix(apply(rank_data, MARGIN = 2, FUN = sum), 
                        ncol = dim(sum_in)[2], 
                        nrow = dim(sum_in)[1])
  predicts    <- sum_in/sum_all
  
  cell_NV     <- matrix(0, ncol=length(celltypes), nrow=length(celltypes))
  colnames(cell_NV) <- colnames(cell_labels)
  rownames(cell_NV) <- colnames(cell_labels)
  
  for(i in seq_len(dim(cell_labels)[2])){
    predicts_temp <- predicts
    
    matching_celltype <- match(pheno$StudyID_CT, colnames(cell_labels)[i])
    unique_studyID    <- unique(pheno[!is.na(matching_celltype),"study_id"])
    matching_studyID  <- match(pheno$study_id, unique_studyID)
    pheno2            <- pheno[!is.na(matching_studyID),]
    predicts_temp     <- predicts_temp[!is.na(matching_studyID),]
    predicts_temp     <- apply(abs(predicts_temp), 
                               MARGIN = 2, 
                               FUN = rank, 
                               na.last= "keep", 
                               ties.method="average")
    
    
    filter  <- matrix(0, ncol=length(celltypes), nrow=dim(pheno2)[1])
    matches <- match(pheno2$StudyID_CT, colnames(cell_labels)[i])
    filter[!is.na(matches),seq_along(celltypes)] <- 1
    
    negatives = which(filter == 0, arr.ind = TRUE)
    positives = which(filter == 1, arr.ind = TRUE)
    
    predicts_temp[negatives] <- 0
    
    np <- colSums(filter, na.rm = TRUE)
    nn <- apply(filter, MARGIN = 2, FUN = function(x) sum(x==0,na.rm=TRUE))
    p  <- apply(predicts_temp, MARGIN = 2, FUN = sum, na.rm = TRUE)
    
    cell_NV[i,]= (p/np - (np+1)/2)/nn
  }
  
  cell_NV <- (cell_NV+t(cell_NV))/2
  return(cell_NV)
}



FindAllMarkers_parallel <- function(
  object,
  genes.use = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  return.thresh = 1e-2,
  do.print = FALSE,
  random.seed = 1,
  min.cells.gene = 3,
  min.cells.group = 3,
  latent.vars = NULL,
  assay.type = "RNA",
  idents.use=NULL,
  ...
) {
  
  require(doParallel)
  require(foreach)
  
  data.1 <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- Seurat:::SetIfNull(x = genes.use, default = rownames(x = data.1))
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh = 0.7
  }
  idents.all <- sort(x = unique(x = object@ident))
  
  if(!is.null(idents.use)){
    idents.all = intersect(idents.all,idents.use)
  }
  
  genes.de <- list()
  #if (max.cells.per.ident < Inf) {
  #  object <- SubsetData(
  #    object = object,
  #    max.cells.per.ident = max.cells.per.ident,
  #    random.seed = random.seed
  #  )
  #}
  
  registerDoParallel(4)
  
  genes.de = foreach(i =1:length(x = idents.all),.packages="Seurat") %dopar% {
    tryCatch(
      {
        FindMarkers(
          object = object,
          assay.type = assay.type,
          ident.1 = idents.all[i],
          ident.2 = NULL,
          genes.use = genes.use,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          print.bar = print.bar,
          min.cells.gene = min.cells.gene,
          min.cells.group = min.cells.group,
          latent.vars = latent.vars,
          max.cells.per.ident = max.cells.per.ident,
          ...
        )
      },
      error = function(cond){
        return(NULL)
      }
    )	
  }  
  
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      # if (test.use == "roc") {
      #   gde <- subset(
      #     x = gde,
      #     subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
      #   )
      # } else {
      #   gde <- gde[order(gde$p_val, -gde$avg_logFC), ]
      #   gde <- subset(x = gde, subset = p_val < return.thresh)
      # }
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if ((only.pos) && nrow(gde.all) > 0) {
    return(subset(x = gde.all, subset = avg_logFC > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(gde.all) == 0) {
    warning("No DE genes identified.")
  }
  return(gde.all)
}




#Use a similar strategy as Hongkui zen's paper
FindSimilarClusters <- function(sobj1, sobj2, q.cutoff = 0.05, training_genes =NULL){
  
  require(dplyr)
  
  if(is.null(training_genes)){
    sobj1 = FindVariableGenes(sobj1)
    sobj2 = FindVariableGenes(sobj2)
    
    training_genes = intersect(sobj1@var.genes, sobj2@var.genes)
    warning("The common variable genes are used for training")
  }
  
  
  # Train a random classifier to assign-cell type
  
  
  
  
  sobj1_classifier  = BuildRFClassifier(sobj1,
                                        training.genes = training_genes,
                                        training.classes =sobj1@ident)
  
  preIdent <- predict(sobj1_classifier, as.matrix(t(sobj2@data[training_genes,])))
  sobj2@meta.data$Predicted_ident = preIdent$predictions
  
  
  # Get Average expression of the reference dataset (don't consider outliers)
  sobj1_avgExpr <- RobustAverageExpression(sobj1)	
  
  common_genes = intersect(rownames(sobj1@data), rownames(sobj2@data))
  
  # Calculate correlation between each cell and the best cluster
  sobjs_cells_cor <- matrix(0,nrow=ncol(sobj2@data), ncol=1)
  rownames(sobjs_cells_cor) = colnames(sobj2@data)
  colnames(sobjs_cells_cor) = "Correlations"
  
  sobjs_cells_cor = as.data.frame(sobjs_cells_cor)
  
  sobjs_cells_cor$Correlations <- sapply(colnames(sobj2@data), function(cellname){
    
    predId = obj2@meta.data[cellname,"Predicted_ident"]
    cr = cor(as.numeric(sobj2@data[common_genes,cellname]), as.numeric(sobj1_avgExpr[common_genes,predId]))
    cr		
  })
  
  # Calculate correlations z-score
  sobjs_cells_cor$zscore = scale(sobjs_cells_cor$Correlations, scale=T, center =T)
  
  sobjs_cells_cor$Cluster = sobj2@meta.data$Predicted_ident
  
  # Find clusters that don't show a high mean correlation
  
  lower_cutoff = qnorm(q.cutoff)
  
  
  # Calculate the mean z-score for each cluster
  
  cluster_mean = sobjs_cells_cor %>% 
    group_by(Cluster) %>%
    dplyr::summarize(Mean_zscore = mean(zscore, na.rm=TRUE))
  
  
  cluster_mean$Pass = "Yes"
  cluster_mean$Pass[cluster_mean$Mean_zscore < lower_cutoff] = "No"
  
  return(cluster_mean)
}



# Adaptation from scrattch.hicat

#' Map a dataset to a reference, and compare existing cluster calls to the reference comparison
#' 
#' @param ref.dat Training data matrix, usually log-transformed CPM
#' @param ref.cl Training cluster factor object
#' @param map.dat Data for cells to map to the training set. Should have the same genes as train.dat.
#' @param map.cl Cluster assignments for the training set to compare to results of mapping.
#' 
#' @return a list object with two objects:  
#' \itemize{
#' \item map.df: A data.frame with the mapping results for each sample in map.dat to the reference
#' \item cl.map.df: A data.frame with cluster-level frequency of mapping for each cluster in map.cl to ref.cl
#' }
map_cl_summary.my <- function(ref.dat, 
                           ref.cl, 
                           map.dat, 
                           map.cl)
{
  # Map the training set to the reference
  map.result <- map_by_cor.my(ref.dat, ref.cl, map.dat)
  cor.matrix <- map.result$cor.matrix
  
  
  map.df <- map.result$pred.df
  colnames(map.df)[1] <- "map.cl"
  map.df$org.cl <- map.cl[row.names(map.df)]
  
  # Compute the fraction of times each sample was mapped to each cluster
  cl.size <- table(map.cl)
  cl.map.df <- as.data.frame(with(map.df, table(org.cl, map.cl)))
  cl.map.df$Prob <- round(cl.map.df$Freq / cl.size[as.character(cl.map.df$org.cl)], digits = 2)
  
  # Compute the mean fraction of mapping for all samples in the training set clusters
  # to the training set clusters.
  cl.map.df$pred.score <- 0
  for(i in 1:nrow(cl.map.df)){
    select <- names(map.cl)[map.cl == as.character(cl.map.df$org.cl[i])]
    cl.map.df$pred.score[i] <- mean(cor.matrix[select, as.character(cl.map.df$map.cl[i])])
  }
  cl.map.df$pred.score <- round(cl.map.df$pred.score, digits = 2)
  
  # Remove comparisons with no mapping
  cl.map.df <- cl.map.df[cl.map.df$Freq > 0, ]
  
  # Return output
  out_list <- list(map.df = map.df,
                   cl.map.df = cl.map.df)
  
  return(out_list)
}



map_by_cor.my <- function(train.dat, 
                       train.cl, 
                       test.dat,
                       method = "mean")
{
  # Get medians or means for each cluster
  if(method == "median"){
    library(matrixStats)
    cl.meds <- tapply(names(train.cl), 
                      train.cl, 
                      function(x) {
                        train.mat <- train.dat[, x, drop = F]
                        train.mat <- as.matrix(train.mat)
                        rowMedians(train.mat)
                      }
    )
    
    cl.dat <- do.call("cbind", cl.meds)
  } else {
    cl.dat <- get_cl_means(train.dat, train.cl)
  }
  row.names(cl.dat) <- row.names(train.dat)
  
  # Perform correlations
  test.cl.cor <- cor(as.matrix(test.dat), cl.dat, method = "spearman")
  test.cl.cor[is.na(test.cl.cor)] <- 0
  
  # Find maximum correlation
  max.cl.cor <- apply(test.cl.cor, 1, which.max)
  pred.cl <- colnames(test.cl.cor)[max.cl.cor]
  pred.cl <- setNames(pred.cl, row.names(test.cl.cor))
  
  # Get maximum correlation values
  pred.score <- apply(test.cl.cor, 1, max)
  
  # Convert to factor if train.cl was a factor and match levels.
  if(is.factor(train.cl)){
    pred.cl <- setNames(factor(pred.cl, levels = levels(train.cl)), names(pred.cl))
  }
  
  # Output results
  pred.df <- data.frame(pred.cl = pred.cl,
                        pred.score = pred.score)
  
  out_list <- list(pred.df = pred.df,
                   cor.matrix = test.cl.cor)
  
  return(out_list)    
}


map_sampling.my <- function(train.dat, 
                         train.cl, 
                         test.dat, 
                         markers, 
                         markers.perc = 0.8, 
                         iter = 100, 
                         method = "mean",
                         verbose = TRUE)
{
  # Perform mapping iter times.
  map.result <- sapply(1:iter, 
                       function(i){
                         if(verbose) {
                           cat("\r", paste0("Running iteration ",i," of ",iter,".        "))
                           flush.console()
                         }
                         tmp.markers <- sample(markers, round(length(markers) * markers.perc))
                         map_by_cor.my(train.dat[tmp.markers,], 
                                    train.cl, 
                                    test.dat[tmp.markers,], 
                                    method = method)
                       }, simplify = F)
  
  # Extract predicted cluster assignments from each iteration
  map.cl <- sapply(map.result, 
                   function(x) {
                     x$pred.df$pred.cl
                   }
  )
  # Compute fraction of times each sample mapped to each cluster
  row.names(map.cl) <- colnames(test.dat)
  map <- as.data.frame(as.table(as.matrix(map.cl)))
  map.freq <- table(map$Var1, map$Freq)
  
  # Find the most frequently mapped cluster for each sample
  max.freq <- apply(map.freq, 1, which.max)
  pred.cl <- colnames(map.freq)[max.freq]
  pred.cl <- setNames(pred.cl, row.names(map.freq))
  
  # Gather results
  map.df <- data.frame(pred.cl = pred.cl, 
                       prob = rowMaxs(map.freq) / iter)
  
  # output results
  out_list <- list(map.df = map.df,
                   map.freq = map.freq)
  
  return(out_list)
}

map_cv <- function(norm.dat, 
                   cl, 
                   markers, 
                   n.bin = 5,
                   g.perc = 1) {
  
  bins <- tapply(names(cl), 
                 cl, 
                 function(x){
                   if(length(x) > n.bin){
                     tmp <- rep_len(1:n.bin, length(x))
                   }else{
                     tmp <- sample(1:n.bin, length(x))
                   }
                   setNames(tmp[sample(length(tmp))], x)
                 })
  names(bins)=NULL
  bins <- unlist(bins)
  bins <- bins[names(cl)]
  pred.cl <- setNames(rep(NA, length(cl)), names(cl))
  
  for(i in 1:n.bin) {
    print(i)
    train.cells <- names(cl)[bins != i]
    test.cells <- names(cl)[bins == i]
    select.markers <- sample(markers, round(length(markers) * g.perc))
    
    map.result <- map_by_cor(norm.dat[select.markers,], 
                             cl[train.cells], 
                             norm.dat[select.markers, 
                                      test.cells])$pred.df
    
    pred.cl[test.cells] <- as.character(map.result[test.cells, "pred.cl"])
  }
  
  return(pred.cl)
}



plotAnnoStats <- function(peaks_anno, ttl, simplify=F, category_order=NULL) {
  anno_df <- data.frame()
  
  for(i in 1:length(peaks_anno)){
    tmp <- peaks_anno[[i]]@annoStat
    tmp$Cell <- names(peaks_anno)[i]
    tmp$Feature <- as.character(tmp$Feature)
    
    if(simplify){
      pos <- grep("Exon",tmp$Feature)
      tmp$Feature[pos] <- "Exon"
      
      pos <- grep("Intron",tmp$Feature)
      tmp$Feature[pos] <- "Intron"
      
      pos <- grep("Downstream",tmp$Feature)
      tmp$Feature[pos] <- "Distal Intergenic"
      
      pos <- grep("Promoter",tmp$Feature)
      tmp$Feature[pos] <- "Promoter"
      
      tmp$Feature <- gsub("5' UTR","Exon",tmp$Feature)
      tmp$Feature <- gsub("3' UTR","Exon",tmp$Feature)
      tmp$Feature <- gsub("Distal Intergenic","Intergenic",tmp$Feature)
      tmp$Feature <- factor(tmp$Feature, levels = rev(c("Promoter","Exon","Intron","Intergenic")))
    }
    
    anno_df <- rbind(anno_df, tmp)
  }
  
  anno_df = anno_df %>% group_by(Feature,Cell) %>% summarize(Frequency = sum(Frequency))
  
  anno_df$Feature <- factor(anno_df$Feature)  
  
  if(is.null(category_order)){
    anno_df$Cell <- factor(anno_df$Cell)  
  }else{
    anno_df$Cell <- factor(anno_df$Cell, levels = category_order)
  }
  
  
  g <- ggplot(anno_df, aes(x=Cell, y=Frequency, fill=Feature)) + 
    coord_flip() + 
    geom_bar(stat="identity") + theme_base() + 
    #scale_fill_brewer(palette = "Set1") +
    scale_fill_locuszoom(breaks=c("Promoter","Exon","Intron","Intergenic")) +
    ggtitle(ttl)+
    xlab("Cell") + ylab("Percentage (%)") + theme(plot.title = element_text(hjust = 0.5), 
                                                  axis.text = element_text(color="black",face = "bold")) #+ boldText
  g
}

KNNSmooth <- function(sobj, k=10){
  
  if(is.null(sobj@snn) | nrow(sobj@snn)==0){
    sobj = BuildSNN(sobj)
  }
  
  exprs  = expm1(sobj@data)
  smoothed = exprs * 0
  
  for(i in 1:nrow(sobj@snn)){
    neighbors <- which(sobj@snn[i,]>0)
    
    knn = order(sobj@snn[i,neighbors])[1:k]
    knn = names(neighbors)[knn]
    
    smoothed[,i] = rowMeans(sobj@data[,knn])
  }
  
  sobj = SetAssayData(sobj,  assay.type ="smoothed", "knn", smoothed)
  sobj
}



FindAllMarkers_parallel <- function(
  object,
  genes.use = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  return.thresh = 1e-2,
  do.print = FALSE,
  random.seed = 1,
  min.cells.gene = 3,
  min.cells.group = 3,
  latent.vars = NULL,
  assay.type = "RNA",
  ...
) {
  
  require(doParallel)
  require(foreach)
  
  data.1 <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  genes.use <- Seurat:::SetIfNull(x = genes.use, default = rownames(x = data.1))
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh = 0.7
  }
  idents.all <- sort(x = unique(x = object@ident))
  genes.de <- list()
  #if (max.cells.per.ident < Inf) {
  #  object <- SubsetData(
  #    object = object,
  #    max.cells.per.ident = max.cells.per.ident,
  #    random.seed = random.seed
  #  )
  #}
  
  registerDoParallel(4)
  
  genes.de = foreach(i =1:length(x = idents.all),.packages="Seurat") %dopar% {
    tryCatch(
      {
        FindMarkers(
          object = object,
          assay.type = assay.type,
          ident.1 = idents.all[i],
          ident.2 = NULL,
          genes.use = genes.use,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          print.bar = print.bar,
          min.cells.gene = min.cells.gene,
          min.cells.group = min.cells.group,
          latent.vars = latent.vars,
          max.cells.per.ident = max.cells.per.ident,
          ...
        )
      },
      error = function(cond){
        return(NULL)
      }
    )	
  }  
  
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      } else {
        gde <- gde[order(gde$p_val, -gde$avg_logFC), ]
        gde <- subset(x = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if ((only.pos) && nrow(gde.all) > 0) {
    return(subset(x = gde.all, subset = avg_logFC > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(gde.all) == 0) {
    warning("No DE genes identified.")
  }
  return(gde.all)
}




#https://www.r-bloggers.com/bigcor-large-correlation-matrices-in-r/
# Not very optimized but it should be much faster.

bigcor <- function(x, nblocks = 10, verbose = TRUE, method ="s", nb.cores=4)
{
  library(bigmemory, quietly = TRUE)
  library(bigstatsr, quietly = TRUE)
  
  
  NCOL <- ncol(x)
  
  ## test if ncol(x) %% nblocks gives remainder 0
  #if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
  
  
  blockIds  = cut(1:NCOL,nblocks)
  blockIds=as.numeric(blockIds)
  
  MAT = as_FBM(x)
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  corMAT <- FBM(NCOL, NCOL)
  
  ## split column numbers into 'nblocks' groups
  #SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(sort(unique(blockIds)), sort(unique(blockIds)))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  
  fct <- function(cMat,ind, mat, combs,splits, method){                  
    COMB <- combs[ind, ]        
    G1 <- which(splits == COMB[1])
    G2 <- which(splits == COMB[2])
    
    cMat[G1,G2] <- cor(mat[, G1], mat[, G2], method = method)    
    cMat[G2, G1] <- t(cMat[G1,G2])
    NULL    
  }
  
  big_apply(corMAT,a.FUN=fct,a.combine = 'c', 
            ind= rows_along(COMBS),
            mat = MAT,
            combs=COMBS, 
            splits=blockIds,        
            method=method,
            block.size=1, # We process a combination at a time (don't change this value)
            ncores=nb.cores)
  
  gc()
  
  res = corMAT[1:nrow(corMAT),1:ncol(corMAT)]
  colnames(res) = colnames(x)
  rownames(res) = colnames(x)
  return(res)
}



MetaNeighborUS_parallel <- function(var_genes, dat, i = 1, study_id, cell_type){
  library(bigmemory)
  library(bigstatsr)
  
  #dat    <- SummarizedExperiment::assay(dat, i = i)
  
  if(class(dat) !="matrix"){
    stop("dat should be a matrix")
  }
  
  samples <- colnames(dat)
  
  #check obj contains study_id
  if(length(study_id)!=length(samples)){
    stop('study_id length does not match number of samples')
  }
  
  #check obj contains cell_type
  if(length(cell_type)!=length(samples)){
    stop('cell_type length does not match number of samples')
  } 
  
  pheno <- as.data.frame(cbind(study_id,cell_type), stringsAsFactors = FALSE)
  pheno$StudyID_CT <- paste(pheno$study_id, pheno$cell_type, sep = "|")
  celltypes   <- unique(pheno$StudyID_CT)
  cell_labels <- matrix(0, ncol=length(celltypes), nrow=dim(pheno)[1])
  rownames(cell_labels) <-colnames(dat)
  colnames(cell_labels) <- celltypes
  
  for(i in seq_along(celltypes)){
    type <- celltypes[i]
    matching_celltype <- match(pheno$StudyID_CT, type)
    cell_labels[!is.na(matching_celltype),i]  <- 1
  }
  
  matching_vargenes <- match(rownames(dat), var_genes)
  matching_vargenes_count   <- sum(!is.na(matching_vargenes))
  
  if(matching_vargenes_count < 2){
    stop("matching_vargenes should have more than 1 matching genes!",
         call. = TRUE)
  } else if(matching_vargenes_count < 5) {
    warning("matching_vargenes should have more matching genes!", 
            immediate. = TRUE)
  }
  
  print("calculating correlation")
  cor_data    <- bigcor(dat[!is.na(matching_vargenes),],nblocks=30, method="s")
  rm(dat)
  gc()
  rank_data   <- cor_data*0
  rank_data[] <- rank(cor_data, ties.method = "average", na.last = "keep")
  rank_data[is.na(rank_data)] <- 0
  rank_data   <- rank_data/max(rank_data)
  sum_in      <- (rank_data) %*% cell_labels
  sum_all     <- matrix(apply(rank_data, MARGIN = 2, FUN = sum), 
                        ncol = dim(sum_in)[2], 
                        nrow = dim(sum_in)[1])
  predicts    <- sum_in/sum_all # nGenes * nClusters
  
  cell_NV     <- matrix(0, ncol=length(celltypes), nrow=length(celltypes))
  colnames(cell_NV) <- colnames(cell_labels)
  rownames(cell_NV) <- colnames(cell_labels)
  
  for(i in seq_len(dim(cell_labels)[2])){
    predicts_temp <- predicts
    
    matching_celltype <- match(pheno$StudyID_CT, colnames(cell_labels)[i])
    unique_studyID    <- unique(pheno[!is.na(matching_celltype),"study_id"])
    matching_studyID  <- match(pheno$study_id, unique_studyID)
    pheno2            <- pheno[!is.na(matching_studyID),]
    predicts_temp     <- predicts_temp[!is.na(matching_studyID),]
    predicts_temp     <- apply(abs(predicts_temp), 
                               MARGIN = 2, 
                               FUN = rank, 
                               na.last= "keep", 
                               ties.method="average")
    
    
    filter  <- matrix(0, ncol=length(celltypes), nrow=dim(pheno2)[1])
    matches <- match(pheno2$StudyID_CT, colnames(cell_labels)[i])
    filter[!is.na(matches),seq_along(celltypes)] <- 1
    
    negatives = which(filter == 0, arr.ind = TRUE)
    positives = which(filter == 1, arr.ind = TRUE)
    
    predicts_temp[negatives] <- 0
    
    np <- colSums(filter, na.rm = TRUE)
    nn <- apply(filter, MARGIN = 2, FUN = function(x) sum(x==0,na.rm=TRUE))
    p  <- apply(predicts_temp, MARGIN = 2, FUN = sum, na.rm = TRUE)
    
    cell_NV[i,]= (p/np - (np+1)/2)/nn
  }
  
  cell_NV <- (cell_NV+t(cell_NV))/2
  return(cell_NV)
}



MetaNeighborUS_CCA <- function(dat, i = 1, study_id, cell_type){
  library(bigmemory)
  library(bigstatsr)
  
  #dat    <- SummarizedExperiment::assay(dat, i = i)
  
  if(class(dat) !="matrix"){
    stop("dat should be a matrix")
  }
  
  samples <- colnames(dat)
  
  #check obj contains study_id
  if(length(study_id)!=length(samples)){
    stop('study_id length does not match number of samples')
  }
  
  #check obj contains cell_type
  if(length(cell_type)!=length(samples)){
    stop('cell_type length does not match number of samples')
  } 
  
  pheno <- as.data.frame(cbind(study_id,cell_type), stringsAsFactors = FALSE)
  pheno$StudyID_CT <- paste(pheno$study_id, pheno$cell_type, sep = "|")
  celltypes   <- unique(pheno$StudyID_CT)
  cell_labels <- matrix(0, ncol=length(celltypes), nrow=dim(pheno)[1])
  rownames(cell_labels) <-colnames(dat)
  colnames(cell_labels) <- celltypes
  
  for(i in seq_along(celltypes)){
    type <- celltypes[i]
    matching_celltype <- match(pheno$StudyID_CT, type)
    cell_labels[!is.na(matching_celltype),i]  <- 1
  }
  
  #matching_vargenes <- match(rownames(dat), var_genes)
  #matching_vargenes_count   <- sum(!is.na(matching_vargenes))
  
  # if(matching_vargenes_count < 2){
  #   stop("matching_vargenes should have more than 1 matching genes!",
  #        call. = TRUE)
  # } else if(matching_vargenes_count < 5) {
  #   warning("matching_vargenes should have more matching genes!", 
  #           immediate. = TRUE)
  # }
  
  print("calculating correlation")
  cor_data    <- bigcor(dat,nblocks=30, method="s")
  rm(dat)
  gc()
  rank_data   <- cor_data*0
  rank_data[] <- rank(cor_data, ties.method = "average", na.last = "keep")
  rank_data[is.na(rank_data)] <- 0
  rank_data   <- rank_data/max(rank_data)
  sum_in      <- (rank_data) %*% cell_labels
  sum_all     <- matrix(apply(rank_data, MARGIN = 2, FUN = sum), 
                        ncol = dim(sum_in)[2], 
                        nrow = dim(sum_in)[1])
  predicts    <- sum_in/sum_all # nGenes * nClusters
  
  cell_NV     <- matrix(0, ncol=length(celltypes), nrow=length(celltypes))
  colnames(cell_NV) <- colnames(cell_labels)
  rownames(cell_NV) <- colnames(cell_labels)
  
  for(i in seq_len(dim(cell_labels)[2])){
    predicts_temp <- predicts
    
    matching_celltype <- match(pheno$StudyID_CT, colnames(cell_labels)[i])
    unique_studyID    <- unique(pheno[!is.na(matching_celltype),"study_id"])
    matching_studyID  <- match(pheno$study_id, unique_studyID)
    pheno2            <- pheno[!is.na(matching_studyID),]
    predicts_temp     <- predicts_temp[!is.na(matching_studyID),]
    predicts_temp     <- apply(abs(predicts_temp), 
                               MARGIN = 2, 
                               FUN = rank, 
                               na.last= "keep", 
                               ties.method="average")
    
    
    filter  <- matrix(0, ncol=length(celltypes), nrow=dim(pheno2)[1])
    matches <- match(pheno2$StudyID_CT, colnames(cell_labels)[i])
    filter[!is.na(matches),seq_along(celltypes)] <- 1
    
    negatives = which(filter == 0, arr.ind = TRUE)
    positives = which(filter == 1, arr.ind = TRUE)
    
    predicts_temp[negatives] <- 0
    
    np <- colSums(filter, na.rm = TRUE)
    nn <- apply(filter, MARGIN = 2, FUN = function(x) sum(x==0,na.rm=TRUE))
    p  <- apply(predicts_temp, MARGIN = 2, FUN = sum, na.rm = TRUE)
    
    cell_NV[i,]= (p/np - (np+1)/2)/nn
  }
  
  cell_NV <- (cell_NV+t(cell_NV))/2
  return(cell_NV)
}



readIPAResults <- function(path, pathway_pval=0.05, disease_pval =0.05){
  require(readxl)
  IPA_results <- list.files(path,pattern = "*.xls",full.names = T)
  
  
  pathway_enrichment <- data.frame()
  disease_enrichemt <- data.frame()
  
  for(f in IPA_results){
    
    smp = gsub(".xls","",basename(f))
    print(smp)
    #smp = unlist(strsplit(smp,split = "_"))
    
    tmp = read_excel(f)
    
    colnames(tmp)[1]="info"
    
    headers_pos = grep("My Project", tmp$info)
    
    # Get pathways
    pathways_hdr = grep("Canonical Pathways",tmp$info[headers_pos])
    pathway_tbl = tmp[(headers_pos[pathways_hdr]+1):(headers_pos[pathways_hdr+1]-1),]
    colnames(pathway_tbl) = pathway_tbl[1,]
    pathway_tbl = pathway_tbl[-1,]
    pathway_tbl = pathway_tbl[,!is.na(colnames(pathway_tbl))]
    pathway_tbl = pathway_tbl[,colnames(pathway_tbl)!="NA"]
    
    pathway_tbl$cluster = smp
    colnames(pathway_tbl)[1:2] = c("Pathway","logPval")
    pathway_tbl$logPval = as.numeric(pathway_tbl$logPval)
    
    pathway_tbl <- subset(pathway_tbl, logPval > -log10(pathway_pval))
    pathway_tbl <- pathway_tbl[order(pathway_tbl$logPval,decreasing = T),]
    #pathway_enrichment <- rbind(pathway_enrichment,head(pathway_tbl,5))
    pathway_enrichment <- rbind(pathway_enrichment,pathway_tbl)
    
    # Get diseases
    
    diseases_hdr = grep("Diseases and Bio Functions",tmp$info[headers_pos])
    disease_tbl = tmp[(headers_pos[diseases_hdr]+1):(headers_pos[diseases_hdr+1]-1),]
    colnames(disease_tbl) = disease_tbl[1,]
    disease_tbl = disease_tbl[-1,]
    disease_tbl = disease_tbl[,!is.na(colnames(disease_tbl))]
    disease_tbl = disease_tbl[,colnames(disease_tbl) != "NA"]
    disease_tbl$cluster = smp
    
    colnames(disease_tbl)[c(3:4,10)] <- c("Disease","pvalue","nbMolecules")
    disease_tbl$pvalue = as.numeric(disease_tbl$pvalue)
    disease_tbl <- subset(disease_tbl, pvalue < disease_pval)
    
    disease_tbl <- disease_tbl[order(disease_tbl$pvalue),]
    
    #disease_enrichemt <- rbind(disease_enrichemt, head(disease_tbl,5))
    disease_enrichemt <- rbind(disease_enrichemt, disease_tbl)
  }
  
  
  res = list(disease = disease_enrichemt, pathway = pathway_enrichment)
  
  res
}





my_eset2Phase <- function(eset, low.prob = 0.99) {
  Y <- round(exprs(eset))
  Cell0 = colMeans(Y == 0)
  par1 = apply(Y, 2, function(yy) {
    yy = yy[yy <= 15]
    RobustPoi0(yy)
  })
  pi0.hat = Cell0/(par1[1, ] + (1 - par1[1, ]) * dpois(0, 
                                                       par1[2, ]))
  if (any((pi0.hat > 1))) {
    warning("Zero proportion is greater than estimation.")
  }
  pi0.hat <- pmin(pi0.hat, 1)
  prob0 = pi0.hat * par1[1, ] + pi0.hat * (1 - par1[1, ]) * 
    dpois(0, par1[2, ])
  x0 = qpois(pmax(1 - (1 - low.prob)/(1 - par1[1, ]), 0), 
             par1[2, ])
  Z = sweep(Y, 2, x0) > 0
  L = colSums(Y * Z)/1e+06
  mu.g1 = log2(rowSums(Z * Y)/rowSums(sweep(Z, 2, L, FUN = "*")))
  mu.g1[is.na(mu.g1)] = 0
  n.g1 = rowSums(Z)
  y1 = log2(sweep(Y, 2, L, FUN = "/") + 1)
  s.g1 = sqrt(rowSums(Z * sweep(y1, 1, mu.g1)^2)/(n.g1 - 1))
  mu.g2 = shrink.mu(mu.g1, s.g1, n.g1)
  res.g1 = log2(sweep(Y, 2, L, FUN = "/") + 1) - mu.g1
  tmp = array(0, dim = c(dim(res.g1), 2))
  tmp[, , 1] = res.g1
  tmp[, , 2] = Z
  sd.g1 = apply(tmp, 1, function(xx) my.mad(xx[xx[, 2] == 
                                                 1, 1]))
  sd.g1[is.na(sd.g1)] = 0
  sd.prior = squeezeVar(sd.g1^2, n.g1 - 1)
  sd.g2 = sqrt(sd.prior$var.post)
  den.fg = den.bg = NA * Y
  for (i in 1:ncol(Y)) {
    den.bg[, i] = dZinf.pois(Y[, i], par1[1, i], par1[2, 
                                                      i])
    den.fg[, i] = dLNP2(x = Y[, i], mu = mu.g1, sigma = sd.g2, 
                        l = L[i])
  }
  Z.fg = sweep(den.fg, 2, 1 - pi0.hat, FUN = "*")
  Z.bg = sweep(den.bg, 2, pi0.hat, FUN = "*")
  post.Z = Z.fg/(Z.fg + Z.bg)
  post.Z[is.na(post.Z)] <- 1
  den.fg2 = NA * Y
  for (i in 1:ncol(Y)) {
    den.fg2[, i] = dLNP2(x = Y[, i], mu = mu.g2, sigma = sd.g2, 
                         l = L[i])
  }
  Z.fg2 = sweep(den.fg2, 2, 1 - pi0.hat, FUN = "*")
  post.Z2 = Z.fg2/(Z.fg2 + Z.bg)
  post.Z2[is.na(post.Z2)] <- 1
  Offset = Y * 0
  Ylim = range(log2(1 + Y) - mu.g1)
  Xlim = range(mu.g1)
  for (i in 1:ncol(Y)) {
    tmp.y = log2(1 + Y[, i]) - mu.g2
    subset = post.Z2[, i] > 0.99
    lm1 <- loess(tmp.y ~ mu.g1, weights = post.Z2[, i] * 
                   mu.g2, subset = subset, degree = 1, span = 0.3)
    Offset[subset, i] = lm1$fitted
  }
  fdata <- fData(eset)
  fdata2 <- as.data.frame(cbind(fdata, mu.g2, sd.g2))
  colnames(fdata2) <- c(colnames(fdata), "mean", "sd")
  fvar <- rbind(fvarMetadata(eset), mean = "shrinkage estimated foreground mean", 
                sd = "shrinkage estimated foreground standard deviation")
  featureData <- new("AnnotatedDataFrame", data = fdata2, 
                     varMetadata = fvar)
  pdata <- pData(eset)
  pdata2 <- as.data.frame(cbind(pdata, par1[1, ], par1[2, 
                                                       ], L))
  colnames(pdata2) <- c(colnames(pdata), "p0", "lambda", "L")
  pvar <- rbind(varMetadata(eset), p0 = "proportion of zero inflation", 
                lambda = "mean of background poisson", L = "foreground library size")
  phenoData <- new("AnnotatedDataFrame", data = pdata2, varMetadata = pvar)
  out <- new("sc2pSet", exprs = Y, Z = post.Z2, Offset = Offset, 
             phenoData = phenoData, featureData = featureData, experimentData = experimentData(eset), 
             annotation = annotation(eset))
  out
}


classNNMap <- function(sobj, feature.name = NULL, k.param=500, min.k = 10){
  
  assay <- DefaultAssay(object = sobj[["iNMF"]])
  data.use <- Embeddings(object = sobj[["iNMF"]])  
  
  data.use <- data.use[,1:20]    
  
  data.use = as.matrix(data.use)
  
  
  pos_merfish <- which(sobj@meta.data$Method == "MERFISH")
  pos_rnaseq <- which(sobj@meta.data$Method == "scRNASeq")
  
  nn.ranked <- Seurat:::NNHelper(
    data = data.use,
    k = k.param,
    method = 'rann',
    searchtype = "standard",
    eps = 0,
    metric = "euclidean")
  
  nn.ranked <- nn.ranked$nn.idx
  
  ## check if each merfish cell has at least min.k scRNASeq neighbors
  min.neighbors = k.param    
  
  
  res = pblapply(pos_merfish, function(i) {
    #for(i in pos_merfish){
    nn = intersect(nn.ranked[i,],pos_rnaseq)
    
    ## if there are less scRNA-seq neighbors calculate individually
    if(length(nn) < min.k){
      cells.use = c(i,pos_rnaseq)
      
      data.use2 = data.use[cells.use,]
      
      nn.ranked2 <- Seurat:::NNHelper(
        data = data.use2,
        k = min.k+1,
        method = 'rann',
        searchtype = "standard",
        eps = 0,
        metric = "euclidean")
      
      nn = nn.ranked2$nn.idx[1,-1]
      nn = pos_rnaseq[nn-1]
    }
    return( nn[1:min.k])
  })
  #  nn.scrnaseq[p,] = nn[1:min.k]
  #  p = p+1
  #}
  
  nn.scrnaseq = do.call("rbind",res)
  rownames(nn.scrnaseq) = rownames(sobj@meta.data)[pos_merfish]
  return(nn.scrnaseq)
}




rotate <- function(x,y,degree){
  
  theta = (degree * pi) / (180)
  
  res = list( new_x =  x * cos(theta) - y * sin(theta),
              new_y =  x * sin(theta) + y * cos(theta)
  )
  return(res)
}

fixMERFISH_coords <- function(sobj, rotations = NULL, flips= NULL, slices = 1:12, ncol=4){
  
  for(slice_id in slices){
    
    slice_pos <- which(sobj@meta.data$sliceNumber == slice_id)
    
    new_x = sobj@meta.data$centroid_x[slice_pos]
    new_y = sobj@meta.data$centroid_y[slice_pos]
    
    if(flips[slice_id]){
      new_y = new_y + 2 * ((1500 - 5500* ((slice_id-1) %/% ncol)) - new_y )    
    }
    
    rot = rotate(new_x, new_y, rotations[slice_id])
    new_x = rot[['new_x']]
    new_y = rot[['new_y']]
    
    xcenter = (max(new_x) + min(new_x))/2 
    ycenter = (max(new_y) + min(new_y))/2  
    
    # Slide to the new position
    new_x = (new_x - xcenter + 1500) + 5500* ((slice_id-1) %%  ncol)
    new_y = (new_y - ycenter + 1500) - 5500* ((slice_id-1) %/% ncol)   
    
    
    sobj@meta.data$centroid_x[slice_pos] = new_x
    sobj@meta.data$centroid_y[slice_pos] = new_y        
  }
  
  return(sobj)
}

RobustAverageExpression <- function(sobj, assay ="RNA",slot="data"){
  
  
  if(assay != "RNA"){
    sobj_expr <- GetAssayData(sobj, assay.type = assay, slot = slot)
    sobj_expr = expm1(sobj_expr)
  }else{
    sobj_expr = expm1(sobj@data)  
  }
  
  AE_expression <- matrix(0,nrow = nrow(sobj_expr), ncol = length(levels(sobj@ident)))
  rownames(AE_expression)= rownames(sobj_expr)
  colnames(AE_expression) <- levels(sobj@ident)
  
  for(clus in levels(sobj@ident)){
    
    pos <- which(sobj@ident == clus)
    
    cls_mean  <- apply(sobj_expr[,pos],1,function(x){
      q99= as.numeric(quantile(x,probs = 0.99))
      x = x[x<=q99]
      if(length(x)>0){
        mean(x)  
      }else{
        0
      }
    })
    
    AE_expression[,clus]= cls_mean
    print(paste0("Finished processing cluster: ",clus))
  }
  
  AE_expression = log2(AE_expression+1)
  AE_expression
}


RobustAverageDetection <- function(sobj, assay ="RNA",slot="data"){
  
  
  if(assay != "RNA"){
    sobj_expr <- GetAssayData(sobj, assay.type = assay, slot = slot)
    sobj_expr = expm1(sobj_expr)
  }else{
    sobj_expr = expm1(sobj@data)  
  }
  
  AE_detection <- matrix(0,nrow = nrow(sobj_expr), ncol = length(levels(sobj@ident)))
  rownames(AE_detection)= rownames(sobj_expr)
  colnames(AE_detection) <- levels(sobj@ident)
  
  for(clus in levels(sobj@ident)){
    
    pos <- which(sobj@ident == clus)
    
    cls_percent  <- apply(sobj_expr[,pos],1,function(x){
      
      q99= as.numeric(quantile(x,probs = 0.99))
      x = x[x<q99]
      if(length(x)>0){
        sum(x>=1)/length(x)        
      }else{
        0
      }
      
    })
    
    AE_detection[,clus]= cls_percent
    print(paste0("Finished processing cluster: ",clus))
  }
  
  AE_detection
}



getPercentExpressed <- function(sobj){
  
  pctExpr = levels(Idents(sobj)) %>%
    map(.f = function(cls){
      print(cls)
      cells.1 = WhichCells(object = sobj, idents = cls)
      
      pct.1 <- round(
        x = Matrix::rowSums(x = sobj[["RNA"]]@data[, cells.1, drop = FALSE] > 0) /
          length(x = cells.1),
        digits = 3
      )
      
      df = tibble( pct = pct.1,.name_repair = ~cls)
      return(df)
    }) %>%
    do.call(what = "cbind") %>%
    data.matrix()
  
  rownames(pctExpr) = rownames(sobj[["RNA"]]@data)
  
  return(pctExpr)
}
