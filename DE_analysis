DE_analysis <- function(threshold = 0.01,
                        geo_accession,
                        gpl_id,
                        gsms, 
                        outdir = "data/"){
  
  # Required libraries
  library(GEOquery)
  library(limma)
  library(data.table)
  library(dplyr)
  
  dir.create(outdir)
  # Load series and platform data from GEO
  gset <- getGEO(geo_accession, GSEMatrix = TRUE, AnnotGPL = TRUE)
  if (length(gset) > 1) idx <- grep(gpl_id, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  
  # Make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  
  # Group membership for all samples
  sml <- strsplit(gsms, split="")[[1]]
  # sml <- paste0(c(rep("0", 3),rep("1", 3)))
  
  # Log2 transformation
  ex <- exprs(gset)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (LogC) { 
    ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) 
  }
  exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # Normalize data
  
  # Assign samples to groups and set up design matrix
  gs <- factor(sml)
  groups <- make.names(c("Control", "Case"))
  levels(gs) <- groups
  gset$group <- gs
  design <- model.matrix(~group + 0, gset)
  colnames(design) <- levels(gs)
  
  gset <- gset[complete.cases(exprs(gset)), ] # Skip missing values
  
  
  # calculate precision weights and show plot of mean-variance trend
  v <- vooma(gset, design, plot=T)
  # OR weights by group
  # v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
  v$genes <- fData(gset) # attach gene annotations
  
  # fit linear model
  fit  <- lmFit(v)
  
  # set up contrasts of interest and recalculate model coefficients
  cts <- paste(groups[1], groups[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
  
  # tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","Gene.symbol","Gene.title"))
  tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.ID","Gene.symbol", "Gene.title"))
  # write.table(tT, file=stdout(), row.names=F, sep="\t")
  
  
  fwrite(tT, file=paste0(outdir, "DE_genes_limma.csv"))
  return(tT)
  
}
