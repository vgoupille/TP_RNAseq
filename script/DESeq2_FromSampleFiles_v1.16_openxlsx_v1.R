##########################################################################################
#
# DESeq2_package_v1.16.R
#
# OBJ : differential analysis of raw counts with DESEq2
#
# input : one file (.txt) = raw counts for one sample (integer format)
#
# version fonctionnant avec DESeq2 version 1.16.1
#
# parameters
#   table     : sample info with association sample names & Condition
#   dir       : sample files directory (one file = raw counts of one sample)
#   pval      : pvalue for filteringDE genes
#   FC        : |FC| for filtering DE genes
#   analysis  : analysis type = "All" | "Ctrl"
#   condCtrl  : control condition (NULL by default)
#   PreFilt   : pre-filtering : NULL | integer : keep only rows that have at least
#                      PreFilt reads total (NULL per default : no pre-filtering)
#   Filt      : HTS Filtering = TRUE | FALSE (FALSE per default)
#   NbCountMin : min number of normalized counts in 100% samples (NULL per default)
#   paired    : paired analysis = TRUE | FALSE (FALSE per default)
#   correction : multiple testing correction = "bonferroni" | "fdr"
#   record    : enregistrer les resultats DESeq2 sans filtrage par p ou FC = TRUE | FALSE (FALSE per default)
#
# initial   2017 12 19 Fabrice Chatonnet
# modified  2018 04 18 DR
#   -ajout moderated FC et plotMA
#   -ajout variable 'record' = TRUE : sauvegarde des r�sultats DE sans filtrage
#     Attention: cette option TRUE peut creer une erreur java d'enregistrement dans fichier .xlsx
#
##########################################################################################

# ---- function for samples PCA analysis
plot_PCAplot = function(dds_pca, title_pca, nb=NULL){
  library(ggplot2)
  if (is.null(nb)) {nb = length(rownames(dds_pca))}
  vsd = varianceStabilizingTransformation(dds_pca)
  PCA_plot = plotPCA(vsd, intgroup=colnames(colData(dds_pca)[1]), ntop=nb,
                     returnData=T)
  percentVar = round(100 * attr(PCA_plot, "percentVar"))
  print(ggplot(PCA_plot, aes(PC1, PC2, color=PCA_plot[,4], label=PCA_plot[,3])) +
          geom_point(size=3) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")) + ggtitle(title_pca) +
          scale_color_discrete(name=colnames(PCA_plot)[4])
        + geom_text(aes(label = rownames(colData(dds_pca))), hjust = 0, nudge_x = 0.01*(max(PCA_plot$PC1) - min(PCA_plot$PC1)))
        )
}


############################################################################################
# ---- function for samples clustering analysis
plot_heatmap = function(dds_hm, title_hm){
  library(RColorBrewer)
  library(gplots)
  hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
  rld = rlogTransformation(dds_hm, blind=TRUE)
  distsRL = dist(t(assay(rld)))
  mat = as.matrix(distsRL)
  rownames(mat) = colnames(mat) = paste(rownames(colData(dds_hm)), colData(dds_hm)[,1], sep="_")
  heatmap.2(mat, col = rev(hmcol), main=title_hm, key=F,
            density.info="none", trace="none", margins = c(10,10))
}


############################################################################################
# ---- function for comparisons between conditions and counts table generation
DESeq2_DE_res = function(Exp, dds_res, pval, FC, i, j, correction, record){
  
  logFC = log(FC,2)
  Nb_i = length(colnames(dds_res[,colData(dds_res)[,1]==levels(colData(dds_res)[,1])[i]]))
  Nb_j = length(colnames(dds_res[,colData(dds_res)[,1]==levels(colData(dds_res)[,1])[j]]))
  if (i>j)
  {
    name_ij=paste(colnames(Exp)[3], levels(colData(dds_res)[,1])[i], "vs",
                  levels(colData(dds_res)[,1])[j], sep="_")
  } else {
    name_ij=paste(colnames(Exp)[3], levels(colData(dds_res)[,1])[j], "vs",
                  levels(colData(dds_res)[,1])[i], sep="_")  
  }

  #print(paste("contrast : ", name_ij, sep=" "))
  res_ij = results(dds_res, contrast=c(colnames(Exp)[3], levels(colData(dds_res)[,1])[j],
                                       levels(colData(dds_res)[,1])[i]))

  # --- plot MA - no moderated FC
  # print("--- plot MA - no moderated FC")
  plotMA(res_ij, main=paste(name_ij, " - no moderated FC", "\nDESeq2 (test de Wald)", sep=""))
  
  # --- plot MA - moderated FC
  # print("--- calculation of moderated FC")
  coeff <- which(resultsNames(dds_res) == name_ij)
  # print(paste("coeff : ", coeff, sep=" "))
  res.lfc = lfcShrink(dds_res, coef=coeff, res=res_ij)
  #print(head(res.lfc))
  
  # print("--- plot MA - with moderated FC")
  plotMA(res.lfc, main=paste(name_ij, " - with moderated FC", "\nDESeq2 (test de Wald)", sep=""))
  
  
  # --- calcul FC lin�aire
  # print("--- calculation of linear FC")
  res.lfc <- as.data.frame(res.lfc)
  res.lfc$FC <- 2^res.lfc$log2FoldChange
  res.lfc$FC[which(res.lfc$log2FoldChange<0)] <- -2^(-res.lfc$log2FoldChange[which(res.lfc$log2FoldChange<0)])
  
  
  if(Nb_i==1 || Nb_j==1){DE = res.lfc$pvalue < pval & !is.na(res.lfc$pvalue)}
  else {res.lfc$padj = p.adjust(res.lfc$pvalue, method = correction)
  DE = res.lfc$padj < pval & !is.na(res.lfc$padj)}
  res_DE = res.lfc[DE,]
  dds_DE = dds_res[DE,]

  # ---- Filtering by fold change
  FCsel = abs(res_DE$log2FoldChange) > logFC
  res_FC = res_DE[FCsel,]
  dds_resC = dds_DE[FCsel,]
  if(length(row.names(res_DE))==0) {res_DE="none"}
  if(length(row.names(res_FC))==0) {res_FC="none"}

  # ---- Calculation of count means by condition
  if(Nb_i==1){bM_i = counts(dds_res, normalized = TRUE)[,colData(dds_res)[,1]==
                                                          levels(colData(dds_res)[,1])[i]]}
  else {bM_i = rowMeans(counts(dds_res, normalized=T)[,colData(dds_res)[,1]==
                                                        levels(colData(dds_res)[,1])[i]])}
  if (Nb_j==1){bM_j = counts(dds_res, normalized = TRUE)[,colData(dds_res)[,1]==
                                                           levels(colData(dds_res)[,1])[j]]}
  else{bM_j = rowMeans(counts(dds_res, normalized=T)[,colData(dds_res)[,1]==
                                                       levels(colData(dds_res)[,1])[j]])}
  bM_i_DE = bM_i[DE]
  bM_j_DE = bM_j[DE]
  bM_i_FC = bM_i_DE[FCsel]
  bM_j_FC = bM_j_DE[FCsel]

  # ---- Correlation graph of all genes (black), DE (red), FoldChange > FC (blue)
  # print("--- correlation graph")
  plot(bM_i, bM_j, pch=20, col = 1, log="xy",
       xlim=c(10,max(counts(dds_res, normalized=T))),
       ylim=c(10,max(counts(dds_res, normalized=T))),
       xlab=levels(colData(dds_res)[,1])[i], ylab=levels(colData(dds_res)[,1])[j],
       main=paste("Genes differentially expressed between", levels(colData(dds_res)[,1])[i],"and",
                  levels(colData(dds_res)[,1])[j], sep=" ")
  )
  points(bM_i_DE, bM_j_DE, pch=19, col = "red")
  points(bM_i_FC, bM_j_FC, pch=19, col = rgb(0,150,130, maxColorValue=255))
  legend("bottomright", c(paste("Genes not differentially expressed (padj > ", pval, ")", sep=""),
                          paste("Genes differentially expressed (padj < ", pval, ")", sep= ""),
                          paste("Genes differentially expressed (padj < ", pval, ")",  
                                " with a fold change > ", FC, sep= "")),
        pch=20, col=c(1, "red", rgb(0,150,130, maxColorValue=255)))

  # ---- Count files (xls format for DE and FC genes)
  # creation of tables with means for each level and count values for each sample of this level
  counts_Filt = cbind(as.data.frame(bM_i), as.data.frame(bM_j),
                      counts(dds_res, normalized = TRUE)[,colData(dds_res)[,1]==levels(colData(dds_res)[,1])[i]],
                      counts(dds_res, normalized = TRUE)[,colData(dds_res)[,1]==levels(colData(dds_res)[,1])[j]])
  colnames(counts_Filt) = c(levels(colData(dds_res)[,1])[i], levels(colData(dds_res)[,1])[j],
                            colnames(counts(dds_res))[colData(dds_res)[,1]==levels(colData(dds_res)[,1])[i]],
                            colnames(counts(dds_res))[colData(dds_res)[,1]==levels(colData(dds_res)[,1])[j]])
  counts_DE = counts_Filt[DE,]
  if(length(row.names(counts_DE))==0) {counts_DE="none"; counts_FC="none"
  cat("\n --- No DE genes found")}
  else {counts_FC = counts_DE[FCsel,]
        cat(paste("\n --- DESEq2 analysis retrieved", length(row.names(counts_DE)), "genes with a padj of",
                    pval, "of which", length(row.names(counts_FC)), "had a fold change of at least", FC))}

  
  wb <- createWorkbook()
  
  # ---- file creation
  if (record==TRUE)
  {
    addWorksheet(wb, "counts")
    writeDataTable(wb, sheet = "counts", x=as.data.frame(counts_Filt), rowNames = TRUE)

    addWorksheet(wb, "pval")
    writeDataTable(wb, sheet = "pval", x=as.data.frame(res.lfc), rowNames = TRUE)
    
    addWorksheet(wb, paste("counts", "DE", sep=" "))
    writeDataTable(wb, sheet = paste("counts", "DE", sep=" "), x=as.data.frame(counts_DE), rowNames = TRUE)
    

  } else {
    addWorksheet(wb, paste("counts", "DE", sep=" "))
    writeDataTable(wb, sheet = paste("counts", "DE", sep=" "), x=as.data.frame(counts_DE), rowNames = TRUE)
    
  }
  
  addWorksheet(wb, paste("pval", "DE", sep=" "))
  writeDataTable(wb, sheet = paste("pval", "DE", sep=" "), x=as.data.frame(res_DE), rowNames = TRUE)
  
  addWorksheet(wb, paste("counts", "FC", sep=" "))
  writeDataTable(wb, sheet = paste("counts", "FC", sep=" "), x=as.data.frame(counts_FC), rowNames = TRUE)
  
  addWorksheet(wb, paste("pval", "FC", sep=" "))
  writeDataTable(wb, sheet = paste("pval", "FC", sep=" "), x=as.data.frame(res_FC), rowNames = TRUE)
  
  Info=c(paste("contains for differentially expressed (DE) genes (padj <", pval, "): ",
               "counts means for the compared conditions of the ", colnames(Exp)[3], " factor: ",
               levels(colData(dds_res)[,1])[i], " (column A) and ", levels(colData(dds_res)[,1])[j], " (column B)",
               " and normalized counts for each sample ordered by condition", sep =""),
         paste("contains for differentially expressed (DE) genes (padj <", pval, "): ",
               "statistical data (baseMean, log2foldChange, pval, padj...) for the compared conditions of the ",
               colnames(Exp)[3], " factor: ", levels(colData(dds_res)[,1])[i], " and ", levels(colData(dds_res)[,1])[j],
               " as given by the DESeq2 results() function", sep =""),
         paste("contains the same info as Sheet counts DE ", name_ij,
               " but only for genes with a fold change (FC) > ", FC, sep =""),
         paste("contains the same info as Sheet pval DE ", name_ij,
               " but only for genes with a fold change (FC) > ", FC, sep ="")
        )
  
  addWorksheet(wb, "Contents information")
  writeDataTable(wb, sheet = "Contents information", 
                 x=as.data.frame(Info, row.names=c(paste("Sheet counts DE ", name_ij, sep=""),
                                                   paste("Sheet pval DE ", name_ij, sep=""),
                                                   paste("Sheet counts FC ", name_ij, sep=""),
                                                   paste("Sheet pval FC ", name_ij, sep=""))), rowNames = TRUE)

  saveWorkbook(wb, paste(name_ij,".xlsx", sep=""), overwrite=TRUE)
  
    }


############################################################################################
# ---- main function
DESeq2_FromSampleFiles = function(table, dir, pval, FC, analysis=c("All", "Ctrl"), condCtrl=NULL, PreFilt=NULL, 
                                  Filt=T, NbCountMin=NULL, paired=F, correction=c("fdr", "BH"), record=F) {
  suppressMessages(library(DESeq2))
  library(openxlsx)
  library(lattice)

  cat("Parameters")
  cat(paste("\n --- Table :", table, sep=""))
  cat(paste("\n --- Dir :", dir, sep=" "))
  cat(paste("\n --- Pvalue :", pval, sep=" "))
  cat(paste("\n --- FC :", FC, sep=" "))
  cat(paste("\n --- Analysis Type :", analysis, sep=" "))
  cat(paste("\n --- Control condition :", condCtrl, sep=" "))
  cat(paste("\n --- Pre-filtering : keep rows with at least ", PreFilt, " reads total", sep=" "))
  cat(paste("\n --- HTS Filtering :", Filt, sep=" "))
  cat(paste("\n --- Manual Filtering - Number of min counts :", NbCountMin, sep=" "))
  cat(paste("\n --- Paired analysis :", paired, sep=" "))
  cat(paste("\n --- Multiple Test Correction :", correction, sep=" "))
  cat(paste("\n --- Record even no filtered DESeq2 results :", record, sep=" "))
  
  

  # ---- Pour travailler sur les gene_symbol et non sur les ENSG 
  #
  # correspondance avec GRCh38 Ensembl version 94 => fichier d'annotations grch38_ensembl94.xlsx
  #
  # contacter Vonnick ou Delphine pour utiliser une autre version d'Ensembl
  #
  # noter qu'avec cette transformation des variables, certains ENSG sont supprim�s car ils n'ont pas de gene_symbol
  # ce qui explique une PCA l�g�rement diff�rente (en contributions)
  
  grch38 <- read.xlsx("data/grch38_ensembl94.xlsx", sheet=1)
  
  
  # ---- import conditions table
  Exp = read.csv(table, sep = " ", header=T) # Conditions table
  


  if (paired==F)
    {
      d = formula(paste("~",colnames(Exp)[3]))
      cat("\n Unpaired analysis")
    } else {
      d = formula(paste("~", colnames(Exp)[3], "+", colnames(Exp)[4]))
      cat("\n Paired analysis")
    }
  
  



  
  pdf(file=paste("Graph_report_",colnames(Exp)[3],".pdf", sep=""), width = 11.69, height = 8.27) # creates a graph report pdf file


  # ---- creation of the counts table
  cat("\n Creation of the counts table")
  dds = DESeqDataSetFromHTSeqCount(Exp, directory=dir, design=d)



  # --- pre-filtering
  if (!is.null(PreFilt))
  {
    cat("\n Pre-filtrage : Nombre de lignes initial : ", length(rowSums(counts(dds))))
    keepPreFilt <- rowSums(counts(dds)) >= PreFilt
    dds <- dds[keepPreFilt, ]
    cat("\n Pre-filtrage : Nombre de lignes conservees : ", length(which(keepPreFilt==TRUE)))
      
  }
  
  # ---- nb of levels in the comparison factor
  NbLvl = length(levels(colData(dds)[,1]))
  NbSamples = length(colData(dds)[,1])
  logFC = log(FC,2)

  # ---- DE global analysis and outputs
  cat("\n DESeq - test global LRT")
  dds_a <- DESeq(dds, test = "LRT", reduced= ~1, quiet=TRUE)
  res = results(dds_a)
  
  plotDispEsts(dds_a, main = "Global Dispersion") # dispersion graph
  plot_PCAplot(dds_a, "PCA Global Gene Expression") # PCA plot
  plot_heatmap(dds_a, "Global Gene Expression") # samples heatmap
  
  wb <- createWorkbook()
  addWorksheet(wb, "Counts_GlobalAnalysis")
  writeDataTable(wb, sheet = "Counts_GlobalAnalysis", x=as.data.frame(counts(dds_a, normalized = TRUE)), rowNames = TRUE)
  addWorksheet(wb, "Pval_GlobalAnalysis")
  writeDataTable(wb, sheet= "Pval_GlobalAnalysis", x= as.data.frame(res[,c(1,5,6)]), rowNames=TRUE)
               
  # ---- Filtering
  cat("\n Filtering")
  dds_f <- dds_a
  
  if(Filt==F&is.null(NbCountMin))
  {
    cat("\n No filtering of data")
  } else {
    if(Filt==T&!is.null(NbCountMin))
    {
      cat("\n No filtering of data : you haven't chosen one filtering method, either HTSFilter or manual")
    } else {
      if(Filt==T)
      {
        library(HTSFilter)
        cat("\n Filtering by HTSFilter")
        colData(dds_f)=colData(dds_f)[1]
        dds_f=HTSFilter(dds_f, normalization="DESeq", plot=T)$filteredData
        title("Normalized read threshold computed by HTSFilter")
        dds_f = dds_a[rownames(dds_f),]
      } else {
        cat("\n Manual filtering")
        dds_f = dds_a[!rowSums(counts(dds_a, normalized=T)<=NbCountMin)==NbSamples,]
      }
      plot_PCAplot(dds_f, "PCA Filtered Gene Expression")
      plot_heatmap(dds_f, "Filtered Gene Expression")
      
      addWorksheet(wb, "Counts_GlobalAnalysis_Filt")
      writeDataTable(wb, sheet = "Counts_GlobalAnalysis_Filt", x=as.data.frame(counts(dds_f, normalized = TRUE)), rowNames = TRUE)

      res_f = res[rownames(dds_f),]
      addWorksheet(wb, "Pval_GlobalAnalysis_Filt")
      writeDataTable(wb, sheet = "Pval_GlobalAnalysis_Filt", x=as.data.frame(res_f[,c(1,5,6)]), rowNames = TRUE)

    }
  }

  saveWorkbook(wb, paste(colnames(Exp)[3],"_Counts_global.xlsx"), overwrite=TRUE)
  
  
  # ---- DESeq analysis
  cat("\n DESeq analysis - Wald tests")
  
  if (!((analysis == "Ctrl") || (analysis == "All"))){cat("\n Analysis should be one of Ctrl or All")}

  # ---- DESeq analysis each level against control level

  if (analysis == "Ctrl")
    {
      if((is.null(condCtrl)) || (length(grep(condCtrl, colData(dds)[,1]))==0))
        {
          cat("\n Please provide a valid control condition")
        } else {
          cat(paste("\n Comparison of all conditions against", condCtrl))
          
          colData(dds)[,'Condition'] = relevel(colData(dds)[,'Condition'], condCtrl) 
          #print(levels(colData(dds)[,1]))
          
          dds_ana = DESeq(dds, quiet=TRUE)
          dds_ana = dds_ana[rownames(dds_f),]
        
          
          for (k in 2:NbLvl){
            DESeq2_DE_res(Exp, dds_ana, pval, FC, 1, k, correction, record)
            }
        }
  }

  # ---- DESeq analysis each level against each other

  if (analysis == "All")
    {
      if (!is.null(condCtrl))
      {
        cat("\n No control condition needed, all comparisons performed")
      }
      cat("\n Comparison of all conditions two-by-two")
      
      dds = DESeqDataSetFromHTSeqCount(Exp, directory=dir, design=d)
      
      Exp1 <- Exp
      for (i in 1:(NbLvl-1))
        {
          if (i>1)
          {
            Exp1 <- Exp1[which(!(Exp1$Condition == levels(colData(dds)[,1])[1])),]
            Exp1$Condition <- droplevels(Exp1$Condition)

            dds = DESeqDataSetFromHTSeqCount(Exp1, directory=dir, design=d)
          }
        
          j=2
          #colData(dds)[,'Condition'] = relevel(colData(dds)[,'Condition'], levels(colData(dds)[,1])[1]) 
          #print(levels(colData(dds)[,1]))
          while(j<=length(levels(colData(dds)[,1])))
            {
            cat(paste("\n condition", levels(colData(dds)[,1])[j], "vs", levels(colData(dds)[,1])[1], sep="_"))
            
            dds_ana = DESeq(dds, quiet=TRUE)
            dds_ana = dds_ana[rownames(dds_f),]

            DESeq2_DE_res(Exp1, dds_ana, pval, FC, 1, j, correction, record)
            j=j+1
          }
      }
      
      
  }

  gc()
  graphics.off()

}
