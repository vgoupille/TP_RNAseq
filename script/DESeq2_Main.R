####################################################################
#
# Module_6_DESeq2_Main.R
#
# 2018 05 03 DR
#
# objectif : DESeq2 sur les donn�es du Module_6_TP
#   input : counts bruts produits par featureCount.sh sur genouest
#   output : counts, moderated FC, Pvalue, Padj dans fichier .xlsx
#
# DESeq2  version 1.16 nrequired
#
####################################################################


# --- r�pertoire de travail
#setwd("d:/home/drossill/bureau/ENSEIGNEMENT/2018 Formation R/MODULE 6/Module_6_TP/5_DESeq2")
dir <- getwd()



############################################################################
#################### DESeq2 from sample files ##############################
############################################################################


###### chargement des fonctions DESeq2 pour sample files
source("script/DESeq2_FromSampleFiles_v1.16_openxlsx_v1.R")


###### creer la table des informations sur les echantillons
# a faire avec blocnotes ou excel -> .txt
# exemple: sampleinfo_ForSampleFiles.txt


###### cas : toutes les comparaisons possibles
# --- setting parameters - All

table="script/Table_Annotations.csv"

#dir <- getwd()
dir  <- "/Users/valentingoupille/Library/Mobile Documents/com~apple~CloudDocs/University/Master_Bioinfo/GNF/Sibut/TP_RNAseq/data"

pval <- 0.01                # seuil de filtrage sur p (si adjP actif alors p=adjp)
FC <- 2                     # seuil de filtrage sur FC (FC lineaire et non logFC)
analysis = "All"
condCtrl= "DN"              # condition de reference : non prise en compte si analysis = All
PreFilt = 10                # pre-filtrage : conservation des genes ayant au min PreFilt reads au total
Filt = T                    # HTSFilter actif uniquement si NbCountMin = NULL
NbCountMin <- NULL          # min nb counts normalises par gene pour un echantillon 
paired <- F                 # analyse non appariee
correction="bonferroni"     # correction bonferroni pour adjp
record = FALSE              # pas d'enregistrement ds resultats DESeq2 sans filtrage

# --- un fichier log sera genere
log <- paste("log_Module_6_TP_", format(Sys.time(), '%Y%m%d_%Hh%M'), ".txt", sep="") 
sink(log, type="output")
cat(log, "\n")
cat(date(), "\n")



DESeq2_FromSampleFiles(table, dir, pval, FC, analysis=analysis,condCtrl=condCtrl, 
                       PreFilt=PreFilt, Filt=Filt, NbCountMin=NbCountMin, 
                       paired=paired, correction=correction, record=record)
cat("\n")
cat(date(), "\n")
sink()


