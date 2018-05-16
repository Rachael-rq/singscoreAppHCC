##------------ Performance evaluation 

## score TGFb data using 5 different methods against MES an EPI gene signature

library(GSEABase)
library(GSVA)
##----- read Epi and Mes gene signatures
epi_mes = read.table("./data/Thiery_EMTsignature_both_tumour_cellLine_EntrezIDs.txt", header=T, sep = "\t")
epi_cll = GeneSet( as.character(na.omit(epi_mes$EntrezGene.ID[epi_mes$epiMes_cellLine == "epi"])))
mes_cll = GeneSet( as.character(na.omit(epi_mes$EntrezGene.ID[epi_mes$epiMes_cellLine == "mes"])))

##----- read comBat_corrected_tfgb_dataset

tgfb_expr = read.table("./data/comBat_corrected_Foroutan_et_al_2017.txt", header=T, sep="\t")
rn = tgfb_expr$EntrezID
tgfb_expr = as.matrix( tgfb_expr[, 2:ncol(tgfb_expr)] )

#rownames as ENTREZID
row.names(tgfb_expr) = rn

tgfbAnnot = data.frame(SampleID = colnames(tgfb_expr),
                       Batch = sapply(colnames(tgfb_expr), function(x) unlist(strsplit(x, "[_]"))[1]), 
                       Type = NA)

tgfbAnnot$Type[grepl("Ctrl", tgfbAnnot$SampleID)] = "Control"
tgfbAnnot$Type[grepl("TGFb", tgfbAnnot$SampleID)] = "TGFb"




##==================================== Scoring tgfb data against EPI =======================================
#install.packages("~/Documents/davisLab/singscore_0.99.9.tar.gz", repos = NULL, type = "source")
library(singscore)

# TGFB samples
# tgfbSamples = tgfb_expr[,tgfbAnnot$Type == "TGFb"]
# dim(tgfbSamples)
# 28 treated samples
##----------------- simple scoring 
tgfb_expr_ranked = rankGenes(expr = tgfb_expr)
cor_simple = simpleScore(rankData = tgfb_expr_ranked, upSet = epi_cll)
cor_simple


##------------------ ssGSEA

cor_ssgsea = gsva(as.matrix(tgfb_expr), list(epi = geneIds(epi_cll)),  
                  method="ssgsea",mx.diff = TRUE,
                  verbose = T, abs.ranking = F, kcdf = "Gaussian") 
cor_ssgsea

##------------------ GSVA


cor_gsva = gsva(as.matrix(tgfb_expr), list(epi = geneIds(epi_cll)),  
                method="gsva", mx.diff = TRUE,
                verbose = T, abs.ranking = F, kcdf = "Gaussian")

cor_gsva


##------------------ z-score



cor_zscore = gsva(as.matrix(tgfb_expr), list(epi = geneIds(epi_cll)), 
                  method="zscore", mx.diff = TRUE,
                  verbose = T, abs.ranking = F, kcdf = "Gaussian")
cor_zscore

##------------------ PLAGE

cor_plage = gsva(as.matrix(tgfb_expr), list(epi = geneIds(epi_cll)), 
                 method="plage", mx.diff = TRUE,
                 verbose = T, abs.ranking = F, kcdf = "Gaussian")
#cor_plage = -cor_plage

cor_plage

##----------------- merge into one data.frame 
allMethScoresTGFb_epi = data.frame(rbind(cor_simple$TotalScore, 
                                         cor_ssgsea, 
                                         cor_gsva, 
                                         cor_zscore, 
                                         cor_plage))
rownames(allMethScoresTGFb_epi) = c("simpleScore", "ssgsea", "gsva", "zscore",
                                    "plage")


control = grepl("Ctrl", colnames(allMethScoresTGFb_epi))
treated = grepl("TGFb", colnames(allMethScoresTGFb_epi))

results = t(allMethScoresTGFb_epi)

test_result_epi = apply(results,2,function(x){
  #browser()
  ctrl = x[control]
  treatment = x[treated]
  data.frame("pval" = t.test(ctrl,treatment, alternative = "two.sided")$p.value,
             "t-statistic" = t.test(ctrl,treatment, alternative = "two.sided")$statistic)
})

test_result_epi


# ------------- plotting results EPI for tgfb control and cases
library(ggplot2)
results = as.data.frame(results)
results$Annotation[control] = "control"
results$Annotation[treated] = "TGFb"
results$Samples = rownames(results)
rnames = colnames(allMethScoresTGFb_cor)


tras = sapply(1:5, function(x){
  return(data.frame("method" = colnames(results)[x], "samples" = rnames, "scores_epi" = results[,x]))
}, simplify = FALSE)

library(plyr)
df <- ldply(tras, data.frame)

df$group = "group"
df$group[grepl("Ctrl", df$samples)] = "Ctrl"
df$group[grepl("TGFb", df$samples)] = "TGFb"


tgfb_plot_epi = ggplot(data = df)+
  geom_histogram(aes(x = scores_epi, fill = group), bins = 20)+
  facet_grid(~ method, scales="free")
tgfb_plot_epi


df_test_result = ldply(test_result_epi, data.frame)
colnames(df_test_result) = c("methods", "pval", "t.statistic")

# df for t test results
df_test_result_epi = data.frame(df_test_result, "GS" = "epi")

##==================================== Scoring MES GS=======================================
##----------------- simple scoring 


cor_simple = simpleScore(rankData = tgfb_expr_ranked, upSet = mes_cll )
cor_simple


##------------------ ssGSEA

cor_ssgsea = gsva(as.matrix(tgfb_expr), list(epi = geneIds(mes_cll)),  
                  method="ssgsea",mx.diff = TRUE,
                  verbose = T, abs.ranking = F, kcdf = "Gaussian") 
cor_ssgsea

##------------------ GSVA


cor_gsva = gsva(as.matrix(tgfb_expr), list(epi = geneIds(mes_cll)),  
                method="gsva", mx.diff = TRUE,
                verbose = T, abs.ranking = F, kcdf = "Gaussian")

cor_gsva


##------------------ z-score



cor_zscore = gsva(as.matrix(tgfb_expr), list(epi = geneIds(mes_cll)), 
                  method="zscore", mx.diff = TRUE,
                  verbose = T, abs.ranking = F, kcdf = "Gaussian")
cor_zscore

##------------------ PLAGE

cor_plage = gsva(as.matrix(tgfb_expr), list(epi = geneIds(mes_cll)), 
                 method="plage", mx.diff = TRUE,
                 verbose = T, abs.ranking = F, kcdf = "Gaussian")
#cor_plage = -cor_plage

cor_plage
##------------------ Combine scoring results for GS mes
allMethScoresTGFb_mes = data.frame(rbind(cor_simple$TotalScore, 
                                         cor_ssgsea, 
                                         cor_gsva, 
                                         cor_zscore, 
                                         cor_plage))
rownames(allMethScoresTGFb_mes) = c("simpleScore", "ssgsea", "gsva", "zscore",
                                    "plage")


control = grepl("Ctrl", colnames(allMethScoresTGFb_mes))
treated = grepl("TGFb", colnames(allMethScoresTGFb_mes))

results = t(allMethScoresTGFb_mes)

test_result_mes = apply(results,2,function(x){
  #browser()
  ctrl = x[control]
  treatment = x[treated]
  data.frame("pval" = t.test(ctrl,treatment, alternative = "two.sided")$p.value,
             "t-statistic" = t.test(ctrl,treatment, alternative = "two.sided")$statistic)
})

test_result_mes


# ------------- plotting results 
library(ggplot2)
results = as.data.frame(results)
results$Annotation[control] = "control"
results$Annotation[treated] = "TGFb"
results$Samples = rownames(results)
rnames = colnames(allMethScoresTGFb_mes)


tras = sapply(1:5, function(x){
  return(data.frame("method" = colnames(results)[x], "samples" = rnames, "scores_mes" = results[,x]))
}, simplify = FALSE)

library(plyr)
df <- ldply(tras, data.frame)

df$group = "group"
df$group[grepl("Ctrl", df$samples)] = "Ctrl"
df$group[grepl("TGFb", df$samples)] = "TGFb"


tgfb_plot_mes = ggplot(data = df)+
  geom_histogram(aes(x = scores_mes, fill = group), bins = 20)+
  facet_grid(~ method, scales="free")

tgfb_plot_mes

df_test_result = ldply(test_result, data.frame)
colnames(df_test_result) = c("methods", "pval", "t.statistic")

df_test_result_mes = data.frame(df_test_result, "GS" = "mes")

##-------------- p value / t test statistic / plot : control vs. case

df_test_result_mes[df_test_result_mes$methods == "plage","t.statistic"] = -df_test_result_mes[df_test_result_mes$methods == "plage","t.statistic"]
#df_test_result_epi[df_test_result_epi$methods == "plage","t.statistic"] = -df_test_result_epi[df_test_result_epi$methods == "plage","t.statistic"]

ggplot(data = rbind(df_test_result_mes,df_test_result_epi),
       aes(x = methods, y = pval, group = GS, color = GS))+
  geom_point()+
  geom_line()

ggplot(data = rbind(df_test_result_mes,df_test_result_epi),aes(x = methods, y = t.statistic, group = GS, color = GS))+
  geom_point()+
  geom_line()


#---------- Score new samples and see what the new samples lay

# Download new samples from GEO

library(GEOquery)
# Note that GSEMatrix=TRUE is the default
GSE79235 <- getGEO('GSE79235',GSEMatrix=TRUE)
show(GSE79235)

# Get the 4 samples we need

eset = GSE79235[[1]]

summary(exprs(eset)[,1:5])

boxplot(log2(exprs(eset)))

pData(eset)

## subset
subset <- eset[,1:4]
dim(exprs(subset))
dim(pData(subset))
dim(fData(subset))

#------------ ge matrix
# for two samples 
# new_ge = exprs(subset)[,c(1,3)]
new_ge = exprs(subset)


library(RColorBrewer)
nsamples <- ncol(new_ge)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(log2(new_ge[,1])), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm") 
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(log2(new_ge[,i]))
  lines(den$x, den$y, col=col[i], lwd=2)
}

## Filter out lowly expressed genes
log2_trans = log2(new_ge)
# thres = 5
Homo_sapiens <- read.delim("~/Documents/davisLab/HCC/Homo_sapiens.gene_info", quote="")

genes_in_tgfb = merge(data.frame("Symbol"= Homo_sapiens$Symbol, "GeneID" = Homo_sapiens$GeneID),
                      data.frame("GeneID" = rownames(tgfb_expr)), by.x ="GeneID", by.y = "GeneID")

sum(duplicated(genes_in_tgfb$Symbol))
gene_annot = fData(subset)
symbs = gene_annot$GENE_SYMBOL
symbs[symbs == ""] = "None"
rownames(log2_trans) = symbs

# remove none rows

log2_trans = log2_trans[!rownames(log2_trans) == "None",]

# find the duplicates
dup_gene = rownames(log2_trans)[duplicated(rownames(log2_trans))]
dup_gene = unique(dup_gene)

dup_log2 = log2_trans[rownames(log2_trans) %in% dup_gene,]

dups = data.frame(dup_log2, genes = rownames(dup_log2))

unqs = unique(dups$genes)

median_i = sapply(unqs, function(x){
  temp = dups[dups$genes == x,]
  re = apply(temp[,c(1:4)],2,median)
  re
})

median_i = t(median_i)
rownames(median_i) = unqs

dim(median_i)
# remove the duplicate, and merge the median for duplicate

re_dup_log2 = log2_trans[!(rownames(log2_trans) %in% dup_gene),]
dim(re_dup_log2)
# 

cleaned_log2 = rbind(re_dup_log2,median_i)
sum(duplicated(rownames(cleaned_log2)))
dim(cleaned_log2)

# keep the genes in tgfb gene set

filtered_expr <- cleaned_log2[rownames(cleaned_log2) %in% genes_in_tgfb$Symbol,]
dim(filtered_expr)

plot(density(filtered_expr[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")

# sum(!gene_annot$GENE_SYMBOL == rownames(filtered_expr))


# Rename the rows with gene symbol
# rownames(filtered_expr) = gene_annot$GENE_SYMBOL
# Remove those rows without symbols
# filtered_expr = filtered_expr[!rownames(filtered_expr) == "",]

#remove duplicated rows 
# dup = duplicated(rownames(filtered_expr))
# 
# mat <- match(unique(rownames(filtered_expr)), rownames(filtered_expr))

# remove duplication 
#filtered_expr = filtered_expr[mat,]

dim(filtered_expr)

sum(duplicated(rownames(filtered_expr)))


# mapping = merge(data.frame("Symbol" = rownames(filtered_expr)), 
#                 data.frame("Symbol" = Homo_sapiens$Symbol, 
#       "EntrezID" = Homo_sapiens$GeneID), sort = FALSE)

# uniq = match(unique(mapping$Symbol), mapping$Symbol)
# mapping = mapping[uniq,]

# sum(duplicated(mapping$Symbol))

# filtered_expr = filtered_expr[mapping$Symbol,]

# filtered_expr = data.frame(Symbol = rownames(filtered_expr), filtered_expr)

# filtered_expr = merge(filtered_expr, mapping, by.x = "Symbol", by.y = "Symbol", 
#                       sort = FALSE)
# rownames(filtered_expr) = filtered_expr$EntrezID

# filtered_expr = filtered_expr[, c(2:5)]


#---------- Gene signatures
# The EPI and MES gene signatures represented by gene symbol
epi_cll_syb = GeneSet( unique(as.character(na.omit(epi_mes$officialSymbol[epi_mes$epiMes_cellLine == "epi"]))))
mes_cll_syb = GeneSet( unique(as.character(na.omit(epi_mes$officialSymbol[epi_mes$epiMes_cellLine == "mes"]))))




#----------------- simple score 4 new samples mes
filter_rank = rankGenes(filtered_expr)

new_simple_mes = simpleScore(rankData = filter_rank, upSet = mes_cll_syb )

new_simple_mes


##------------------ ssGSEA

new_ssgsea_mes = gsva(as.matrix(filter_rank), list(epi = geneIds(mes_cll_syb)),  
                      method="ssgsea",mx.diff = TRUE,
                      verbose = T, abs.ranking = F, kcdf = "Gaussian") 
new_ssgsea_mes

##------------------ GSVA


new_gsva_mes = gsva(as.matrix(filter_rank), list(epi = geneIds(mes_cll_syb)),  
                    method="gsva", mx.diff = TRUE,
                    verbose = T, abs.ranking = F, kcdf = "Gaussian")

new_gsva_mes


##------------------ z-score



new_zscore_mes = gsva(as.matrix(filter_rank), list(epi = geneIds(mes_cll_syb)), 
                      method="zscore", mx.diff = TRUE,
                      verbose = T, abs.ranking = F, kcdf = "Gaussian")
new_zscore_mes

##------------------ PLAGE

new_plage_mes = gsva(as.matrix(filter_rank), list(epi = geneIds(mes_cll_syb)), 
                     method="plage", mx.diff = TRUE,
                     verbose = T, abs.ranking = F, kcdf = "Gaussian")
#new_plage_mes = -new_plage_mes

new_plage_mes


## ------------ combine new samples' scores mes

new_data_df = rbind(new_simple_mes$TotalScore, new_ssgsea_mes, new_gsva_mes, new_zscore_mes, new_plage_mes)

rownames(new_data_df) = c("simpleScore", "ssgsea", "gsva", "zscore", "plage")
temp = sapply(1:5, function(x){
  return(data.frame("method" = rownames(new_data_df)[x], 
                    "samples" = colnames(new_data_df),
                    "sampleType" = c("Ctrl","Ctrl","TGFb","TGFb"),
                    "scores_mes" = new_data_df[x,]))
}, simplify = FALSE)

new_data_df_mes = ldply(temp,data.frame)

tgfb_plot_mes+
  geom_vline(data = new_data_df_mes, aes(xintercept=scores_mes, linetype = samples, color = sampleType))


## new sample scores EPI

#----------------- simple score 4 new samples EPI

new_simple_epi = simpleScore(rankData = filter_rank, upSet = epi_cll_syb )

new_simple_epi


##------------------ ssGSEA

new_ssgsea_epi = gsva(as.matrix(filter_rank), list(epi = geneIds(epi_cll_syb)),  
                      method="ssgsea",mx.diff = TRUE,
                      verbose = T, abs.ranking = F, kcdf = "Gaussian") 
new_ssgsea_epi

##------------------ GSVA


new_gsva_epi = gsva(as.matrix(filter_rank), list(epi = geneIds(epi_cll_syb)),  
                    method="gsva", mx.diff = TRUE,
                    verbose = T, abs.ranking = F, kcdf = "Gaussian")

new_gsva_epi


##------------------ z-score


new_zscore_epi = gsva(as.matrix(filter_rank), list(epi = geneIds(epi_cll_syb)), 
                      method="zscore", mx.diff = TRUE,
                      verbose = T, abs.ranking = F, kcdf = "Gaussian")
new_zscore_epi

##------------------ PLAGE

new_plage_epi = gsva(as.matrix(filter_rank), list(epi = geneIds(epi_cll_syb)), 
                     method="plage", mx.diff = TRUE,
                     verbose = T, abs.ranking = F, kcdf = "Gaussian")
#new_plage_epi = -new_plage_epi

new_plage_epi


## ------------ combine new samples' scores epi


new_data_df = rbind(new_simple_epi$TotalScore, new_ssgsea_epi, new_gsva_epi, new_zscore_epi, new_plage_epi)

rownames(new_data_df) = c("simpleScore", "ssgsea", "gsva", "zscore", "plage")


temp = sapply(1:5, function(x){
  return(data.frame("method" = rownames(new_data_df)[x], 
                    "samples" = colnames(new_data_df),
                    "sampleType" = c("Ctrl","Ctrl","TGFb","TGFb"),
                    "scores_epi" = new_data_df[x,]))
}, simplify = FALSE)

new_data_df_epi = ldply(temp,data.frame)

tgfb_plot_epi+
  geom_vline(data = new_data_df_epi, aes(xintercept=scores_epi, linetype = samples, color = sampleType))


