cat("###################################################\n\n")
cat("      Welcome to trinity2go                        \n\n")
cat(" v0.4.0                  2016 NTOU-MVIL Lafudoci   \n")
cat("###################################################\n\n")
expFile <- readline("STEP-1/4: Please enter input filename ( e.g. edgeR.DE_results ): \n步驟1/4:請輸入欲分析的完整目標檔案名稱（範例：edgeR.DE_results ）\n")

expDE <- read.delim(expFile)
expDE <- cbind(row.names(expDE),expDE)
row.names(expDE) <- NULL
names(expDE)[1] <- "Name"

cat(c("\n>>>>", expFile, "載入成功..\n\n"))

expName <- readline("STEP-2/4: Please enter experment shortname ( e.g. BRnnv ): \n步驟2/4:請簡短命名實驗名稱，不可包含空格（範例：BRnnv ）\n")

oppsiteSign <- readline("STEP-3/4: Oppsite sign ? ( Y/N ): \n步驟3/4:是否反轉FC sign?（ Y/N ）\n")

pathLv <- readline("STEP-4/4: Please enter KEGG path expression level render scale ( e.g. 2 ): \n步驟4/4:請輸入KEGG表現量級距渲染上限logFC值（範例：2 ）\n")

cat(c("\n>>>>", expName, "工作開始執行..\n\n"))

library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)


#Process clc expression table

refAnnot <- read.csv("~/Documents/2016denovo/trinity_out_dir_all6/uniList_transcripts.txt")
#視需要修改annotation的來源文件

if( oppsiteSign == "Y"){
  expDE[c(2,3)] <- c( -expDE[2] , -expDE[3])
} else if ( oppsiteSign == "y") {
  expDE[c(2,3)] <- c( -expDE[2] , -expDE[3])
} else if ( oppsiteSign == "N") {
} else if ( oppsiteSign == "n") {
} else { stop("錯誤:無法辨識反轉FC sign指令") }

expAnnot <- merge(expDE, refAnnot, by= "Name")
write.csv(expAnnot, paste0(expName,"Annot.csv"), row.names=F)

expTable <- expAnnot[c(1,2,5,6,7)]
expTableFilted <- subset(expTable, expTable$FDR<= 0.05 )
expTableVar <- subset(expTableFilted, abs(expTableFilted$logFC) >= 0.5 )


#write.csv(expTableFilted, paste0(expName,"Table.csv"), row.names=F)
write.csv(expTableVar, paste0(expName,"VarTable.csv"), row.names=F)

spFCMean <-aggregate(cbind(expTableVar[2],expTableVar[3]), by=list(UNIPROT=expTableVar[[4]]), FUN=mean)
write.csv(spFCMean, paste0(expName,"spFCMean.csv"), row.names=F)

Sp2Enz <- bitr(spFCMean[,1], fromType=c("UNIPROT"), toType=c("ENTREZID"), OrgDb ="org.Hs.eg.db")
Sp2Path<- bitr(spFCMean[,1], fromType=c("UNIPROT"), toType=c("PATH"), OrgDb="org.Hs.eg.db")

ezFC <- merge(spFCMean[1:2],Sp2Enz, all.y=T)

uniezFC <- aggregate(ezFC[2], by=list(ENTREZID=ezFC[[3]]), FUN=mean)
uniezFC.t <- uniezFC[2]
row.names(uniezFC.t) <- uniezFC[[1]]
names(uniezFC.t) <-"logFC"

write.csv(uniezFC.t, paste0(expName,"UniezFC.t.csv"))

#GO enrichment
cat("\nGO enrichment...\nGO 富集化分析中...\n\n")

GoEBP <- enrichGO(gene = uniezFC[,1], OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff  = 0.05, readable = T)
GoEMF <- enrichGO(gene = uniezFC[,1], OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff  = 0.05, readable = T)
GoECC <- enrichGO(gene = uniezFC[,1], OrgDb = org.Hs.eg.db, ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff  = 0.05, readable = T)

bpbar <- barplot(GoEBP, showCategory=15, font.size=15)
mfbar <- barplot(GoEMF, showCategory=15, font.size=15)
ccbar <- barplot(GoECC, showCategory=15, font.size=15)

bpdot <- dotplot(GoEBP, showCategory=15, font.size=15)
mfdot <- dotplot(GoEMF, showCategory=15, font.size=15)
ccdot <- dotplot(GoECC, showCategory=15, font.size=15)

write.csv(summary(GoEBP), paste0(expName,"GoEBP.csv"))
write.csv(summary(GoEMF), paste0(expName,"GoEMF.csv"))
write.csv(summary(GoECC), paste0(expName,"GoECC.csv"))

png(paste0(expName,"GoEBP.bar.png"), width=750, height=550)
print(bpbar)
dev.off()
png(paste0(expName,"GoEMF.bar.png"), width=750, height=550)
print(mfbar)
dev.off()
png(paste0(expName,"GoECC.bar.png"), width=750, height=550)
print(ccbar)
dev.off()

png(paste0(expName,"GoEBP.dot.png"), width=600, height=550)
print(bpdot)
dev.off()
png(paste0(expName,"GoEMF.dot.png"), width=600, height=550)
print(mfdot)
dev.off()
png(paste0(expName,"GoECC.dot.png"), width=600, height=550)
print(ccdot)
dev.off()


#KEGG pathway mapping

unipath <- unique(Sp2Path[,2])
pathlist <- unipath[unipath!= "00511" & unipath!= "00514" & unipath!= "01100" & unipath!= "00533" ]
write.csv(pathlist, paste0(expName,"Pathlist.csv"), row.names=F)

uniezFC.m <-data.matrix(uniezFC.t)


dir.create (paste0(expName,"KEGG-supdata"))

pv.out <- pathview(gene.data = uniezFC.m[,1], pathway.id = pathlist, species = "hsa", out.suffix = expName, kegg.native = T, kegg.dir = paste0(expName,"KEGG-supdata"),  limit = as.integer(pathLv))

cat("\nclc2go run finish!\nclc2go執行完成！\n\n")
