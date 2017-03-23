## Oncopanel CNV analysis

## looks genes most frequently suffering high-level gains or complete loss
all.cnvs <- merge(all.cnvs, master.sheet[, c("Cancer_Type_Broad", "SAMPLE_ACCESSION_NBR")], "SAMPLE_ACCESSION_NBR")
all.cnv.high <- subset(all.cnvs, CNV_TYPE_CD %in% c("2DEL", "HA"))
temp <- all.cnv.high
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
pdf("all high-level CNVs.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:50]),
       aes(x = GENE, fill = CNV_TYPE_CD)) + geom_bar() + rameen_theme
dev.off()


## look just at deletions, colored by tumor type
temp <- subset(all.cnv.high, CNV_TYPE_CD == "2DEL")
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
pdf("all high-level CNVs deletions.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:50]),
       aes(x = GENE, fill = Cancer_Type_Broad)) + geom_bar() + rameen_theme + ggtitle(label = "2 Copy deletions")
dev.off()


## look just at amplifications, colored by tumor type
temp <- subset(all.cnv.high, CNV_TYPE_CD == "HA")
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
temp <- merge(temp, master.sheet[, c("Cancer_Type_Broad", "SAMPLE_ACCESSION_NBR")], "SAMPLE_ACCESSION_NBR")
pdf("all high-level CNVs amplifications.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:50]),
       aes(x = GENE, fill = Cancer_Type_Broad)) + geom_bar() + rameen_theme +ggtitle(label = "High level amplifications")
dev.off()

temp <- all.cnv.high
temp$Cancer_Type_Broad <- factor(temp$Cancer_Type_Broad, levels = names(sort(table(temp$Cancer_Type_Broad), decreasing = T)))
pdf("all high-level CNVs by tumor.pdf", width = 14)
ggplot(data = temp, aes(x = Cancer_Type_Broad, fill = CNV_TYPE_CD)) + geom_bar(position = "dodge") + rameen_theme +
    ggtitle(label = "high level CNVs stratified by tumor type")
dev.off()




## look at low level amplifications and deletions
all.cnv.low <- subset(all.cnvs, CNV_TYPE_CD %in% c("1DEL", "LA"))

temp <- all.cnv.low
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
pdf("all low-level CNVs.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:100]),
       aes(x = GENE, fill = CNV_TYPE_CD)) + geom_bar() + rameen_theme
dev.off()


## look just at deletions, colored by tumor type
temp <- subset(all.cnv.low, CNV_TYPE_CD == "1DEL")
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
pdf("all low-level CNVs deletions.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:50]),
       aes(x = GENE, fill = Cancer_Type_Broad)) + geom_bar() + rameen_theme + ggtitle(label = "1 Copy deletions")
dev.off()


## look just at amplifications, colored by tumor type
temp <- subset(all.cnv.low, CNV_TYPE_CD == "LA")
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
pdf("all low-level CNVs amplifications.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:50]),
       aes(x = GENE, fill = Cancer_Type_Broad)) + geom_bar() + rameen_theme +ggtitle(label = "Low level amplifications")
dev.off()

temp <- all.cnv.low
temp$Cancer_Type_Broad <- factor(temp$Cancer_Type_Broad, levels = names(sort(table(temp$Cancer_Type_Broad), decreasing = T)))
pdf("all low-level CNVs by tumor.pdf", width = 14)
ggplot(data = temp, aes(x = Cancer_Type_Broad, fill = CNV_TYPE_CD)) + geom_bar(position = "dodge") + rameen_theme +
    ggtitle(label = "low level CNVs stratified by tumor type")
dev.off()