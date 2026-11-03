
suppressPackageStartupMessages({library(DESeq2); library(ggplot2)})
dir.create('results', showWarnings=FALSE, recursive=TRUE)
counts <- read.csv('data/counts.csv', check.names=FALSE, stringsAsFactors=FALSE)
stopifnot('gene_id' %in% names(counts), !anyDuplicated(counts$gene_id))
rownames(counts) <- counts$gene_id; counts$gene_id <- NULL
meta <- read.csv('data/sample_metadata.csv', stringsAsFactors=FALSE)
stopifnot(all(c('sample_id','group') %in% names(meta)), !anyDuplicated(meta$sample_id))
stopifnot(setequal(colnames(counts),meta$sample_id))
meta <- meta[match(colnames(counts),meta$sample_id),,drop=FALSE]
stopifnot(all(meta$group %in% c('shNC_TGFb1','shMARCKS_TGFb1')))
stopifnot(all(vapply(counts,is.numeric,logical(1))), all(as.matrix(counts) >= 0), all(as.matrix(counts) == round(as.matrix(counts))))
meta$group <- factor(meta$group, levels=c('shNC_TGFb1','shMARCKS_TGFb1'))
dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(counts)), colData=meta, design=~group)
# Prefilter is a declared example choice: confirm against the actual original analysis.
dds <- dds[rowSums(counts(dds)) >= 10,]
dds <- DESeq(dds)
res <- results(dds, contrast=c('group','shMARCKS_TGFb1','shNC_TGFb1'), alpha=0.05)
df <- as.data.frame(res); df$gene_id <- rownames(df)
df <- df[,c('gene_id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
df$change <- 'NOT_SIGNIFICANT'
df$change[!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange > 1.5] <- 'UP'
df$change[!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange < -1.5] <- 'DOWN'
write.csv(df,'results/DESeq2_all_genes.csv',row.names=FALSE)
write.csv(subset(df,change!='NOT_SIGNIFICANT'),'results/MARCKS_DEG.csv',row.names=FALSE)
write.csv(data.frame(change=c('UP','DOWN'), count=c(sum(df$change=='UP'),sum(df$change=='DOWN'))), 'results/DEG_counts.csv',row.names=FALSE)
p <- ggplot(df, aes(log2FoldChange,-log10(pmax(padj,.Machine$double.xmin)),color=change)) +
  geom_point(alpha=0.65,size=0.8,na.rm=TRUE) +
  geom_vline(xintercept=c(-1.5,1.5),linetype='dashed') +
  geom_hline(yintercept=-log10(0.05),linetype='dashed') +
  labs(x='log2 fold change (shMARCKS / shNC)',y='-log10(BH-adjusted P)',title='MARCKS knockdown differential expression') + theme_classic()
ggsave('results/volcano.pdf',p,width=7,height=5)
saveRDS(dds,'results/deseq2_object.rds')

suppressPackageStartupMessages({library(clusterProfiler);library(org.Mm.eg.db);library(ggplot2);library(enrichplot)})
dir.create('results',showWarnings=FALSE,recursive=TRUE)
f <- 'results/DESeq2_all_genes.csv'
if (!file.exists(f)) stop('Please first create results/DESeq2_all_genes.csv using 01_deseq2.R or your validated complete DESeq2 output.')
x <- read.csv(f,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(all(c('gene_id','padj','log2FoldChange') %in% names(x)))

all_symbols <- unique(x$gene_id[!is.na(x$pvalue)])
degs <- x$gene_id[!is.na(x$padj) & x$padj < .05 & abs(x$log2FoldChange) > 1.5]
map_all <- suppressMessages(bitr(all_symbols, fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Mm.eg.db))
map_deg <- suppressMessages(bitr(unique(degs),fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Mm.eg.db))
universe <- unique(map_all$ENTREZID)
ids <- intersect(unique(map_deg$ENTREZID),universe)
stopifnot(length(ids)>0,length(universe)>length(ids))
ego <- enrichGO(gene=ids, universe=universe, OrgDb=org.Mm.eg.db,
                keyType='ENTREZID',ont='BP',pAdjustMethod='BH',
                pvalueCutoff=1,qvalueCutoff=1,readable=TRUE)
go <- as.data.frame(ego)
write.csv(go,'results/GO_BP_complete.csv',row.names=FALSE)
go_sig <- subset(go,!is.na(p.adjust) & p.adjust < 0.05)
write.csv(go_sig,'results/GO_BP_FDR_lt_0.05.csv',row.names=FALSE)
wnt <- go[grepl('regulation of Wnt signaling pathway',go$Description,ignore.case=TRUE,fixed=TRUE),,drop=FALSE]
write.csv(wnt,'results/GO_WNT_terms.csv',row.names=FALSE)
if(nrow(wnt)>0) {
  genes <- unique(unlist(strsplit(wnt$geneID,'/',fixed=TRUE)))
  write.csv(data.frame(gene_symbol=genes),'results/GO_WNT_contributing_genes.csv',row.names=FALSE)
}
if(nrow(go_sig)>0){
  ego_sig <- ego; ego_sig@result <- go_sig
  p <- dotplot(ego_sig,showCategory=min(15,nrow(go_sig)),color='p.adjust',label_format=60) + ggtitle('GO Biological Process (BH FDR < 0.05)')
  ggsave('results/GO_BP_dotplot.pdf',p,width=10,height=7)
}
saveRDS(ego,'results/GO_BP_enrichResult.rds')

suppressPackageStartupMessages({library(clusterProfiler);library(org.Mm.eg.db)})
x <- read.csv('results/DESeq2_all_genes.csv',stringsAsFactors=FALSE)
stopifnot(all(c('gene_id','pvalue','padj','log2FoldChange') %in% names(x)))
all_sym <- unique(x$gene_id[!is.na(x$pvalue)])
deg_sym <- unique(x$gene_id[!is.na(x$padj)&x$padj<0.05&abs(x$log2FoldChange)>1.5])
bg <- suppressMessages(bitr(all_sym,'SYMBOL','ENTREZID',org.Mm.eg.db))
fg <- suppressMessages(bitr(deg_sym,'SYMBOL','ENTREZID',org.Mm.eg.db))
universe <- unique(bg$ENTREZID)
genes <- intersect(unique(fg$ENTREZID),universe)
kk <- enrichKEGG(gene=genes,universe=universe,organism='mmu',pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1)
write.csv(as.data.frame(kk),'results/KEGG_complete.csv',row.names=FALSE)
