#Author: A.T
#Script to process a single-cell 10x dataset with Seurat
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  #cat("\n",args[1])
  stop("Usage: Rscript script <directory(with genenames,.mtx file,barcodes)>
       <name>", call.=FALSE)
}

suppressPackageStartupMessages({
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(Rtsne)
library(cowplot)
})


mainDir = getwd()
subDir = args[2]
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir = paste0(args[2],"/otherplots")
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir = paste0(args[2],"/featureplots")
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir = paste0(args[2],"/savedobjects")
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir = paste0(args[2],"/markergenes")
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)

outputname = args[2]
cat("\nOutputname:",outputname,"\n")
myobj.data = Read10X(data.dir = args[1])
cat("\nDimension:",myobj.data@Dim,"\n")
myobj = CreateSeuratObject(raw.data = myobj.data,
                           project = outputname,display.progress = FALSE,
                           min.genes = 10)
treatment.status = data.frame(treatment = rep(outputname,
                              length(myobj@cell.names)),
                              row.names = myobj@cell.names)
myobj = AddMetaData(object = myobj, metadata = treatment.status)
cat("\nDimension:",dim(myobj@raw.data),"\n")
num.genes = Matrix::colSums(myobj@raw.data > 0)
num.cells = Matrix::colSums(myobj@raw.data)/1e5

#Distribution of number of genes (nGene) in each library having non-zero counts.
cat("\nsummary of non-zero rows:\n",summary(num.genes))

#Distribution of per-million library sizes of cells
cat("\nsumary of lib size:\n",summary(num.cells),"\n")

get_stats = function(obj){
get_stats = function(obj){
reads.per.cell = Matrix::colSums(obj@raw.data)
cat("mean reads per cell:",mean(reads.per.cell))
cat("\nsds per cell:",sd(reads.per.cell))
genes.per.cell = Matrix::colSums(obj@raw.data > 0)
cat("\nmean genes per cell:",mean(genes.per.cell))
cat("\nsd genes per cell:",sd(genes.per.cell))
}
}

cat("\nStats:\n")
get_stats(myobj)
cat("\n=======================\n")

myobj = NormalizeData(object = myobj, normalization.method = "LogNormalize",
                       scale.factor = 1e4,display.progress=FALSE)
myobj = FindVariableGenes(object = myobj,do.plot=FALSE)
cat("\nNo of variable genes:",length(x = myobj@var.genes),"\n")
myobj = ScaleData(object = myobj,display.progress=FALSE)
#First find appropriate total number of  dimension to use
cat("Dimension of scale.data:",ncol(x=myobj@scale.data),"\n")
noc=50
noc=ifelse(noc > ncol(x=myobj@scale.data),40,noc)
noc=ifelse(noc > ncol(x=myobj@scale.data),30,noc)
noc=ifelse(noc > ncol(x=myobj@scale.data),20,noc)
noc=ifelse(noc > ncol(x=myobj@scale.data),10,noc)
noc=ifelse(noc > ncol(x=myobj@scale.data),5,noc)
cat("Dimension to use:",noc,"\n")
myobj = RunPCA(object = myobj, pc.genes = myobj@var.genes, pcs.compute = noc,
                do.print = FALSE)
ggcol=c("#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3",
"#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3",
"#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3",
"#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3",
"#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3",
"#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3",
"#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3",
"#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3", "#984EA3",
"#984EA3", "#984EA3")
DATA  =data.frame(PC=paste0("PC",c(1:noc)),
                  sd=round((myobj@dr$pca@sdev)^2*100/sum((myobj@dr$pca@sdev)^2),2),
                  c=rep("same",noc))
DATA$PC = reorder(DATA$PC,DATA$sd)
ggplot(DATA, aes(x=PC,y=sd,fill = c)) +
       geom_bar(stat="identity",position=position_dodge(0.8)) +
       coord_flip() + scale_fill_manual(values=ggcol) + theme_bw() +
       theme(#panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.border = element_blank(),
       panel.background = element_blank()) +
       ylab(paste0("% Variance in first ",noc," 50 PCs")) +
       xlab("PCS") + theme(legend.position="none")
ggsave(filename =paste0(args[2],"/otherplots/",outputname,
                        "_PC_variance.pdf"),dpi=300)

p12 = PCAPlot(object = myobj, dim.1 = 1, dim.2 = 2,do.return = TRUE,
              no.legend = TRUE)
p13 = PCAPlot(object = myobj, dim.1 = 1, dim.2 = 3,do.return = TRUE,
              no.legend = TRUE)
p23 = PCAPlot(object = myobj, dim.1 = 2, dim.2 = 3,do.return = TRUE,
              no.legend = TRUE)
p14 = PCAPlot(object = myobj, dim.1 = 1, dim.2 = 4,do.return = TRUE,
              no.legend = TRUE)

plot_grid(p12,p13,p14,p23,ncol=2)
ggsave(filename =paste0(args[2],"/otherplots/",outputname,"_PC1-4.pdf"),dpi=300)

cat("\nPlot PCA\n")
myobj = ProjectPCA(object = myobj, do.print = FALSE)

pdf(file=paste0(args[2],"/otherplots/",outputname,"_PCelbowplot.pdf"))
PCElbowPlot(object = myobj,num.pc = noc)
dev.off()

cat("\nFind Cluster\n")
myobj = FindClusters(object = myobj, reduction.type = "pca", dims.use = 1:noc,
                      resolution = 0.6, print.output = 0, save.SNN = TRUE)
myruntsne = function(obj,noc){
    tryCatch({

    tempobj = RunTSNE(object = obj, dims.use = 1:noc, do.fast = TRUE,
                     check_duplicates = FALSE)
    },
     error = function(e){
         tryCatch({
         cat("\nDecreasing perplexity to 20\n")
         tempobj = RunTSNE(object = obj, dims.use = 1:noc, do.fast = TRUE,
                     check_duplicates = FALSE,perplexity=20)
         return(tempobj)

         },
         error = function(e2){
         cat("\nDecreasing perplexity to 1\n")
         tempobj  = RunTSNE(object = obj, dims.use = 1:noc, do.fast = TRUE,
                     check_duplicates = FALSE,perplexity=1)
         return(tempobj)

         })
     }
    )
}
myobj = myruntsne(myobj,noc)
tryCatch({
pdf(paste0(args[2],"/",outputname,"_cluster_dendogram.pdf"))
myobj = BuildClusterTree(myobj, do.reorder = TRUE, reorder.numeric = TRUE)
dev.off()

},
error = function(e){
    cat("\nOnly one cluster!\n")
}
)

pdf(paste0(args[2],"/",outputname,"_tSNE.pdf"))
TSNEPlot(object = myobj,do.label = TRUE, plot.title=outputname)
dev.off()
save(myobj, file = paste0(args[2],"/","savedobjects/",outputname,".Robj"))

cat("\nFind all Markers\n")
myobj.markers = FindAllMarkers(object = myobj, min.pct = 0.25,
                                return.thresh = 0.05, test.use = "wilcox",
                                print.bar = FALSE)
save(myobj.markers, file = paste0(args[2],"/","savedobjects/",
                                  outputname,".markers.Robj"))

plotfeatureplots = function(obj,df, status,cluster, no.genes=20){
  no.genes = ifelse(nrow(df)>10,no.genes,nrow(df))
  df = head(df,n=no.genes)
  for( genename in df$gene){
    #print(genename)
    if(status=="up"){
    pdf(file = paste0(args[2],"/","featureplots/cluster",cluster,"/",
                      genename,"_UP.pdf"))
    FeaturePlot(object = obj, features.plot = genename,
                cols.use =c("lightgray","red"), reduction.use = "tsne",
                no.legend = FALSE)
    dev.off()
    }else{
    pdf(file = paste0(args[2],"/","featureplots/cluster",cluster,"/",
                      genename,"_DOWN.pdf"))
    FeaturePlot(object = obj, features.plot = genename,
                cols.use =c("lightgray","red"), reduction.use = "tsne",
                no.legend = FALSE)
    dev.off()
    }

  }
}

cat("\nDE\n")
for(cluster in unique(myobj.markers$cluster)){
  df = myobj.markers[myobj.markers$cluster==cluster,]
  df = df[order(df$p_val_adj),]
  df = df[,-c(3,4)]
  df = df[,c(2,1,3,4,5)]
  df.up = df[df$avg_logFC > 0,]
  df.down = df[df$avg_logFC < 0,]
  df.significant.up = df[df$p_val_adj < 0.01 & df$avg_logFC > 0,]
  df.significant.up = df.significant.up[order(-df.significant.up$avg_logFC),]
  df.significant.down = df[df$p_val_adj < 0.01 & df$avg_logFC < 0,]
  df.significant.down = df.significant.down[order(-df.significant.down$avg_logFC),]
  mainDir = getwd()
  subDir = paste0(args[2],"/","featureplots/cluster",cluster)
  dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
  plotfeatureplots(myobj, df.significant.up,"up",cluster,20)
  plotfeatureplots(myobj, df.significant.down,"down",cluster,20)
  write.table(df.up,file=paste0(args[2],"/","markergenes/",outputname,
              "_cluster",cluster,"_UP.tsv"),sep="\t", quote=F,col.names = NA)
  write.table(df.down,file=paste0(args[2],"/","markergenes/",
              outputname,"_cluster",cluster,"_DOWN.tsv"),sep="\t",
              quote=F,col.names = NA)

  myobj.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
# setting slim.col.label to TRUE will print
# just the cluster IDS instead of every cell name
pdf(paste0(args[2],"/otherplots/",outputname,
           "_heatmapbasedon_top10markergenes.pdf"))
DoHeatmap(object = myobj, genes.use = top10$gene, slim.col.label = TRUE,
          cex.row = 5,rotate.key = TRUE,group.cex = 8,group.spacing = 0.3)
dev.off()
print("\nSuccess!!!\n")
}
