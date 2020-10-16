library(ggplot2)
library(edgeR)
library(dplyr)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(data.table)
library(limma)
library(Biobase)
library(DESeq2)
library(biomaRt)
font_import()

theme_science <- function (base_size = 12, base_family = "Arial Black") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour= "black", size=1),  axis.line.y = element_line(colour= "black", size=1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}

#this function stacks ggplot graphs, downloaded from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#used with scale_axis_continuous to only show two decimal points, downloaded from https://code-examples.net/en/q/24eda9a
scaleFUN <- function(x) sprintf("%.2f", x)

#color blind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("/Users/cottr/Box Sync/BRCA/Essential Genes/")

#cell line annotations from Marcotte et al., 2016, https://github.com/neellab/bfg/tree/gh-pages/data/annotations
cell_line_subtypes <- read.delim("cell_line_subtypes.txt")

#read in DEMETER2 scores from DepMap portal, get just ADAR scores
depmap <- fread("D2_combined_gene_dep_scores.csv")
depmap <- subset(depmap, V1 == "ADAR (103)")
#make row.names = ADAR (103)
depmap <- data.frame(depmap, row.names = depmap$V1)
#get rid of column 1
depmap$V1 <- NULL
#transpose, rows are now cells and the one column is ADAR DEMETER2
depmap <- as.data.frame(t(depmap))
#change column name
colnames(depmap) <- "DEMETER2_ADAR"

#split row.names by '_' and make new columns with cell_line and Site
out <- strsplit(as.character(row.names(depmap)), "_", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

depmap$cell_line <- tolower(out$V1)
depmap$Site <- out$V2


#Choose only Breast cell lines and get rid of Site
depmap <- subset(depmap, Site == "BREAST")
depmap$Site <- NULL

#read in CERES data and annotations for mapping
CERES <- fread("Achilles_gene_effect.csv")
Cell_annotations <- fread("cell_lines_annotations_20181226.txt")

#New data.frame for CERES with just ADAR scores and cell lines, merge with annotations
CERES <- data.frame(CERES = CERES$`ADAR (103)`, cell_line = CERES$V1)

CERES <- merge(CERES, Cell_annotations, by.x = "cell_line", by.y = "depMapID")

#split row.names by '_' and make new columns with cell_line and Site
out <- strsplit(as.character(CERES$CCLE_ID), "_", fixed = TRUE)
out <- do.call(rbind, out)
out <- as.data.frame(out)

CERES <- data.frame(CERES = as.numeric(as.character(CERES$CERES)), cell_line = tolower(out$V1), Site = out$V2)

#Choose only Breast cell lines
CERES <- subset(CERES, Site == "BREAST")

#merge demeter2 and ceres scores
Dependency <- merge(depmap, CERES, by = "cell_line", all = TRUE)

#merge with cell line annotations
Dependency <- merge(cell_line_subtypes, Dependency, by = "cell_line", all = TRUE)

#Make plots of interest
ggplot(subset(Dependency, !is.na(DEMETER2_ADAR) & !is.na(subtype_three_receptor)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = subtype_three_receptor)) + 
  geom_bar(stat = "identity") + scale_x_discrete(position = "top") +
  scale_fill_manual(values = cbPalette) + theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "Subtype")
ggsave("ADAR_BCCL_Depmap.tiff")

ggplot(subset(Dependency, !is.na(CERES) & !is.na(subtype_three_receptor)), aes(reorder(cell_line, -CERES), CERES, fill = subtype_three_receptor)) + 
  geom_bar(stat = "identity") + scale_x_discrete(position = "top") +
  scale_fill_manual(values = cbPalette) + theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "Subtype")
ggsave("ADAR_BCCL_CERES.tiff")


ggplot(subset(Dependency, !is.na(DEMETER2_ADAR) & !is.na(subtype_neve)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = subtype_neve)) + 
  geom_bar(stat = "identity") + scale_x_discrete(position = "top") +
  scale_fill_manual(values = cbPalette) + theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "Subtype")
ggsave("ADAR_BCCL_Depmap_neve.tiff")

ggplot(subset(Dependency, !is.na(CERES) & !is.na(subtype_neve)), aes(reorder(cell_line, -CERES), CERES, fill = subtype_neve)) + 
  geom_bar(stat = "identity") + scale_x_discrete(position = "top") +
  scale_fill_manual(values = cbPalette) + theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "Subtype")
ggsave("ADAR_BCCL_CERES_neve.tiff")


#read in CCLE rnaseq data
CCLE <- fread("CCLE_RNAseq_genes_counts_20180929.gct")

#make a data.frame of identities for gene names and ensembl gene id
CCLE_idents <- data.frame(Name = CCLE$Name, Gene = CCLE$Description)

#make row.names ensembl ids
CCLE <- data.frame(CCLE, row.names = CCLE$Name)
#keep only breast cell lines
CCLE <- dplyr::select(CCLE, contains("BREAST"))

#split cell_line identifiers by '_' and use to name columns
out <- strsplit(as.character(colnames(CCLE)), "_", fixed = TRUE)
out <- as.data.frame(do.call(rbind, out))

#now assin the column names
colnames(CCLE) <- as.character(out$V1)

#remove genes with no counts in all samples
CCLE <- CCLE[rowSums(CCLE == 0) != ncol(CCLE),]

#normalize
logCPM <- cpm(CCLE, prior.count=2, log=TRUE)

#calculate z-scores 
logCPM <- scale(t(logCPM))

#transpose 
CCLE_Z <- as.data.frame(t(logCPM))

#make a column with ensembl ids
CCLE_Z$ensembl <- row.names(CCLE_Z)

#merge with ccle_idents by ensembl ids
CCLE_Z <- merge(CCLE_Z, CCLE_idents, by.x = "ensembl", by.y = "Name")

#subset to remove duplicated gene names
CCLE_Z <- subset(CCLE_Z, !duplicated(Gene))

#Assign row names based on gene name
row.names(CCLE_Z) <- CCLE_Z$Gene

#remove the gene name and ensembl id
CCLE_Z$ensembl <- NULL
CCLE_Z$Gene <- NULL

#transpose
CCLE_Z <- as.data.frame(t(as.matrix(CCLE_Z)))

#make a new column with cell line names
CCLE_Z$Cell_Line <- tolower(row.names(CCLE_Z))


#merge with dependency data
Dependency_CCLE <- merge(Dependency, CCLE_Z, by.x = "cell_line", by.y = "Cell_Line")

#make row.names cell_line names
row.names(Dependency_CCLE) <- Dependency_CCLE$cell_line

#read core_isgs 
Core_ISGs <- read.delim("Core_ISGs.txt")

ISGS <- as.character(Core_ISGs$ISGs)

#get core_isg_score, select only core_isgs from rnaseq data
CCLE_Z_ISG_Core <- dplyr::select(Dependency_CCLE, one_of(ISGS)) 

#get median z-score
Core_ISG_scores <- data.frame(Core_ISG_median = rowMedians(as.matrix(CCLE_Z_ISG_Core)), Cell_Line = Dependency_CCLE$cell_line)

#merge core isg scores with dependency data
Dependency_CCLE <- merge(Dependency_CCLE, Core_ISG_scores, by.x = "cell_line", by.y = "Cell_Line")


#make plots of interest
ggplot(subset(Dependency_CCLE, !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = Core_ISG_median)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "Core ISG Score")
ggsave("ADAR_BCCL_Depmap_ISG.tiff")

ggplot(subset(Dependency_CCLE, !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = Core_ISG_median)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "Core ISG Score")
ggsave("ADAR_BCCL_CERES_ISG.tiff")


ggplot(subset(Dependency_CCLE, subtype_three_receptor == "TNBC"), aes(DEMETER2_ADAR, Core_ISG_median)) + geom_point() + theme_science() + 
  geom_smooth(method = "lm", colour = "#56B4E9") + labs(y = "Core ISG Score", x="ADAR DEMETER2 Score", colour = "Subtype") + 
  geom_hline(yintercept = 0, colour = "grey") + geom_vline(xintercept = -0.5, colour = "grey")
ggsave("Depmap_ISG_point_TNBC.tiff", width = 4, height = 4)

cor.test(Dependency_CCLE$DEMETER2_ADAR[Dependency_CCLE$subtype_three_receptor == "TNBC"], Dependency_CCLE$Core_ISG_median[Dependency_CCLE$subtype_three_receptor == "TNBC"])

ggplot(Dependency_CCLE, aes(DEMETER2_ADAR, Core_ISG_median)) + geom_point() + theme_science() + geom_smooth(method = "lm") + 
  geom_hline(yintercept = 0.20, colour = "grey") + geom_vline(xintercept = -0.5, colour = "grey") +labs(y = "Core ISG Score", x="ADAR DEMETER2 Score")
ggsave("DEMETER2_ISG_point.tiff")

cor.test(Dependency_CCLE$DEMETER2_ADAR, Dependency_CCLE$Core_ISG_median)


ggplot(Dependency_CCLE, aes(CERES, Core_ISG_median)) + geom_point() + theme_science() + geom_smooth(method = "lm") + geom_hline(yintercept = 0.25)
ggsave("CERES_ISG_point.tiff")

cor.test(Dependency_CCLE$CERES, Dependency_CCLE$Core_ISG_median)


ggplot(subset(Dependency_CCLE, !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = ADAR)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "ADAR z-score")
ggsave("ADAR_BCCL_Depmap_ADAR.tiff")

ggplot(subset(Dependency_CCLE, !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = ADAR)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "ADAR z-score")
ggsave("ADAR_BCCL_CERES_ADAR.tiff")

ggplot(Dependency_CCLE, aes(DEMETER2_ADAR, ADAR)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("Depmap_ADAR_point.tiff")

cor.test(Dependency_CCLE$DEMETER2_ADAR, Dependency_CCLE$ADAR)

ggplot(Dependency_CCLE, aes(CERES, ADAR)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("CERES_ADAR_point.tiff")

cor.test(Dependency_CCLE$CERES, Dependency_CCLE$ADAR)


ggplot(subset(Dependency_CCLE, !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = EIF2AK2)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "PKR z-score")
ggsave("ADAR_BCCL_Depmap_PKR.tiff")

ggplot(subset(Dependency_CCLE, !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = EIF2AK2)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "PKR z-score")
ggsave("ADAR_BCCL_CERES_PKR.tiff")


ggplot(Dependency_CCLE, aes(DEMETER2_ADAR, EIF2AK2)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("Depmap_EIF2AK2_point.tiff", height = 4, width = 4, units = "in")

cor.test(Dependency_CCLE$DEMETER2_ADAR, Dependency_CCLE$EIF2AK2)

ggplot(Dependency_CCLE, aes(CERES, EIF2AK2)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("CERES_EIF2AK2_point.tiff", height = 4, width = 4, units = "in")

cor.test(Dependency_CCLE$CERES, Dependency_CCLE$EIF2AK2)


ggplot(subset(Dependency_CCLE, !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = ATF4)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "ATF4 z-score")
ggsave("ADAR_BCCL_Depmap_ATF4.tiff")

ggplot(subset(Dependency_CCLE, !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = ATF4)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "ATF4 z-score")
ggsave("ADAR_BCCL_CERES_ATF4.tiff")

ggplot(Dependency_CCLE, aes(DEMETER2_ADAR, ATF4)) + geom_point() + theme_science() + geom_smooth(method = "lm") + labs(x = "ADAR DEMETER2 Score", y = "ATF4 z-score")
ggsave("Depmap_ATF4_point.tiff")

cor.test(Dependency_CCLE$DEMETER2_ADAR, Dependency_CCLE$ATF4)

ggplot(Dependency_CCLE, aes(CERES, ATF4)) + geom_point() + theme_science() + geom_smooth(method = "lm") + labs(x = "ADAR CERES Score", y = "ATF4 z-score")
ggsave("CERES_ATF4_point.tiff")

cor.test(Dependency_CCLE$CERES, Dependency_CCLE$ATF4)


plot2 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = subtype_three_receptor)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = cbPalette) + theme_science(base_size = 16)  +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + labs(x = "", y = "ADAR DEMETER2 Score", fill = "Subtype")  + 
  theme(legend.position="bottom", legend.title = element_blank())
plot1 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), EIF2AK2)) + 
  geom_bar(stat = "identity", fill = "grey") + scale_fill_manual(values = cbPalette) + theme_science(base_size = 16) + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  xlab("BRCA Cell Lines") + ylab("PKR Z-Score") + scale_y_continuous(labels=scaleFUN)
multiplot(plot1, plot2, cols = 1)

tiff("PKR_Type_ADAR_DEMETER2.tiff", height = 8, width = 7, units = "in", res = 300)
multiplot(plot1, plot2, cols = 1)
dev.off()

plot2 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = subtype_three_receptor)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = cbPalette) + theme_science(base_size = 16)  +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + labs(x = "", y = "ADAR CERES Score", fill = "Subtype")  + 
  theme(legend.position="bottom", legend.title = element_blank())
plot1 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(CERES)), aes(reorder(cell_line, -CERES), EIF2AK2)) + 
  geom_bar(stat = "identity", fill = "grey") + scale_fill_manual(values = cbPalette) + theme_science(base_size = 16) + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  xlab("BRCA Cell Lines") + ylab("PKR Z-Score") + scale_y_continuous(labels=scaleFUN)
multiplot(plot1, plot2, cols = 1)

tiff("PKR_Type_ADAR_CERES.tiff", height = 8, width = 7, units = "in", res = 300)
multiplot(plot1, plot2, cols = 1)
dev.off()


plot2 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = subtype_three_receptor)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = cbPalette) + theme_science(base_size = 16)  +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + labs(x = "", y = "ADAR DEMETER2 Score", fill = "Subtype")  + 
  theme(legend.position="bottom", legend.title = element_blank())
plot1 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), Core_ISG_median)) + 
  geom_bar(stat = "identity", fill = "grey") + scale_fill_manual(values = cbPalette) + theme_science(base_size = 16) + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  xlab("BRCA Cell Lines") + ylab("Core ISG Score") + scale_y_continuous(labels=scaleFUN)
multiplot(plot1, plot2, cols = 1)

tiff("ISG_Type_ADAR_DEMETER2.tiff", height = 8, width = 7, units = "in", res = 300)
multiplot(plot1, plot2, cols = 1)
dev.off()

plot2 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = subtype_three_receptor)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = cbPalette) + theme_science(base_size = 16)  +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + labs(x = "", y = "ADAR CERES Score", fill = "Subtype")  + 
  theme(legend.position="bottom", legend.title = element_blank())
plot1 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(CERES)), aes(reorder(cell_line, -CERES), Core_ISG_median)) + 
  geom_bar(stat = "identity", fill = "grey") + scale_fill_manual(values = cbPalette) + theme_science(base_size = 16) + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  xlab("BRCA Cell Lines") + ylab("Core ISG Score") + scale_y_continuous(labels=scaleFUN)
multiplot(plot1, plot2, cols = 1)

tiff("ISG_Type_ADAR_CERES.tiff", height = 8, width = 7, units = "in", res = 300)
multiplot(plot1, plot2, cols = 1)
dev.off()


ggplot(subset(Dependency_CCLE, !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = RNASEL)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "RNASEL z-score")
ggsave("ADAR_BCCL_Depmap_RNASEL.tiff")
ggplot(subset(Dependency_CCLE, !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = OAS1)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "OAS1 z-score")
ggsave("ADAR_BCCL_Depmap_OAS1.tiff")
ggplot(subset(Dependency_CCLE, !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = OAS2)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "OAS2 z-score")
ggsave("ADAR_BCCL_Depmap_OAS2.tiff")
ggplot(subset(Dependency_CCLE, !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = OAS3)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "OAS3 z-score")
ggsave("ADAR_BCCL_Depmap_OAS3.tiff")

ggplot(Dependency_CCLE, aes(DEMETER2_ADAR, OAS1)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("Depmap_OAS1_point.tiff")

ggplot(Dependency_CCLE, aes(DEMETER2_ADAR, OAS2)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("Depmap_OAS2_point.tiff")

ggplot(Dependency_CCLE, aes(DEMETER2_ADAR, OAS3)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("Depmap_OAS3_point.tiff")

ggplot(subset(Dependency_CCLE, !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = RNASEL)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "RNASEL z-score")
ggsave("ADAR_BCCL_CERES_RNASEL.tiff")
ggplot(subset(Dependency_CCLE, !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = OAS1)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "OAS1 z-score")
ggsave("ADAR_BCCL_CERES_OAS1.tiff")
ggplot(subset(Dependency_CCLE, !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = OAS2)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "OAS2 z-score")
ggsave("ADAR_BCCL_CERES_OAS2.tiff")
ggplot(subset(Dependency_CCLE, !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = OAS3)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "OAS3 z-score")
ggsave("ADAR_BCCL_CERES_OAS3.tiff")

ggplot(Dependency_CCLE, aes(CERES, OAS1)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("CERES_OAS1_point.tiff")

ggplot(Dependency_CCLE, aes(CERES, OAS2)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("CERES_OAS2_point.tiff")

ggplot(Dependency_CCLE, aes(CERES, OAS3)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("Depmap_OAS3_point.tiff")


#get OAS score like ISG previously
CCLE_Z_OAS <- dplyr::select(Dependency_CCLE, one_of("OAS1", "OAS2", "OAS3")) 

OAS_scores <- data.frame(OAS_median = rowMedians(as.matrix(CCLE_Z_OAS)), Cell_Line = Dependency_CCLE$cell_line)

Dependency_CCLE <- merge(Dependency_CCLE, OAS_scores, by.x = "cell_line", by.y = "Cell_Line")

ggplot(subset(Dependency_CCLE, !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = OAS_median)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR DEMETER2 Score", fill = "Median OAS \nz-score")
ggsave("ADAR_BCCL_Depmap_OAS_median.tiff")

ggplot(Dependency_CCLE, aes(DEMETER2_ADAR, OAS_median)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("Depmap_OAS_point.tiff")

ggplot(subset(Dependency_CCLE, !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = OAS_median)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(position = "top") + geom_hline(yintercept = 0) + 
  scale_fill_gradient2(midpoint = 0, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") + 
  theme_science() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  labs(x = "BRCA Cell Lines", y="ADAR CERES Score", fill = "Median OAS \nz-score")
ggsave("ADAR_BCCL_CERES_OAS_median.tiff")

ggplot(Dependency_CCLE, aes(CERES, OAS_median)) + geom_point() + theme_science() + geom_smooth(method = "lm")
ggsave("CERES_OAS_point.tiff")


plot2 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), DEMETER2_ADAR, fill = subtype_three_receptor)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = cbPalette) + theme_science()  +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + labs(x = "", y = "ADAR DEMETER2 Score", fill = "Subtype")  + 
  theme(legend.position="bottom", legend.title = element_blank())
plot1 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(DEMETER2_ADAR)), aes(reorder(cell_line, -DEMETER2_ADAR), OAS_median)) + 
  geom_bar(stat = "identity", fill = "grey") + scale_fill_manual(values = cbPalette) + theme_science() + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  xlab("BRCA Cell Lines") + ylab("Median OAS z-score") + scale_y_continuous(labels=scaleFUN)
multiplot(plot1, plot2, cols = 1)

tiff("OAS_Type_ADAR_DEMETER2.tiff", height = 8, width = 7, units = "in", res = 300)
multiplot(plot1, plot2, cols = 1)
dev.off()

plot2 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(CERES)), aes(reorder(cell_line, -CERES), CERES, fill = subtype_three_receptor)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = cbPalette) + theme_science()  +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + labs(x = "", y = "ADAR CERES Score", fill = "Subtype")  + 
  theme(legend.position="bottom", legend.title = element_blank()) + scale_y_continuous(labels=scaleFUN)
plot1 <- ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor) & !is.na(CERES)), aes(reorder(cell_line, -CERES), OAS_median)) + 
  geom_bar(stat = "identity", fill = "grey") + scale_fill_manual(values = cbPalette) + theme_science() + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  xlab("BRCA Cell Lines") + ylab("Median OAS z-score") + scale_y_continuous(labels=scaleFUN)
multiplot(plot1, plot2, cols = 1)

tiff("OAS_Type_ADAR_CERES.tiff", height = 8, width = 7, units = "in", res = 300)
multiplot(plot1, plot2, cols = 1)
dev.off()


ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor)), aes(as.factor(subtype_three_receptor), ADAR)) + 
  geom_boxplot(width = 1, notch = TRUE, outlier.alpha = 0) + geom_jitter(width=0.2) + theme_science(base_size = 16) + 
  ylab("ADAR Expression\nZ-Score") + xlab("")
ggsave("ADAR_Receptor_Subtype.tiff")
ggplot(subset(Dependency_CCLE, !is.na(subtype_neve)), aes(as.factor(subtype_neve), ADAR)) +
  geom_boxplot(width = 1, notch = TRUE, outlier.alpha = 0) + geom_jitter(width=0.2) + theme_science() + ylab("ADAR Expression Z-Score") + xlab("")
ggsave("ADAR_Neve_Subtype.tiff")

ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor)), aes(as.factor(subtype_three_receptor), EIF2AK2)) + 
  geom_boxplot(width = 1, notch = TRUE, outlier.alpha = 0) + geom_jitter(width=0.2) + theme_science(base_size = 16) + ylab("PKR Expression\nz-score") + xlab("")
ggsave("PKR_Receptor_Subtype.tiff")
ggplot(subset(Dependency_CCLE, !is.na(subtype_neve)), aes(as.factor(subtype_neve), EIF2AK2)) +
  geom_boxplot(width = 1, notch = TRUE, outlier.alpha = 0) + geom_jitter(width=0.2) + theme_science() + ylab("PKR Expression z-score") + xlab("")
ggsave("PKR_Neve_Subtype.tiff")

ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor)), aes(as.factor(subtype_three_receptor), EIF2AK2, colour = DEMETER2_ADAR)) + 
  geom_boxplot(width = 1, notch = TRUE, outlier.alpha = 0) + geom_jitter(width=0.2) + scale_colour_gradient2(midpoint = -0.5, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab")+
  theme_science(base_size = 16) + labs(y = "PKR Expression\nz-score", x = "", colour = "ADAR \nDEMETER2\nScore")
ggsave("PKR_Receptor_Subtype_demeter.tiff")
ggplot(subset(Dependency_CCLE, !is.na(subtype_neve)), aes(as.factor(subtype_neve), EIF2AK2, colour = DEMETER2_ADAR)) +
  geom_boxplot(width = 1, notch = TRUE, outlier.alpha = 0) + geom_jitter(width=0.2) + scale_colour_gradient2(midpoint = -0.5, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") +
  theme_science() + labs(y = "PKR Expression z-score", x = "", colour = "ADAR \nDEMETER2\nScore")
ggsave("PKR_Neve_Subtype_demeter.tiff")

ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor)), aes(as.factor(subtype_three_receptor), Core_ISG_median)) +
  geom_boxplot(width = 1, notch = TRUE, outlier.alpha = 0) + geom_jitter(width=0.2) + theme_science(base_size = 16) + ylab("Core ISG Score") + xlab("")
ggsave("ISG_Receptor_Subtype.tiff")
ggplot(subset(Dependency_CCLE, !is.na(subtype_neve)), aes(as.factor(subtype_neve), Core_ISG_median)) + 
  geom_boxplot(width = 1, notch = TRUE, outlier.alpha = 0) + geom_jitter(width=0.2) + theme_science() + ylab("Core ISG Score") + xlab("")
ggsave("ISG_Neve_Subtype.tiff")

ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor)), aes(as.factor(subtype_three_receptor), Core_ISG_median, colour = DEMETER2_ADAR)) +
  geom_boxplot(width = 1, notch = TRUE, outlier.alpha = 0) + geom_jitter(width=0.2) + scale_colour_gradient2(midpoint = -0.5, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab") +
  theme_science() + labs(y = "Core ISG Score", x = "", colour = "ADAR \nDEMETER2\nScore")
ggsave("ISG_Receptor_Subtype_demeter.tiff")
ggplot(subset(Dependency_CCLE, !is.na(subtype_neve)), aes(as.factor(subtype_neve), Core_ISG_median, colour = DEMETER2_ADAR)) + 
  geom_boxplot(width = 1, notch = TRUE, outlier.alpha = 0) + geom_jitter(width=0.2) + scale_colour_gradient2(midpoint = -0.5, low = "#E69F00", mid = "white", high = "#56B4E9", space = "Lab")+ 
  theme_science() + labs(y = "Core ISG Score", x = "", colour = "ADAR \nDEMETER2\nScore")
ggsave("ISG_Neve_Subtype_demeter.tiff")


ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor)), aes(as.factor(subtype_three_receptor), OAS_median)) + 
  geom_boxplot(width = 1, notch = TRUE) + geom_jitter(width=0.2) + theme_science() + ylab("OAS Median Expression Z-Score") + xlab("")
ggsave("OAS_Receptor_Subtype.tiff")
ggplot(subset(Dependency_CCLE, !is.na(subtype_neve)), aes(as.factor(subtype_neve), OAS_median)) +
  geom_boxplot(width = 1, notch = TRUE) + geom_jitter(width=0.2) + theme_science() + ylab("OAS Median Expression Z-Score") + xlab("")
ggsave("OAS_Neve_Subtype.tiff")

ggplot(subset(Dependency_CCLE, !is.na(subtype_three_receptor)), aes(as.factor(subtype_three_receptor), OAS_median, colour = DEMETER2_ADAR)) + 
  geom_boxplot(width = 1, notch = TRUE) + geom_jitter(width=0.2) + theme_science() + ylab("OAS Median Expression Z-Score") + xlab("")
ggsave("OAS_Receptor_Subtype_demeter.tiff")
ggplot(subset(Dependency_CCLE, !is.na(subtype_neve)), aes(as.factor(subtype_neve), OAS_median, colour = DEMETER2_ADAR)) +
  geom_boxplot(width = 1, notch = TRUE) + geom_jitter(width=0.2) + theme_science() + ylab("OAS Median Expression Z-Score") + xlab("")
ggsave("OAS_Neve_Subtype_demeter.tiff")


                           

                            
#read CCLE rnaseq transcript data - pre-normalized
CCLE_transcripts <- fread("CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt")

#get only ADAR transcripts
CCLE_transcripts_ADAR <- subset(CCLE_transcripts, grepl("ENSG00000160710", CCLE_transcripts$gene_id))

#make row.names = transcript_id
CCLE_transcripts_ADAR <- data.frame(CCLE_transcripts_ADAR, row.names = CCLE_transcripts_ADAR$transcript_id)

#remove transcrpt_ids and gene_id
CCLE_transcripts_ADAR$transcript_id <- NULL
CCLE_transcripts_ADAR$gene_id <- NULL

#transpose
CCLE_transcripts_ADAR <- as.data.frame(t(CCLE_transcripts_ADAR))

#split apart cell_line identifiers and add new columns for cell line and site
out <- strsplit(as.character(row.names(CCLE_transcripts_ADAR)), "_", fixed = TRUE)

out <- do.call(rbind, out)

out <- as.data.frame(out)

CCLE_transcripts_ADAR$cell_line <- tolower(as.character(out$V1))
CCLE_transcripts_ADAR$Site <- as.character(out$V2)

#keep only breast cell lines
CCLE_transcripts_ADAR <- subset(CCLE_transcripts_ADAR, Site == "BREAST")

#rename columns for p150 and p110 trnascripts
colnames(CCLE_transcripts_ADAR)[colnames(CCLE_transcripts_ADAR) == "ENST00000368474.4"] <- "p150"
colnames(CCLE_transcripts_ADAR)[colnames(CCLE_transcripts_ADAR) == "ENST00000368471.3"] <- "p110"

#merge transcript data with cell line annotations
CCLE_transcripts_ADAR <- merge(CCLE_transcripts_ADAR, cell_line_subtypes, by.x = "cell_line", by.y = "cell_line")

#make plots of interest
plot1 <- ggplot(CCLE_transcripts_ADAR, aes(reorder(cell_line, -as.numeric(subtype_three_receptor)), p150, fill = subtype_three_receptor)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cbPalette) + theme_science(base_size = 16) + scale_y_continuous(limits = c(0,160)) + scale_x_discrete(position = "top") +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank()) + labs(x ="Breast Cancer Cell Lines", y ="p150 Expression", fill ="Subtype") 

plot2 <-ggplot(CCLE_transcripts_ADAR, aes(reorder(cell_line, -as.numeric(subtype_three_receptor)), p110, fill = subtype_three_receptor)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cbPalette) + theme_science(base_size = 16) + scale_y_continuous(limits = c(0,160)) +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank()) + labs(x ="", y ="p110 Expression", fill = "Subtype") 

multiplot(plot1, plot2, cols = 1)

tiff("ADAR_variant_TN.tiff", height = 8, width = 9, units = "in", res = 300)
multiplot(plot1, plot2, cols = 1)
dev.off()


plot1 <- ggplot(CCLE_transcripts_ADAR, aes(reorder(cell_line, -as.numeric(subtype_three_receptor)), p150, fill = subtype_neve)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cbPalette) + theme_science() + scale_y_continuous(limits = c(0,160)) + scale_x_discrete(position = "top") +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank()) + labs(x ="Breast Cancer Cell Lines", y ="p150 Expression", fill ="Subtype") 

plot2 <-ggplot(CCLE_transcripts_ADAR, aes(reorder(cell_line, -as.numeric(subtype_three_receptor)), p110, fill = subtype_neve)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cbPalette) + theme_science() + scale_y_continuous(limits = c(0,160)) +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank()) + labs(x ="", y ="p110 Expression", fill = "Subtype") 

multiplot(plot1, plot2, cols = 1)

tiff("ADAR_variant_neve.tiff", height = 8, width = 9, units = "in", res = 300)
multiplot(plot1, plot2, cols = 1)
dev.off()

CCLE_transcripts_ADAR <- merge(CCLE_transcripts_ADAR, Dependency, by = "cell_line")

ggplot(CCLE_transcripts_ADAR, aes(p110, DEMETER2_ADAR)) + geom_point() +
  theme_science() + labs(x ="p110 Expression (RSEM TPM)", y ="DEMETER2 ADAR") + geom_smooth(method = "lm") + 
  annotate("text", x = 120, y = -.8, label = paste("Pearson r =", round(cor.test(CCLE_transcripts_ADAR$p110, CCLE_transcripts_ADAR$DEMETER2_ADAR)$estimate, digits = 2))) +
  annotate("text", x = 120, y = -.9, label = paste("p-value =", round(cor.test(CCLE_transcripts_ADAR$p110, CCLE_transcripts_ADAR$DEMETER2_ADAR)$p.value, digits = 2)))
ggsave("p110_v_DEMETER2ADAR.tiff", height = 4, width = 4, unit = "in")

ggplot(CCLE_transcripts_ADAR, aes(p150, DEMETER2_ADAR)) + geom_point() +
  theme_science() + labs(x ="p150 Expression (RSEM TPM)", y ="DEMETER2 ADAR") + geom_smooth(method = "lm") + 
  annotate("text", x = 50, y = 0, label = paste("Pearson r =", round(cor.test(CCLE_transcripts_ADAR$p150, CCLE_transcripts_ADAR$DEMETER2_ADAR)$estimate, digits = 2))) +
  annotate("text", x = 50, y = -0.1, label = paste("p-value =", round(cor.test(CCLE_transcripts_ADAR$p150, CCLE_transcripts_ADAR$DEMETER2_ADAR)$p.value, digits = 2)))
ggsave("p150_v_DEMETER2ADAR.tiff", height = 4, width = 4, unit = "in")


#get RPPA data
CCLE_RPPA <- fread("CCLE_RPPA_20181003.csv")
#split apart cell_line identifiers and add new columns for cell line and site
out <- strsplit(as.character(CCLE_RPPA$V1), "_", fixed = TRUE)

out <- do.call(rbind, out)

out <- as.data.frame(out)

CCLE_RPPA$cell_line <- tolower(as.character(out$V1))

CCLE_RPPA <- merge(CCLE_RPPA, cell_line_subtypes, by.x = "cell_line", by.y = "cell_line")


#make plots of interest
ggplot(CCLE_RPPA, aes(reorder(cell_line, -as.numeric(ADAR1)), ADAR1, fill = subtype_three_receptor)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cbPalette) + theme_science() + 
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank()) + labs(x="BRCA Cell Lines", y = "ADAR Expression", fill = "Subtype")
ggsave("RPPA_ADAR_TN.tiff")


ggplot(CCLE_RPPA, aes(reorder(cell_line, -as.numeric(ADAR1)), ADAR1, fill = subtype_neve)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cbPalette) + theme_science() + 
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank()) + labs(x="BRCA Cell Lines", y = "ADAR Expression", fill = "Subtype")
ggsave("RPPA_ADAR_neve.tiff")



