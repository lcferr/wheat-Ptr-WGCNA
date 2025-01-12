###Load packages ####
pacman::p_load(tidyverse, WGCNA, data.table, ggdendro,magrittr, DESeq2)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

###Create input for WGCNA####
dds<- readRDS("dds.rds")
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
input_mat = t(expr_normalized)

write.csv(input_mat, "inputWGCNA.csv", row.names = T, quote = F)

###Read input files####
datExpr <- read.csv("inputWGCNA.csv", row.names = 1)
annot = read.csv("FunctionalAnnotation.csv");
traitData = read.csv("traitData.csv");

dim(traitData)
names(traitData)
rnaSamples = rownames(datExpr);
traitRows = match(rnaSamples, traitData$Name);
datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage();

###Check for genes and samples with too many missing values####
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

###Dendograms####
sampleTree = hclust(dist(datExpr), method = "average")
par(mfrow = c(1,2));
cex1 = 0.6;
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
abline(h = 85, col = "red");

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 85, minSize = 10)
table(clust)

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE);

png("plotDendroAndColors.png", units = "in",  width=10, height=5, res=300)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "")


dev.off()

###Network construction####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  datExpr,             
  powerVector = powers,
  verbose = 5)

png("scalefreeR2_meanK.png", units = "in",  width=10, height=5, res=300)
par(mfrow = c(1,2));
par(mar = c(4, 5, 2, 1))
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = expression("Scale Free Topology Model Fit, signed " ~R^2),
     main = "")

text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]-0.02) * sft$fitIndices[, 2]-0.03,
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
abline(h = 20, col = "red")

dev.off()

softPower = 16;
adjacency = adjacency(datExpr, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

###Module detection####
picked_power = 16
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
net <- blockwiseModules(datExpr,                # <= input here
                        
                        # == Adjacency Function ==
                        power = picked_power,                # <= power here
                        networkType = "signed",
                        
                        # == Tree and Block Options ==
                        deepSplit = 2,
                        pamRespectsDendro = F,
                        # detectCutHeight = 0.75,
                        minModuleSize = 30,
                        maxBlockSize = 3000,
                        
                        # == Module Adjustments ==
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        
                        # == TOM == Archive the run results in TOM file (saves time)
                        saveTOMs = T,
                        saveTOMFileBase = "ER",
                        
                        # == Output Options
                        numericLabels = T,
                        verbose = 3)

cor <- temp_cor     # Return cor function to original namespace

table(net$colors) #see how many modules were identified and what the module sizes are

#The label 0 is reserved for genes outside of all modules

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

png("plotDendroAndColors.png", units = "in",  width=10, height=5, res=300)
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

###Corr heatmap####
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
png("labeledHeatmap_detail.png", units = "in",  width=8, height=8, res=600)
par(mar = c(6, 8, 8, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = "")
dev.off()

###Plot Expression Profiles####
modules_of_interest = c("black","pink","blue")

# Pull out list of genes in that module

module_df <- data.frame(
  gene_id = names(net$colors),
  colors = labels2colors(net$colors)
)
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
subexpr = t(datExpr)[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df$modules_f = factor(submod_df$module, levels = c("black","pink","blue"))

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none", axis.title = element_text(face ="plain")) +
  facet_grid(rows = vars(modules_f)) +
  labs(x = "Samples",
       y = "Normalised expression") +
  scale_color_manual(values = c("#000000", "#0000FF","#FFC0CB"))

ggsave("normalised_expression2.png", width = 6, height = 5, device = "png", dpi = 600)

###Gene Significance and Module Membership####
# Define variable weight containing the weight column of datTrait

#Blue module-HerI_RobI
weight = as.data.frame(datTraits$HerI_RobI);
names(weight) = "HerI_RobI"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
png("geneModuleMembership_blue.png", units = "in",  width=5, height=5, res=300)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Ptr infection",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#Black module-Rob
weight = as.data.frame(datTraits$Rob);
names(weight) = "Rob"
geneTraitSignificanceRob = as.data.frame(cor(datExpr, weight, use = "p"));
names(geneTraitSignificanceRob) = paste("GS.", names(weight), sep="");
GSPvalueRob = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceRob), nSamples));
names(geneTraitSignificanceRob) = paste("GS.", names(weight), sep="");
names(GSPvalueRob) = paste("p.GS.", names(weight), sep="")

module = "black"
column = match(module, modNames);
moduleGenes = moduleColors==module;
png("geneModuleMembership_black.png", units = "in",  width=5, height=5, res=300)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceRob[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#Pink module-Her
weight = as.data.frame(datTraits$Her);
names(weight) = "Her"
geneTraitSignificanceHer = as.data.frame(cor(datExpr, weight, use = "p"));
names(geneTraitSignificanceHer) = paste("GS.", names(weight), sep="");
GSPvalueHer = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceHer), nSamples));
names(geneTraitSignificanceHer) = paste("GS.", names(weight), sep="");
names(GSPvalueHer) = paste("p.GS.", names(weight), sep="")

module = "pink"
column = match(module, modNames);
moduleGenes = moduleColors==module;
png("geneModuleMembership_pink.png", units = "in",  width=5, height=5, res=300)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceHer[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

###Hub genes####

hubs    = chooseTopHubInEachModule(datExpr, mergedColors)
hubs

write.csv(hubs, "hubs.csv")

x <- cbind(geneModuleMembership, geneTraitSignificance, MMPvalue, GSPvalue, geneTraitSignificanceRob, geneTraitSignificanceHer)

topgenes <- subset(x, abs(GS.HerI_RobI) > 0.6 & MMblue > 0.8)
topgenes <- topgenes[c("GS.HerI_RobI", "p.GS.HerI_RobI", "MMblue", "p.MMblue")]
nrow(topgenes)

write.csv(topgenes, "bluehubs.csv")

###Exporting to Cytoscape####

TOM = TOMsimilarityFromExpr(datExpr, power = 16);

# Select modules

modules= "black"
# Select module probes
probes = names(as.data.frame(datExpr))
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$Human.Readable.Description[match(modProbes, annot$Gene.ID)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               altNodeNames=modGenes,
                               nodeAttr = moduleColors[inModule])

