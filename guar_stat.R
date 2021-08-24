#packages <- c("Rcpp", "RcppEigen", "RcppParallel", "parallel", "doParallel", 
#    "foreach", "bigmemory", "devtools")
#install.packages(packages)
#devtools::install_github(repo = "amkusmec/FarmCPUpp")

setwd("/media/array/guar_proj/production_callset")
library(ggplot2)
library(reshape2)
library(ggsci)
library(VennDiagram)
library(CMplot)
library(FarmCPUpp)
require(bigmemory)
library(lattice)

mycol1 = rgb(110, 224, 255, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol3 = rgb(255, 177, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)


stat_table = read.table('ALL_variant_stats.tsv', sep='\t', header=F)
colnames(stat_table) = c('locus', 'tool', 'AF', 'HET_gt', 'ALT_gt', 
                         'AB_ref', 'AB_het', 'AB_hom', 
                         'DP_ref', 'DP_het', 'DP_hom')
head(stat_table)
stat_table = stat_table[!is.na(stat_table$AF), ]
table(stat_table$tool)

ggplot(stat_table, aes(x=AF, fill=tool)) + geom_density(col='black') +
  theme_bw() + facet_wrap(~tool, scales='free')

cnts <- ggplot(stat_table, aes(x=tool, fill=tool)) +
  geom_histogram(col='black', stat='count') + theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x=element_blank()) +
  scale_fill_manual(values=c(mycol2, mycol1, mycol3)) + guides(fill=F) +
  xlab('') + ylab('Variant count') + guides(fill=F)
print(cnts)

a <- ggplot(stat_table, aes(x=tool, y=AF, fill=tool)) + 
  geom_violin(col='black', scale = "width") +
  geom_boxplot(outlier.shape=NA, col='black', width=0.25, fill='white') +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x=element_blank()) +
  scale_fill_manual(values=c(mycol2, mycol1, mycol3)) + guides(fill=F) + 
  xlab('') + ylab('Non-reference AF')

b <- ggplot(stat_table[as.numeric(stat_table$AF) > 0.05, ], aes(x=tool, y=HET_gt, fill=tool)) + 
  geom_violin(col='black', scale = "width") +
  geom_boxplot(outlier.shape=NA, col='black', width=0.25, fill='white') +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x=element_blank()) +
  scale_fill_manual(values=c(mycol2, mycol1, mycol3)) + guides(fill=F) +
  xlab('') + ylab('Number of HETs')

c <- ggplot(stat_table[as.numeric(stat_table$AF) > 0.05, ], aes(x=tool, y=AB_het, fill=tool)) + 
  geom_violin(col='black', scale = "width") +
  geom_boxplot(outlier.shape=NA, col='black', width=0.25, fill='white') +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x=element_blank()) +
  scale_fill_manual(values=c(mycol2, mycol1, mycol3)) + guides(fill=F) +
  xlab('') + ylab('HET allele balance')

d <- ggplot(stat_table[as.numeric(stat_table$AF) > 0.05, ], aes(x=tool, y=DP_het, fill=tool)) + 
  geom_violin(col='black', scale = "width") +
  geom_boxplot(outlier.shape=NA, col='black', width=0.25, fill='white') +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x=element_blank()) +
  scale_fill_manual(values=c(mycol2, mycol1, mycol3)) + guides(fill=F) +
  scale_y_continuous(limits=c(0, 50)) +
  xlab('') + ylab('Depth at variant site')

top = plot_grid(a, b, c, d, nrow=1, labels=c('a' , '', ''),
                rel_widths = c(1, 1, 1, 1))
print(top)

overlaps = c()
for (i in c('GATK', 'NGSEP', 'TASSEL')){
  for (j in c('GATK', 'NGSEP', 'TASSEL')){
    vars_a = stat_table[stat_table$tool == i, 'locus']
    vars_b = stat_table[stat_table$tool == j, 'locus']
    denom = min(length(vars_a), length(vars_b))
    overlaps = c(overlaps, length(intersect(vars_a, vars_b))/denom)
  }
}

comparisons = matrix(overlaps, nrow=3, ncol=3)
rownames(comparisons) = c('GATK', 'NGSEP', 'TASSEL')
colnames(comparisons) = c('GATK', 'NGSEP', 'TASSEL')

drawCoolHM = function(df, p_df){
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y,  p_df[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet, 
                   at=seq(min(df), max(df), length.out=100), 
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=45)), xlab=list(label=''), 
                   ylab=list(label=''), panel=myPanel_a))
}

#mypal = colorRampPalette(c('white', '#ecad2f'))
mypal = colorRampPalette(c('white', '#d50f0dff'))
jet = mypal(100)
#jet = mypal(100)

e <- drawCoolHM(comparisons, round(comparisons, digits=3))
pholder <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text = element_blank(), 
                     axis.title = element_blank(),
                     axis.ticks = element_blank(),
                     panel.border = element_blank())

# Venn Diagram
isec <- get.venn.partitions(list(gatk=stat_table[stat_table$tool == 'GATK', 'locus'],
          ngsep=stat_table[stat_table$tool == 'NGSEP', 'locus'],
          tassel=stat_table[stat_table$tool == 'TASSEL', 'locus']))[, c(1, 2, 3, 6)]
vd <- venn.diagram(list(GATK=stat_table[stat_table$tool == 'GATK', 'locus'],
        NGSEP=stat_table[stat_table$tool == 'NGSEP', 'locus'],
        TASSEL=stat_table[stat_table$tool == 'TASSEL', 'locus']), filename = NULL,
        fill=c(mycol2, mycol1, mycol3), euler.d=T, margin=0.1)

bottom_left <- plot_grid(e, vd, labels=c('b', 'c'), ncol=1, 
                         rel_heights = c(1, 1))
bottom <- plot_grid(bottom_left, pholder, rel_widths = c(0.7, 1))
plot_grid(top, bottom, nrow=2, rel_heights = c(0.45, 1))


# Validation diagram - AF vs HET genotypes

final_vars = read.table('/media/array/guar_proj/production_callset/final_variants.list',
                        header=F)$V1

gatk_vars = stat_table[stat_table$tool == 'GATK' & 
      as.character(stat_table$locus) %in% final_vars, ]
gatk_vars$bin = ceiling(gatk_vars$AF/0.1)

agg_data = aggregate(HET_gt~bin, gatk_vars, function(x) mean(x)/192)
agg_data$mean_AF = aggregate(AF~bin, gatk_vars, mean)$AF
agg_data$expected = 2*agg_data$mean_AF*(1-agg_data$mean_AF)
het_freqs = melt(agg_data, id.vars=c('bin', 'mean_AF'))
#ggplot(het_freqs, aes(x=mean_AF, y=value, col=variable)) + 
#  geom_line() + geom_point() + theme_bw()

gatk_vars$H_e = sapply(gatk_vars$AF, function(x) 
               rbinom(1, 192,  2 * x * (1 - x)))
colnames(gatk_vars)[4] = 'H_o'
gatk_vars = gatk_vars[gatk_vars$H_o <= 100, ]
hets = melt(gatk_vars, id.vars=c('locus', 'AF'), 
            measure.vars=c('H_o', 'H_e'))
het_gt_plot <- ggplot(hets, aes(x=AF, y=value, col=variable)) + 
  geom_point(alpha=0.01) + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) + 
  theme_bw() + 
  theme(panel.grid=element_blank(), 
        legend.position = c(0.85, 0.85),
        legend.title = element_blank()) +
  xlab('Non-reference allele frequency') + 
  ylab('Number of heterozygotes')
print(het_gt_plot)
#  guides(col=F)


# PCA 
pc_scores = read.table('/media/array/guar_proj/production_callset/out_pc_scores/guar_filtered.pc_scores.csv',
                       sep=',', header=T, row.names=1)
head(pc_scores)
ancestry = read.table('/media/array/guar_proj/production_callset/guar_ancestry_v3.tsv',
                      sep='\t', header=T)

pc_scores$country = sapply(pc_scores$genome_id, function(x)
  ancestry[as.character(ancestry$plant) == as.character(x), 'country'])


pca_plot_country <- ggplot(pc_scores, aes(x=PC1, y=PC2)) + 
  geom_point(aes(fill=country), size=2, pch=21, col='black') + theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'bottom') +
  xlab('PC1 - 22.8% of variance') + ylab('PC2 - 15.4% of variance') +
  scale_fill_npg()

# InbreedingCoeff

head(gatk_vars)
gatk_vars$Inbreeding = 1 - gatk_vars$H_o/gatk_vars$H_e
inbr <- ggplot(gatk_vars[gatk_vars$AF > 0.05, ], aes(x=1, y=Inbreeding)) + 
  geom_violin() + 
  geom_boxplot(outlier.shape=NA, width=0.4) +
  coord_flip() + theme_bw() + theme(panel.grid = element_blank())



# LD figure
ld_data = read.table('/media/array/guar_proj/production_callset/ld_data.tsv',
                     sep='\t')
colnames(ld_data) = c('loc_A', 'pos_A', 'AF_A', 'loc_B', 'pos_B', 'AF_B', 'D')
ld_data$pos_diff = abs(ld_data$pos_A - ld_data$pos_B)/1000000
hist(ld_data$pos_diff)
ld_data$D = abs(ld_data$D)
ld_data = ld_data[ld_data$pos_diff <= 1, ]

ggplot(ld_data, aes(x=pos_diff, y=D)) + 
  geom_point(alpha=0.01) +
  geom_smooth(method='lm', formula=(y~exp(-x))) + 
  theme_bw() + 
  scale_x_continuous(limits=c(0, 1)) +
  xlab('Physical distance (Mbp)') +
  ylab('Pairwise r2 score')

ld_data$binned = as.factor(floor(ld_data$pos_diff/0.05)*0.05)
bin_ld <- ggplot(ld_data, aes(x=binned, y=D)) + 
  geom_boxplot(outlier.shape=NA) +
  theme_bw() + theme(axis.text.x=element_blank(), panel.grid = element_blank())


cd <- plot_grid(inbr, bin_ld, labels=c('c', 'd'), nrow=2, rel_heights = c(0.5, 1)) 
top_2 <- plot_grid(pca_plot_country, cd, ncol=2, labels=c('a', ''))
print(top_2)

# ADMIXTURE

plot_admix <- function(admix_file) {
  admix_3 = read.table(admix_file,
                     sep=' ', header=F)
  colnames(admix_3) = paste0('P', 1:ncol(admix_3))
  admix_3$origin = pc_scores$country
  admix_3$pops <- apply(admix_3, 1, which.max)
  admix_3 = admix_3[order(admix_3$pops), ]
  admix_3$fake_id = pc_scores$genome_id

  adm3 = melt(admix_3, id.vars=c('fake_id', 'pops', 'origin'))
  p <- ggplot(adm3, aes(x=fake_id, y=value, fill=variable)) + 
    geom_bar(stat='identity', width=1) + theme_bw() + 
    theme(panel.grid=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank()) +
    facet_grid(~origin, scales = 'free', space = 'free') +
    scale_fill_simpsons() + guides(fill=F)
  return(p)
}

d1 <- plot_admix('/media/array/guar_proj/production_callset/GUAR_2of3_AF5pct_hail.3.Q')
d2 <- plot_admix('/media/array/guar_proj/production_callset/GUAR_2of3_AF5pct_hail.4.Q')
d3 <- plot_admix('/media/array/guar_proj/production_callset/GUAR_2of3_AF5pct_hail.5.Q')

bottom_2 <- plot_grid(d1, d2, d3, nrow=3)
print(bottom_2)

plot_grid(top_2, bottom_2, nrow=2, rel_heights = c(1, 0.8), labels=c('', 'b'))

plot_admix('/media/array/guar_proj/production_callset/GUAR_2of3_AF5pct_hail.2.Q')
plot_admix('/media/array/guar_proj/production_callset/GUAR_2of3_AF5pct_hail.6.Q')


# Plotting GWAS  results

gwas_tab = read.table('/media/array/guar_proj/production_callset/glm_pvals/branch_height.tsv',
                      sep='\t', header=T)
head(gwas_tab)

previous_chr = ''
chr_index = 0
numbered_chr = c()
for (i in gwas_tab$locus){
  chrom = strsplit(i,':')[[1]][1]
  if (chrom != previous_chr) {
    chr_index = chr_index + 1
    previous_chr = chrom
  }
  numbered_chr = c(numbered_chr, chr_index)
}
gwas_tab$chr = numbered_chr
gwas_tab$pos = sapply(gwas_tab$locus, function(x) strsplit(x, ':')[[1]][2])
gwas_tab$SNP = gwas_tab$locus
gwas_tab = gwas_tab[, c('SNP', 'chr', 'pos')]

for (i in dir('/media/array/guar_proj/production_callset/glm_pvals/')) {
  this_trait = read.table(paste0('/media/array/guar_proj/production_callset/glm_pvals/', i), 
          sep='\t', header=T)
  trait_name = strsplit(i, '\\.')[[1]][1]
  print(trait_name)
  gwas_tab[, trait_name] = this_trait$p_value
}

CMplot(gwas_tab,plot.type="q",threshold=1e-5,cex=2.5,
       ylab.pos=2,signal.cex=3,signal.col="red",
       conf.int=TRUE,box=FALSE,multracks=TRUE,
       cex.axis=2,cex.lab=2,file="jpg",memo="",dpi=300,
       file.output=TRUE,verbose=TRUE,ylim=c(0,6),width=4,height=5)

CMplot(gwas_tab, plot.type="m",multracks=TRUE,
       threshold=c(1e-5,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), 
       amplify=TRUE,bin.size=1e6,
       signal.cex=1, file="jpg",memo="",dpi=500,
       file.output=TRUE,verbose=TRUE)

CMplot(gwas_tab,type="p",plot.type="c",r=0.4,cir.legend=TRUE,
       outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,
       chr.den.col="black",file="jpg",
       memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)


# Making FarmCPU and plotting

myY <- read.table('../gwas_farmcpu/FarmCPU_Y_revamp.tsv',
                  header = TRUE, stringsAsFactors = FALSE)
myGM <- read.table('../gwas_farmcpu/FarmCPU_GM.revamp.tsv',
                   header = TRUE, stringsAsFactors = FALSE)
myGM$Chromosome = sapply(as.character(myGM$Chromosome), 
                         function(x) which(unique(as.character(myGM$Chromosome)) == x))
head(myGM)

myGD <- read.big.matrix('../gwas_farmcpu/FarmCPU_GD.revamp.tsv',
                        type = "double", sep = "\t", header = TRUE,
                        col.names = myGM$SNP, ignore.row.names = FALSE,
                        has.row.names = TRUE)
rownames(myGD)
rownames(myY) = myY$taxa
myY = myY[rownames(myGD), ]

good_traits = c(2:6)
gwas_fc_tab = myGM
for (i in good_traits){
  myResults <- farmcpu(Y = myY[, c(1, i)], GD = myGD, GM = myGM, maxLoop = 100)
  gwasRes <- myResults[[colnames(myY)[i]]]$GWAS[, 1:4]
  colnames(gwasRes)[4] = colnames(myY)[i]
  write.table(gwasRes, 
              file=paste0('farmcpu_pvals/', colnames(myY)[i], '.tsv'), 
              sep='\t', row.names=F, quote=F)
  gwas_fc_tab[colnames(myY)[i]] = gwasRes[colnames(myY)[i]]
  CMplot(gwasRes,plot.type="q",box=FALSE,file="jpg",memo="",dpi=300,
         conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
         file.output=TRUE,verbose=TRUE,width=5,height=5)
  
}

CMplot(gwas_fc_tab,plot.type="q",threshold=1e-4,cex=2.5,
       ylab.pos=2,signal.cex=2.5,signal.col="red",
       conf.int=TRUE,box=FALSE,multracks=TRUE,
       cex.axis=2,cex.lab=2,file="jpg",memo="",dpi=300,
       file.output=TRUE,verbose=TRUE,ylim=c(0,8),width=4,height=5)

CMplot(gwas_fc_tab, plot.type="m",multracks=TRUE,
       threshold=c(1e-5,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), 
       amplify=TRUE,bin.size=1e6,
       signal.cex=1, file="jpg",memo="",dpi=500,
       file.output=TRUE,verbose=TRUE)

hits = gwas_fc_tab[gwas_fc_tab$maturation_pct < 1e-4, ]

write.table(hits, file='bean_maturation_hits.tsv', sep='\t', row.names=F, quote=F)

myY_tocorr <- read.table('/media/array/guar_proj//gwas_farmcpu/FarmCPU_Y_revamp.tsv',
                  header = TRUE, stringsAsFactors = FALSE)
dim(myY_tocorr)

corr_matrix <- cor(as.matrix(myY_tocorr[, 2:6]), use = "pairwise.complete.obs")
image(corr_matrix)

drawCorrHM = function(df, p_df){
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y,  p_df[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet, 
                   at=seq(-1, 1, length.out=100), 
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=45)), xlab=list(label=''), 
                   ylab=list(label=''), panel=myPanel_a))
}
jet = colorspace::diverge_hsv(100)
drawCorrHM(corr_matrix, round(corr_matrix, digits=3))


# Top hits box-plots

top_hits <- read.table('hits_maturation_stats.tsv', sep='\t', header=F)
colnames(top_hits) = c('variant_id', 'gt', 'maturation_pct')

ggplot(top_hits, aes(x=variant_id, y=maturation_pct, fill=gt)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

top_hits$ref_nonref = ifelse(top_hits$gt == 'REF', 'REF', 'NONREF')

for (i in unique(top_hits$variant_id)){
  data_piece = top_hits[top_hits$variant_id == i, ]
  print(i)
  print(wilcox.test(maturation_pct~ref_nonref, data_piece)$p.value)
}
