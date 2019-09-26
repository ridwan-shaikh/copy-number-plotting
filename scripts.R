library(knitr)
library(kableExtra)
library(XML)
library(methods)
library(dplyr)
library(tidyr)
library(reshape2)
library(jsonlite)
library(readr)
library(ggplot2)
suppressMessages(library(karyoploteR))
copy_number_file <- list.files(path = ".", pattern = ".cnr")
copy_number <- read.delim(copy_number_file)
#subset copy number for only the antitarget regions
copy_number_plot <- copy_number[which(copy_number$gene == "Antitarget"),]

#obtain patient report for amplificiations and deletions
reported_copy_number_file <- list.files(path = ".", pattern = "patient.report")
reported_copy_number <- read.delim(reported_copy_number_file, skip = 1, sep = "\t")
Amplificaition <- reported_copy_number[grep(pattern = "Amplification",
                                            ignore.case = T,
                                            x = reported_copy_number$Gene_Result),]
Deletions <- reported_copy_number[grep(pattern = "Deletion",
                                       ignore.case = T,
                                       x = reported_copy_number$Gene_Result),]
Amps_dels <-rbind(Amplificaition, Deletions)
Amps_dels <- subset(Amps_dels,
                    select = c("Gene",
                               "Gene_ROI",
                               "Gene_Result"))
names(Amps_dels)[names(Amps_dels) == 'Gene_ROI'] <-'ROI'
names(Amps_dels)[names(Amps_dels) == 'Gene_Result'] <-'Result'

#find locations of the amplifications and deletions and seperate the gene name out of the gene column
locations_amps_dels <- copy_number[grep(pattern = paste0(Amps_dels$Gene,collapse = "|","_"), x = copy_number$gene),]
locations_amps_dels <- locations_amps_dels %>% separate(gene, into = c("gene","info"), sep = "_")

#convert to granges and make the gene a factor
grange_locations_amps_dels <- toGRanges(locations_amps_dels)
grange_locations_amps_dels$gene <- as.factor(grange_locations_amps_dels$gene)
#condense the granges into a range for the whole gene to plot as markers
markers <-unlist(range(split(grange_locations_amps_dels, ~gene)))
split_markers <- split(markers, names(markers))
length(split_markers)
length(markers)
names(markers)

as.data.frame(split(markers, as.factor(markers)))

#attempt at using the base plot function seen as though zoom in plotKaryotypes doesnt allow more than one argument
find_FGFR <- markers[1] + 1e8
FGFR <- subsetByOverlaps(grange_copy_number, find_FGFR)
plot(FGFR$log2, pch = 20)


#https://github.com/bernatgel/CopyNumberPlots.git
list_regions <-unlist(levels(as.factor(markers)))
list_regions <-unlist(levels(as.factor(markers + 5e6)))

for (i in seq_along(markers)) {
  plotKaryotype(genome = "hg19", plot.type = 4, zoom = (markers[i] + 5e6))
}


#zoom into gains and losses
zoomkp <- plotKaryotype(genome = "hg19", plot.type = 4, zoom = (markers[1] + 5e6))
kpAddBaseNumbers(zoomkp, minor.tick.dist = 1000000)
#this needs context... plot all the regions... not sure how i will make this work if you have multiple amplifications and losses
grange_copy_number <- toGRanges(copy_number)
kpAxis(zoomkp, ymin = -5, ymax = 5, cex=0.5, tick.pos = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5))
kpPoints(zoomkp, data = grange_copy_number, y = grange_copy_number$log2, ymin = -10, ymax = 10, col = "red")
copy_number_plot <- toGRanges(copy_number_plot)
kpPoints(zoomkp, data = copy_number_plot, y=copy_number_plot$log2, ymin = -10, ymax = 10)
#kpPlotRegions(zoomkp, data = markers, col = "#D55E00")
kpPlotMarkers(zoomkp, data = markers, labels = markers@ranges@NAMES, text.orientation = "horizontal", r1=1.29, line.color = "#D55E00", lty = 3)


#Gviz version
library(Gviz)
chroms <- c("chr1", "chr2", "chr3", "chr4")
maTrack <- AnnotationTrack(range = GRanges(seqnames = chroms,
                                           ranges = IRanges(start = 1, width = c(100, 40, 200, 1000)),
                                           strand = c("+", "+", "-", "+")),
                           genome = "mm9", chromosome = "chr1", name = "foo")
mdTrack <- DataTrack(range = GRanges(seqnames = rep(chroms, c(10, 40, 20 ,100)),
                                     ranges = IRanges(start = c(seq(1, 100, len = 10),
                                                                seq(1, 400, len = 40),
                                                                seq(1, 200, len = 20),
                                                                seq(1, 1000, len = 100)),
                                                      width = 9),
                                     values = runif(170)),
                     data = "values",
                     chromosome = "chr1",
                     genome = "mm9", name = "bar")
mgTrack <- GenomeAxisTrack(scale = 50, labelPos = "below", exponent = 3)
itrack <- IdeogramTrack(genome = "hg19")
chromosome(itrack) <- "chr1"
ncols <- 2
nrows <- length(chroms)%/%ncols
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrows, ncols)))
for (i in seq_along(chroms)) {
  pushViewport(viewport(layout.pos.col = ((i - 1)%%ncols) + 1, layout.pos.row = (((i) - 1) %/% ncols) + 1))
  plotTracks(list(itrack, maTrack, mdTrack, mgTrack), chromosome = chroms[i], add = TRUE)
  popViewport(1)
}

pushViewport(viewport(layout = grid.layout(nrows, ncols)))
for (i in seq_along(markers)) {
  pushViewport(viewport(layout.pos.col = ((i - 1) %% ncols) + 1, layout.pos.row = (((i) -1 ) %/% ncols) + 1))
  plotKaryotype(genome = "hg19", plot.type = 4, zoom = markers[i])
  popViewport(1)
}

library(gridExtra)
library(grid)
library(lattice)
library(karyoploteR)
library(magrittr)
#this works!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for (i in seq_along(markers)) {
pdf(file = paste((i),'.pdf', sep=""))
plotKaryotype(genome = "hg19", plot.type = 4, zoom = (markers[i] + 5e6)) %>%
  kpAddBaseNumbers(minor.tick.dist = 1000000) %>%
  kpAxis(ymin = -5, ymax = 5, cex=0.5, tick.pos = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)) %>%
  kpPoints(data = grange_copy_number, y = grange_copy_number$log2, ymin = -10, ymax = 10, col = "red") %>%
  kpPoints(data = copy_number_plot, y=copy_number_plot$log2, ymin = -10, ymax = 10) %>%
  kpPlotMarkers(data = markers, labels = markers@ranges@NAMES, text.orientation = "horizontal", r1=1.29, line.color = "#D55E00", lty = 3)
dev.off()
}

for (i in seq_along(markers)) {
  png(file = paste((i),'.png', sep=""))
  plotKaryotype(genome = "hg19", plot.type = 4, zoom = (markers[i] + 5e6)) %>%
    kpAddBaseNumbers(minor.tick.dist = 1000000) %>%
    kpAxis(ymin = -5, ymax = 5, cex=0.5, tick.pos = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)) %>%
    kpPoints(data = grange_copy_number, y = grange_copy_number$log2, ymin = -10, ymax = 10, col = "red") %>%
    kpPoints(data = copy_number_plot, y=copy_number_plot$log2, ymin = -10, ymax = 10) %>%
    kpPlotMarkers(data = markers, labels = markers@ranges@NAMES, text.orientation = "horizontal", r1=1.29, line.color = "#D55E00", lty = 3)
  dev.off()
}

plotKaryotype(genome = "hg19", plot.type = 4, zoom = (markers[1] + 5e6)) %>%
  kpAddBaseNumbers(minor.tick.dist = 1000000) %>%
  kpAxis(ymin = -5, ymax = 5, cex=0.5, tick.pos = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)) %>%
  kpRect(data = markers, y0 = 0, y1 = 1, col="#FFDDDD", border=NA, r0 = 0.05, r1 = 1) %>%
  kpPoints(data = grange_copy_number, y = grange_copy_number$log2, ymin = -10, ymax = 10, col = "red") %>%
  kpPoints(data = copy_number_plot, y=copy_number_plot$log2, ymin = -10, ymax = 10) %>%
  kpPlotMarkers(data = markers, labels = markers@ranges@NAMES, text.orientation = "horizontal", r1=0, line.color = "#D55E00", lty = 0)






#####IchorCNA output in karyoploteR######
library(karyoploteR)

input <- read.delim("1601104-NMB866-T.cna.seg", header = T)#[,c(-4:-7)]
input$chr <- paste0("chr", input$chr)
colnames(input)[4] <- "copy.number"
colnames(input)[5] <- "event"
colnames(input)[6] <- "logR"
gr_input <- toGRanges(input)
plot.params <- getDefaultPlotParams(plot.type=5)
plot.params$leftmargin <- 0.01
plot.params$rightmargin <- 0.01
kp <- plotKaryotype(genome = "hg19", plot.type = 5, labels.plotter = NULL, plot.params = plot.params)
kpAddChromosomeNames(kp, cex=0.7)
kpDataBackground(kp, data.panel = 1, color = c("grey","lightgrey"))
#kpDataBackground(kp, data.panel = 1, color = "white")
kpPoints(kp, data = gr_input, y = log(gr_input$logR), data.panel = 1, r0 = 0.95, r1 = 0, ymin = ceiling(min(log(gr_input$logR), na.rm = T)), ymax = ceiling(max(log(gr_input$logR), na.rm = T)), cex = 0.3, pch = 21)
kpAxis(kp, ymin = -2, ymax = 2, cex=0.5, tick.pos = c(-2, -1, 0, 1, 2), r0=1, r1=0)
kpAbline(kp, h = 0.5, col ="black")
kpPoints(kp, data = gr_input[gr_input$event == "HOMD"], y = gr_input[gr_input$event == "HOMD"]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#FEE0D2", pch = 21)
kpPoints(kp, data = gr_input[gr_input$event == "HETD"], y = gr_input[gr_input$event == "HETD"]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#FC9272", pch = 21)
kpPoints(kp, data = gr_input[gr_input$event == "NEUT"], y = gr_input[gr_input$event == "NEUT"]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="grey90", pch = 21)
kpPoints(kp, data = gr_input[gr_input$event == "GAIN"], y = gr_input[gr_input$event == "GAIN"]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#EFF3FF", pch = 21)
kpPoints(kp, data = gr_input[gr_input$event == "AMP"], y = gr_input[gr_input$event == "AMP"]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#C6DBEF", pch = 21)
kpPoints(kp, data = gr_input[gr_input$event == "HLAMP"], y = gr_input[gr_input$event == "HLAMP"]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#9ECAE1", pch = 21)
kpPoints(kp, data = gr_input[gr_input$event == "HLAMP1"], y = gr_input[gr_input$event == "HLAMP1"]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#6BAED6", pch = 21)
kpPoints(kp, data = gr_input[gr_input$event == "HLAMP2"], y = gr_input[gr_input$event == "HLAMP2"]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col= "#3182BD", pch = 21)
kpPoints(kp, data = gr_input[gr_input$event == "HLAMP3"], y = gr_input[gr_input$event == "HLAMP3"]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#08519C", pch = 21)
kpPlotMarkers(kp, data = paed_markers, labels = paed_markers@ranges@NAMES, text.orientation = "vertical", r1=1.35, line.color = "#D55E00", lty = 0, adjust.label.position = F)
kpRect(kp, data = markers, y0 = 0, y1 = 1, col="grey", border="grey", r0 = 1, r1 = 0)
markers <- toGRanges(paed_bed)
paed_markers <-unlist(range(split(markers, ~V5)))


length(levels(gr_input$event))


cfc_output <- read.delim("1601104-NMB866-T.sort.bam_ratio.txt")
cfc_output$Chromosome <- paste0("chr",cfc_output$Chromosome)
length <- as.numeric(cfc_output$Start[2] - cfc_output$Start[1])
cfc_output$End <- as.numeric((cfc_output$Start + length) - 1)
cfc_output <- cfc_output[,c(1,2,6,3:5)]
gr_cfcout <- toGRanges(cfc_output)

plot.params <- getDefaultPlotParams(plot.type=5)
plot.params$leftmargin <- 0.01
plot.params$rightmargin <- 0.01
kp <- plotKaryotype(genome = "hg19", plot.type = 5, labels.plotter = NULL, plot.params = plot.params)
kpAddChromosomeNames(kp, cex=0.7)
kpDataBackground(kp, data.panel = 1, color = c("grey", "lightgrey"))
kpAxis(kp, ymin = -10, ymax = 10, cex=0.5, tick.pos = c(-10, -9, -8, -7, -6, -5, -4 , -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), r0=1, r1=0)
kpAbline(kp, h = 0.5, col ="black")
kpPoints(kp, data = gr_cfcout[gr_cfcout$Ratio != -1], y = gr_cfcout[gr_cfcout$Ratio != -1]$Ratio, data.panel = 1, r0 = 1, r1 = 0, ymin = -10, ymax = 10, pch = 1)
kpLines(kp, data = gr_cfcout[gr_cfcout$Ratio != -1], y = gr_cfcout[gr_cfcout$Ratio != -1]$MedianRatio, data.panel = 1, r0 = 1, r1 = 0, ymin = -10, ymax = 10, col = "red")

library(RColorBrewer)
brewer.pal(length(levels(gr_input$event)), "RdBu")

colours <- data.frame(
  levels = factor(unique(gr_input$event), levels = c("HOMD", "HETD", "NEUT", "GAIN", "AMP", "HLAMP", "HLAMP1", "HLAMP2", "HLAMP3"))
)
colours$levels <- sort(colours$levels)
colours$colour <- likertColorBrewer(nc = length(unique(colours$levels)), ReferenceZero = as.numeric(grep("NEUT", colours$levels)), BrewerPaletteName = "RdBu")


#this works!
library(HH)
likertColorBrewer(nc = 9, ReferenceZero = 3, BrewerPaletteName = "RdBu")
likertColorBrewer(nc = length(unique(colours$levels)), ReferenceZero = as.numeric(grep("NEUT", colours$levels)), BrewerPaletteName = "RdBu")



copy_number <- factor(unique(gr_input$event), levels = c("HOMD", "HETD", "NEUT", "GAIN", "AMP", "HLAMP", "HLAMP1", "HLAMP2", "HLAMP3"))
copy_number <- reorder(colours$levels, c("HOMD", "HETD", "NEUT", "GAIN", "AMP", "HLAMP", "HLAMP1", "HLAMP2", "HLAMP3"))

for (i in brewer.pal(length(levels(gr_input$event)), "RdBu")) {
  print(i)
}

kpPoints(kp, data = gr_GAIN, y = gr_GAIN$X1601104.NMB866.T.sort.bam.logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#bae4bc")
kpPoints(kp, data = gr_HETD, y = gr_HETD$X1601104.NMB866.T.sort.bam.logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#43a2ca")
kpPoints(kp, data = gr_HLAMP3, y = gr_HLAMP3$X1601104.NMB866.T.sort.bam.logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#0868ac")
kpPoints(kp, data = gr_NEUT, y = gr_NEUT$X1601104.NMB866.T.sort.bam.logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, col="#f0f9e8")
kpAxis(kp, ymin = -2, ymax = 2, cex=0.5, tick.pos = c(-2, -1, 0, 1, 2), r0=1, r1=0)
kpAbline(kp, h = 0.5, col ="black")

kpPoints(kp, data = gr_input, y = gr_input$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4)


gr_input <- toGRanges(input)
list_df = do.call("list", mget(grep("gr", ls(), value = T)))
#use the assign function to create seperate dataframes for all the different possible events
for (i in levels(input$X1601104.NMB866.T.sort.bam.event)) {
  assign(paste0("gr_",i), toGRanges(input[which(input$X1601104.NMB866.T.sort.bam.event == (i)),]))
}

#create a dataframe of types and their colours
types <- data.frame(
  types = levels(gr_input$event),
  colours = brewer.pal(length(levels(gr_input$event)), "RdBu")
)

for (i in levels(gr_input$event)) {
  kpPoints(kp, data = gr_input[gr_input$event == (i)], y = gr_input[gr_input$event == (i)]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4, alpha = gr_input[gr_input$event == (i)]$logR)
}


for (t in types$types) {
  for (c in types$colours) {
    kpPoints(kp, data = gr_input[gr_input$event == (t)], y = gr_input[gr_input$event == (t)]$logR, data.panel = 1, r0 = 1, r1 = 0, ymin = -2, ymax = 2, cex = 0.4)
  }
}



##to combine geneCN and lcWGS and plot together with markers

library(karyoploteR)

gene_cn <- read.delim("1808530-UK0542WC-T_log2ratios.txt", sep = ",")
gr_gene_cn <- toGRanges(gene_cn)
gene_cn_markers <- gr_gene_cn[gr_gene_cn$feature != "background"]
gr_gene_cn_markers <- unlist(range(split(gene_cn_markers, ~feature)))


lcwgs_cn <- read.delim("1808530-SIOPEN-T.cna.seg")
lcwgs_cn$chr <- paste0("chr",lcwgs_cn$chr)
gr_lcwgs_cn <- toGRanges(lcwgs_cn)


plot.params <- getDefaultPlotParams(plot.type=3)
plot.params$leftmargin <- 0.02
plot.params$rightmargin <- 0.01

kp <- plotKaryotype(genome = "hg19", plot.type = 3, labels.plotter = NULL, plot.params = plot.params, chromosomes = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))

kp <- plotKaryotype(genome = "hg19", plot.type = 3, labels.plotter = NULL, plot.params = plot.params, chromosomes = "chr3")

kpDataBackground(kp, data.panel = 1, color = c("grey", "lightgrey"))
kpDataBackground(kp, data.panel = 2, color = c("grey", "lightgrey"))
kpAddChromosomeNames(kp, cex=0.7, xoffset = -212.5)
kpAxis(kp, data.panel = 1, ymin = -5, ymax = 5, cex=0.5, tick.pos = c( -5, -4 , -3, -2, -1, 0, 1, 2, 3, 4, 5))
kpAddLabels(kp, labels = "Log fold change", srt=90, pos = 3, data.panel = 1, cex = 0.7)
kpAbline(kp, h = 0.5, col ="darkgrey", data.panel = 1)

kpPoints(kp, data = gr_gene_cn[gr_gene_cn$feature == "background"], y = gr_gene_cn[gr_gene_cn$feature == "background"]$log2ratio, data.panel = 1, ymin = -5, ymax = 5, pch = 21, cex = 0.3, col="#00000050")
kpPoints(kp, data = gr_gene_cn[gr_gene_cn$feature == "other"], y = gr_gene_cn[gr_gene_cn$feature == "other"]$log2ratio, data.panel = 1, ymin = -5, ymax = 5, pch = 21, cex = 0.3, col="#00000050")
#kpPoints(kp, data = gr_gene_cn[gr_gene_cn$feature != "background"], y = gr_gene_cn[gr_gene_cn$feature != "background"]$log2ratio, data.panel = 1, ymin = -5, ymax = 5, pch = 21, cex = 0.3, col="#ff000050")

#plot DNAcopy
#kpLines(kp, chr = DNAcopy_gene.cn$chrom, x = DNAcopy_gene.cn$loc.start, y = DNAcopy_gene.cn$seg.mean, data.panel = 1, col="red", ymin = -5, ymax = 5)
#kpPoints(kp, chr = DNAcopy_gene.cn$chrom, x = DNAcopy_gene.cn$loc.start, y = DNAcopy_gene.cn$seg.mean, data.panel = 1, col="red", ymin = -5, ymax = 5)
#kpPoints(kp, data = gr_DNAcopy_gene.cn, y = gr_DNAcopy_gene.cn$seg.mean, data.panel = 1, col="red", ymin = -5, ymax = 5)
#kpLines(kp, data = gr_DNAcopy_gene.cn, y = gr_DNAcopy_gene.cn$seg.mean, data.panel = 1, col = "red", ymin = -5, ymax = 5)
#kpRect(kp, data = gr_DNAcopy_gene.cn, y0 = 0, y1 = gr_DNAcopy_gene.cn$seg.mean, data.panel = 1, ymin = -5, ymax = 5)

#Base segments colour
#kpSegments(kp, data = gr_DNAcopy_gene.cn, y0 = gr_DNAcopy_gene.cn$seg.mean, y1 = gr_DNAcopy_gene.cn$seg.mean, data.panel = 1, ymin = -5, ymax = 5, col="red", lwd=1.5)
kpSegments(kp, data = gr_DNAcopy_gene.cn[gr_DNAcopy_gene.cn$seg.mean >-0.5 & gr_DNAcopy_gene.cn$seg.mean < 2.0], y0 = gr_DNAcopy_gene.cn[gr_DNAcopy_gene.cn$seg.mean >-0.5 & gr_DNAcopy_gene.cn$seg.mean < 2.0]$seg.mean, y1 = gr_DNAcopy_gene.cn[gr_DNAcopy_gene.cn$seg.mean >-0.5 & gr_DNAcopy_gene.cn$seg.mean < 2.0]$seg.mean, data.panel = 1, ymin = -5, ymax = 5, col="#009e74", lwd=1.5)

#colour segments depending on gain/loss
kpSegments(kp, data = gr_DNAcopy_gene.cn[gr_DNAcopy_gene.cn$seg.mean < -0.5], y0 = gr_DNAcopy_gene.cn[gr_DNAcopy_gene.cn$seg.mean < -0.5]$seg.mean, y1 = gr_DNAcopy_gene.cn[gr_DNAcopy_gene.cn$seg.mean < -0.5]$seg.mean, data.panel = 1, ymin = -5, ymax = 5, col="red", lwd=1.5)
kpSegments(kp, data = gr_DNAcopy_gene.cn[gr_DNAcopy_gene.cn$seg.mean > 1], y0 = gr_DNAcopy_gene.cn[gr_DNAcopy_gene.cn$seg.mean > 1]$seg.mean, y1 = gr_DNAcopy_gene.cn[gr_DNAcopy_gene.cn$seg.mean > 1]$seg.mean, data.panel = 1, ymin = -5, ymax = 5, col="green", lwd=1.5)

overlaps <- findOverlaps(gr_gene_cn_markers, gr_DNAcopy_gene.cn)
gr_gene_cn_overlaps <- gr_gene_cn_markers[queryHits(overlaps)]
mcols(gr_gene_cn_overlaps) <- cbind.data.frame(
  mcols(gr_gene_cn_overlaps),
  mcols(gr_DNAcopy_gene.cn[subjectHits(overlaps)])
  )

signif_loss <- gr_gene_cn_overlaps[gr_gene_cn_overlaps$seg.mean < -0.5]
signif_gain <- gr_gene_cn_overlaps[gr_gene_cn_overlaps$seg.mean > 2.0]

kpPlotMarkers(kp, data = gr_gene_cn_markers, labels = gr_gene_cn_markers@ranges@NAMES, adjust.label.position = F, text.orientation = "vertical", data.panel = 1, cex = 0.7, r1 = 1.35, lty=0)
#kpSegments(kp, chr = as.character(seqnames(gr_gene_cn_markers)), x0 = start(gr_gene_cn_markers), x1 = start(gr_gene_cn_markers), y0 = gr_DNAcopy_gene.cn$seg.mean, y1 = 5, ymin = -5, ymax = 5)
kpSegments(kp, chr = as.character(seqnames(gr_gene_cn_overlaps)), x0 = start(gr_gene_cn_overlaps), x1 = start(gr_gene_cn_overlaps), y0 = gr_gene_cn_overlaps$seg.mean, y1 = 5, ymin = -5, ymax = 5, lty=2, col="#0000ff50")

#lines to just the deletions
#kpSegments(kp, chr = as.character(seqnames(gr_gene_cn_overlaps[gr_gene_cn_overlaps$seg.mean < -0.5])), x0 = start(gr_gene_cn_overlaps[gr_gene_cn_overlaps$seg.mean < -0.5]), x1 = start(gr_gene_cn_overlaps[gr_gene_cn_overlaps$seg.mean < -0.5]), y0 = gr_gene_cn_overlaps[gr_gene_cn_overlaps$seg.mean < -0.5]$seg.mean, y1 = 5, ymin = -5, ymax = 5, lty=2, col="#0000ff50")



kpPlotMarkers(kp, data = gr_gene_cn_overlaps[gr_gene_cn_overlaps$seg.mean < -0.5], labels = names(gr_gene_cn_overlaps[gr_gene_cn_overlaps$seg.mean < -0.5]), text.orientation = "vertical", data.panel = 1, cex = 0.7, r1 = 1.35, line.color="#ff000050", lty=1, label.color = "#ff0000")
kpPlotMarkers(kp, data = gr_gene_cn_overlaps[gr_gene_cn_overlaps$seg.mean > 1], labels = names(gr_gene_cn_overlaps[gr_gene_cn_overlaps$seg.mean > 1]), text.orientation = "vertical", data.panel = 1, cex = 0.7, r1 = 1.35, line.color="#0072b250", lty=1, label.color = "#0072b2")

kpPlotMarkers(kp, data = gr_gene_cn_markers, labels = gr_gene_cn_markers@ranges@NAMES, adjust.label.position = T, text.orientation = "vertical", data.panel = 2, cex = 0.7, label.dist = 0.0000000001, max.iter = 1000000, r1 = 1.35, lty=3)

kpAxis(kp, data.panel = 2, ymin = -5, ymax = 5, cex=0.5, tick.pos = c( -5, -4 , -3, -2, -1, 0, 1, 2, 3, 4, 5))
kpAddLabels(kp, labels = "Log fold change", srt=90, pos = 3, data.panel = 2, cex = 0.7)
kpAbline(kp, h = 0.5, col ="darkgrey", data.panel = 2)
#kpPoints(kp, data = gr_lcwgs_cn, y = lcwgs_cn$X1808530.SIOPEN.T.logR, data.panel = 2, ymin = -5, ymax = 5, pch = 21, cex = 0.3, r0 = 1, r1 = 0)
kpPoints(kp, data = gr_lcwgs_cn[gr_lcwgs_cn$X1808530.SIOPEN.T.Corrected_Call == "NEUT"], y = gr_lcwgs_cn[gr_lcwgs_cn$X1808530.SIOPEN.T.Corrected_Call == "NEUT"]$X1808530.SIOPEN.T.logR, data.panel = 2, ymin = -5, ymax = 5, pch = 21, cex = 0.3, r0 = 1, r1 = 0, col="#70707050")
kpPoints(kp, data = gr_lcwgs_cn[gr_lcwgs_cn$X1808530.SIOPEN.T.Corrected_Call == "HETD"], y = gr_lcwgs_cn[gr_lcwgs_cn$X1808530.SIOPEN.T.Corrected_Call == "HETD"]$X1808530.SIOPEN.T.logR, data.panel = 2, ymin = -5, ymax = 5, pch = 21, cex = 0.3, r0 = 1, r1 = 0, col="#ff0000")
kpPoints(kp, data = gr_lcwgs_cn[gr_lcwgs_cn$X1808530.SIOPEN.T.Corrected_Call == "AMP"], y = gr_lcwgs_cn[gr_lcwgs_cn$X1808530.SIOPEN.T.Corrected_Call == "AMP"]$X1808530.SIOPEN.T.logR, data.panel = 2, ymin = -5, ymax = 5, pch = 21, cex = 0.3, r0 = 1, r1 = 0, col="#009e74")
kpPoints(kp, data = gr_lcwgs_cn[gr_lcwgs_cn$X1808530.SIOPEN.T.Corrected_Call == "GAIN"], y = gr_lcwgs_cn[gr_lcwgs_cn$X1808530.SIOPEN.T.Corrected_Call == "GAIN"]$X1808530.SIOPEN.T.logR, data.panel = 2, ymin = -5, ymax = 5, pch = 21, cex = 0.3, r0 = 1, r1 = 0, col="#0072b2")




#KP gene expression results follow

BiocManager::install("pasilla")
BiocManager::install("DESeq2")
BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library("pasilla")
library("DESeq2")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")

pasCts <- system.file("extdata", "pasilla_gene_counts.tsv", package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv", package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]

rownames(coldata) <- sub("fb", "", rownames(coldata))
cts <- cts[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "untreated")
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, coef = 2, res = res)

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
dm.genes <- genes(txdb)

mcols(dm.genes) <- res[names(dm.genes), c("log2FoldChange", "stat", "pvalue", "padj")]

library(karyoploteR)
ordered <- dm.genes[order(dm.genes$padj, na.last = T),]
kp <- plotKaryotype(genome = "dm6")
kp <- kpPlotMarkers(kp, ordered[1:10], labels = names(ordered[1:10]), text.orientation = "horizontal")


filtered.dm.genes <- dm.genes[!is.na(dm.genes$padj)]
log.pval <- -log10(filtered.dm.genes$padj)
mcols(filtered.dm.genes)$logpval <- log.pval

sign.genes <- filtered.dm.genes[filtered.dm.genes$padj < 0.05,]

kp <- plotKaryotype(genome = "dm6")
kpPoints(kp, data = sign.genes, y = sign.genes$logpval, ymax = max(sign.genes$logpval))

#calculation of the ceiling of the range of log2foldchange
fc.ymax <- ceiling(max(abs(range(sign.genes$log2FoldChange))))
fc.ymin <- -fc.ymax

kp <- plotKaryotype(genome="dm6")
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin)

cex.val <- sqrt(sign.genes$logpval)/3

kp <- plotKaryotype(genome="dm6")
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex = cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin)

top.genes <- ordered[1:20]
kp <- plotKaryotype(genome="dm6")
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex = cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal")


points.top <- 0.8
kp <- plotKaryotype(genome="dm6")
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top)


kp <- plotKaryotype(genome="dm6")
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top)

gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2

kpSegments(kp, chr = as.character(seqnames(top.genes)), x0 = gene.mean, x1 = gene.mean, y0 = top.genes$log2FoldChange, y1 = fc.ymax, ymax = fc.ymax, ymin = fc.ymin, r1 = points.top)


#use only the list of genes mike gave:
#MET, PDGFRA/KIT, CDKN2A/B, CCND1, CDK4, MYCN, TERT, ATRX, TSC1, TSC2, EGFR, TP53, ATM, NF1, RB1, CCND2, MDM2, MDM4, KRAS, PIK3CA, IGF1R, FGFR1, FGFR2, AKT1, PTEN, EZH2, CIC, SMARCA4, SMARCB1, IGF2, DROSHA, CDK4, CDK6, BCOR

list.of.genes <- c(c("MET",
                     "PDGFRA",
                     "KIT",
                     "CDKN2A",
                     "CDKN2B",
                     "CCND1",
                     "CDK4",
                     "MYCN",
                     "TERT",
                     "ATRX",
                     "TSC1",
                     "TSC2",
                     "EGFR",
                     "TP53",
                     "ATM",
                     "NF1",
                     "RB1", 
                     "CCND2",
                     "MDM2",
                     "MDM4",
                     "KRAS",
                     "PIK3CA",
                     "IGF1R",
                     "FGFR1",
                     "FGFR2",
                     "AKT1",
                     "PTEN",
                     "EZH2",
                     "CIC",
                     "SMARCA4",
                     "SMARCB1",
                     "IGF2",
                     "DROSHA",
                     "CDK4",
                     "CDK6",
                     "BCOR"))


#this works!!
subset(gr.ichor.overlaps, names(gr.ichor.overlaps) %in% list.of.genes)















