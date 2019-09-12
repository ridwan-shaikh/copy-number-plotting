##Script to plot both gene.cn and ichorCNA data together on the same plot using lcWGS

library(DNAcopy)
library(karyoploteR)

#load in results from gene.cn for circular binary segmentation

gene.cn.file <- list.files(path = ".", pattern = "log2ratios.txt")
gene.cn.output <- read.delim(gene.cn.file, sep = ",")

#get sampleID from the filename
sampleID <- unlist(strsplit(x = gene.cn.file,split = "_"))[1]

#create a DNAcopy object from the output file
gene.cn.object <- CNA(cbind(gene.cn.output$log2ratio), gene.cn.output$chr, gene.cn.output$start, data.type = "logratio", sampleid = sampleID)
segment.gene.cn.object <- segment(gene.cn.object)
#plot(segment.gene.cn.object, plot.type="w")

#convert to dataframe and then genomicranges object for karyoploteR
df.segment.gene.cn.ojbect <- as.data.frame(segment.gene.cn.object$output)
#drop the ID row
df.segment.gene.cn.ojbect[,1] <- NULL
gr.segment.gene.cn <- toGRanges(df.segment.gene.cn.ojbect)
#sort by genomic locations
gr.segment.gene.cn <- sortSeqlevels(gr.segment.gene.cn)
gr.segment.gene.cn <- sort(gr.segment.gene.cn)

#convert the gene.cn output to a genomic ranges to get the markers
gr.gene.cn <- toGRanges(gene.cn.output)
gr.gene.cn.markers <- gr.gene.cn[gr.gene.cn$feature != "background"]
gr.gene.cn.markers <- unlist(range(split(gr.gene.cn.markers, ~feature)))

#set the plotting parameters for karyoploteR
plot.params <- getDefaultPlotParams(plot.type = 3)
plot.params$leftmargin <- 0.015
plot.params$rightmargin < 0.001

#plot all chromosomes apart from chrY
kp <- plotKaryotype(genome = "hg19", plot.type = 3, labels.plotter = NULL, plot.params = plot.params, chromosomes = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))
#Add background colours
kpDataBackground(kp, data.panel = 1, color = c("grey", "lightgrey"))
kpDataBackground(kp, data.panel = 2, color = c("grey", "lightgrey"))
#Add chromosome names
kpAddChromosomeNames(kp, cex=0.7, xoffset = -212.5)
#Addd axis and labels
kpAxis(kp, data.panel = 1, ymin = -5, ymax = 5, cex=0.5, tick.pos = c(-5,-4,-3,-2,-1,0,1,2,3,4,5))
kpAxis(kp, data.panel = 2, ymin = -5, ymax = 5, cex=0.5, tick.pos = c(-5,-4,-3,-2,-1,0,1,2,3,4,5), r0 = 1, r1 = 0)
kpAddLabels(kp, labels = "log fold change", srt = 90, pos = 3, data.panel = 1, cex = 0.7)
kpAddLabels(kp, labels = "log fold change", srt = 90, pos = 3, data.panel = 2, cex = 0.7)
#Add abline
kpAbline(kp, h = 0.5, col="darkgrey", data.panel = 1)
kpAbline(kp, h = 0.5, col="darkgrey", data.panel = 2)


#Plot gene.cn data
#Plot only the probes labelled as background and other (other = chrX) as the data is too messy with the gene probes ontop in the WGS view
kpPoints(kp,
         data = gr.gene.cn[gr.gene.cn$feature == "background"],
         y = gr.gene.cn[gr.gene.cn$feature == "background"]$log2ratio,
         data.panel = 1,
         ymin = -5,
         ymax = 5,
         pch = 21,
         cex = 0.4,
         col= "#00000050")

kpPoints(kp,
         data = gr.gene.cn[gr.gene.cn$feature == "other"],
         y = gr.gene.cn[gr.gene.cn$feature == "other"]$log2ratio,
         data.panel = 1,
         ymin = -5,
         ymax = 5,
         pch = 21,
         cex = 0.4,
         col= "#00000050")

#Plot the segmention lines from CBS based on the following cut off: losses = < 0.5, gains= > 2.4
kpSegments(kp,
           data = gr.segment.gene.cn[gr.segment.gene.cn$seg.mean > -0.5 & gr.segment.gene.cn$seg.mean < 2.4],
           y0 = gr.segment.gene.cn[gr.segment.gene.cn$seg.mean > -0.5 & gr.segment.gene.cn$seg.mean < 2.4]$seg.mean,
           y1 = gr.segment.gene.cn[gr.segment.gene.cn$seg.mean > -0.5 & gr.segment.gene.cn$seg.mean < 2.4]$seg.mean,
           data.panel = 1,
           ymin = -5,
           ymax = 5,
           col="red",
           lwd=1.5)


#Plot losses
kpSegments(kp, 
           data = gr.segment.gene.cn[gr.segment.gene.cn$seg.mean < -0.5],
           y0 = gr.segment.gene.cn[gr.segment.gene.cn$seg.mean < -0.5]$seg.mean,
           y1 = gr.segment.gene.cn[gr.segment.gene.cn$seg.mean < -0.5]$seg.mean,
           data.panel = 1,
           ymin = -5,
           ymax = 5,
           col = "#0071b2",
           lwd = 1.5) #"#0071b2" = blue

#Plot gains
kpSegments(kp, 
           data = gr.segment.gene.cn[gr.segment.gene.cn$seg.mean > 2.4],
           y0 = gr.segment.gene.cn[gr.segment.gene.cn$seg.mean > 2.4]$seg.mean,
           y1 = gr.segment.gene.cn[gr.segment.gene.cn$seg.mean > 2.4]$seg.mean,
           data.panel = 1,
           ymin = -5,
           ymax = 5,
           col = "#009e74",
           lwd = 1.5) ##009e74 = green


#find ichor files
ichor.cna.file <- list.files(".", pattern = "cna.seg")
ichor.cna <- read.delim(ichor.cna.file)

#add chr character to convert to genomic ranges
ichor.cna$chr <- paste0("chr", ichor.cna$chr)
gr.ichor.cna <- toGRanges(ichor.cna)


#kpPoints(kp, data = gr.ichor.cna, y = lcwgs_cn$X1808530.SIOPEN.T.logR, data.panel = 2, ymin = -5, ymax = 5, pch = 21, cex = 0.3, r0 = 1, r1 = 0)
kpPoints(kp, data = gr.ichor.cna[gr.ichor.cna$X1808530.SIOPEN.T.Corrected_Call == "NEUT"], y = gr.ichor.cna[gr.ichor.cna$X1808530.SIOPEN.T.Corrected_Call == "NEUT"]$X1808530.SIOPEN.T.logR, data.panel = 2, ymin = -5, ymax = 5, pch = 21, cex = 0.4, r0 = 1, r1 = 0, col="#70707050")
kpPoints(kp, data = gr.ichor.cna[gr.ichor.cna$X1808530.SIOPEN.T.Corrected_Call == "HETD"], y = gr.ichor.cna[gr.ichor.cna$X1808530.SIOPEN.T.Corrected_Call == "HETD"]$X1808530.SIOPEN.T.logR, data.panel = 2, ymin = -5, ymax = 5, pch = 21, cex = 0.4, r0 = 1, r1 = 0, col="#0071b2")
kpPoints(kp, data = gr.ichor.cna[gr.ichor.cna$X1808530.SIOPEN.T.Corrected_Call == "AMP"], y = gr.ichor.cna[gr.ichor.cna$X1808530.SIOPEN.T.Corrected_Call == "AMP"]$X1808530.SIOPEN.T.logR, data.panel = 2, ymin = -5, ymax = 5, pch = 21, cex = 0.4, r0 = 1, r1 = 0, col="#009e74")
kpPoints(kp, data = gr.ichor.cna[gr.ichor.cna$X1808530.SIOPEN.T.Corrected_Call == "GAIN"], y = gr.ichor.cna[gr.ichor.cna$X1808530.SIOPEN.T.Corrected_Call == "GAIN"]$X1808530.SIOPEN.T.logR, data.panel = 2, ymin = -5, ymax = 5, pch = 21, cex = 0.4, r0 = 1, r1 = 0, col="#e7298a")

