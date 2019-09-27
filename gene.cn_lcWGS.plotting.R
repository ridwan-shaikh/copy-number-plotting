##Script to plot both geneCN and ichorCNA data together on the same plot using lcWGS

library(DNAcopy)
library(karyoploteR)

#load in results from geneCN for circular binary segmentation

geneCN.file <- list.files(path = ".", pattern = "log2ratios.txt")
if (length(geneCN.file) == 0) {}

geneCN.output <- read.delim(geneCN.file, sep = ",")

#get sampleID from the filename
sampleID <- unlist(strsplit(x = geneCN.file,split = "_"))[1]

#create a DNAcopy object from the output file
geneCN.object <- CNA(cbind(geneCN.output$log2ratio), geneCN.output$chr, geneCN.output$start, data.type = "logratio", sampleid = sampleID)
segment.geneCN.object <- segment(geneCN.object)
#plot(segment.geneCN.object, plot.type="w")

#convert to dataframe and then genomicranges object for karyoploteR
df.segment.geneCN.ojbect <- as.data.frame(segment.geneCN.object$output)
#drop the ID row
df.segment.geneCN.ojbect[,1] <- NULL
gr.segment.geneCN <- toGRanges(df.segment.geneCN.ojbect)
#sort by genomic locations
gr.segment.geneCN <- sortSeqlevels(gr.segment.geneCN)
gr.segment.geneCN <- sort(gr.segment.geneCN)

#convert the geneCN output to a genomic ranges to get the markers
gr.geneCN <- toGRanges(geneCN.output)
gr.geneCN.markers <- gr.geneCN[gr.geneCN$feature != "background"]
gr.geneCN.markers <- unlist(range(split(gr.geneCN.markers, ~feature)))

#intersect to get seg.mean values onto the markers genomicrange object
overlaps <- findOverlaps(gr.geneCN.markers, gr.segment.geneCN)
gr.geneCN.overlaps <- gr.geneCN.markers[queryHits(overlaps)]
mcols(gr.geneCN.overlaps) <- cbind.data.frame(
  mcols(gr.geneCN.overlaps),
  mcols(gr.segment.geneCN[subjectHits(overlaps)])
)

#get gains and losses and a genomicrange object
gr.geneCN.markers.losses <- gr.geneCN.overlaps[gr.geneCN.overlaps$seg.mean < -0.5 & names(gr.geneCN.overlaps) != "other"]
gr.geneCN.markers.gains <- gr.geneCN.overlaps[gr.geneCN.overlaps$seg.mean > 2.4 & names(gr.geneCN.overlaps) != "other"]

#compress to have one one gene name
gr.geneCN.markers.losses <- unlist(range(split(gr.geneCN.markers.losses, ~names(gr.geneCN.markers.losses))))
gr.geneCN.markers.gains <- unlist(range(split(gr.geneCN.markers.gains, ~names(gr.geneCN.markers.gains))))

##Create the ichor dataframes

ichorCNA.file <- list.files(".", pattern = "cna.seg")
ichorCNA.output <- read.delim(ichorCNA.file[1])

#Add chr character to convert to genomicranges object
ichorCNA.output$chr <- paste0("chr", ichorCNA.output$chr)

#Strip sample names out of the headers by renaming them
colnames(ichorCNA.output)[4] <- "copy.number"
colnames(ichorCNA.output)[5] <- "event"
colnames(ichorCNA.output)[6] <- "logR"
colnames(ichorCNA.output)[7] <- "subclone.status"
colnames(ichorCNA.output)[8] <- "corrected.copy.number"
colnames(ichorCNA.output)[9] <- "corrected.call"
colnames(ichorCNA.output)[10] <- "copy.number.est"

#convert to genomicranges object
gr.ichorCNA <- toGRanges(ichorCNA.output)

#intersect the markers with the ichor dataset to plot the genes
ichor.overlaps <- findOverlaps(gr.geneCN.markers, gr.ichorCNA)
gr.ichor.overlaps <- gr.geneCN.markers[queryHits(ichor.overlaps)]
mcols(gr.ichor.overlaps) <- cbind.data.frame(
  mcols(gr.ichor.overlaps),
  mcols(gr.ichorCNA[subjectHits(ichor.overlaps)])
)

#list of known oncogenes for the lcWGS plots
list.of.genes <- c("MET",
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
                   "BCOR")

gr.ichor.markers <- subset(gr.ichor.overlaps, names(gr.ichor.overlaps) %in% list.of.genes)

##Plotting
#set the plotting parameters for karyoploteR
plot.params <- getDefaultPlotParams(plot.type = 3)
plot.params$leftmargin <- 0.015
plot.params$rightmargin <- 0.01

#plot all chromosomes apart from chrY
kp <- plotKaryotype(genome = "hg19", plot.type = 3, labels.plotter = NULL, plot.params = plot.params, chromosomes = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))
#Plot single chromosome
#kp <- plotKaryotype(genome = "hg19", plot.type = 3, labels.plotter = NULL, plot.params = plot.params, chromosomes = "chr17")
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


#Plot geneCN data
#Plot only the probes labelled as background and other (other = chrX) as the data is too messy with the gene probes ontop in the WGS view
kpPoints(kp,
         data = gr.geneCN[gr.geneCN$feature == "background"],
         y = gr.geneCN[gr.geneCN$feature == "background"]$log2ratio,
         data.panel = 1,
         ymin = -5,
         ymax = 5,
         pch = 21,
         cex = 0.4,
         col= "#70707050")

kpPoints(kp,
         data = gr.geneCN[gr.geneCN$feature == "other"],
         y = gr.geneCN[gr.geneCN$feature == "other"]$log2ratio,
         data.panel = 1,
         ymin = -5,
         ymax = 5,
         pch = 21,
         cex = 0.4,
         col= "#70707050")

#Plot the segmention lines from CBS based on the following cut off: losses = < 0.5, gains= > 2.4
kpSegments(kp,
           data = gr.segment.geneCN[gr.segment.geneCN$seg.mean > -0.5 & gr.segment.geneCN$seg.mean < 2.4],
           y0 = gr.segment.geneCN[gr.segment.geneCN$seg.mean > -0.5 & gr.segment.geneCN$seg.mean < 2.4]$seg.mean,
           y1 = gr.segment.geneCN[gr.segment.geneCN$seg.mean > -0.5 & gr.segment.geneCN$seg.mean < 2.4]$seg.mean,
           data.panel = 1,
           ymin = -5,
           ymax = 5,
           col="#00000090",
           lwd=2) #"#0071b2" = blue


#Plot losses
if (length(gr.segment.geneCN[gr.segment.geneCN$seg.mean < -0.5]) > 0) {
  kpSegments(kp,
             data = gr.segment.geneCN[gr.segment.geneCN$seg.mean < -0.5],
             y0 = gr.segment.geneCN[gr.segment.geneCN$seg.mean < -0.5]$seg.mean,
             y1 = gr.segment.geneCN[gr.segment.geneCN$seg.mean < -0.5]$seg.mean,
             data.panel = 1,
             ymin = -5,
             ymax = 5,
             col = "red",
             lwd = 2)
}

#Plot gains
if (length(gr.segment.geneCN[gr.segment.geneCN$seg.mean > 2.4]) > 0) {
  kpSegments(kp, 
             data = gr.segment.geneCN[gr.segment.geneCN$seg.mean > 2.4],
             y0 = gr.segment.geneCN[gr.segment.geneCN$seg.mean > 2.4]$seg.mean,
             y1 = gr.segment.geneCN[gr.segment.geneCN$seg.mean > 2.4]$seg.mean,
             data.panel = 1,
             ymin = -5,
             ymax = 5,
             col = "#009e74",
             lwd = 2) ##009e74 = green
}


#Plot markers
if (length(gr.geneCN.markers.losses) > 0) {
  suppressWarnings(
    kpPlotMarkers(kp,
                  data = gr.geneCN.markers.losses,
                  labels = names(gr.geneCN.markers.losses),
                  text.orientation = "vertical",
                  data.panel = 1,
                  cex = 0.5,
                  r1 = 1.35,
                  line.color = "#ff000025",
                  label.color = "#ff0000")
  )
}

if (length(gr.geneCN.markers.gains) > 0) {
  suppressWarnings(
    kpPlotMarkers(kp,
                  data = gr.geneCN.markers.gains,
                  labels = names(gr.geneCN.markers.gains),
                  text.orientation = "vertical",
                  data.panel = 1,
                  cex = 0.5,
                  r1 = 1.35,
                  line.color = "#009e7425",
                  label.color = "#009e74")
  )
}

#Plot ichorCNA data
if (length(gr.ichorCNA[gr.ichorCNA$corrected.call == "NEUT"]) > 0) {
  kpPoints(kp,
           data = gr.ichorCNA[gr.ichorCNA$corrected.call == "NEUT"],
           y = gr.ichorCNA[gr.ichorCNA$corrected.call == "NEUT"]$logR,
           data.panel = 2,
           ymin = -5,
           ymax = 5,
           pch = 21,
           cex = 0.4,
           r0 = 1,
           r1 = 0,
           col="#70707050") #70707050 = grey 50% opacity
}

if (length(gr.ichorCNA[gr.ichorCNA$corrected.call == "HETD"]) > 0) {
  kpPoints(kp,
           data = gr.ichorCNA[gr.ichorCNA$corrected.call == "HETD"],
           y = gr.ichorCNA[gr.ichorCNA$corrected.call == "HETD"]$logR,
           data.panel = 2,
           ymin = -5,
           ymax = 5,
           pch = 21,
           cex = 0.4,
           r0 = 1,
           r1 = 0,
           col="#FF0000") #FF0000 = red
}

if (length(gr.ichorCNA[gr.ichorCNA$corrected.call == "HOMD"]) > 0) {
  kpPoints(kp,
           data = gr.ichorCNA[gr.ichorCNA$corrected.call == "HOMD"],
           y = gr.ichorCNA[gr.ichorCNA$corrected.call == "HOMD"]$logR,
           data.panel = 2,
           ymin = -5,
           ymax = 5,
           pch = 21,
           cex = 0.4,
           r0 = 1,
           r1 = 0,
           col="#FF0000") #FF0000 = red
}

if (length(gr.ichorCNA[gr.ichorCNA$corrected.call == "AMP"]) > 0) {
  kpPoints(kp,
           data = gr.ichorCNA[gr.ichorCNA$corrected.call == "AMP"],
           y = gr.ichorCNA[gr.ichorCNA$corrected.call == "AMP"]$logR,
           data.panel = 2,
           ymin = -5,
           ymax = 5,
           pch = 21,
           cex = 0.4,
           r0 = 1,
           r1 = 0,
           col="#009e74") ##009e74 = green
}

if (length(gr.ichorCNA[gr.ichorCNA$corrected.call == "GAIN"]) > 0) {
  kpPoints(kp,
           data = gr.ichorCNA[gr.ichorCNA$corrected.call == "GAIN"],
           y = gr.ichorCNA[gr.ichorCNA$corrected.call == "GAIN"]$logR,
           data.panel = 2,
           ymin = -5,
           ymax = 5,
           pch = 21,
           cex = 0.4,
           r0 = 1,
           r1 = 0,
           col="#009e74") ##009e74 = green
}

if (length(unique(gr.ichor.markers[gr.ichor.markers$corrected.call == "HETD"])) > 0) {
  kpPlotMarkers(kp,
                data = unique(gr.ichor.markers[gr.ichor.markers$corrected.call == "HETD"]),
                labels = unique(names(gr.ichor.markers[gr.ichor.markers$corrected.call == "HETD"])),
                text.orientation = "vertical",
                data.panel = 2,
                cex = 0.5,
                r0 = 0,
                r1 = 1.35,
                line.color = "#ff000025",
                label.color = "#ff0000")
}

if (length(unique(gr.ichor.markers[gr.ichor.markers$corrected.call == "HOMD"])) > 0) {
  kpPlotMarkers(kp,
                data = unique(gr.ichor.markers[gr.ichor.markers$corrected.call == "HOMD"]),
                labels = unique(names(gr.ichor.markers[gr.ichor.markers$corrected.call == "HOMD"])),
                text.orientation = "vertical",
                data.panel = 2,
                cex = 0.5,
                r0 = 0,
                r1 = 1.35,
                line.color = "#ff000025",
                label.color = "#ff0000")
}

if (length(unique(gr.ichor.markers[gr.ichor.markers$corrected.call == "AMP"])) > 0) {
  kpPlotMarkers(kp,
                data = unique(gr.ichor.markers[gr.ichor.markers$corrected.call == "AMP"]),
                labels = unique(names(gr.ichor.markers[gr.ichor.markers$corrected.call == "AMP"])),
                text.orientation = "vertical",
                data.panel = 2,
                cex = 0.5,
                r0 = 0,
                r1 = 1.35,
                line.color = "#009e7425",
                label.color = "#009e74")
}

if (length(unique(gr.ichor.markers[gr.ichor.markers$corrected.call == "GAIN"])) > 0) {
  kpPlotMarkers(kp,
                data = unique(gr.ichor.markers[gr.ichor.markers$corrected.call == "GAIN"]),
                labels = names(gr.ichor.markers[gr.ichor.markers$corrected.call == "GAIN"]),
                text.orientation = "vertical",
                data.panel = 2,
                cex = 0.5,
                r0 = 0,
                r1 = 1.35,
                line.color = "#009e7425",
                label.color = "#009e74",
                label.dist = 0.00001)
}

kpAddMainTitle(kp, main = as.character(sampleID))






