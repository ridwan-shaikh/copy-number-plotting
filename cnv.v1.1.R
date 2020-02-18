## script to analysis panel CNV data
'
Description - CNV Pipeline to plot both genome-wide and focal CNV plots for targeted 
  panels. In order to detect larger chromosomal abberations as well as focal changes. 
  Both plots originate from the same normalised data using different plotting logic 
  to extrapolate valuable data for either focal or genome wide calling.

Author: Sabri Jamal, Lina Yuan, Ridwan Shaikh
Date: 14/10/2019
updated by Lina 19/10/2019

'

## normalized_coverage from Picard: The ratio of "coverage" in this GC bin vs. the mean coverage of all GC bins.
## A number of 1 represents mean coverage, a number less than one represents lower than mean coverage (e.g. 0.5 means half as much coverage as average)
### while a number greater than one represents higher than mean coverage (e.g. 3.1 means this GC bin has 3.1 times more reads per window than average).


library(optparse)
library(DNAcopy)
library(affy)
library(GenomicRanges)
library(gtools)
library(ggplot2)
library(plotly)
library(tidyverse)


option_list <- list(
  make_option(c("--samplecoverage"), type = "character", help = "Path to tumour probe/bait coverage file. Required."),
  make_option(c("--refcoverage"), type = "character", default=NULL, help = "Path to panel reference coverage file. Required"),
  make_option(c("--id"), type = "character", default="test", help = "Patient ID. Default: [%default]"),
  make_option(c("--cnvdir"), type = "character", help = "Path to output dir")
 )
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

options(stringsAsFactors=FALSE)
options(bitmapType='cairo')


sampid <- opt$id
tumour_file <- opt$samplecoverage
ref_file <- opt$refcoverage
cnvdir <-opt$cnvdir
outputfile=paste('/', cnvdir, '/', sampid, '.cnv.tsv', sep='')
outputhtml=paste('/', cnvdir, '/', sampid, '.cnv.html', sep='')


sampid <- "2000520-2076608-T"
tumour_file <-"2000520-2076608-T.qc.coverage"
ref_file<-"RMH200.male.cnv.ref"


tumour_depth=read.table(tumour_file,header=T,sep='\t', stringsAsFactors =F, check.names=F, na.strings = "NA")
ref_depth=read.table(ref_file,header=T,sep=',', stringsAsFactors =F, check.names=F, na.strings = "NA")
ref_depth$order=c(1:nrow(ref_depth))
data=merge(tumour_depth[,c("chrom", "start", "end", "name", "%gc", "normalized_coverage")], ref_depth, by="name")
data=data[order(data$order),]
colnames(data)[5]='gc'

## filling in zero
data$median[data$median == 0]=1
data$nratio=data$normalized_coverage/data$median
fit=loess(nratio ~ gc, data=data)
bias=predict(fit, data.frame(gc=data$gc))
data$gcnratio=data$nratio-bias
#data$probe <- stringr::str_remove(data$name, "_[a-z]+")

#write data for webserver display
deletions <- data %>%
  select(name, chrom, start, end, gcnratio) %>%
  separate_rows(name, sep = "\\|") %>%
  filter(str_detect(name, "_mut_")) %>%
  mutate(gene = str_remove(name, "_mut")) %>%
  separate(gene, into = c("gene", "probe_no"), sep = "_") %>%
  filter(gcnratio <= -0.5) %>%
  arrange(as.numeric(probe_no)) %>%
  group_by(gene) %>%
  summarise(exons_deleted = toString(probe_no)) %>%
  ungroup()

amplifications <- data %>%
  select(name, chrom, start, end, gcnratio) %>%
  separate_rows(name, sep = "\\|") %>%
  filter(str_detect(name, "_mut_")) %>%
  mutate(gene = str_remove(name, "_mut")) %>%
  separate(gene, into = c("gene", "probe_no"), sep = "_") %>%
  filter(gcnratio >= 1) %>%
  arrange(as.integer(probe_no)) %>%
  group_by(gene) %>%
  summarise(exons_amplified = toString(probe_no)) %>%
  ungroup()

genes <- data %>%
  select(name) %>%
  separate_rows(name, sep = "\\|") %>%
  filter(str_detect(name, "_mut_")) %>%
  separate(name, into = c("gene", "probe_no"), sep = "_mut") %>%
  distinct(gene)

genes %>%
  left_join(deletions, "gene") %>%
  left_join(amplifications, "gene") %>%
  arrange(gene) %>%
  rename(GENE = gene, DELEXON = exons_deleted, AMPEXON = exons_amplified) %>%
  write.table(outputfile, quote = F, sep = "\t", row.names = F)

#Use DNAcopy to segment the log2 ratios
DNAcopy.obj <- DNAcopy::CNA(cbind(data$gcnratio), data$chrom, data$start, data.type = "logratio")
DNAcopy.obj.seg <- DNAcopy::segment(DNAcopy.obj)

#retrieve the segmentation values as a dataframe and drop 1st column of sampleID
df.DNAcopy.seg <- base::as.data.frame(DNAcopy.obj.seg$output)
df.DNAcopy.seg$ID <- NULL

#covert the input file into a genomicranges object and sort to ensure it is in the correct order
gr.input <- GRanges(data)
gr.input <- sortSeqlevels(gr.input)
gr.input <- sort(gr.input)

#covert the DNAcopy df into genomicranges object and sort to ensure it is in the correct order
gr.DNAcopy.seg <- GRanges(df.DNAcopy.seg)
gr.DNAcopy.seg <- sortSeqlevels(gr.DNAcopy.seg)
gr.DNAcopy.seg <- sort(gr.DNAcopy.seg)

#annotate segmentation data with cytobands and output for webserver display

#read in cytoband bedfile
cytobands <- read.delim("/mnt/scratch/DMP/DUDMP/TRANSGEN/transgen-mdx/ngs/7.resources/hg19/bedfiles/cytoband.bed")

#covert to GenomicRanges object
gr.cytobands <- GRanges(cytobands)

#find overlaps between the segmentation data and the cytoband information
overlaps <- findOverlaps(gr.DNAcopy.seg, gr.cytobands)
gr.cytobands.overlaps <- gr.DNAcopy.seg[queryHits(overlaps)]
mcols(gr.cytobands.overlaps) <- cbind.data.frame(
  mcols(gr.cytobands.overlaps),
  mcols(gr.cytobands[subjectHits(overlaps)])
)

#covert to dataframe
df.cytobands.overlaps <- as.data.frame(gr.cytobands.overlaps)

#collapse the cytoband information into a range
group.cytobands <- df.cytobands.overlaps %>%
  group_by(seqnames, start, end, num.mark, seg.mean) %>%
  summarise(cytobands = toString(name)) %>%
  ungroup()
group.cytobands$cytoband1 <- as.character(lapply(strsplit(as.character(group.cytobands$cytobands), split = ","), "[", 1))
group.cytobands$cytoband2 <- as.character(lapply(strsplit(as.character(group.cytobands$cytobands), split = ","), tail, n=1))
group.cytobands$cytoband2 <- trimws(group.cytobands$cytoband2, which = "left")
group.cytobands$cytoband <- ifelse(
  as.character(group.cytobands$cytoband1) == as.character(group.cytobands$cytoband2),
  paste0(group.cytobands$cytoband1), ifelse(
    as.character(group.cytobands$cytoband1) < as.character(group.cytobands$cytoband2),
    paste0(group.cytobands$cytoband1, "-", group.cytobands$cytoband2),
    paste0(group.cytobands$cytoband2, "-", group.cytobands$cytoband1)
  ))

#write the annotated segmentation data
group.cytobands %>%
  select(seqnames, start, end, cytoband, num.mark, seg.mean) %>%
  rename("Chromsome" = seqnames, "Start" = start, "End" = end, "Number_of_Probes" = num.mark, "Segmentation_Mean" = seg.mean) %>%
  write.table(paste0(outputfile, ".seg"), row.names = F, quote = F)

#intersect the segmentation and input file to get a segmented mean value per probe
overlaps <- findOverlaps(gr.input, gr.DNAcopy.seg)
gr.overlaps <- gr.input[queryHits(overlaps)]
mcols(gr.overlaps) <- cbind.data.frame(
  mcols(gr.overlaps),
  mcols(gr.DNAcopy.seg[subjectHits(overlaps)])
)

#convert to resulting intersect to a dataframe and prepare for plotting
df.plot <- as.data.frame(gr.overlaps)

#change num.mark column to a sequence of numbers to ensure the plot it set out correctly
df.plot$num.mark <- seq(1:nrow(df.plot))

#add colours to the fields depending on chromosome
df.plot$chr <- stringr::str_extract(as.character(df.plot$seqnames), "[0-9]+")
df.plot$colours <- ifelse(gtools::odd(as.numeric(df.plot$chr)) | df.plot$seqnames == "chrX", "#253494", "#74add1")

#ggplot removes data points if they are outside ylim
#to prevent this: new column named scale and set max log2ratio to 4.
df.plot$scale <- ifelse(df.plot$gcnratio >= 4, 4, df.plot$gcnratio)
df.plot$scale.seg.mean <- ifelse(df.plot$seg.mean >= 4, 4, df.plot$seg.mean)


#plot!
g <- ggplot(df.plot, aes(x = num.mark, y = scale, text = paste0(df.plot$name,
                                                                '<br>', df.plot$seqnames,
                                                                '<br>log2 ratio: ', as.numeric(round(gcnratio, 3)),
								                                                '<br>Segmented mean: ', as.numeric(round(seg.mean, 3))))) +
  geom_hline(yintercept = c(-1, 0, 1), colour = c("darkgrey", "lightgrey", "darkgrey"), linetype = c("longdash", "solid", "longdash")) +
  geom_point(colour = df.plot$colours, size = 0.8, show.legend = F) +
  scale_y_continuous(breaks = seq(-4,4,2), limits = c(-4,4), labels = c("Min", -2, 0, 2, "Max")) +
  ylab('Log2 Ratio of Normalised Depths') +
  geom_crossbar(aes(y = df.plot$scale.seg.mean, ymin = df.plot$scale.seg.mean, ymax = df.plot$scale.seg.mean), colour = "darkorange") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

#covert to plotly graphic
gp <- ggplotly(g, tooltip = c("text"), width = 2000, height = 500)

#export as html
library(htmlwidgets)
htmlwidgets::saveWidget(gp, outputhtml)

