library(ggdendro)
library(data.table)
library(PeakSegPipeline)
library(ggplot2)
macs.pvalues <- fread("macs.pvalues.csv")
macs.pvalues[, peakBases := end - start]
source("getMeanMat.R")
source("getSizeMat.R")

exp.stats.dt <- data.table(joint_peaks_mat.csv=Sys.glob("labels/H*_TDH_immune/samples/monocyte/McGill0001/joint_peaks_mat.csv"))[, {
  mean.mat <- readMeanMat(joint_peaks_mat.csv)
  mean.tall <- data.table(melt(mean.mat))
  setnames(mean.tall, c("peak", "pos", "mean"))
  i.stats <- mean.tall[, list(
    median=median(mean),
    q25=quantile(mean, 0.25),
    q75=quantile(mean, 0.75),
    peaks=.N
    ), by=list(pos)]
  sample.dir <- dirname(joint_peaks_mat.csv)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  data.table(experiment, i.stats)
}, by=list(joint_peaks_mat.csv)]

joint.peak.counts <- exp.stats.dt[pos==1]
setkey(joint.peak.counts, experiment)

## This is for creating macs2.broad_peaks.bedGraph files, for input to
## getMeanMat above.
macs.pvalues[sample.id=="McGill0001", {
  joint.peak.count <- joint.peak.counts[experiment, peaks]
  set.name <- paste0(experiment, "_TDH_immune")
  sample.dir <- file.path("labels", set.name, "samples/monocyte/McGill0001")
  peaks.bedGraph <- file.path(sample.dir, paste0(model, ".topN_peaks.bedGraph"))
  top.peaks <- .SD[order(-`-log10(pvalue)`)][1:joint.peak.count, list(
    chr, start, end, mean=NA)]
  cat(sprintf("writing top %d peaks to %s\n", nrow(top.peaks), peaks.bedGraph))
  fwrite(top.peaks,
         peaks.bedGraph,
         sep="\t",
         col.names=FALSE)
}, by=list(model, experiment)]

## Parallel version.
peaks.bedGraph.vec <- Sys.glob("labels/H*_TDH_immune/samples/monocyte/McGill0001/*_peaks.bedGraph")
library(parallel)
options(mc.cores=4)
exp.stats.dt.list <- mclapply(peaks.bedGraph.vec, function(peaks.bedGraph){
  print(peaks.bedGraph)
  mean.mat <- getMeanMat(peaks.bedGraph)
  mean.tall <- data.table(melt(mean.mat))
  setnames(mean.tall, c("peak", "pos", "mean"))
  i.stats <- mean.tall[, list(
    median=median(mean),
    q25=quantile(mean, 0.25),
    q75=quantile(mean, 0.75),
    count=.N
    ), by=list(pos)]
  sample.dir <- dirname(peaks.bedGraph)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  model <- sub("_.*", "", basename(peaks.bedGraph))
  data.table(experiment, model, i.stats)
})
exp.stats.dt <- do.call(rbind, exp.stats.dt.list)

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ model, scales="free")+
  geom_ribbon(aes(pos, ymin=q25, ymax=q75), alpha=0.5, data=exp.stats.dt)+
  geom_line(aes(pos, median), data=exp.stats.dt)+
  scale_x_continuous(
    "position relative to peak",
    breaks=c(5.5, 15.5),
    labels=c("peakStart", "peakEnd"))+
  scale_y_continuous("aligned read coverage (median and quartiles)")+
  geom_text(aes(
    5.5, -Inf, label=paste(count, "peaks")), data=exp.stats.dt[pos==5],
            hjust=0,
            vjust=-0.5)
print(gg)

png("figure-macs-pvalues-metapeak-topN.png", 1200, 500, res=100)
print(gg)
dev.off()

exp.stats.dt.list <- mclapply(peaks.bedGraph.vec, function(peaks.bedGraph){
  print(peaks.bedGraph)
  mean.mat <- getMeanMat(peaks.bedGraph)
  pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
  size.dt <- data.table(str_match_named(rownames(mean.mat), pattern, list(
    peakStart=as.integer,
    peakEnd=as.integer)))
  size.dt[, peakBases := peakEnd - peakStart]
  size.dt[, peakHeight := rowMeans(mean.mat[, -c(1:5, 16:20)]) ]
  sample.dir <- dirname(peaks.bedGraph)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  model <- sub("_.*", "", basename(peaks.bedGraph))
  data.table(experiment, model, size.dt)
})
exp.stats.dt <- do.call(rbind, exp.stats.dt.list)
mean.dt <- exp.stats.dt[, list(
  mean.peakBases=mean(peakBases),
  mean.peakHeight=mean(peakHeight),
  median.peakBases=median(peakBases),
  median.peakHeight=median(peakHeight)
  ), by=list(experiment, model)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ model)+
  geom_bin2d(aes(log10(peakBases), log10(peakHeight), fill=..density..), data=exp.stats.dt)+
  scale_fill_gradient(low="white", high="black")+
  ## geom_point(aes(
  ##   log10(mean.peakBases), log10(mean.peakHeight)
  ##   ), data=mean.dt,
  ##            color="red")+
  geom_point(aes(
    log10(median.peakBases), log10(median.peakHeight)
    ), data=mean.dt,
             shape=1,
             color="red")+
  geom_text(aes(
    4, 0,
    label=paste0(
      "medians:
height=", round(median.peakHeight, 1), " reads
width=", round(median.peakBases), " bases")),
            data=mean.dt, color="red")
print(gg)

png("figure-macs-pvalues-compare-size-height.png", 1200, 500, res=100)
print(gg)
dev.off()

## macs 
sample.peaks.dt <- macs.pvalues[experiment=="H3K36me3" & model=="macs.broad" & sample.id=="McGill0001"][order(-`-log10(pvalue)`)]
N <- 100
sample.peaks.dt[, peak.name := sprintf("%s:%d-%d", chr, start, end)]
topN.lik.dt <- sample.peaks.dt[1:N, list(`-log10(pvalue)`, peak.name)]
setnames(topN.lik.dt, c("lik", "peak.name"))
topN.peak.vec <- topN.lik.dt$peak.name
pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
topN.peak.mat <- str_match_named(topN.peak.vec, pattern)
topN.peak.dt <- data.table(topN.peak.mat, count=0)
peaks.bedGraph <- paste0(
  "labels/H3K36me3_TDH_immune/samples/monocyte/McGill0001/macs2.broad.top",
  N, "_peaks.bedGraph")
fwrite(topN.peak.dt, peaks.bedGraph, quote=FALSE, sep="\t", col.names=FALSE)
mean.mat <- getMeanMat(peaks.bedGraph)
d.mat <- dist(mean.mat)
h.fit <- hclust(d.mat, "average")
h.list <- dendro_data(h.fit)
rownames(h.list$labels) <- h.list$labels$label
mean.dt <- data.table(
  peak.name=rownames(mean.mat)[row(mean.mat)],
  position=as.integer(col(mean.mat)),
  mean=as.numeric(mean.mat))
mean.dt[, peak.x := h.list$labels[peak.name,]$x ]
topN.lik.dt[, peak.x := h.list$labels[peak.name,]$x ]
mean.dt[, pos.norm := position / max(position) * max(h.list$segments$y)]

ggplot()+
  geom_segment(aes(y, x, xend=yend, yend=xend), data=h.list$segments)+
  geom_tile(aes(-pos.norm, peak.x, fill=mean), data=mean.dt, color=NA)+
  scale_fill_gradient(low="white", high="blue")+
  geom_point(aes(0, peak.x, color=log10(lik)), data=topN.lik.dt)+
  scale_color_gradient(low="white", high="red")+
  xlab("position relative to peakStart/end | euclidean distance (average linkage)")+
  ylab("peak")

## macs top11 peaks absolute x axis.
sample.peaks.dt <- macs.pvalues[experiment=="H3K36me3" & model=="macs.broad" & sample.id=="McGill0001"][order(-`-log10(pvalue)`)]
N <- 100
sample.peaks.dt[, peak.name := sprintf("%s:%d-%d", chr, start, end)]
topN.lik.dt <- sample.peaks.dt[1:N, list(`-log10(pvalue)`, peak.name)]
setnames(topN.lik.dt, c("lik", "peak.name"))
topN.peak.vec <- topN.lik.dt$peak.name
pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
topN.peak.mat <- str_match_named(topN.peak.vec, pattern)
topN.peak.dt <- data.table(topN.peak.mat, count=0)
peaks.bedGraph <- paste0(
  "labels/H3K36me3_TDH_immune/samples/monocyte/McGill0001/macs2.broad.top",
  N, "_peaks.bedGraph")
fwrite(topN.peak.dt, peaks.bedGraph, quote=FALSE, sep="\t", col.names=FALSE)
mean.mat <- getSizeMat(peaks.bedGraph)
d.mat <- dist(mean.mat)
h.fit <- hclust(d.mat, "average")
h.list <- dendro_data(h.fit)
rownames(h.list$labels) <- h.list$labels$label
mean.dt <- data.table(
  peak.name=rownames(mean.mat)[row(mean.mat)],
  position=as.integer(col(mean.mat)),
  mean=as.numeric(mean.mat))
mean.dt[, peak.x := h.list$labels[peak.name,]$x ]
topN.lik.dt[, peak.x := h.list$labels[peak.name,]$x ]
mean.dt[, pos.norm := position / max(position) * max(h.list$segments$y)]

ggplot()+
  geom_segment(aes(y, x, xend=yend, yend=xend), data=h.list$segments)+
  geom_tile(aes(-pos.norm, peak.x, fill=mean), data=mean.dt, color=NA)+
  scale_fill_gradient(low="white", high="blue")+
  geom_point(aes(0, peak.x, color=log10(lik)), data=topN.lik.dt)+
  scale_color_gradient(low="white", high="red")+
  xlab("position relative to peakStart/end | euclidean distance (average linkage)")+
  ylab("peak")

## Clustering via first principal component, not as interesting as
## hierarchical.
pc.fit <- princomp(mean.mat)
pc.ord.vec <- structure(order(pc.fit$scores[,1]), names=rownames(mean.mat))
mean.dt[, peak.pc.ord := pc.ord.vec[peak.name] ]
ggplot()+
  geom_tile(aes(-pos.norm, peak.pc.ord, fill=mean), data=mean.dt, color=NA)+
  scale_fill_gradient(low="white", high="blue")+
  xlab("position relative to peakStart/end")+
  ylab("peak")

## Experiment: heatmap cluster peak mean matrix.
lik.dt <- fread("zcat labels/H3K36me3_TDH_immune/peaks_matrix_likelihood.tsv.gz")
lik.ord.dt <- lik.dt[order(-`monocyte/McGill0001`), list(`monocyte/McGill0001`, peak.name)]
N <- 100
topN.lik.dt <- lik.ord.dt[1:N]
setnames(topN.lik.dt, c("lik", "peak.name"))
setkey(topN.lik.dt, peak.name)
topN.peak.vec <- topN.lik.dt$peak.name
pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
topN.peak.mat <- str_match_named(topN.peak.vec, pattern)
topN.peak.dt <- data.table(topN.peak.mat, count=0)
peaks.bedGraph <- paste0(
  "labels/H3K36me3_TDH_immune/samples/monocyte/McGill0001/joint.top",
  N, "_peaks.bedGraph")
fwrite(topN.peak.dt, peaks.bedGraph, quote=FALSE, sep="\t", col.names=FALSE)
mean.mat <- getMeanMat(peaks.bedGraph)
d.mat <- dist(mean.mat)
h.fit <- hclust(d.mat, "average")
h.list <- dendro_data(h.fit)
rownames(h.list$labels) <- h.list$labels$label
mean.dt <- data.table(
  peak.name=rownames(mean.mat)[row(mean.mat)],
  position=as.integer(col(mean.mat)),
  mean=as.numeric(mean.mat))
mean.dt[, peak.x := h.list$labels[peak.name,]$x ]
topN.lik.dt[, peak.x := h.list$labels[peak.name,]$x ]
max.dist <- max(h.list$segments$y)
mean.dt[, pos.norm := position / max(position) * max.dist]
match.df <- str_match_named(mean.dt$peak.name, pattern, list(
  peakStart=as.integer,
  peakEnd=as.integer))
mean.dt[, peakBases := with(match.df, peakEnd-peakStart)]

ggplot()+
  geom_text(aes(
    -max.dist, peak.x,
    label=paste0(as.integer(peakBases/1e3), "kb  ")),
            data=mean.dt[position==1], hjust=1)+
  geom_segment(aes(y, x, xend=yend, yend=xend), data=h.list$segments)+
  geom_tile(aes(-pos.norm, peak.x, fill=mean), data=mean.dt, color=NA)+
  scale_fill_gradient(low="white", high="blue")+
  geom_point(aes(0, peak.x, color=log10(lik)), data=topN.lik.dt)+
  scale_color_gradient(low="white", high="red")+
  xlab("position relative to peakStart/end | euclidean distance (average linkage)")+
  ylab("peak")

## combined plot
sample.peaks.dt <- macs.pvalues[experiment=="H3K36me3" & model=="macs.broad" & sample.id=="McGill0001"][order(-`-log10(pvalue)`)]
N <- 100
sample.peaks.dt[, peak.name := sprintf("%s:%d-%d", chr, start, end)]
lik.dt <- fread("zcat labels/H3K36me3_TDH_immune/peaks_matrix_likelihood.tsv.gz")
lik.ord.dt <- lik.dt[order(-`monocyte/McGill0001`), list(`monocyte/McGill0001`, peak.name)]
topN.list <- list(
  joint=lik.ord.dt[1:N],
  macs2.broad=sample.peaks.dt[1:N, list(`-log10(pvalue)`, peak.name)])

heatmap.dt.list <- list()
size.dt.list <- list()
for(model.name in names(topN.list)){
  topN.lik.dt <- topN.list[[model.name]]
  setnames(topN.lik.dt, c("lik", "peak.name"))
  topN.peak.vec <- topN.lik.dt$peak.name
  pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
  topN.peak.mat <- str_match_named(topN.peak.vec, pattern, list(
    peakStart=as.integer,
    peakEnd=as.integer))
  topN.peak.dt <- data.table(topN.peak.mat, count=0)
  topN.peak.dt[, peakBases := peakEnd - peakStart]
  peaks.bedGraph <- paste0(
    "labels/H3K36me3_TDH_immune/samples/monocyte/McGill0001/",
    model.name, ".top",
    N, "_peaks.bedGraph")
  fwrite(topN.peak.dt, peaks.bedGraph, quote=FALSE, sep="\t", col.names=FALSE)
  mean.mat <- getMeanMat(peaks.bedGraph)
  d.mat <- dist(mean.mat)
  h.fit <- hclust(d.mat, "average")
  h.list <- dendro_data(h.fit)
  rownames(h.list$labels) <- h.list$labels$label
  mean.dt <- data.table(
    peak.name=rownames(mean.mat)[row(mean.mat)],
    position=as.integer(col(mean.mat)),
    mean=as.numeric(mean.mat))
  mean.dt[, peak.x := h.list$labels[peak.name,]$x ]
  topN.peak.dt[, peak.x := h.list$labels[topN.lik.dt$peak.name,]$x ]
  heatmap.dt.list[[model.name]] <- data.table(model.name, mean.dt)
  size.dt.list[[model.name]] <- data.table(model.name, topN.peak.dt)
}
heatmap.dt <- do.call(rbind, heatmap.dt.list)
size.dt <- do.call(rbind, size.dt.list)

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.name)+
  geom_tile(aes(position-21, peak.x, fill=mean), data=heatmap.dt, color=NA)+
  geom_segment(aes(0, peak.x, xend=peakBases/1e4, yend=peak.x), data=size.dt)+
  scale_fill_gradient(low="white", high="blue")+
  scale_x_continuous(
    "Left: position relative to peakStart/end, Right: peak size",
    breaks=c(-15.5, -5.5, 10, 20),
    labels=c("peakStart", "peakEnd", "100kb", "200kb"))+
  ylab("peak")
png("labels/H3K36me3_TDH_immune/figure-macs-pvalues-top100-heatmap-compare.png", 1000, 500, res=100)
print(gg)
dev.off()

## combined plot on abs scale
sample.peaks.dt <- macs.pvalues[experiment=="H3K36me3" & model=="macs.broad" & sample.id=="McGill0001"][order(-`-log10(pvalue)`)]
N <- 100
sample.peaks.dt[, peak.name := sprintf("%s:%d-%d", chr, start, end)]
lik.dt <- fread("zcat labels/H3K36me3_TDH_immune/peaks_matrix_likelihood.tsv.gz")
lik.ord.dt <- lik.dt[order(-`monocyte/McGill0001`), list(`monocyte/McGill0001`, peak.name)]
dfilter.dt <- fread("sed 's/ \\+/\t/g' labels/H3K36me3_TDH_immune/samples/monocyte/McGill0001/DFilter_peaks.tsv")[order(-maxScore)]
dfilter.dt[, peak.name := sprintf("%s:%d-%d", chromosome, `peak-start`, `peak-end`)]
topN.list <- list(
  DFilter=dfilter.dt[1:N, list(maxScore, peak.name)],
  joint=lik.ord.dt[1:N],
  macs2.broad=sample.peaks.dt[1:N, list(`-log10(pvalue)`, peak.name)])

heatmap.dt.list <- list()
size.dt.list <- list()
bases.per.bin <- 20000
for(model.name in names(topN.list)){
  topN.lik.dt <- topN.list[[model.name]]
  setnames(topN.lik.dt, c("lik", "peak.name"))
  topN.peak.vec <- topN.lik.dt$peak.name
  pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
  topN.peak.mat <- str_match_named(topN.peak.vec, pattern, list(
    peakStart=as.integer,
    peakEnd=as.integer))
  topN.peak.dt <- data.table(topN.peak.mat, count=0)
  peaks.bedGraph <- paste0(
    "labels/H3K36me3_TDH_immune/samples/monocyte/McGill0001/",
    model.name, ".top",
    N, "_peaks.bedGraph")
  fwrite(topN.peak.dt, peaks.bedGraph, quote=FALSE, sep="\t", col.names=FALSE)
  topN.peak.dt[, peakBases := peakEnd - peakStart]
  ord.peak.dt <- topN.peak.dt[order(peakBases)]
  mean.mat <- getSizeMat(peaks.bedGraph, bases.per.bin)
  d.mat <- dist(mean.mat)
  h.fit <- hclust(d.mat, "average")
  h.list <- dendro_data(h.fit)
  rownames(h.list$labels) <- h.list$labels$label
  mean.dt <- data.table(
    peak.name=rownames(mean.mat)[row(mean.mat)],
    position=as.integer(col(mean.mat)),
    mean=as.numeric(mean.mat))
  mean.dt[, peak.x := h.list$labels[peak.name,]$x ]
  ord.peak.dt[, peak.x := h.list$labels[topN.lik.dt$peak.name,]$x ]
  ord.peak.dt[, size.x := 1:.N]
  ord.peak.dt[, peak.name := sprintf("%s:%d-%d", chrom, peakStart, peakEnd)]
  mean.dt[, size.x := ord.peak.dt[peak.name, on=list(peak.name)]$size.x]
  heatmap.dt.list[[model.name]] <- data.table(model.name, mean.dt)
  size.dt.list[[model.name]] <- data.table(model.name, ord.peak.dt)
}
heatmap.dt <- do.call(rbind, heatmap.dt.list)
size.dt <- do.call(rbind, size.dt.list)

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.name)+
  geom_tile(aes(position-21, peak.x, fill=mean), data=heatmap.dt, color=NA)+
  geom_segment(aes(0, peak.x, xend=peakBases/1e4, yend=peak.x), data=size.dt)+
  scale_fill_gradient(low="white", high="blue")+
  scale_x_continuous(
    "Left: position relative to peakStart/end, Right: peak size",
    breaks=c(-15.5, -5.5, 10, 20),
    labels=c("peakStart", "peakEnd", "100kb", "200kb"))+
  ylab("peak")

size.dt[, halfNorm := peakBases/bases.per.bin/2]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.name)+
  geom_tile(aes(position-10.5, size.x, fill=mean), data=heatmap.dt, color=NA)+
  ## geom_segment(aes(
  ##   -halfNorm, size.x, xend=halfNorm, yend=size.x), alpha=0.5, data=size.dt)+
  geom_point(aes(-halfNorm, size.x), data=size.dt, alpha=0.5, shape=1)+
  geom_point(aes(halfNorm, size.x), data=size.dt, alpha=0.5, shape=1)+
  scale_fill_gradient("mean\ncoverage", low="white", high="blue")+
  scale_x_continuous(
    "position relative to peak center (kb = kilo bases)",
    breaks=c(0, -5.5, 5.5),
    labels=c("0", "-100", "100"))+
  ylab("top 100 most likely/significant peaks ordered by size")+
  coord_equal()
print(gg)

png("labels/H3K36me3_TDH_immune/figure-macs-pvalues-top100-heatmap-compare-absolute.png", 500, 700, res=100)
print(gg)
dev.off()

## combined plot on abs scale, for H3K4me3
sample.peaks.dt <- macs.pvalues[experiment=="H3K4me3" & model=="macs.default" & sample.id=="McGill0001"][order(-`-log10(pvalue)`)]
N <- 100
sample.peaks.dt[, peak.name := sprintf("%s:%d-%d", chr, start, end)]
lik.dt <- fread("zcat labels/H3K4me3_TDH_immune/peaks_matrix_likelihood.tsv.gz")
lik.ord.dt <- lik.dt[order(-`monocyte/McGill0001`), list(`monocyte/McGill0001`, peak.name)]
dfilter.dt <- fread("sed 's/ \\+/\t/g' labels/H3K4me3_TDH_immune/samples/monocyte/McGill0001/DFilter_peaks.tsv")[order(-maxScore)]
dfilter.dt[, peak.name := sprintf("%s:%d-%d", chromosome, `peak-start`, `peak-end`)]
topN.list <- list(
  DFilter=dfilter.dt[1:N, list(maxScore, peak.name)],
  joint=lik.ord.dt[1:N],
  macs2.default=sample.peaks.dt[1:N, list(`-log10(pvalue)`, peak.name)])

heatmap.dt.list <- list()
size.dt.list <- list()
bases.per.bin <- 1000
for(model.name in names(topN.list)){
  topN.lik.dt <- topN.list[[model.name]]
  setnames(topN.lik.dt, c("lik", "peak.name"))
  topN.peak.vec <- topN.lik.dt$peak.name
  pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
  topN.peak.mat <- str_match_named(topN.peak.vec, pattern, list(
    peakStart=as.integer,
    peakEnd=as.integer))
  topN.peak.dt <- data.table(topN.peak.mat, count=0)
  peaks.bedGraph <- paste0(
    "labels/H3K4me3_TDH_immune/samples/monocyte/McGill0001/",
    model.name, ".top",
    N, "_peaks.bedGraph")
  fwrite(topN.peak.dt, peaks.bedGraph, quote=FALSE, sep="\t", col.names=FALSE)
  topN.peak.dt[, peakBases := peakEnd - peakStart]
  ord.peak.dt <- topN.peak.dt[order(peakBases)]
  mean.mat <- getSizeMat(peaks.bedGraph, bases.per.bin)
  d.mat <- dist(mean.mat)
  h.fit <- hclust(d.mat, "average")
  h.list <- dendro_data(h.fit)
  rownames(h.list$labels) <- h.list$labels$label
  mean.dt <- data.table(
    peak.name=rownames(mean.mat)[row(mean.mat)],
    position=as.integer(col(mean.mat)),
    mean=as.numeric(mean.mat))
  mean.dt[, peak.x := h.list$labels[peak.name,]$x ]
  ord.peak.dt[, peak.x := h.list$labels[topN.lik.dt$peak.name,]$x ]
  ord.peak.dt[, size.x := 1:.N]
  ord.peak.dt[, peak.name := sprintf("%s:%d-%d", chrom, peakStart, peakEnd)]
  mean.dt[, size.x := ord.peak.dt[peak.name, on=list(peak.name)]$size.x]
  heatmap.dt.list[[model.name]] <- data.table(model.name, mean.dt)
  size.dt.list[[model.name]] <- data.table(model.name, ord.peak.dt)
}
heatmap.dt <- do.call(rbind, heatmap.dt.list)
size.dt <- do.call(rbind, size.dt.list)

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.name)+
  geom_tile(aes(position-21, peak.x, fill=mean), data=heatmap.dt, color=NA)+
  geom_segment(aes(0, peak.x, xend=peakBases/1e4, yend=peak.x), data=size.dt)+
  scale_fill_gradient(low="white", high="blue")+
  scale_x_continuous(
    "Left: position relative to peakStart/end, Right: peak size",
    breaks=c(-15.5, -5.5, 10, 20),
    labels=c("peakStart", "peakEnd", "100kb", "200kb"))+
  ylab("peak")

size.dt[, halfNorm := peakBases/bases.per.bin/2]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.name)+
  geom_tile(aes(position-10.5, size.x, fill=mean), data=heatmap.dt, color=NA)+
  geom_point(aes(-halfNorm, size.x), data=size.dt, alpha=0.5, shape=1)+
  geom_point(aes(halfNorm, size.x), data=size.dt, alpha=0.5, shape=1)+
  scale_fill_gradient("mean\ncoverage", low="white", high="blue")+
  scale_x_continuous(
    "position relative to peak center (kb = kilo bases)",
    breaks=c(0, -5.5, 5.5),
    labels=c("0", c(-1,1)*bases.per.bin*5/1e3))+
  ylab("top 100 most likely/significant peaks ordered by size")+
  coord_equal()
print(gg)

png("labels/H3K4me3_TDH_immune/figure-macs-pvalues-top100-heatmap-compare-absolute.png", 500, 700, res=100)
print(gg)
dev.off()

## combined plot on abs scale, for H3K4me3, with JAMM
sample.dir <- "labels/H3K4me3_TDH_immune/samples/tcell/McGill0107"
sample.dir <- "labels/H3K4me3_TDH_immune/samples/monocyte/McGill0001"
sample.id <- basename(sample.dir)
group.dir <- dirname(sample.dir)
samples.dir <- dirname(group.dir)
set.dir <- dirname(samples.dir)
set.name <- basename(set.dir)
experiment <- sub("_.*", "", set.name)
select.dt <- data.table(
  experiment,
  model=ifelse(experiment=="H3K36me3", "macs.broad", "macs.default"),
  sample.id)
sample.peaks.dt <- macs.pvalues[select.dt, on=list(experiment, model, sample.id)][order(-`-log10(pvalue)`)]
N <- 100
sample.peaks.dt[, peak.name := sprintf("%s:%d-%d", chr, start, end)]
lik.gz <- file.path("labels", set.name, "peaks_matrix_likelihood.tsv.gz")
lik.dt <- fread(paste("zcat", lik.gz))
cell.type <- basename(group.dir)
sample.path <- paste0(cell.type, "/", sample.id)
lik.vec <- lik.dt[[sample.path]]
lik.ord.dt <- data.table(lik=lik.vec, peak=lik.dt$peak.name)[order(-lik)]
dfilter.dt <- fread(paste0("sed 's/ \\+/\t/g' ", sample.dir, "/DFilter_peaks.tsv"))[order(-maxScore)]
dfilter.dt[, peak.name := sprintf("%s:%d-%d", chromosome, `peak-start`, `peak-end`)]
jamm.narrowPeak <- file.path(
  "~",
  "JAMM",
  experiment,
  cell.type,
  "peaks",
  "filtered.peaks.narrowPeak")
jamm <- fread(jamm.narrowPeak, select=c(1:3, 7))
setnames(jamm, c("chrom", "chromStart", "chromEnd", "score"))
jamm[, peak.name := sprintf("%s:%d-%d", chrom, chromStart, chromEnd)]
topN.list <- list(
  PeakSegPipeline=lik.ord.dt[1:N]
  ,JAMM=jamm[1:N, list(score, peak.name)]
  ,DFilter=dfilter.dt[1:N, list(maxScore, peak.name)]
  ,MACS2=sample.peaks.dt[1:N, list(`-log10(pvalue)`, peak.name)]
  )

heatmap.dt.list <- list()
size.dt.list <- list()
bases.per.bin <- 1000
for(model.name in names(topN.list)){
  topN.lik.dt <- topN.list[[model.name]]
  setnames(topN.lik.dt, c("lik", "peak.name"))
  topN.peak.vec <- topN.lik.dt$peak.name
  pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
  topN.peak.mat <- str_match_named(topN.peak.vec, pattern, list(
    peakStart=as.integer,
    peakEnd=as.integer))
  topN.peak.dt <- data.table(topN.peak.mat, count=0)
  peaks.bedGraph <- paste0(
    sample.dir,
    "/",
    model.name, ".top",
    N, "_peaks.bedGraph")
  fwrite(topN.peak.dt, peaks.bedGraph, quote=FALSE, sep="\t", col.names=FALSE)
  topN.peak.dt[, peakBases := peakEnd - peakStart]
  ord.peak.dt <- topN.peak.dt[order(peakBases)]
  mean.mat <- getSizeMat(peaks.bedGraph, bases.per.bin)
  d.mat <- dist(mean.mat)
  mean.dt <- data.table(
    peak.name=rownames(mean.mat)[row(mean.mat)],
    position=as.integer(col(mean.mat)),
    mean=as.numeric(mean.mat))
  ord.peak.dt[, size.x := 1:.N]
  ord.peak.dt[, peak.name := sprintf("%s:%d-%d", chrom, peakStart, peakEnd)]
  mean.dt[, size.x := ord.peak.dt[peak.name, on=list(peak.name)]$size.x]
  model.fac <- factor(model.name, names(topN.list))
  heatmap.dt.list[[model.name]] <- data.table(model.name, model.fac, mean.dt)
  size.dt.list[[model.name]] <- data.table(model.name, model.fac, ord.peak.dt)
}
heatmap.dt <- do.call(rbind, heatmap.dt.list)
size.dt <- do.call(rbind, size.dt.list)

size.dt[, halfNorm := peakBases/bases.per.bin/2]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.fac)+
  geom_tile(aes(position-10.5, size.x, fill=mean), data=heatmap.dt, color=NA)+
  geom_point(aes(-halfNorm, size.x), data=size.dt, alpha=0.5, shape=1)+
  geom_point(aes(halfNorm, size.x), data=size.dt, alpha=0.5, shape=1)+
  scale_fill_gradient("mean\ncoverage", low="white", high="blue")+
  scale_x_continuous(
    "position relative to peak center (kb = kilo bases)",
    breaks=c(0, -5.5, 5.5),
    labels=c("0", c(-1,1)*bases.per.bin*5/1e3))+
  ylab("top 100 most likely/significant peaks ordered by size")+
  coord_equal()
print(gg)

png.name <- file.path(
  "labels",
  set.name,
  "figure-macs-pvalues-top100-heatmap-compare-absolute.png")
print(png.name)
png(png.name, 500, 700, res=100)
print(gg)
dev.off()

## combined plot on abs scale, for H3K36me3, with JAMM
sample.dir <- "labels/H3K36me3_TDH_immune/samples/monocyte/McGill0001"
sample.id <- basename(sample.dir)
group.dir <- dirname(sample.dir)
samples.dir <- dirname(group.dir)
set.dir <- dirname(samples.dir)
set.name <- basename(set.dir)
experiment <- sub("_.*", "", set.name)
select.dt <- data.table(
  experiment,
  model=ifelse(experiment=="H3K36me3", "macs.broad", "macs.default"),
  sample.id)
sample.peaks.dt <- macs.pvalues[select.dt, on=list(experiment, model, sample.id)][order(-`-log10(pvalue)`)]
N <- 100
sample.peaks.dt[, peak.name := sprintf("%s:%d-%d", chr, start, end)]
lik.gz <- file.path("labels", set.name, "peaks_matrix_likelihood.tsv.gz")
lik.dt <- fread(paste("zcat", lik.gz))
cell.type <- basename(group.dir)
sample.path <- paste0(cell.type, "/", sample.id)
lik.vec <- lik.dt[[sample.path]]
lik.ord.dt <- data.table(lik=lik.vec, peak=lik.dt$peak.name)[order(-lik)]
dfilter.dt <- fread(paste0("sed 's/ \\+/\t/g' ", sample.dir, "/DFilter_peaks.tsv"))[order(-maxScore)]
dfilter.dt[, peak.name := sprintf("%s:%d-%d", chromosome, `peak-start`, `peak-end`)]
jamm.narrowPeak <- file.path(
  "~",
  "JAMM",
  experiment,
  cell.type,
  "peaks",
  "filtered.peaks.narrowPeak")
jamm <- fread(jamm.narrowPeak, select=c(1:3, 7))
setnames(jamm, c("chrom", "chromStart", "chromEnd", "score"))
jamm[, peak.name := sprintf("%s:%d-%d", chrom, chromStart, chromEnd)]
topN.list <- list(
  PeakSegPipeline=lik.ord.dt[1:N]
  ,JAMM=jamm[1:N, list(score, peak.name)]
  ,DFilter=dfilter.dt[1:N, list(maxScore, peak.name)]
  ,MACS2=sample.peaks.dt[1:N, list(`-log10(pvalue)`, peak.name)]
  )

heatmap.dt.list <- list()
size.dt.list <- list()
bases.per.bin <- 20000
for(model.name in names(topN.list)){
  topN.lik.dt <- topN.list[[model.name]]
  setnames(topN.lik.dt, c("lik", "peak.name"))
  topN.peak.vec <- topN.lik.dt$peak.name
  pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
  topN.peak.mat <- str_match_named(topN.peak.vec, pattern, list(
    peakStart=as.integer,
    peakEnd=as.integer))
  topN.peak.dt <- data.table(topN.peak.mat, count=0)
  peaks.bedGraph <- paste0(
    sample.dir,
    "/",
    model.name, ".top",
    N, "_peaks.bedGraph")
  fwrite(topN.peak.dt, peaks.bedGraph, quote=FALSE, sep="\t", col.names=FALSE)
  topN.peak.dt[, peakBases := peakEnd - peakStart]
  ord.peak.dt <- topN.peak.dt[order(peakBases)]
  mean.mat <- getSizeMat(peaks.bedGraph, bases.per.bin)
  d.mat <- dist(mean.mat)
  mean.dt <- data.table(
    peak.name=rownames(mean.mat)[row(mean.mat)],
    position=as.integer(col(mean.mat)),
    mean=as.numeric(mean.mat))
  ord.peak.dt[, size.x := 1:.N]
  ord.peak.dt[, peak.name := sprintf("%s:%d-%d", chrom, peakStart, peakEnd)]
  mean.dt[, size.x := ord.peak.dt[peak.name, on=list(peak.name)]$size.x]
  model.fac <- factor(model.name, names(topN.list))
  heatmap.dt.list[[model.name]] <- data.table(model.name, model.fac, mean.dt)
  size.dt.list[[model.name]] <- data.table(model.name, model.fac, ord.peak.dt)
}
heatmap.dt <- do.call(rbind, heatmap.dt.list)
size.dt <- do.call(rbind, size.dt.list)

size.dt[, halfNorm := peakBases/bases.per.bin/2]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.fac)+
  geom_tile(aes(position-10.5, size.x, fill=mean), data=heatmap.dt, color=NA)+
  geom_point(aes(-halfNorm, size.x), data=size.dt, alpha=0.5, shape=1)+
  geom_point(aes(halfNorm, size.x), data=size.dt, alpha=0.5, shape=1)+
  scale_fill_gradient("mean\ncoverage", low="white", high="blue")+
  scale_x_continuous(
    "position relative to peak center (kb = kilo bases)",
    breaks=c(0, -5.5, 5.5),
    labels=c("0", c(-1,1)*bases.per.bin*5/1e3))+
  ylab("top 100 most likely/significant peaks ordered by size")+
  coord_equal()
print(gg)

png.name <- file.path(
  "labels",
  set.name,
  "figure-macs-pvalues-top100-heatmap-compare-absolute.png")
print(png.name)
png(png.name, 500, 700, res=100)
print(gg)
dev.off()

## combined plot on abs scale, for both experiments, with JAMM
sample.dir.vec <- c(
  "labels/H3K4me3_TDH_immune/samples/monocyte/McGill0001",
  "labels/H3K36me3_TDH_immune/samples/monocyte/McGill0001")
heatmap.dt.list <- list()
size.dt.list <- list()
for(sample.dir in sample.dir.vec){
  sample.id <- basename(sample.dir)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  select.dt <- data.table(
    experiment,
    model=ifelse(experiment=="H3K36me3", "macs.broad", "macs.default"),
    sample.id)
  sample.peaks.dt <- macs.pvalues[select.dt, on=list(experiment, model, sample.id)][order(-`-log10(pvalue)`)]
  N <- 100
  sample.peaks.dt[, peak.name := sprintf("%s:%d-%d", chr, start, end)]
  lik.gz <- file.path("labels", set.name, "peaks_matrix_likelihood.tsv.gz")
  lik.dt <- fread(paste("zcat", lik.gz))
  cell.type <- basename(group.dir)
  sample.path <- paste0(cell.type, "/", sample.id)
  lik.vec <- lik.dt[[sample.path]]
  lik.ord.dt <- data.table(lik=lik.vec, peak=lik.dt$peak.name)[order(-lik)]
  dfilter.dt <- fread(paste0("sed 's/ \\+/\t/g' ", sample.dir, "/DFilter_peaks.tsv"))[order(-maxScore)]
  dfilter.dt[, peak.name := sprintf("%s:%d-%d", chromosome, `peak-start`, `peak-end`)]
  jamm.narrowPeak <- file.path(
    "~",
    "JAMM",
    experiment,
    cell.type,
    "peaks",
    "filtered.peaks.narrowPeak")
  jamm <- fread(jamm.narrowPeak, select=c(1:3, 7))
  setnames(jamm, c("chrom", "chromStart", "chromEnd", "score"))
  jamm[, peak.name := sprintf("%s:%d-%d", chrom, chromStart, chromEnd)]
  topN.list <- list(
    PeakSegPipeline=lik.ord.dt[1:N]
    ,JAMM=jamm[1:N, list(score, peak.name)]
    ,DFilter=dfilter.dt[1:N, list(maxScore, peak.name)]
    ,MACS2=sample.peaks.dt[1:N, list(`-log10(pvalue)`, peak.name)]
    )
  bases.per.bin <- ifelse(experiment=="H3K36me3", 20000, 1000)
  for(model.name in names(topN.list)){
    topN.lik.dt <- topN.list[[model.name]]
    setnames(topN.lik.dt, c("lik", "peak.name"))
    topN.peak.vec <- topN.lik.dt$peak.name
    pattern <- paste0(
      "(?<chrom>chr[^:]+)",
      ":",
      "(?<peakStart>[^-]+)",
      "-",
      "(?<peakEnd>.*)")
    topN.peak.mat <- str_match_named(topN.peak.vec, pattern, list(
      peakStart=as.integer,
      peakEnd=as.integer))
    topN.peak.dt <- data.table(topN.peak.mat, count=0)
    peaks.bedGraph <- paste0(
      sample.dir,
      "/",
      model.name, ".top",
      N, "_peaks.bedGraph")
    fwrite(topN.peak.dt, peaks.bedGraph, quote=FALSE, sep="\t", col.names=FALSE)
    topN.peak.dt[, peakBases := peakEnd - peakStart]
    ord.peak.dt <- topN.peak.dt[order(peakBases)]
    mean.mat <- getSizeMat(peaks.bedGraph, bases.per.bin)
    d.mat <- dist(mean.mat)
    mean.dt <- data.table(
      peak.name=rownames(mean.mat)[row(mean.mat)],
      position=as.integer(col(mean.mat)),
      mean=as.numeric(mean.mat))
    mean.dt[, pos.bases := (position-10.5) * bases.per.bin]
    mean.dt[, min.bases := pos.bases -bases.per.bin/2]
    mean.dt[, max.bases := pos.bases +bases.per.bin/2]
    ord.peak.dt[, size.x := 1:.N]
    ord.peak.dt[, peak.name := sprintf("%s:%d-%d", chrom, peakStart, peakEnd)]
    mean.dt[, size.x := ord.peak.dt[peak.name, on=list(peak.name)]$size.x]
    model.fac <- factor(model.name, names(topN.list))
    heatmap.dt.list[[paste(experiment, model.name)]] <- data.table(
      experiment, model.name, model.fac, mean.dt)
    size.dt.list[[paste(experiment, model.name)]] <- data.table(
      experiment, model.name, model.fac, ord.peak.dt)
  }
}
heatmap.dt <- do.call(rbind, heatmap.dt.list)
size.dt <- do.call(rbind, size.dt.list)

psize <- 0.5
size.dt[, halfNorm := peakBases/bases.per.bin/2]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ model.fac, scales="free")+
  geom_rect(aes(
    xmin=size.x-0.5, xmax=size.x+0.5,
    ymin=min.bases/1e3, ymax=max.bases/1e3, fill=mean),
            data=heatmap.dt, color=NA)+
  geom_point(aes(size.x, -peakBases/2e3), data=size.dt, alpha=0.5, shape=1, size=psize)+
  geom_point(aes(size.x, peakBases/2e3), data=size.dt, alpha=0.5, shape=1, size=psize)+
  scale_fill_gradient("mean\ncoverage", low="white", high="blue")+
  scale_y_continuous(
    "position relative to peak center (kb)")+
  xlab("top 100 most likely/significant peaks ordered by size")
##print(gg)
png.name <- file.path(
  "labels",
  set.name,
  "figure-top100-absolute.png")
print(png.name)
png(png.name, 1200, 300, res=100)
print(gg)
dev.off()
pdf(sub("png$", "pdf", png.name), 12, 3)
print(gg)
dev.off()
