library(data.table)
library(flexclust)
library(ggplot2)
set.vec <- c(
  "H3K4me3_TDH_immune",
  "H3K36me3_TDH_immune",
  "H3K4me1_TDH_BP",
  "H3K9me3_TDH_BP",
  "H3K27ac_TDH_some",
  "H3K27me3_TDH_some")

## experiment 1: plot distribution of peak sizes in two samples per data set.
bg.dt <- data.table(set.name=set.vec)[, {
  data.table(joint_peaks.bedGraph=Sys.glob(file.path("labels", set.name, "samples/*/*/joint_peaks.bedGraph")))
}, by=set.name]
some.bg <- bg.dt[, {
  data.table(experiment=sub("_.*", "", set.name), .SD[1:2])
}, by=set.name]
stats.dt <- some.bg[, {
  peaks <- fread(joint_peaks.bedGraph)
  setnames(peaks, c("chrom", "peakStart", "peakEnd", "mean"))
  peaks[, peakBases := peakEnd-peakStart]
  peaks[, list(
    q25=quantile(peakBases, 0.25),
    median=as.numeric(median(peakBases)),
    q75=quantile(peakBases, 0.75)
  )]
}, by=list(experiment, joint_peaks.bedGraph)]
stats.dt[, id := sub(".*samples/", "", sub("/joint_peaks.bedGraph", "", joint_peaks.bedGraph))]
exp.sizes <- stats.dt[, list(mean=mean(median)), by=experiment][order(mean)]
stats.dt[, exp.fac := factor(experiment, exp.sizes$experiment)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(exp.fac ~ ., scales="free", space="free")+
  geom_point(aes(median, id), data=stats.dt)+
  geom_segment(aes(
    q25, id,
    xend=q75, yend=id), data=stats.dt)+
  scale_x_log10("peak size (bases)")+
  ylab("sample ID")

## experiment 2: plot distribution of peak sizes in all joint peak calls.
summary.tsv.vec <- Sys.glob("labels/*/peaks_summary.tsv")
some.vec <- grep("-", summary.tsv.vec, invert=TRUE, value=TRUE)
stats.dt <- data.table(summary.tsv=some.vec)[, {
  summary.dt <- fread(summary.tsv)
  summary.dt[, list(
    set.name=basename(dirname(summary.tsv)),
    q25=quantile(peakBases, 0.25),
    median=as.numeric(median(peakBases)),
    q75=quantile(peakBases, 0.75)
  )]
}, by=summary.tsv]
set.sizes <- stats.dt[, list(mean=mean(median)), by=set.name][order(mean)]
stats.dt[, set.fac := factor(set.name, set.sizes$set.name)]
stats.dt[, samples := sub(".*_", "", set.name)]
stats.dt[, center := ifelse(
  samples=="BP", "BluePrint", ifelse(
  samples=="adipose", "unknown", ifelse(
  samples=="ENCODE", "ENCODE", "CEEHRC")))]

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_point(aes(median, set.fac, color=center), data=stats.dt)+
  geom_segment(aes(
    q25, set.fac,
    color=center,
    xend=q75, yend=set.fac), data=stats.dt)+
  scale_x_log10("peak size (bases)", limits=c(1e2, NA), breaks=10^seq(2, 5, by=1))+
  ylab("data set")
print(gg)

png("figure-compare-size-centers.png", 600, 200, res=100)
print(gg)
dev.off()
system("cp figure-compare-size-centers.png labels/H3K36me3_TDH_ENCODE")

stats.dt[, experiment := sub("_.*", "", set.name)]
stats.dt[, kb := median/1e3]
stats.dt[, show.bases := ifelse(kb < 1, paste(round(median), "b"), sprintf("%.1f kb", kb))]
stats.dt[, samples := sub(".*_", "", set.name)]
stats.dt[, show.samples := ifelse(experiment %in% c("H3K36me3", "H3K4me3"), samples, "")]
stats.dt[, exp.fac := factor(experiment, c("H3K36me3", "H3K9me3", "H3K27me3", "H3K4me1", "H3K27ac", "H3K4me3", "ATAC", "CTCF"))]
stats.dt[, side := ! kb < 1]
gg <- ggplot()+
  theme_bw()+
  theme(
    panel.margin=grid::unit(0, "lines"),
    axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  geom_text(aes(
    show.samples, ifelse(
      side, 0, Inf),
    hjust=ifelse(side, 0, 1),
    label=ifelse(side, paste0(" ", show.bases), paste0(show.bases, " "))),
            angle=90,
            size=3,
            data=stats.dt)+
  facet_grid(. ~ exp.fac, scales="free")+
  geom_point(aes(show.samples, median), data=stats.dt, shape=1)+
  geom_segment(aes(
    show.samples, q25, 
    yend=q75, xend=show.samples), data=stats.dt)+
  scale_y_log10("peak size (bases)", limits=c(1e2, NA), breaks=10^seq(2, 5, by=1))+
  xlab("data set")
print(gg)

png("labels/H3K36me3_TDH_ENCODE/figure-compare-size-panels.png", 600, 300, res=100)
print(gg)
dev.off()

pdf("labels/H3K36me3_TDH_ENCODE/figure-compare-size-panels.pdf", 6, 3)
print(gg)
dev.off()


## experiment 3: plot distribution of peak sizes of peaks present in
## all samples.
summary.tsv.vec <- Sys.glob("labels/*/peaks_summary.tsv")
some.vec <- grep("-", summary.tsv.vec, invert=TRUE, value=TRUE)
stats.dt <- data.table(summary.tsv=some.vec)[, {
  summary.dt <- fread(summary.tsv)
  no.input <- if("n.Input" %in% names(summary.dt)){
    summary.dt[n.Input==0]
  }else{
    summary.dt[, n.samples := n.samples.up]
    summary.dt[Input.up==FALSE]
  }
  browser(expr=grepl("H3K4me3_TDH_ENC", summary.tsv))
  max.samples <- max(no.input$n.samples)
  no.input[n.samples==max.samples, list(
    set.name=basename(dirname(summary.tsv)),
             max.samples,
    q25=quantile(peakBases, 0.25),
    median=as.numeric(median(peakBases)),
    q75=quantile(peakBases, 0.75),
             min=min(peakBases),
             max=max(peakBases),
               peaks=.N
  )]
}, by=summary.tsv]
set.sizes <- stats.dt[, list(mean=mean(median)), by=set.name][order(mean)]
stats.dt[, set.fac := factor(set.name, set.sizes$set.name)]
stats.dt[, samples := sub(".*_", "", set.name)]
stats.dt[, center := ifelse(samples=="BP", "BluePrint", ifelse(samples=="adipose", "unknown", "CEEHRC"))]
text.dt <- melt(stats.dt, measure.vars=c("q25", "q75"), id.vars=c("set.fac", "center"))

gg <- ggplot()+
  ggtitle("Inter-quartile range (25-75%) of sizes of peaks up in all samples")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_text(aes(
    value, set.fac, color=center,
    hjust=ifelse(variable=="q75", 0, 1),
    label=paste0(
      ifelse(variable=="q75", " ", ""),
      round(value),
      ifelse(variable=="q75", "", " "))),
    vjust=-0.2,
    data=text.dt)+
  geom_text(aes(
    median, set.fac, color=center,
    label=sprintf("N=%d peaks up in %d samples", as.integer(peaks/2), max.samples)),
            vjust=1.3,
            data=stats.dt)+
  geom_point(aes(
    median, set.fac, color=center),
             shape=1,
             data=stats.dt)+
  geom_segment(aes(
    q25, set.fac,
    color=center,
    xend=q75, yend=set.fac),
               size=1,
               data=stats.dt)+
  ## geom_segment(aes(
  ##   min, set.fac,
  ##   color=center,
  ##   xend=max, yend=set.fac),
  ##              size=0.5,
  ##              data=stats.dt)+
  scale_x_log10("peak size (bases)", limits=c(1e1, 1e6), breaks=10^seq(2, 5, by=1))+
  ylab("data set")
print(gg)

png("figure-compare-size-centers-allup.png", 1000, 500, res=100)
print(gg)
dev.off()
system("cp figure-compare-size-centers-allup.png labels/CTCF_TDH_ENCODE/")
  

## experiment 3B: plot distribution of peak sizes of top 1000 peaks.
summary.tsv.vec <- Sys.glob("labels/*/peaks_summary.tsv")
some.vec <- grep("-", summary.tsv.vec, invert=TRUE, value=TRUE)
stats.dt <- data.table(summary.tsv=some.vec)[, {
  summary.dt <- fread(summary.tsv)
  no.input <- if("n.Input" %in% names(summary.dt)){
    summary.dt[n.Input==0]
  }else{
    summary.dt[, n.samples := n.samples.up]
    summary.dt[, loss.diff := group.loss.diff]
    summary.dt[!grepl("Input", str.groups.up)]
  }
  no.input[order(-loss.diff)][1:1000, list(
    set.name=basename(dirname(summary.tsv)),
    q25=quantile(peakBases, 0.25),
    median=as.numeric(median(peakBases)),
    q75=quantile(peakBases, 0.75),
             min=min(peakBases),
             max=max(peakBases),
               peaks=.N
  )]
}, by=summary.tsv]
set.sizes <- stats.dt[, list(mean=mean(median)), by=set.name][order(mean)]
stats.dt[, set.fac := factor(set.name, set.sizes$set.name)]
stats.dt[, samples := sub(".*_", "", set.name)]
stats.dt[, center := ifelse(samples=="BP", "BluePrint", ifelse(samples=="adipose", "unknown", "CEEHRC"))]
text.dt <- melt(stats.dt, measure.vars=c("q25", "q75"), id.vars=c("set.fac", "center"))

gg <- ggplot()+
  ggtitle("Inter-quartile range (25-75%) of sizes of top 1000 peaks")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_text(aes(
    value, set.fac, color=center,
    hjust=ifelse(variable=="q75", 0, 1),
    label=paste0(
      ifelse(variable=="q75", " ", ""),
      round(value),
      ifelse(variable=="q75", "", " "))),
    vjust=-0.2,
    data=text.dt)+
  geom_point(aes(
    median, set.fac, color=center),
             shape=1,
             data=stats.dt)+
  geom_segment(aes(
    q25, set.fac,
    color=center,
    xend=q75, yend=set.fac),
               size=1,
               data=stats.dt)+
  scale_x_log10("peak size (bases)", limits=c(1e1, 1e6), breaks=10^seq(2, 5, by=1))+
  ylab("data set")
print(gg)

png("figure-compare-size-centers-top1000.png", 1000, 500, res=100)
print(gg)
dev.off()
system("cp figure-compare-size-centers-top1000.png labels/CTCF_TDH_ENCODE/")
  

## peak size for all four
labels.bed.vec <- grep(
  "immune|other|ENCODE",
  Sys.glob("labels/H*_TDH_*/samples/*/*/joint_peaks.bedGraph"),
  value=TRUE)
sample.count.vec <- table(basename(dirname(dirname(dirname(dirname(labels.bed.vec))))))

sample.count.dt <- data.table(
  samples=as.integer(sample.count.vec),
  experiment=sub("_.*", "", names(sample.count.vec)),
  test.samples=sub(".*_", "", names(sample.count.vec)))
error.dt <- data.table(labels.bed=labels.bed.vec)[, {
  print(labels.bed)
  labels.dt <- fread(labels.bed)
  setnames(labels.dt, c("chrom", "chromStart", "chromEnd", "annotation"))
  sample.dir <- dirname(labels.bed)
  sample.id <- basename(sample.dir)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  proj.dir <- dirname(samples.dir)
  set.name <- basename(proj.dir)
  experiment <- sub("_.*", "", set.name)
  test.samples <- sub(".*_", "", set.name)
  train.samples <- ifelse(test.samples=="immune", "other", "immune")
  train.proj.dir <- paste0(experiment, "_TDH_", train.samples)
  cell.type <- basename(group.dir)
  jamm.narrowPeak <- file.path(
    "~",
    "JAMM",
    experiment,
    cell.type,
    "peaks",
    "filtered.peaks.narrowPeak")
  jamm <- fread(jamm.narrowPeak, select=1:3)
  setnames(jamm, c("chrom", "chromStart", "chromEnd"))
  joint_peaks.bedGraph <- file.path(
    "labels", train.proj.dir, "samples",
    cell.type, sample.id,
    "joint_peaks.bedGraph")
  DFilter_peaks.tsv <- file.path(
    "labels", train.proj.dir, "samples",
    cell.type, sample.id,
    "DFilter_peaks.tsv")
  DFilter.dt <- fread(paste("sed 's/ \\+/\t/g'", DFilter_peaks.tsv))
  PeakSeg.dt <- fread(joint_peaks.bedGraph)
  setnames(PeakSeg.dt, c("chrom", "chromStart", "chromEnd", "mean"))
  select.dt <- data.table(
    experiment, sample.id,
    model=paste0("macs2.", ifelse(experiment=="H3K4me3", "default", "broad")))
  pred.dt <- rbind(
    data.table(model="JAMM", jamm),
    DFilter.dt[, data.table(
      model="DFilter",
      chrom=chromosome,
      chromStart=`peak-start`,
      chromEnd=`peak-end`)],
    PeakSeg.dt[, data.table(model="PeakSegPipeline", chrom, chromStart, chromEnd)],
    macs.default.peaks[select.dt, list(
      model="MACS2", chrom, chromStart=peakStart, chromEnd=peakEnd)])
  error.dt <- pred.dt[, {
    PeakError::PeakError(.SD, labels.dt)
  }, by=list(model)]
  data.table(experiment, cell.type, test.samples, train.samples, sample.id, error.dt)
}, by=list(labels.bed)]

sample.totals <- error.dt[, {
  fp <- sum(fp)
  fn <- sum(fn)
  errors <- fp+fn
  data.table(
    fp, fn, errors,
    labels=.N,
    percent.errors=errors/.N,
    possible.fp=sum(possible.fp),
    possible.fn=sum(possible.tp))
}, by=list(experiment, test.samples, train.samples, model, sample.id)]
set.stats.dt <- sample.totals[, list(
  min.labels=min(labels), # over models, to check.
  max.labels=max(labels),
  total.labels=sum(labels),
  median=median(percent.errors),
  lo=quantile(percent.errors, 0.25),
  hi=quantile(percent.errors, 0.75),
  mean=mean(percent.errors),
  sd=sd(percent.errors)
  ), by=list(model, experiment, test.samples, train.samples)]
set.stats.dt[, stopifnot(min.labels==max.labels)]
text.dt <- set.stats.dt[model=="PeakSegPipeline"][sample.count.dt, on=list(experiment, test.samples)]
levs <- c("PeakSegPipeline", "MACS2", "JAMM", "DFilter", "")
text.dt[, model.fac := factor("", levs)]
set.stats.dt[, model.fac := factor(model, levs)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ train.samples + test.samples, labeller=label_both)+
  geom_point(aes(
    mean, model.fac),
             data=set.stats.dt)+
  geom_segment(aes(
    mean-sd, model.fac, xend=mean+sd, yend=model), data=set.stats.dt)+
  xlab("Percent incorrect labels (mean ± sd over labeled samples in test set)")+
  geom_text(aes(
    0.35, model.fac,
    label=sprintf("%d samples, %d labels/sample", samples, min.labels)
    ), data=text.dt, size=3)+
  ylab("model")
print(gg)
png("labels/H3K36me3_TDH_immune/figure-test-error.png", 600, 400, res=100)
print(gg)
dev.off()
pdf("labels/H3K36me3_TDH_immune/figure-test-error.pdf", 6, 4)
print(gg)
dev.off()
