library(data.table)
library(ggplot2)
macs.default.peaks <- fread("macs.default.peaks.csv")
macs.default.peaks[, peakBases := peakEnd - peakStart]
setkey(macs.default.peaks, experiment, sample.id, model, chrom, peakStart, peakEnd)

## test error for both default and broad models for each data set.
labels.bed.vec <- grep(
  "immune|other",
  Sys.glob("labels/H*_TDH_*/samples/*/*/labels.bed"),
  value=TRUE)
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
  joint_peaks.bedGraph <- file.path(
    "labels", train.proj.dir, "samples",
    basename(group.dir), sample.id,
    "joint_peaks.bedGraph")
  if(file.exists(joint_peaks.bedGraph)){
    PeakSeg.dt <- fread(joint_peaks.bedGraph)
    setnames(PeakSeg.dt, c("chrom", "chromStart", "chromEnd", "mean"))
    select.dt <- data.table(experiment, sample.id)
    pred.dt <- rbind(
      PeakSeg.dt[, data.table(model="PeakSeg", chrom, chromStart, chromEnd)],
      macs.default.peaks[select.dt, list(
        model, chrom, chromStart=peakStart, chromEnd=peakEnd)])
    error.dt <- pred.dt[, {
      PeakError::PeakError(.SD, labels.dt)
    }, by=list(model)]
    data.table(experiment, test.samples, train.samples, sample.id, error.dt)
  }
}, by=list(labels.bed)]

error.totals <- error.dt[, {
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

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ .)+
  geom_point(aes(
    percent.errors, model, color=test.samples),
             data=error.totals)


ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ test.samples)+
  geom_point(aes(
    percent.errors, model),
             data=error.totals)

stats.dt <- error.totals[, list(
  median=median(percent.errors),
  lo=quantile(percent.errors, 0.25),
  hi=quantile(percent.errors, 0.75),
  mean=mean(percent.errors),
  sd=sd(percent.errors)
  ), by=list(model, experiment, test.samples, train.samples)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ test.samples)+
  geom_point(aes(
    median, model),
             data=stats.dt)+
  geom_segment(aes(
    lo, model, xend=hi, yend=model), data=stats.dt)

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ train.samples + test.samples, labeller=label_both)+
  geom_point(aes(
    mean, model),
             data=stats.dt)+
  geom_segment(aes(
    mean-sd, model, xend=mean+sd, yend=model), data=stats.dt)+
  xlab("Percent incorrect labels (mean ± sd over labeled samples in test set)")
png("labels/H3K36me3_TDH_immune/figure-test-error-macs-default-broad.png", 800, res=100)
print(gg)
dev.off()

## test error for only one of default or broad model for each data set.
labels.bed.vec <- grep(
  "immune|other",
  Sys.glob("labels/H*_TDH_*/samples/*/*/labels.bed"),
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
