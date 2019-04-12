### Compute data set to give to Isaac.
library(data.table)
macs.default.peaks <- fread("macs.default.peaks.csv")
joint_peaks.bedGraph.vec <- grep(
  "immune|other",
  Sys.glob("labels/H*_TDH_*/samples/*/*/joint_peaks.bedGraph"),
  value=TRUE)

joint.dt <- data.table(joint_peaks.bedGraph=joint_peaks.bedGraph.vec)[, {
  print(joint_peaks.bedGraph)
  sample.dir <- dirname(joint_peaks.bedGraph)
  sample.id <- basename(sample.dir)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  proj.dir <- dirname(samples.dir)
  set.name <- basename(proj.dir)
  experiment <- sub("_.*", "", set.name)
  PeakSeg.dt <- fread(joint_peaks.bedGraph, select=1:3)
  setnames(PeakSeg.dt, c("chrom", "peakStart", "peakEnd"))
  data.table(experiment, model="PeakSegPipeline", sample.id, PeakSeg.dt)
}, by=list(joint_peaks.bedGraph)]

all.peaks <- rbind(
  macs.default.peaks,
  joint.dt[, names(macs.default.peaks), with=FALSE])

hg19.contigs <- fread("hg19_contigs.txt", col.names=c("chrom", "contigStart", "contigEnd"))
setkey(hg19.contigs, chrom, contigStart, contigEnd)
setkey(all.peaks, chrom, peakStart, peakEnd)
(over.dt <- foverlaps(all.peaks, hg19.contigs, nomatch=0L))
rbind(nrow(over.dt), nrow(all.peaks))
peak.contigs <- over.dt[, list(
  peaks=.N,
  bases=contigEnd-contigStart
), by=list(model, sample.id, experiment, chrom, contigStart, contigEnd)]

fwrite(peak.contigs, "peak.contigs.csv")
