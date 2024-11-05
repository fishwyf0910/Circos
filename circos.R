#conda activate r4.3

library(stringr)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)

rm(list = ls())

setwd("/data01/wangyf/project2/CyprinusCarpio/15.pop/10.circos/")

seq_stat <-  read.delim("SPL01.genome.chromesome",header = T,sep = "\t",stringsAsFactors = FALSE)

out_dir = 'output'
if (!file.exists(out_dir)) dir.create(out_dir)

sample_name="Cyprinus_carpio"

pdf(str_c(out_dir, '/', sample_name, '.4.15.pdf'), width = 8, height = 8)

circos.clear()
circle_size = unit(1.2, 'snpc')
circos.par(#gap.degree =0,
  start.degree = 84.2,gap.after = c("A01" = 0.5, "A02" = 0.5, "A03" = 0.5, "A04" = 0.5, "A05" = 0.5, "A06" = 0.5,
                                    "A07" = 0.5, "A08" = 0.5, "A09" = 0.5, "A10" = 0.5, "A11" = 0.5, "A12" = 0.5,
                                    "A13" = 0.5, "A14" = 0.5, "A15" = 0.5, "A16" = 0.5, "A17" = 0.5, "A18" = 0.5,
                                    "A19" = 0.5, "A20" = 0.5, "A21" = 0.5, "A22" = 0.5, "A23" = 0.5, "A24" = 0.5,
                                    "A25" = 0.5,"B01" = 12, "B02" = 0.5, "B03" = 0.5, "B04" = 0.5, "B05" = 0.5, #此处调整亚基因组间gap
                                    "B06" = 0.5, "B07" = 0.5, "B08" = 0.5, "B09" = 0.5, "B10" = 0.5, "B11" = 0.5,
                                    "B12" = 0.5, "B13" = 0.5, "B14" = 0.5, "B15" = 0.5, "B16" = 0.5, "B17" = 0.5,
                                    "B18" = 0.5, "B19" = 0.5, "B20" = 0.5, "B21" = 0.5, "B22" = 0.5, "B23" = 0.5,
                                    "B24" = 0.5, "B25" = 0.5))

#1,Chr
chromosome.index = c(paste0("A0", c(1:9)),paste0("A", c(10:25)),
                     rev(paste0("B", c(10:25))),rev(paste0("B0", c(1:9))))

circos.initializeWithIdeogram(seq_stat, chromosome.index = chromosome.index,
                              plotType = c("labels", "axis"), track.height = 0.01,
                              major.by = 20000000)

set_track_gap(mm_h(0.01))
circos.genomicTrackPlotRegion(
  seq_stat, track.height = 0.05, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col ='#979797', border = '#979797', ...)
  } )

#2,special windows
window <- read.table("window.txt", header = T,sep = "\t")
window=na.omit(window)

circos.genomicTrack(window,ylim=c(0,1),panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, ytop= 0.95, ybottom = 0.05, col = '#FF0000',border='#DB7093',lwd=0.4)
}, stack = F,track.height = 0.05)

#3,SNPs density
snpden <- read.table("non-overlap-window.snpden.0.4-2.65", header = T,sep = "\t")
snpden=na.omit(snpden)

col_fun = colorRamp2(c(0,3), c("blue", "red"))
col=rainbow(5)

circos.genomicHeatmap(snpden, col = col_fun, border = "white",heatmap_height = 0.05)

:<<!
circos.genomicTrack(data=snpden,panel.fun=function(region,value,...) {
  circos.genomicHeatmap(snpden, col = col_fun, border = NA)
},track.height=0.08,bg.border=F)

circos.genomicTrack(data=snpden,panel.fun=function(region,value,...) {
circos.genomicHeatmap(snpden, col = col_fun, border = NA)
# circos.genomicLines(region,value,type="l",col="#1a237e",lwd=0.6,alpha=1)
},track.height=0.08,bg.border=F)
!
:<<!

#4,PI
PI <- read.table("non-overlap-window.pi", header = T,sep = "\t")
PI=na.omit(PI)

circos.genomicTrack(data=PI,ylim=c(0,0.01),panel.fun=function(region,value,...) {
 circos.genomicLines(region,value,type="l",col="#FF8C00",lwd=0.8,alpha=1)
},track.height=0.08,bg.border=F)

#5,XJ-ref fst
fst <- read.table("non-overlap-window.xj-ref.fst", header = T,sep = "\t")
fst=na.omit(fst)

circos.genomicTrack(data=fst,ylim=c(0,0.12),panel.fun=function(region,value,...) {
  circos.genomicLines(region,value,type="l",col="#008000",lwd=0.6,alpha=1)
},track.height=0.08,bg.border=F)

#6,Tajimas D
tajimaD <- read.table("non-overlap-window.tajimad", header = T,sep = "\t")
tajimaD=na.omit(tajimaD)

circos.genomicTrack(data=tajimaD,ylim=c(-0.003,6.2),panel.fun=function(region,value,...) {
  circos.genomicLines(region,value,type="l",col="#8A2BE2",lwd=0.5,alpha=1)
},track.height=0.08,bg.border=F)
!
################################################################################
#南方群体和额河群体的适应
#XP-CLR
xpclr <- read.table("south-Irtysh/overlap-window.xpclr", header = T,sep = "\t")
xpclr=na.omit(xpclr)

circos.genomicTrack(data=xpclr,ylim=c(-1,20),panel.fun=function(region,value,...) {
  circos.genomicLines(region,value,type="l",col="#8A2BE2",lwd=0.5,alpha=1)
},track.height=0.08,bg.border=F)

#tajima-south
tajimaDsouth <- read.table("south-Irtysh/southpop.window.Tajima.D.snp40-265.3-4", header = T,sep = "\t")
tajimaDsouth=na.omit(tajimaDsouth)

circos.genomicTrack(data=tajimaDsouth,ylim=c(-3,4),panel.fun=function(region,value,...) {
  circos.genomicLines(region,value,type="l",col="#8A2BE2",lwd=0.5,alpha=1)
},track.height=0.08,bg.border=F)

dev.off()
