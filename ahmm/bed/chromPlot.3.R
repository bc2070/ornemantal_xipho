#根据ahmm之后生成的BED结果，使用chromplot作图；

#切换工作目录
setwd("/data2/projects/jgeng/Xiphophorus/ahmm_posteriro/ahmm_bed")

#导入bed目录
bed_list<-c("APS1_ahmm2bed","APS2_ahmm2bed","APS3_ahmm2bed","APS4_ahmm2bed","APS5_ahmm2bed",
            "BMK1_ahmm2bed","BMK2_ahmm2bed","BMK3_ahmm2bed","BMK4_ahmm2bed","BMK5_ahmm2bed","BMK6_ahmm2bed","BMK7_ahmm2bed","BMK8_ahmm2bed","BMK9_ahmm2bed","BMK10_ahmm2bed",
            "BML1_ahmm2bed","BML2_ahmm2bed","BML3_ahmm2bed","BML4_ahmm2bed","BML5_ahmm2bed","BML6_ahmm2bed",
            "HJ1_ahmm2bed","HJ2_ahmm2bed","HJ3_ahmm2bed","HJ4_ahmm2bed",
            "HMK1_ahmm2bed","HMK2_ahmm2bed","HMK3_ahmm2bed","HMK4_ahmm2bed","HMK5_ahmm2bed",
            "HW1_ahmm2bed","HW2_ahmm2bed","HW3_ahmm2bed","HW4_ahmm2bed","HW5_ahmm2bed",
            "MK1_ahmm2bed","MK2_ahmm2bed","MK3_ahmm2bed","MK4_ahmm2bed","MK5_ahmm2bed",
            "QMK1_ahmm2bed","QMK2_ahmm2bed","QMK3_ahmm2bed","QMK4_ahmm2bed","QMK5_ahmm2bed","QMK6_ahmm2bed",
            "SHJ1_ahmm2bed","SHJ2_ahmm2bed","SHJ3_ahmm2bed","SHJ4_ahmm2bed",
            "THJ2_ahmm2bed","THJ3_ahmm2bed","THJ4_ahmm2bed","THJ5_ahmm2bed","THJ6_ahmm2bed",
            "WMK1_ahmm2bed","WMK2_ahmm2bed","WMK3_ahmm2bed","WMK4_ahmm2bed","WMK5_ahmm2bed",
            "YS1_ahmm2bed","YS2_ahmm2bed","YS3_ahmm2bed","YS4_ahmm2bed","YS5_ahmm2bed","YS6_ahmm2bed","YS7_ahmm2bed","YS8_ahmm2bed","YS9_ahmm2bed","YS10_ahmm2bed","YS11_ahmm2bed","YS12_ahmm2bed","YS13_ahmm2bed","YS14_ahmm2bed","YS15_ahmm2bed","YS16_ahmm2bed","YS17_ahmm2bed","YS18_ahmm2bed")


library(chromPlot)
t1=read.table("/data2/projects/jgeng/Xiphophorus/ahmm_posteriro/chrom.txt",header = T,sep = "\t",fileEncoding = "utf-8",stringsAsFactors = F)


#根据bed_list循环操作
for (i in c(1:78)){
  a<-read.table(file=bed_list[i],sep = "\t",header = TRUE)
  n<-length(a$chrom)
  num<-c(1:n)
  col<-character(length = n)
  col[a$ancestry=="2,0,0"]<-"cadetblue3"
  col[a$ancestry=="1,1,0"]<-"chartreuse"
  col[a$ancestry=="1,0,1"]<-"coral3"
  col[a$ancestry=="0,2,0"]<-"darkgoldenrod"
  col[a$ancestry=="0,1,1"]<-"darkmagenta"
  col[a$ancestry=="0,0,2"]<-"forestgreen"
  col[a$ancestry=="N"]<-"grey17"
  col[is.na(a$ancestry)]<-"grey60"
  a$Colors<-col
  names(a)<-c("Chrom","Start","End","Ancestry","Colors")
  group<-character(length = n)
  group[a$Ancestry=="2,0,0"]<-"a.Xiphophorus_hellerii"
  group[a$Ancestry=="1,1,0"]<-"b.Xh+Xm"
  group[a$Ancestry=="1,0,1"]<-"c.Xh+Xv"
  group[a$Ancestry=="0,2,0"]<-"d.Xiphophorus_maculatus"
  group[a$Ancestry=="0,1,1"]<-"e.Xm+Xv"
  group[a$Ancestry=="0,0,2"]<-"f.Xiphophorus_variatus"
  group[a$Ancestry=="N"]<-"N"
  group[is.na(a$Ancestry)]<-"NA"

  a$Group<-group
#  a <- subset(a, Start != End)
  name_pdf<-paste(bed_list[i],"pdf",sep = ".")
  setwd("/data2/projects/jgeng/Xiphophorus/ahmm_posteriro/ahmm_bed/out_1")
  pdf(name_pdf)
  chromPlot(bands=a,chr=c(1:25), gaps=t1, 
            figCols=25)
  dev.off()
  setwd("/data2/projects/jgeng/Xiphophorus/ahmm_posteriro/ahmm_bed")
}


