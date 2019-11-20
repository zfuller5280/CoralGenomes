#Test for differences in highest and lowest quartile of PGS
perc1<-vector()
perc2<-vector()
perc3<-vector()
perc4<-vector()
full_pheno<-read.table("~/visual_score.pheno")
names(full_pheno)<-c("IID","FID","chl_content")
for(n in seq(1,100)){
  setwd("/Users/zachfuller/cv_dosages_vs_85_XD")
  file = paste("cv_jk_",n,"_0.00001_bleaching_dosage_scores.sscore",sep="")
  preds<-read.table(file,header=T)
  names(preds)<-c("IID","FID","n","total","avg","SCORESUM")
  pred_df<-merge(pred_df,full_pheno,by="IID")
  pred_df$percentile<-ntile(pred_df$SCORESUM,4)
  perc1[n]<-mean(pred_df[pred_df$percentile==1,"chl_content"])
  perc2[n]<-mean(pred_df[pred_df$percentile==2,"chl_content"])
  perc3[n]<-mean(pred_df[pred_df$percentile==3,"chl_content"])
  perc4[n]<-mean(pred_df[pred_df$percentile==4,"chl_content"])
}
perc1<-data.frame(percentile=rep(1,length(perc1)),avgs=perc1)
perc2<-data.frame(percentile=rep(2,length(perc2)),avgs=perc2)
perc3<-data.frame(percentile=rep(3,length(perc3)),avgs=perc3)
perc4<-data.frame(percentile=rep(4,length(perc4)),avgs=perc4)
pgs_percentiles<-rbind(perc1,perc2,perc3,perc4)
boxplot(avgs ~ percentile, data=pgs_percentiles)
stripchart(avgs ~ percentile, vertical = TRUE, data = pgs_percentiles, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
wilcox.test(perc1,perc4)
plot(c(1,2),c(perc1[1],perc4[1]),ylim=c(1,4),type="n")
for(i in seq(1:length(perc1))){
  lines(c(1,2),c(perc1[i],perc4[i]),type="b",ylim=c(1,4),add=T,col="gray")
}
boxplot(perc1,perc4,add=T)