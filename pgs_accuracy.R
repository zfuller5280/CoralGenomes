vs_scores<-vector()
pred_tile_list=list()
#Intital plot (just to get axes)
plot(pred_tile_df$percentile_10, pred_tile_df$chl_content,type="n",xlim=c(1,10),ylim=c(0,6),xlab="Predicted Phenotype Percentile",ylab="Mean Visual Score")
for(n in seq(1,100)){
  print(n)
  #Split for 85-15 train/test
  setwd("/Users/zachfuller/cv_dosages_vs_85_XD")
  train_covars<-read.table(paste("train_covars_",n,".txt",sep=""),header=T)
  test_covars<-read.table(paste("test_covars_",n,".txt",sep=""),header=T)
  train_pheno<-read.table(paste("train_pheno_",n,".txt",sep=""),header=F)
  names(train_pheno)<-c("IID","FID","chl_content")
  file = paste("cv_jk_",n,"_0.00001_bleaching_dosage_scores.sscore",sep="")
  preds<-read.table(file,header=T)
  names(preds)<-c("IID","FID","n","total","avg","SCORESUM")
  pred_df<-merge(preds,test_covars,by.x="IID",by.y="IID")
  pred_df<-merge(pred_df,full_pheno,by="IID")
  train_file = paste("cv_jk_",n,"_0.00001_bleaching_gwas_scores.sscore",sep="")
  train_preds<-read.table(train_file,header=T)
  names(train_preds)<-c("IID","FID","n","total","avg","SCORESUM")
  train_df<-merge(train_covars,train_pheno,by="IID")
  train_df<-merge(train_df,train_preds,by="IID",by.y="IID")

  pgs_mdl<-lm(chl_content ~  X.D + PC1 + PC2 + PC3 + PC4 + Date + Depth
              + (HC_batch2) + (LC_batch) + EV1 + EV2 + SCORESUM,data=train_df)
  pred_df$percentile<-ntile(pred_df$SCORESUM,4)
  pred_df$predicted_vs<-predict(pgs_mdl, pred_df)

  pred_df$predicted_quartile<-ntile(pred_df$predicted_vs,4)
  pred_df$predicted_cdf<-ecdf(pred_df$predicted_vs)(pred_df$predicted_vs)

  pred_df$vs_tile<-ecdf(pred_df$chl_content)(pred_df$chl_content)

  pred_df$percentile_10<-ntile(pred_df$predicted_vs,10)
  pred_tile_df<-aggregate(chl_content ~ percentile_10, pred_df, mean)
  points(pred_tile_df$percentile_10, pred_tile_df$chl_content,col="gray",type="b")
  pred_tile_list[[n]]<-pred_tile_df

  vs_scores[n]<-list(pred_df$vs_tile)

  perc1[n]<-mean(pred_df[pred_df$percentile==1,"chl_content"])
  perc2[n]<-mean(pred_df[pred_df$percentile==2,"chl_content"])
  perc3[n]<-mean(pred_df[pred_df$percentile==3,"chl_content"])
  perc4[n]<-mean(pred_df[pred_df$percentile==4,"chl_content"])

}
pred_tiles<-do.call(rbind, pred_tile_list)
pred_tile_agg<-aggregate(chl_content ~ percentile_10, pred_tiles, mean)

wilcox.test(pred_tiles[pred_tiles$percentile_10==1,"vs_tile"],pred_tiles[pred_tiles$percentile_10==10,"vs_tile"])
#Plot final mean for each decile
points(pred_tile_agg$percentile_10, pred_tile_agg$chl_content, col="blue", pch=16, cex=2)
