.libPaths("/moto/palab/users/zlf2101/rpackages/")
library(SNPRelate)
require(sp)

##READ INDEX NUMBER FOR JACKKNIFE TRIAL
args <- commandArgs(TRUE)
trial_index <- args[1]

##MAKE COVARIATES FOR BLEACHING GWAS - CV JACKNIFE
##READ IN SAMPLE DATA
#full_pheno<-read.table("~/Amil.chl.missing_rmvd.plink.pheno")
full_pheno<-read.table("amil_nomissing.zoox.pheno")
names(full_pheno)<-c("IID","FID","chl_content")
full_pheno$bleached<-ifelse(full_pheno$chl_content>3,2,1)
full_pheno$Reef<-substr(full_pheno$IID,start=1,stop=2)

##SPLIT DATA INTO TRAIN/TEST
smp_siz = floor(0.90 * nrow(full_pheno))
train_ind = sample(seq_len(nrow(full_pheno)),size=smp_siz)
train = full_pheno[train_ind,]
test = full_pheno[-train_ind,]

train = full_pheno

##MAKE ENVIRONMENTAL PCS FOR COVARIATES
train_reefs<-unique(train$Reef)
test_reefs<-unique(test$Reef)

coords<-read.csv("data/reef_gps_coords.csv",colClasses="character")
coords$Lat<-sub('\xa1', 'd', coords$Lat)
coords$Lat<-as.numeric(sp::char2dms(coords$Lat))
coords$Lon<-sub('\xa1', 'd', coords$Lon)
coords$Lon<-as.numeric(sp::char2dms(coords$Lon))
coords<-coords[order(coords$Reef),]

names(coords)<-c("Reef","Latitude","Longitude")
env_dat<-read.csv("data/Environmental_data.csv")
dhw_dat<-read.csv("data/DHW_data.csv")

datalist = list()
for(i in 1:nrow(coords)){
  dists<-( (coords[i,]$Longitude-env_dat$lon)*cos(90*env_dat$lat ) )^2 + (coords[i,]$Latitude-env_dat$lat)^2
  dhw_dists<-( (coords[i,]$Longitude-dhw_dat$lon)*cos(90*dhw_dat$lat ) )^2 + (coords[i,]$Latitude-dhw_dat$lat)^2
  env_row<-env_dat[ which.min(dists),]
  dhw_row<-dhw_dat[ which.min(dhw_dists),]$annMaxDHW_2017
  env_row$DHW_max<-dhw_row
  datalist[[i]] <- env_row 
}

reef_dat = do.call(rbind, datalist)
reef_dat<-cbind(coords,reef_dat)
reef_dat[,c("PIXEL_ID","lon","lat","REEF_ID","REEF_NAME")]<-NULL
head(reef_dat)

drop_cols<-c("pop","FAMID","SECTOR","CROSS_SHELF","Comp.1","Comp.2","Comp.3","Comp.4","Comp.5","DHW_max","Latitude","Longitude","gcta_fid")
env_covars_names<-c("CRS_NO3_AV","CRS_NO3_SR","CRS_PO4_AV","CRS_PO4_SR","CRS_O2_AV","CRS_O2_SR","CRS_S_AV","CRS_S_SR","CRS_T_AV","CRS_T_SR","CRS_SI_AV","CRS_SI_SR",
                    "GA_BATHY","GA_SLOPE","GA_ASPECT","GBR_BATHY","GA_CRBNT","GA_GRAVEL","GA_SAND","GA_MUD","GMCS_STRESS_TMN","GMCS_STRESS_IQR",
                    "SW_CHLA_AV","SW_CHLA_SR","SW_K490_AV","SW_K490_SR","SW_BIR_AV","SW_BIR_SR","MT_SST_AV","MT_SST_SR","MT_SST_MIN","Primary","Secondary","Tertiary","Latitude","Longitude",
                    "DHW_max","mindistbar","mindistcoa")
rownames(reef_dat)<-reef_dat$Reef

env_pca<-prcomp(reef_dat[,env_covars_names], center = TRUE,scale. = TRUE)
env_pca_loadings<-data.frame(env_pca$x)
env_pca_loadings$Reef<-as.vector(as.character(rownames(env_pca_loadings)))

##READ IN DATA FOR MISSING READ PROPORTION
missing_p<-read.table("data/missing.smiss")
head(missing_p)
missing_p<-missing_p[,c("V1","V4")]
names(missing_p)<-c("IID","prop_missing")

train_covars<-mermage(train, missing_p, by="IID")perm
test_covars<-merge(test, missing_p, by="IID")

train_covars<-merge(train_covars,env_pca_loadings,by="Reef",all.x=T)
test_covars<-merge(test_covars,env_pca_loadings,by="Reef",all.x=T)

locov_read_counts<-read.table("data/locov_merged_read_counts.tsv",header=T)
hicov_read_counts<-read.table("data/merged_read_counts.tsv",header=T)
scaff_lengths<-read.table("data/scaffold_names_lenghts.tsv")
names(scaff_lengths)<-c("scaff","length")
amil_scaffs<-scaff_lengths$scaff[!(scaff_lengths$scaff %in% c("symbiont1","symbiont2","symbiont3","symbiont4","Acropora_millepora_mitochondria"))]
locov_reads<-(locov_read_counts[locov_read_counts$scaff %in% amil_scaffs,2:length(locov_read_counts)])
hicov_reads<-(hicov_read_counts[hicov_read_counts$scaff %in% amil_scaffs,2:length(hicov_read_counts)])
locov_total_read_sums<-colSums(locov_reads)
hicov_total_read_sums<-colSums(hicov_reads)
total_read_counts<-c(hicov_total_read_sums,locov_total_read_sums)
total_read_counts<-data.frame(as.table(total_read_counts))
head(total_read_counts)
names(total_read_counts)<-c("ID","total_read_count")

train_covars<-merge(train_covars,total_read_counts,by.y="ID",by.x="IID")
test_covars<-merge(test_covars,total_read_counts,by.y="ID",by.x="IID")
head(test_covars)
quantNorm =function(x){qnorm(rank(x,ties.method = "average")/(length(x)+1))}
train_covars$norm_read_count<-quantNorm(train_covars$total_read_count)
all_read_counts<-rbind(train_covars[,c("IID","total_read_count")],test_covars[,c("IID","total_read_count")])
all_read_counts$norm_read_counts<-quantNorm(all_read_counts$total_read_count)
test_norm_reads<-all_read_counts[all_read_counts$IID %in% test_covars$IID,]
test_covars$norm_read_count<-test_norm_reads$norm_read_counts

##GENETIC PCS
genofile<-snpgdsOpen("data/Amil.full_samps.vcf.gds", readonly=FALSE, allow.duplicate = T)
train_snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,sample.id = as.factor(train_covars$IID))
train_snpset.id <- unlist(train_snpset)
train_pca <- snpgdsPCA(genofile, snp.id=train_snpset.id, num.thread=2,sample.id = as.factor(train_covars$IID),need.genmat=TRUE)
train_pca_res<-prcomp(train_pca$genmat)


train_tab <- data.frame(sample.id = train_pca$sample.id,
                  EV1 = train_pca$eigenvect[,1],    # the first eigenvector
                  EV2 = train_pca$eigenvect[,2],
                  EV3 = train_pca$eigenvect[,3],
                  EV4 = train_pca$eigenvect[,4],# the second eigenvector
                  stringsAsFactors = FALSE)
head(train_tab)
train_covars<-merge(train_covars,train_tab,by.x="IID",by.y="sample.id")

train_snp_loadings<-snpgdsPCASNPLoading(train_pca, genofile)
test_loadings<-snpgdsPCASampLoading(train_snp_loadings,genofile,sample.id=as.factor(test_covars$IID))
test_tab <- data.frame(sample.id = test_loadings$sample.id,
                       EV1 = test_loadings$eigenvect[,1],    # the first eigenvector
                       EV2 = test_loadings$eigenvect[,2],
                       EV3 = test_loadings$eigenvect[,3],
                       EV4 = test_loadings$eigenvect[,4],# the second eigenvector
                       stringsAsFactors = FALSE)
test_covars<-merge(test_covars,test_tab,by.x="IID",by.y="sample.id")

##ADD DATES

sample_dates<-read.table("data/Amil_sample_dates.txt",sep="\t",header=T)

train_covars<-merge(sample_dates[c("SampleID","Date","Field.ID")],train_covars,by.x = "SampleID",by.y="IID")
train_covars$Date<-as.character(train_covars$Date)
train_covars$Date<-as.numeric(as.Date(train_covars$Date, format = "%m/%d/%y",origin="1970-01-01 00:00:00 +10:00"))
min_date<-min(train_covars$Date)
train_covars$Date<-train_covars$Date - min_date

test_covars<-merge(sample_dates[c("SampleID","Date","Field.ID")],test_covars,by.x = "SampleID",by.y="IID")
test_covars$Date<-as.character(test_covars$Date)
test_covars$Date<-as.numeric(as.Date(test_covars$Date, format = "%m/%d/%y",origin="1970-01-01 00:00:00 +10:00"))
test_covars$Date<-test_covars$Date - min_date

##ADD DEPTH
depth_data<-read.table("data/Amil_depth.txt",sep="\t",header=T)
colnames(train_covars)[colnames(train_covars) == "SampleID"] <- "IID"
train_covars<-merge(depth_data,train_covars,by="IID")

colnames(test_covars)[colnames(test_covars) == "SampleID"] <- "IID"
test_covars<-merge(depth_data,test_covars,by="IID")

##ADD SYMBIONT DATA
chl_reg<-read.table("data/predicted2_chl_vals.tsv",sep="\t",header=T)
colnames(chl_reg)[colnames(chl_reg) == "SampleID"] <- "IID"
train_covars<-merge(train_covars,chl_reg,by="IID")
test_covars<-merge(test_covars,chl_reg,by="IID")

##WRITE FILES
write.table(train_covars[,c("FID","IID","Dom_symbiont","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","EV1","EV2","norm_read_count","prop_missing","Depth","Date")],paste("cv_input/train_covars_",trial_index,".txt",sep=""),row.names=F,quote=F)
write.table(test_covars[,c("FID","IID","Dom_symbiont","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","EV1","EV2","norm_read_count","prop_missing","Depth","Date")],paste("cv_input/test_covars_",trial_index,".txt",sep=""),row.names=F,quote=F)

write.table(train_covars[,c("FID","IID","bleached")],paste("cv_input/train_pheno_",trial_index,".txt",sep=""),quote=F,row.names = F,col.names = F)

write.table(train_covars[,c("FID","IID")],paste("cv_input/plink19_gwas_set_",trial_index,".txt",sep=""),quote=F,row.names = F,col.names = F)
write.table(test_covars[,c("FID","IID")],paste("cv_input/plink19_validation_set_",trial_index,".txt",sep=""),quote=F,row.names = F,col.names = F)

write.table(train_covars[,c("IID")],paste("cv_input/plink2_gwas_set_",trial_index,".txt",sep=""),quote=F,row.names = F,col.names = F)
write.table(test_covars[,c("IID")],paste("cv_input/plink2_validation_set_",trial_index,".txt",sep=""),quote=F,row.names = F,col.names = F)





