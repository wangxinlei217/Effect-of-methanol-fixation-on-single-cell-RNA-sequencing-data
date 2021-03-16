######

#live_ratio and fixed_ratio are files describing mappped reads of each transcript in all cells from either kind
#coverage ratio calculation
#for single cell coverage variation
library(matrixStats)
load("ratio.RData")

live_ratio[is.na(live_ratio)]<-0
live_ratio<-live_ratio[rowSums(live_ratio[,-1])!=0,]
live_ratio<-live_ratio[,colSums(live_ratio[,-1])!=0]
rownames(live_ratio)<-live_ratio$X1
live_ratio<-live_ratio[,-1]

fixed_ratio[is.na(fixed_ratio)]<-0
fixed_ratio<-fixed_ratio[rowSums(fixed_ratio[,-1])!=0,]
fixed_ratio<-fixed_ratio[,colSums(fixed_ratio[,-1])!=0]
rownames(fixed_ratio)<-fixed_ratio$X1
fixed_ratio<-fixed_ratio[,-1]

#####mean without including zeros,
length<-fixed_ratio[,c(1,2)]
colnames(length)<-c("trans","length")


fixed_ratio$sum<-rowSums(fixed_ratio)
fixed_ratiotf<-tf_trans(fixed_ratio,1)
fixed_ratio$freq<-rowSums(fixed_ratiotf)
fixed_ratio$ave<-fixed_ratio$sum/fixed_ratio$freq

fixmap<-as.data.frame(fixed_ratio$ave)
rownames(fixmap)<-rownames(fixed_ratio)




live_ratio$sum<-rowSums(live_ratio)
live_ratiotf<-tf_trans(live_ratio,1)
live_ratio$freq<-rowSums(live_ratiotf)
live_ratio$ave<-live_ratio$sum/live_ratio$freq
livemap<-as.data.frame(live_ratio$ave)
rownames(livemap)<-rownames(live_ratio)


bothmap<-merge(livemap,fixmap,by="row.names")
colnames(bothmap)<-c("trans","live","fixed")
bothmap$ratio<-bothmap$live/bothmap$fixed
bothmap1<-merge(bothmap,length,by="trans")
bothmap1$live1<-bothmap1$live/bothmap1$length
bothmap1$fixed1<-bothmap1$fixed/bothmap1$length
bothmap1$length<-log(bothmap1$length)
bothmap1 <- bothmap1[!is.infinite(bothmap1$live)&!is.infinite(bothmap1$fixed),]
bothmap1$group<-c(rep(1:9,each=13772),rep(10,each=13771))
bothmap1$group<-as.factor(bothmap1$group)



########
#Combining length and GC content in the same matrix
#sorting and indexing all transcripts according to gc and length
#load gc_length file to load in a file with gc and length information of each transcript
load(gc_length.RData)
gc_length<-gc_length[order(gc_length$gcscore,decreasing = F),]
gc_length$gc_ranking<-c(rep(1:99,each=1980),rep(100,time=1982))
gc_length<-gc_length[order(gc_length$em_length,decreasing = F),]
gc_length$length_ranking<-c(rep(1:99,each=1980),rep(100,time=1982))

#making matrix for combining length and gc content

gc_length_mat<-as.data.frame(matrix(nrow=100,ncol=100))

#fill in the matrix
for (i in 1:100) {

  for (j in 1:100) {
    gc_length_mat[i,j]<-mean(bothmap[bothmap$trans%in%gc_length[gc_length$gc_ranking==i&gc_length$length_ranking==j,1],4])

  }
}


is.na(gc_length_mat) <- do.call(cbind,lapply(gc_length_mat, is.infinite))
gc_length_mat[is.na(gc_length_mat)]<-1
gc_length_mat<-log(gc_length_mat)

is.na(gc_length_mat) <- do.call(cbind,lapply(gc_length_mat, is.infinite))
gc_length_mat[is.na(gc_length_mat)]<-0

pheatmap::pheatmap(gc_length_mat,cluster_rows = F,cluster_cols = F,border=NA,color= magma(10),show_rownames = F,show_colnames = F)
