############
#mean and sd for each basepair
 

fixedcov<-read.table("~/hctfixedcov/lengthmat178210",stringsAsFactors = F,fill = T,col.names = c(1:3000),row.names = 1)
fixedcov[is.na(fixedcov)]<-0
fixedcov<-remove_0(fixedcov)

tfixedcov<- apply(t(fixedcov), 1, function(x) x/max(x))

tfixedcov[is.na(tfixedcov)]<-0

fixedcovmat<-as.data.frame(rowSums(tfixedcov))
fixedcovmat<-cbind(c(1:27028),fixedcovmat)
fixedcovmat<-cbind(fixedcovmat,as.data.frame(rowSds(as.matrix(tfixedcov))))
fixedcovmat$type<-c("fixed")
colnames(fixedcovmat)<-c("bp","cov","sd","type")


livecov<-read.table("~/hctlivecov/lengthmat178210",stringsAsFactors = F,fill = T,col.names = c(1:5000),row.names = 1)
livecov[is.na(livecov)]<-0
livecov<-livecov[,colSums(livecov)!=0]

tlivecov<- apply(t(livecov), 1, function(x) x/max(x))

tlivecov[is.na(tlivecov)]<-0

livecovmat<-as.data.frame(rowSums(tlivecov))
livecovmat<-cbind(c(1:28767),livecovmat)

livecovmat<-cbind(livecovmat,as.data.frame(rowSds(as.matrix(tlivecov))))
livecovmat<-cbind(livecovmat,"live")
colnames(livecovmat)<-c("bp","cov","sd","type")


fixedcovmat$cov<-fixedcovmat$cov/max(fixedcovmat$cov)
livecovmat$cov<-livecovmat$cov/max(livecovmat$cov)
twocov<-rbind(fixedcovmat,livecovmat)
colnames(twocov)<-c("bp","cov","sd","type")
twocov$sd<-as.numeric(twocov$sd)
p <- ggplot(twocov, aes(x= bp, y= cov,group=type,colour=type)) +xlim(0,10000)+
  geom_line() + geom_point()+
  geom_errorbar(aes(ymin = cov-sd, ymax= cov+sd),width=0.01,position=position_dodge(0.05),alpha=.01)
p$labels$x<-c("distance from 3' end")
p$labels$y<-c("Percentile of mapped reads")

p+theme_bw() +
  scale_color_manual(values=c("#fcac0c","#2470a0"))+
  scale_alpha_manual(name='',values  = c('visible' = 0.8,'hidden' = 0.5))

p <- ggplot(twocov, aes(x= bp, y= cov,group=type,colour=type)) +
  geom_line()
