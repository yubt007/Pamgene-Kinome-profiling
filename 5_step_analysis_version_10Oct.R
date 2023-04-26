rm(list=ls())
#setwd('~/Desktop/HP')

### ---------------
###
### Create: Bingting Yu
### Date: 2020-Oct-10
### Email: ybting@gamil.com
### 
### Update Log: 2020-Oct-10  version 1.0
### MDL/Erasmus MC
### ---------------

library("tidyr")
library('dplyr')
library('ggplot2')
library(ggplot2)
library(grid)
library(stringr)
library(pheatmap)
library(limma)
library(reshape2)

Raw_data <-read.csv('',header=F)  ##read your Pamgene Kinome assay data
Raw_data <-Raw_data[-1,]. ##remove the chip information.


##Separate the data by different Exposure Time
step0_split_data<-function(data){
  data<-Raw_data ##checking intermediate data
  l_data<-as.numeric(ncol(data))
  r_data<-as.numeric(nrow(data))
  
  Exposure_time<-as.matrix(data[3,5:l_data])
  Exposure_time<-as.numeric(Exposure_time[1,])
  
  gene_name_matrix<-data[,1:4]
  
  data_N<-data[,5:l_data]
  
  data_10s<-cbind(gene_name_matrix,data_N[,which(Exposure_time==10)])
  data_50s<-cbind(gene_name_matrix,data_N[,which(Exposure_time==50)])
  data_200s<-cbind(gene_name_matrix,data_N[,which(Exposure_time==200)])
  
  write.table(data_10s,file = 'data_10s.txt',sep='\t',row.names=FALSE,col.names = FALSE)
  write.table(data_50s,file = 'data_50s.txt',sep='\t',row.names=FALSE,col.names = FALSE)
  write.table(data_200s,file = 'data_200s.txt',sep='\t',row.names=FALSE,col.names = FALSE)
  
}
step0_split_data(data= Raw_data)

#Raw_data<-read.csv('data_200s.txt',sep='\t',header = F)

##selected different Exposure Time to start the whole analysis
start<-function(exposure_time){
  ##based on Exposure Time to create a fold saving results
  exposure_time<-exposure_time
  data_name<-paste('data_',exposure_time,'.txt',sep = '')
  Raw_data<-read.csv(data_name,sep='\t',header = F)
  floder_name<-paste(exposure_time,'_results',sep='')
  dir.create(floder_name)
  setwd(floder_name)
  
  ##Change the unit of circle to during time
  
  step1_change_unit<-function(data){
    #data<-Raw_data ##checking intermediate data
    
    l_data<-as.numeric(ncol(data))
    r_data<-as.numeric(nrow(data))
    sample_name<-as.matrix(data[1,5:l_data])
    
    cycle<-as.matrix(data[2,5:l_data])
    cycle<-gsub(pattern = '\\<32\\>',replacement = '1920', cycle)
    cycle<-gsub(pattern = '\\<37\\>',replacement = '2250', cycle)
    cycle<-gsub(pattern = '\\<42\\>',replacement = '2580', cycle)
    cycle<-gsub(pattern = '\\<47\\>',replacement = '2910', cycle)
    cycle<-gsub(pattern = '\\<52\\>',replacement = '3240', cycle)
    cycle<-gsub(pattern = '\\<57\\>',replacement = '3570', cycle)
    cycle<-gsub(pattern = '\\<62\\>',replacement = '3900', cycle)
    cycle<-gsub(pattern = '\\<67\\>',replacement = '4230', cycle)
    cycle<-gsub(pattern = '\\<72\\>',replacement = '4560', cycle)
    cycle<-gsub(pattern = '\\<77\\>',replacement = '4890', cycle)
    cycle<-gsub(pattern = '\\<82\\>',replacement = '5220', cycle)
    cycle<-gsub(pattern = '\\<87\\>',replacement = '5550', cycle)
    cycle<-gsub(pattern = '\\<92\\>',replacement = '5880', cycle)
    cycle<-gsub(pattern = '\\<94\\>',replacement = '6210', cycle)
    
    cycle_un94<-which(cycle!=6210)
    
    re<-as.matrix(data[4,5:l_data])
    
    Treatment<-as.matrix(data[5,5:l_data])
    names<-paste(sample_name,'_',Treatment,',',re,',',cycle)
    rownames_data<-as.matrix(unlist(Raw_data[7:r_data,1]))[,1]
    
    data<-data[7:r_data,5:l_data]
    #data<-matrix(as.numeric(unlist(data)),nrow=nrow(data)) ##change all of values from character to numeric
    data<-as.matrix(data)
    data[is.na(data)]=0
    
    data<-apply(data,2,as.numeric)
    data[data<0]=0##change all of negative values to 0
    
    data<-(data*2E+12)/6E+23 ##recalculate data in counts to data in moles
    
    colnames(data)<-names
    row.names(data)<-rownames_data
    
    data<-data[,cycle_un94]
    
    write.csv(data,file = 'unit_changed_data.csv')
  }
  step1_change_unit(data = Raw_data)
  
  data_unit_changed<-read.csv(file = 'unit_changed_data.csv',sep=',',header = F)
  ##normalization the data
  step2_normalzation_data<-function(data){
    #data<-data_unit_changed
    colnames(data)<-unlist(data[1,])
    rownames(data)<-data[,1]
    data<-data[-1,-1]
    names_data<-str_split(colnames(data),pattern=',',4,simplify=T)

    
    data_temp<-apply(data,2,as.numeric)
    
    
    rownames(data_temp)<-rownames(data)
    
    data_temp_sum<-apply(data_temp, 2, sum)
    sum_standard<-c()
    for (i in 1:length(table(names_data[,1]))) {
      sum_sample<-0
      for (j in 1:length(table(names_data[,3]))) {
        n<-i*j
        sum_sample<-sum_sample+data_temp_sum[n]
        
      }
      sum_standard<-append(sum_standard,sum_sample)
    }
    
    times<-length(table(names_data[,3]))
    sum_standard_all<-rep(sum_standard,each=times)
    data_temp_standard<-rbind(sum_standard_all,data_temp)
    colnames(data_temp_standard)<-colnames(data_temp)
    write.table(data_temp_standard,file='pre_sum_normalzation.txt',sep='\t')
    data_temp_normalization<-apply(data_temp_standard,2,function(x) x/x[1]*100)
    data_temp_normalization_pre_dooule<-apply(data_temp_standard,2,function(x) x/x[1])
    write.csv(data_temp_normalization,file = 'data_normalization_by_sum.csv')
    
    ##choose the standard kinase to normalization the data
    ART_003<-data_temp[11,]
    sum_ART003<-c(sum(ART_003[1:13]),sum(ART_003[14:26]),sum(ART_003[27:39]),sum(ART_003[40:52]))
    ART_004<-data_temp[12,]
    sum_ART004<-c(sum(ART_004[1:13]),sum(ART_004[14:26]),sum(ART_004[27:39]),sum(ART_004[40:52]))
    sum_ART003_all<-rep(sum_ART003,each=times)
    data_temp_standard<-rbind(sum_ART003_all,data_temp)
    colnames(data_temp_standard)<-colnames(data_temp)
    data_temp_normalization_BYART3<-apply(data_temp_standard,2,function(x) x/x[1]*10)
    write.csv(data_temp_normalization_BYART3,file = 'data_normalization_by_ART003.csv')
  
    sum_ART004_all<-rep(sum_ART004,each=times)
    data_temp_standard<-rbind(sum_ART004_all,data_temp)
    colnames(data_temp_standard)<-colnames(data_temp)
    data_temp_normalization_BYART4<-apply(data_temp_standard,2,function(x) x/x[1]*10)
    write.csv(data_temp_normalization_BYART4,file = 'data_normalization_by_ART004.csv')
    
    data_temp_standard<-rbind(sum_ART003_all,data_temp_normalization_pre_dooule)
    colnames(data_temp_standard)<-colnames(data_temp)
    data_temp_normalization_BYART3<-apply(data_temp_standard,2,function(x) x/x[1]*100)
    write.csv(data_temp_normalization,file = 'data_normalization_by_sum_ART3.csv')
    
    data_orgin<-apply(data_temp,2,function(x) x*100*1000000)
    write.csv(data_orgin,file='data_orgin.csv')
    }
  step2_normalzation_data(data = data_unit_changed)
  
  ##calculated the AUC
  step3_data_analysis<-function(data){
    #data<-data_unit_changed_50s
    #time<-50
    
    colnames(data)<-unlist(data[1,])
    rownames(data)<-data[,1]
    data<-data[-1,-1]
    #names_data<-str_split(colnames(data),pattern=',',4,simplify = T)
    #exprosure_time<-which(names_data[,4]==time)
    
    #data_time<-data[,exprosure_time]
    
    
    data_time<-apply(data,2,as.numeric)
    
    
    rownames(data_time)<-rownames(data)
    
    
    
    data<-rbind(colnames(data_time),data_time)
    data<-as.matrix(cbind(rownames(data),data))
    kinome<-t(data)##transpose the raw data
    
    
    colnames(kinome)<-kinome[1,]
    kinome<-kinome[-1,]##delete the first column  
    result<-as.data.frame(kinome)
    
    result.Clean <- result %>% tidyr::separate(., 'V1', into = c('sample','replicate','time'), sep = ',')##separate the result table name as 3 different factors
    
    
    result.Clean %>% dplyr::group_by(sample, replicate) %>%
      reshape2::melt(.,id.vars=c('sample','replicate','time'))%>% dplyr::mutate(groupID = paste0('sample',':',sample,'\n','replicate',':', replicate,'\n')) ->a##1,change the table form to another format,selected the protein which we need.
    group<-as.matrix(paste(a$sample,a$replicate,sep=','))##make a new factor
    a<-cbind(a,group)
    
    
    
    
    R<-as.numeric(length(table(a$replicate))) ##set a number which is how many replicates we have used
    a$time<-as.numeric(a$time)  
    a$value<-as.numeric(a$value) 
    sample_name_loop<-c()
    protein<-c()
    repulication<-c()
    treatment<-c()
    AUC_value<-c()
    Pmax<-c()
    Pend<-c()
    Phalf<-c()
    slope_32_47<-c()
    slope_52_72<-c()
    slope_77_92<-c()
    for(p in 1:length(table(a$variable))){ 
      #p=8##test number
      pr<-rownames(as.matrix(table(a$variable)))[p]
      protein<-append(pr,protein)
      a_sub<-a[a$variable==pr,]
      for (i in 1:length(table(a_sub$group))) { 
        #i=2
        d=rownames(as.matrix(unlist(table(a_sub$group))))[i]
        sample<-str_split(d,pattern=',',2,simplify=T)[,1]
        rep<-str_split(d,pattern=',',2,simplify=T)[,2]
        sample_name_loop<-append(d,sample_name_loop)
        repulication<-append(rep,repulication)
        ##combine this sample and this replicate as a group in this loop
        group_func<-a_sub[which(a_sub$group==d),] ## then we will check all of data from the replicate of the sample you are working on
        #write.table(group_func,file='check.txt',sep = '\t')
        value_list<-c()
        slope_list<-c()
        for (l in (1:length(group_func$sample))) { 
          #l=1
          
          value_single<-as.numeric(group_func$value[l]) ## get all of the velocities of the replicate of this sample in different cycles
          value_list<-append(value_single,value_list)## save those velocities in a list
          slope_single<-group_func$value/group_func$time## get all of the slopes of the replicate of this sample in different cycles
          slope_list<-append(slope_single,slope_list)## save those slopes in a list
        }
        
        
        slope_32_47_1<-mean(slope_list[1:4]) ##get average slopes from cycle 32 to cycle 47 in the replicate of the samples
        slope_52_72_1<-mean(slope_list[5:9])##get average slopes from cycle 52 to cycle 72 in the replicate of the samples
        slope_77_92_1<-mean(slope_list[10:13])##get average slopes from cycle 77 to cycle 94 in the replicate of the samples
        slope_32_47<-append(slope_32_47_1,slope_32_47)
        slope_52_72<-append(slope_52_72_1,slope_52_72)
        slope_77_92<-append(slope_77_92_1,slope_77_92)
        
        
        
        value_list<-c(as.numeric(value_list))##set all of velocities from the replicate of the sample as a list
        Pmax<-append(max(value_list),Pmax)##find the largest values in the list
        Pend<-append(value_list[1],Pend)##find the last values in this list
        Phalf<-append(value_list[length(value_list)/2+1],Phalf)##find the middle values in this list
        AUC_single<-c()
        AUC<-0## the area under the curve  should start from 0
        
        for (k in 2:(length(value_list))) {
          
          AUC_single<-(value_list[k]+value_list[k-1])*330/2##every AUC from different parts of the whole curve are calculated
          AUC<-AUC+AUC_single##add the next part area to the whole part which we already set
        }
        AUC_value<-append(AUC,AUC_value) ##save all of AUCs in a list
      }
    }
    sample_all<-paste(a$sample,a$replicate,sep = '')
    N_S<-as.numeric(length(table(sample_all)))
    P_table<-cbind(rep(protein,each = N_S),sample_name_loop,repulication,Pend,Phalf,Pmax,AUC_value,slope_32_47,slope_52_72,slope_77_92)## combine all factors
    colnames(P_table)[1]<-'protein'
    write.csv(P_table,file= 'data_analysis.csv')
    
  }
  step3_data_analysis( data = data_unit_changed)
  
  #normalized_data<-read.csv('data_normalization_by_sum.csv',header = F)
  data_temp_normalization_BYART3<-read.csv('data_normalization_by_ART003.csv',header = F)
  #double_normalzation<-read.csv('data_normalization_by_sum_ART3.csv',header=F)
  data_temp_normalization_BYART4<-read.csv('data_normalization_by_ART004.csv',header = F)
  #data_orgin<-read.csv('data_orgin.csv',header = F)
  
  ##Create the Curve graph
  step4_curve_graph_script<-function(data){
    #data<-normalized_data
    #time<-50
    colnames(data)<-unlist(data[1,])
    rownames(data)<-data[,1]
    data<-data[-1,-1]
    #names_data<-str_split(colnames(data),pattern=',',4,simplify = T)
    #exprosure_time<-which(names_data[,4]==time)
    
    #data_time<-data[,exprosure_time]
    
    
    data_time<-apply(data,2,as.numeric)
    
    
    rownames(data_time)<-rownames(data)
    
    
    
    data<-rbind(colnames(data_time),data_time)
    data<-as.matrix(cbind(rownames(data),data))
    kinome<-t(data)##transpose the raw data
    
    
    colnames(kinome)<-kinome[1,]
    kinome<-kinome[-1,]##delete the first column  
    result<-as.data.frame(kinome)
    
    result.Clean <- result %>% tidyr::separate(., 'V1', into = c('sample','replicate','time'), sep = ',')##separate the result table name as 3 different factors
    
    
    result.Clean %>% dplyr::group_by(sample, replicate) %>%
      reshape2::melt(.,id.vars=c('sample','replicate','time'))%>% dplyr::mutate(groupID = paste0('sample',':',sample,'\n','replicate',':', replicate,'\n')) ->a##1,change the table form to another format,selected the protein which we need.
    group<-as.matrix(paste(a$sample,a$replicate,sep=','))##make a new factor
    a<-cbind(a,group)
    
    R<-as.numeric(length(table(a$replicate))) ##set a number which is how many replicates we have used
    a$time<-as.numeric(a$time)  
    a$value<-as.numeric(a$value) 
    
    write.csv(a,file='Long_dataForCurve.csv')
    
    group<-c()
    ##make all of the graphs of the curves, each graph will combine 9 proteins picture 
    for (grap in 1:(length(table(a$variable))/4)) {
      group<-(as.matrix(table(a$variable)))
      group<-rownames(group)[(((grap-1)*4)+1):(((grap-1)*4)+4)]##this and the following lines pick 9 proteins in every loop and produce the graphs
      a_1<-a[which(a$variable==group[1]),]
      a_2<-a[which(a$variable==group[2]),]
      a_3<-a[which(a$variable==group[3]),]
      a_4<-a[which(a$variable==group[4]),]
      #a_5<-a[which(a$variable==group[5]),]
      #a_6<-a[which(a$variable==group[6]),]
      #a_7<-a[which(a$variable==group[7]),]
      #a_8<-a[which(a$variable==group[8]),]
      #a_9<-a[which(a$variable==group[9]),]
      #a_10<-a[which(a$variable==group[10]),]
      #a_11<-a[which(a$variable==group[11]),]
      #a_12<-a[which(a$variable==group[12]),]
      
      ##next few lines help tailor the formed graphs to create the best lay-out
      p_1<-ggplot2::ggplot(a_1, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line() +
        facet_wrap(~variable, scales = 'free')+ ylim(0,max(0.1,max(a_1$value)))
      p_2<-ggplot2::ggplot(a_2, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line()  +
        facet_wrap(~variable, scales = 'free')+ ylim(0,max(0.1,max(a_2$value)))
      p_3<-ggplot2::ggplot(a_3, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line()  +
        facet_wrap(~variable, scales = 'free')+ ylim(0,max(0.1,max(a_3$value)))
      p_4<-ggplot2::ggplot(a_4, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line()  +
        facet_wrap(~variable, scales = 'free')+ ylim(0,max(0.1,max(a_4$value)))
      # p_5<-ggplot2::ggplot(a_5, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line()  +
      #facet_wrap(~variable, scales = 'free')
      # p_6<-ggplot2::ggplot(a_6, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line()  +
      #facet_wrap(~variable, scales = 'free')
      #p_7<-ggplot2::ggplot(a_7, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line()  +
      #facet_wrap(~variable, scales = 'free')
      #p_8<-ggplot2::ggplot(a_8, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line() +
      #facet_wrap(~variable, scales = 'free')
      #p_9<-ggplot2::ggplot(a_9, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line()  +
      #facet_wrap(~variable, scales = 'free')
      #p_10<-ggplot2::ggplot(a_10, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line()  +
      #facet_wrap(~variable, scales = 'free')
      #p_11<-ggplot2::ggplot(a_11, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line()  +
      #facet_wrap(~variable, scales = 'free')
      #p_12<-ggplot2::ggplot(a_12, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line()  +
      #facet_wrap(~variable, scales = 'free')
      pdf(file = paste(grap,'kinome_graph.pdf',sep = '_'),width=30,height = 20)##we will build a blank PDF file to save our graphs in this file
      grid.newpage()  ###we will build a blank area to draw our graphs in this area
      pushViewport(viewport(layout = grid.layout(2,2))) ####split this area as 3 x 3 lay-out
      vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}
      print(p_1, vp = vplayout(1,1))   ###draw the first graph in (1,1)lay-out
      print(p_2, vp = vplayout(1,2))     ###draw the 2nd graph in (1,2)lay-out
      print(p_3 , vp = vplayout(2,1))
      print(p_4 , vp = vplayout(2,2))
      #print(p_5 , vp = vplayout(2,1))
      #print(p_6 , vp = vplayout(2,2))
      #print(p_7 , vp = vplayout(2,3))
      #print(p_8 , vp = vplayout(2,4))
      #print(p_9 , vp = vplayout(3,1))
      #print(p_9 , vp = vplayout(3,2))
      #print(p_9 , vp = vplayout(3,3))
      #print(p_9 , vp = vplayout(3,4))
      
      dev.off()##draw all of the graphs in the pdf file which we build before
    }
  }
  step4_curve_graph_script(data_temp_normalization_BYART3)
  
  data_analysis<-read.csv(file = 'data_analysis.csv',header = T)
  
  ##Create the Heatmap
  step5_Heatmap_script<-function(data,time,cycle_time){
    #data<-data_analysis
    data<-data[,c(2,3,4,8)]
    data_protein<-str_split(data[,'protein'],pattern='_',3,simplify = T)
    ART_protein<-which(data_protein[,1]=='ART')
    data<-data[-c(ART_protein),]
    name<-paste(data[,2],'_',data[,3])
    data<-data[,-2]
    data<-data[,-2]
    data<-cbind(name,data)
    #colnames(data)<-c('protein','name','value')
    data_wide<-dcast(data=data,protein ~ name,mean)
    
    data_time<-apply(data_wide,2,as.numeric)
    
    rownames(data_time)<-data_wide[,1]
    
    data_time<-data_time[,-1]
    
    write.table(data_time,file='AUC_data.txt',sep='\t')
    
    data_LOG<-apply(data_time,2,function(x) log10(x)+10)
    data_LOG[apply(!is.finite(data_LOG), FUN = any, 1), ] = 0
    
    write.table(data_LOG,file='AUC_LOG_data.txt',sep='\t')
    
    data_time_n<-matrix(as.numeric(unlist(data_time)),nrow=nrow(data_time))
    data_t <- t(data_time_n)
    data_temp_sum<-apply(data_t, 2, sum)
    data_temp_standard<-rbind(data_temp_sum,data_t)
    data_temp_normalization<-apply(data_temp_standard,2,function(x) x/x[1])
    data_time_n<-t(data_temp_normalization[-1,])
    data_time_n<-apply(data_time_n,2,as.numeric)
    
    
    
    colnames(data_time_n)<-colnames(data_time)
    rownames(data_time_n)<-rownames(data_time)
    
    write.table(data_time_n,file='AUC_normalization_data.txt',sep='\t')
    
    p2<-pheatmap(data_time,
                 angle_col=45)
    pdf(file = paste('heat_map_graph.pdf'),width=20,height = 80)##we will build a blank PDF file to save our graphs in this file
    
    grid.newpage()  ###we will build a blank area to draw our graphs in this area
    pushViewport(viewport(layout = grid.layout(1,4))) ####split this area as 3 x 3 lay-out
    vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}
    print(p2, vp = vplayout(1,1:4)) 
    dev.off()
    
    p3<-pheatmap(data_LOG,
                 angle_col=45)
    pdf(file = paste('heat_map_graph_LOG10.pdf'),width=20,height = 80)##we will build a blank PDF file to save our graphs in this file
    
    grid.newpage()  ###we will build a blank area to draw our graphs in this area
    pushViewport(viewport(layout = grid.layout(1,4))) ####split this area as 3 x 3 lay-out
    vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}
    print(p3, vp = vplayout(1,1:4)) 
    dev.off()
    
    data_time_n[apply(!is.finite(data_time_n), FUN = any, 1), ] = 0
    
    #data_time_n<-data_time_n[-1,]
    p1<-pheatmap(data_time_n[-1,],
                 angle_col=45)
    pdf(file = paste('heat_map_graph_normalzation.pdf'),width=20,height = 80)##we will build a blank PDF file to save our graphs in this file
    
    grid.newpage()  ###we will build a blank area to draw our graphs in this area
    pushViewport(viewport(layout = grid.layout(1,4))) ####split this area as 3 x 3 lay-out
    vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}
    print(p1, vp = vplayout(1,1:4)) 
    dev.off()##draw all of the graphs in the pdf file which we build before
  }
  step5_Heatmap_script(data = data_analysis)
  
  
  ##choose a specific kinase to create a curve graph
  Specific_protein_graph<-function(data,protein){
    #data<-normalized_data
    #time<-50
    colnames(data)<-unlist(data[1,])
    rownames(data)<-data[,1]
    data<-data[-1,-1]
    #names_data<-str_split(colnames(data),pattern=',',4,simplify = T)
    #exprosure_time<-which(names_data[,4]==time)
    
    #data_time<-data[,exprosure_time]
    
    
    data_time<-apply(data,2,as.numeric)
    
    
    rownames(data_time)<-rownames(data)
    
    
    
    data<-rbind(colnames(data_time),data_time)
    data<-as.matrix(cbind(rownames(data),data))
    kinome<-t(data)##transpose the raw data
    
    
    colnames(kinome)<-kinome[1,]
    kinome<-kinome[-1,]##delete the first column  
    result<-as.data.frame(kinome)
    
    result.Clean <- result %>% tidyr::separate(., 'V1', into = c('sample','replicate','time'), sep = ',')##separate the result table name as 3 different factors
    
    
    result.Clean %>% dplyr::group_by(sample, replicate) %>%
      reshape2::melt(.,id.vars=c('sample','replicate','time'))%>% dplyr::mutate(groupID = paste0('sample',':',sample,'\n','replicate',':', replicate,'\n')) ->a##1,change the table form to another format,selected the protein which we need.
    group<-as.matrix(paste(a$sample,a$replicate,sep=','))##make a new factor
    a<-cbind(a,group)##make a new factor
    
    a$value<-as.numeric(a$value)
    a$time<-as.numeric(a$time)
    #protein<-'41_654_666'
    a_1<-a[which(a$variable==protein),]
    p_1<-ggplot2::ggplot(a_1, aes(x = time, color = groupID, group = groupID, y = value)) + geom_line() +
      facet_wrap(~variable, scales = 'free')
    pdf(file = paste(protein,'_graph.pdf'),width=20,height = 20)
    print(p_1) 
    dev.off()##draw all of the graphs in the pdf file which we build before
  }
  Specific_protein_graph(data = data_temp_normalization_BYART4 ,protein <- c('AKT1_320_332'))
  
  
  ##back to the root manu
  setwd('../')
}
start(exposure_time='10s')





data_unit_changed<-read.csv(file = '200s_results/unit_changed_data.csv',sep=',',header = F)


split_all_data_by_protein<-function(data){
  data<-data_unit_changed
  
  colnames(data)<-unlist(data[1,])
  rownames(data)<-data[,1]
  data<-data[-1,-1]
  
  r_data<-as.numeric(nrow(data))
  dir.create('each_protein')
  
  protein_names<-rownames(data)[3:r_data]
  protein_names<-str_split(protein_names,pattern='\\/',2,simplify=T)[,1]
  
  
  for (i in 1 : length(protein_names)) {
    #i<-1
    single_protein<-protein_names[i]
    data_protein<-data[single_protein,]
    sample_names<-str_split(colnames(data_protein), pattern = ',' , 4 , simplify = T) 
    sample_N<-rownames(as.matrix(table(sample_names[,1])))
    protein_table<-c()
    for (j in 1:length(sample_N)) {
      #j<-1
      sample<-which(sample_names[,1]== sample_N[j])
      data_protein_sample<-data_protein[,sample]
      data_protein_sample<-t(data_protein_sample)
      colnames(data_protein_sample)<-sample_N[j]
      rownames(data_protein_sample)<-rownames(as.matrix(table(sample_names[,3])))
      protein_table<-cbind(protein_table,data_protein_sample)
        }
    write.csv(protein_table,file = paste('each_protein/single_protein_',single_protein,'.csv'))
    
     
  }
  
} ##This function could be use to normalized data and unormalized data                      
split_all_data_by_protein(data = data_unit_changed)


