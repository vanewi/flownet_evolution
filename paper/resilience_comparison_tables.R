source("paper/stat_util.R")

library(gt)
main_folder=getwd()

setwd(main_folder)
base_log=2
input_output_entropy=input_output_entropy_generator(log_base=base_log);
entropy_diff=entropy_diff_generator(log_base=base_log);

node_ensemble=c(30,50,70,100)
orientors=list(TST=tst, AMI=ami,EDiff=entropy_diff)
plotting=c('tst','ami','entropy_diff');
greed_perc=c(0.1,0.3,0.5)
greed_orientor=c('part_stock','proportion_stock','part_finn','proportion_finn')

num_checkpoints=1

main_folder=getwd()
#options(digits = 3)

#----to create resilience table and y_s:----

y_sN=vector(mode='list',length = length(node_ensemble))
for (N in 1:length(node_ensemble)) {
  print(paste0('ensemble: ',node_ensemble[N]))
  pos=1;
  y_s=vector(mode='list',length=length(orientors)*length(greed_orientor)*length(greed_perc))
  to_name=vector(mode='character',length=length(orientors)*length(greed_orientor)*length(greed_perc))
  for (or in 1:length(orientors)) {
    print(paste0('orientor: ',names(orientors[or])))
    
    con=Manipulation1Cp$new(ensemble_size=node_ensemble[N],orientor=names(orientors[or]))
    control=con$open_control()
    cp_list_control=control$results_control_list
    n_nets=length(cp_list_control)
    variables=names(cp_list_control[[1]][[1]])
    
    cp_list_control_scalars=Visualization$new(results_list=cp_list_control,variables=variables)$data
    
    rm(cp_list_control,con,control)
    gc()
    
    y_control=vector(mode='numeric',length=n_nets)
    for (i in 1:n_nets){
      y_control[i]=cp_list_control_scalars[[i]]$b[length(cp_list_control_scalars[[i]]$b)]
    }
    
    inf_error=which(y_control==0)
    if(!length(inf_error)==0) {y_control[inf_error]=(-0.0000000000000001)}
    
    table_or=vector(mode='list',length=length(greed_orientor))
    
    for (greed_or in 1:length(greed_orientor)){
      print(greed_orientor[greed_or])
      table_perc_or=vector(mode='list',length=length(greed_perc))
      for (greed_p in 1:length(greed_perc)){
        print(greed_perc[greed_p])
        
        #greedy_path='cp_greedy_0.3_part_stock'
        greedy_path=paste0('cp_greedy_',greed_perc[greed_p],'_',greed_orientor[greed_or])
        greedy_path
        
        gred=Manipulation1Cp$new(ensemble_size=node_ensemble[N],orientor=names(orientors[or]),greed_orientor=greed_orientor[greed_or],greed_perc=greed_perc[greed_p])
        greedy=gred$open_greedy()
        #cp_list_greedy=greedy$results_greedy_list
        
        cp_list_greedy_scalars=Visualization$new(results_list=greedy$results_greedy_list,variables=variables)$data
        
        y_greedy=vector(mode='numeric',length=n_nets)
        for (i in 1:n_nets){
          y_greedy[i]=cp_list_greedy_scalars[[i]]$b[length(cp_list_greedy_scalars[[i]]$b)]
        }
        
        inf_error2=which(y_greedy==0)
        if(!length(inf_error2)==0) {y_greedy[inf_error2]=-0.0000000000000001}
        
        y=(y_control-y_greedy)/abs(y_control+y_greedy)
        
        print('comparison done')
        
        #resilience table:    
        better=which(y<0)
        GpercY=length(better)/length(y)*100
        GpercY
        if(length(better)>1){
          YG=summary(-1*y[better])
          YC=summary(y[-better])
          
          Yperc=as.data.frame(t(c('',paste0(round(GpercY,2)),'')))
          names(Yperc)=c()
          YrangeG=as.data.frame(t(c(paste0(round(YG[['Min.']],2)),paste0(round(YG[['Median']],2)),paste0(round(YG[['Max.']],2)))))
          names(YrangeG)=c()
          YrangeC=as.data.frame(t(c(paste0(round(YC[['Min.']],2)),paste0(round(YC[['Median']],2)),paste0(round(YC[['Max.']],2)))))
          names(YrangeC)=c()
        }
        if(length(better)==0){
          Yperc=as.data.frame(t(c('','0','')))
          names(Yperc)=c()
          YrangeG=as.data.frame(t(c('0','0','0')))
          names(YrangeG)=c()
          YC=summary(y)
          YrangeC=as.data.frame(t(c(paste0(round(YC[['Min.']],2)),paste0(round(YC[['Median']],2)),paste0(round(YC[['Max.']],2)))))
          names(YrangeC)=c()
        }
        if(length(better)==1){
          YG=summary(-1*y[better])
          YC=summary(y[-better])
          Yperc=as.data.frame(t(c('',paste0(round(GpercY,2)),'')))
          names(Yperc)=c()
          YrangeG=as.data.frame(t(c('',paste0(round(YG[['Min.']],2)),'')))
          names(YrangeG)=c()
          YrangeC=as.data.frame(t(c(paste0(round(YC[['Min.']],2)),paste0(round(YC[['Median']],2)),paste0(round(YC[['Max.']],2)))))
          names(YrangeC)=c()
        }
        
        #experiment=paste0(node_ensemble[N],'_',names(orientors[or]),'_',greedy_path)
        table_perc_or[[greed_p]]=matrix(cbind(Yperc,YrangeG,YrangeC),ncol=3,byrow = TRUE)
        names(table_perc_or[[greed_p]])=c()
        table_perc_or
        print('table_perc_advances')
        y_s[[pos]]=y;
        to_name[pos]=paste0(names(orientors[or]),'_',node_ensemble[N],'_',greed_perc[greed_p],'_',greed_orientor[greed_or])
        pos=pos+1;
      }
      table_or[[greed_or]]=Reduce(rbind,table_perc_or)
      print('table_or_advances')
    }
    final_table=as.data.frame(Reduce(cbind,table_or))
    tag=rep(c('greedy %','greedy range', 'control range'),3)
    final_table=cbind(tag,final_table)
    
    print_table<-final_table%>%gt(rowname_col = 'tag')
    print_table<-print_table %>%  tab_stubhead(label = paste0(names(orientors[or]))
                                               #                        ) %>% tab_row_group(label=c('10%','30%',"50%"),rows=list(c(1:3),c(4:6),c(7,9)))                   
    ) %>% tab_row_group(label = "10%", rows = 1:3
    ) %>% tab_row_group(label = "30%", rows = 4:6
    ) %>% tab_row_group(label = "50%", rows = 7:9
    ) %>% row_group_order(groups = c("10%", "30%","50%")
    ) %>%  tab_spanner(label = 'absolute stock',columns = c(V1,V2,V3)
    ) %>% tab_spanner(label = 'proportion stock',columns = c(V4,V5,V6)
    ) %>% tab_spanner(label = 'absolute finn',columns = c(V7,V8,V9)
    ) %>% tab_spanner(label = 'proportion finn',columns = c(V10,V11,V12),id='prop_finn'
    ) %>% tab_options(row_group.background.color = "grey",column_labels.border.bottom.color ='white',column_labels.border.top.color ='black',table.border.bottom.width = 3,table.border.bottom.color = 'black'
    ) %>% tab_style(locations = cells_body(columns = c(V4,V5,V6)), style = list(cell_fill(color = "lightgrey"))
    ) %>% tab_style(locations = cells_body(columns = c(V10,V11,V12)),style = list(cell_fill(color = "lightgrey"))
    ) %>% tab_style(locations=cells_row_groups(groups = everything()),style=list(cell_text(weight = "bold", size = 24))
    ) %>% tab_style(locations=cells_stubhead(), style=list(cell_text(weight = "bold", size = 24),cell_borders(sides = "bottom", weight = px(3)))
    ) %>% tab_style(locations=cells_column_labels(), style=list(cell_borders(sides = "bottom", weight = px(3)))
    ) %>% tab_style(locations=cells_column_spanners(),style=list(cell_text(weight = "bold", size = 24),cell_borders(sides = "bottom", weight = 0))
    ) %>% cols_label(V1='',V2='',V3='',V4='',V5='',V6='',V7='',V8='',V9='',V10='',V11='',V12='')
    
    #print_table    
    node_folder=paste0(main_folder,'/',node_ensemble[N])
    orientor_folder=paste0(node_folder,'/',names(orientors[or]))
    setwd(orientor_folder)
    saveRDS(list(pretty_table=print_table,data_table=final_table),file = paste0('resilience_tables'))
    print('table_saved')
    setwd(main_folder)
  }
  names(y_s)=to_name;
  y_sN[[N]]=y_s;
}
setwd(main_folder)
saveRDS(y_sN,'y_s')
setwd(main_folder)

#options(digits=7)

xs=readRDS('x_s')


#----to open all tables:----
source("FlowNet.R")
library(gt)
main_folder=getwd()

setwd(main_folder)
node_ensemble=c(30,50,70,100)
orientors=list(EDiff=entropy_diff)#TST=tst,EDiff=entropy_diff,AMI=ami,EDiff=entropy_diff, ASC=ascendency, Finn=finn,b=b)

which_tables='performance_tables'
for (N in 1:length(node_ensemble)) {
  for (or in 1:length(orientors)) {
    node_folder=paste0(main_folder,'/',node_ensemble[N])
    orientor_folder=paste0(node_folder,'/',names(orientors[or]))
    orientor_folder
    setwd(orientor_folder)
    tables=readRDS(which_tables)
    print_table=tables$pretty_table
    print(print_table %>% tab_header(title=paste0(node_ensemble[N],'core nodes')))
    setwd(main_folder)
  }
}



#---------------------
#raster plots -----------------
source("FlowNet.R")
library(gt)
main_folder=getwd()

setwd(main_folder)
input_output_entropy=input_output_entropy_generator(log_base=2);
entropy_diff=entropy_diff_generator(log_base=2);

node_ensemble=c(30,50,70,100)
greed_orientor=c('part_stock','proportion_stock','part_finn','proportion_finn')
orientors=list(TST=tst,AMI=ami,EDiff=entropy_diff)#, ASC=ascendency, AMI=ami,Finn=finn,b=b,TST=tst)

create_raster_per_orientor=function(which_orientor, which_table){
  part_stock_table=matrix(ncol=length(node_ensemble),nrow=3,dimnames = list(c('10%','30%','50%'),node_ensemble))
  proportion_stock_table=matrix(ncol=length(node_ensemble),nrow=3,dimnames = list(c('10%','30%','50%'),node_ensemble))
  part_finn_table=matrix(ncol=length(node_ensemble),nrow=3,dimnames = list(c('10%','30%','50%'),node_ensemble))
  proportion_finn_table=matrix(ncol=length(node_ensemble),nrow=3,dimnames = list(c('10%','30%','50%'),node_ensemble))
  for (N in 1:length(node_ensemble)) {
    setwd(main_folder)
    node_folder=paste0(main_folder,'/',node_ensemble[N])
    orientor_folder=paste0(node_folder,'/',which_orientor)
    setwd(orientor_folder)
    tables=readRDS(which_table)
    #data table for all experiments of an orientor of one ensemble
    
    tab_all=tables$data_table
    part_stock_table[,N]=as.numeric(unlist(tab_all[c(1,4,7),3]))/100;
    proportion_stock_table[,N]=as.numeric(unlist(tab_all[c(1,4,7),6]))/100;
    part_finn_table[,N]=as.numeric(unlist(tab_all[c(1,4,7),9]))/100;
    proportion_finn_table[,N]=as.numeric(unlist(tab_all[c(1,4,7),12]))/100;
  }
  setwd(main_folder)
  #return(list(part_stock=as.data.frame(part_stock_table),proportion_stock=as.data.frame(proportion_stock_table),part_finn=as.data.frame(part_finn_table),proportion_finn=as.data.frame(proportion_finn_table)))
  return(list(part_stock=part_stock_table,proportion_stock=proportion_stock_table,part_finn=part_finn_table,proportion_finn=proportion_finn_table))
}

which_table='correct_performance_tables'
setwd(main_folder)
TST_rasters=create_raster_per_orientor(names(orientors[1]),which_table = which_table)
AMI_rasters=create_raster_per_orientor(names(orientors[2]),which_table = which_table)
Ediff_rasters=create_raster_per_orientor(names(orientors[3]),which_table = which_table)


par(mfrow=c(2,2))

library(corrplot)
corrplot(TST_rasters$part_stock,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='A. absolute stock',mar=c(1,0,1,1),cex.main=1)
corrplot(TST_rasters$proportion_stock,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='B. propotion stock',mar=c(1,1,1,0),cex.main=1)
corrplot(TST_rasters$part_finn,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='C. absolute finn',mar=c(0,0,1,1),cex.main=1)
corrplot(TST_rasters$proportion_finn,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='D. proportion finn',mar=c(0,1,1,0),cex.main=1)
mtext("TST", side = 3, line = -2, outer = TRUE)

par(mfrow=c(2,2))

corrplot(AMI_rasters$part_stock,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='A. absolute stock',mar=c(1,0,1,1),cex.main=1)
corrplot(AMI_rasters$proportion_stock,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='B. propotion stock',mar=c(1,1,1,0),cex.main=1)
corrplot(AMI_rasters$part_finn,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='C. absolute finn',mar=c(0,0,1,1),cex.main=1)
corrplot(AMI_rasters$proportion_finn,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='D. proportion finn',mar=c(0,1,1,0),cex.main=1)
mtext("AMI", side = 3, line = -2, outer = TRUE)

par(mfrow=c(2,2))

corrplot(Ediff_rasters$part_stock,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='A. absolute stock',mar=c(1,0,1,1),cex.main=1)
corrplot(Ediff_rasters$proportion_stock,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='B. propotion stock',mar=c(1,1,1,0),cex.main=1)
corrplot(Ediff_rasters$part_finn,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='C. absolute finn',mar=c(0,0,1,1),cex.main=1)
corrplot(Ediff_rasters$proportion_finn,method = 'color',addCoef.col="black",number.cex = 1,col.lim = c(0,1),title='D. proportion finn',mar=c(0,1,1,0),cex.main=1)
mtext("EDiff", side = 3, line = -2, outer = TRUE)








