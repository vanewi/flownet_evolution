source("paper/stat_util.R")

#PLOTTING FUNCTIONS -----
#RESULTS_CONTROL CONSTRUCTION:
base_log=2
input_output_entropy=input_output_entropy_generator(log_base=base_log);
entropy_diff=entropy_diff_generator(log_base=base_log);

node_ensembles=as.integer(c(30,50,70,100));
orientors=list(TST=tst, ASC=ascendency, AMI=ami,Finn=finn,EDiff=entropy_diff,b=b)
num_checkpoints=1
main_folder=getwd()

for (N in 1:length(node_ensembles)){
  node_folder=paste0(main_folder,'/',node_ensembles[N])
  for (or in 1:length(orientors)){
    orientor_folder=paste0(node_folder,'/',names(orientors[or]))
    setwd(orientor_folder)
    
    results_control=vector(mode='list',length=num_checkpoints);
    for (i in 1:num_checkpoints){
      to_save=readRDS(paste0(orientor_folder,'/','cp_control','/',i))
      results_control[[i]]=to_save
      print(paste0("cp: ",i))
    }
    
    saveRDS(results_control,file='results_control')
    rm(results_control)
    gc()
    print(paste0('orientor:',or))
  }
  print(paste('N: ',N))
}

#PLOTTING FUNCTIONS OF ORIENTORS--
#FOR A FIXED NODE ENSEMBLE
#EXAMPLE WITH N=3 (50 NODES) 100 REPLICATES
#creates corrplot

source("FlowNet.R")
main_folder=getwd()
node_ensembles=as.integer(c(30,50,70,100));
input_output_entropy=input_output_entropy_generator(log_base=2);
entropy_diff=entropy_diff_generator(log_base=2);
N=2
node_folder=paste0(main_folder,'/',node_ensembles[N])
orientors=list(TST=tst,ASC=ascendency,AMI=ami,Finn=finn,EDiff=entropy_diff,b=b);
plotting=c('tst','ascendency','ami','finn','entropy_diff','b');
ylims_max=c(NULL,NULL,NULL,NULL,NULL,NULL)
ylims_min=c(NULL,NULL,NULL,NULL,NULL,NULL)
transforming=list((function(x) return(log(x))),
                  (function(x) return(log(x))),
                  (function(x) return(x)),
                  (function(x) return(x)),
                  (function(x) return(x)),
                  (function(x) return(x))
)

par(mfrow=c(length(orientors),length(orientors)),cex.axis=0.8)#, mai=c(0.1,0.1,0.2,0.1))
for (i0 in 1:length(orientors)){
  print(paste0('orientors_x :',names(orientors[i0])))
  orientor_folder=paste0(node_folder,'/',names(orientors[i0]))
  setwd(orientor_folder)
  
  results_control=readRDS('results_control');
  n_nets=length(results_control[[1]])
  variables=names(results_control[[1]][[1]]$state$extra$history[[1]])
  
  all_control=vector(mode='list',length=n_nets)
  for (cp in 1:length(results_control)){
    for (net in 1:n_nets){
      all_control[[net]]=c(all_control[[net]],results_control[[cp]][[net]]$state$extra$history)
    }
  }
  
  rm(results_control)
  gc()
  
  vis=Visualization$new(results_list=all_control,variables=variables)
  rm(all_control)
  gc()
  
  for(j in 1:length(orientors)){
    if (i0==j){
      par(mar=c(1.9, 1.9, 1, 1))
      plot(x = 0:5, y = 0:5, ann = F,bty = "n",type = "n",
           xaxt = "n", yaxt = "n")
      text(x=2,y=2.5,names(orientors[j]))
    }
    else{
      par(mar=c(1.9, 1.9, 1, 1))
      vis$plot_me2(x_var = plotting[i0],
                   y_var = plotting[j],
                   col_var=NULL,
                   x_func = transforming[[i0]],
                   y_func = transforming[[j]]
                   #y_lim=list(min=ylims_min[j],max=ylims_max[j]), 
      )
      
    }
    print(paste0('orientors_y :',names(orientors[j])))
  }
}


#------------------------------------------------------------------
#PLOTTING FUNCTIONS WITH COLORS-----
source("FlowNet.R")
main_folder=getwd()

setwd(main_folder)
base_log=2
input_output_entropy=input_output_entropy_generator(log_base=base_log);
entropy_diff=entropy_diff_generator(log_base=base_log);

node_ensembles=as.integer(c(30,50,70,100));
orientors=list(TST=tst,ASC=ascendency,AMI=ami,EDiff=entropy_diff,Finn=finn,b=b);
plotting=c('tst','ascendency','ami','ediff_norm','finn','b');
#ylims_max=c(45,45,3.5,4,35,0)
#ylims_min=c(5,5,0.5,-4,0,-0.010)
transforming=list((function(x) return(log(x))),
                  (function(x) return(log(x))),
                  (function(x) return(x)),
                  (function(x) return(x)),
                  (function(x) return(-log(1.00000000000001-x))),
                  (function(x) return(x))
)

N=2

#node_folder=paste0(main_folder,'/',node_ensembles[N])

par(mfrow=c(length(orientors),length(orientors)),cex.axis=0.8)#, mai=c(0.1,0.1,0.2,0.1))
for (i0 in 1:length(orientors)){
  print(paste0('orientors_x :',names(orientors[i0])))
  #orientor_folder=paste0(node_folder,'/',names(orientors[i0]))
  #setwd(orientor_folder)
  
  #ens_evolve_control=readRDS('ens_evolve_control')
  #results_control=ens_evolve_control$get_results()
  #results_control_list=lapply(results_control,function(x) x$state$extra$history)
  
  #variables=names(results_control_list[[1]][[1]])
  #n_nets=length(results_control_list)
  
  #rm(results_control,ens_evolve_control)
  #gc()
  
  #vis=Visualization$new(results_list=results_control_list,variables=variables)
  
  con=Manipulation1Cp$new(ensemble_size=node_ensembles[N],orientor=names(orientors[i0]))
  control=con$open_control()
  vis=Visualization$new(results_list=con$results_control_list,variables=con$variables)
  
  gc()
  
  for(j in 1:length(orientors)){
    if (i0==j){
      par(mar=c(1.9, 1.9, 1, 1))
      plot(x = 0:5, y = 0:5, ann = F,bty = "n",type = "n",
           xaxt = "n", yaxt = "n")
      text(x=2,y=2.5,names(orientors[j]))
    }
    else{
      par(mar=c(1.9, 1.9, 1, 1))
      vis$plot_me2(x_var = plotting[i0],
                   y_var = plotting[j],
                   with_raster = FALSE,
                   col_var='time',
                   x_func = transforming[[i0]],
                   y_func = transforming[[j]],
                 #  y_lim=list(min=ylims_min[j],max=ylims_max[j]) 
      )
      
    }
    print(paste0('orientors_y :',names(orientors[j])))
  }
}

#-------------------------------------------------------------
#--------------------
#CORRELATION PLOT -----
source("FlowNet.R")
library(corrplot)
main_folder=getwd()

setwd(main_folder)
node_ensembles=as.integer(c(30,50,70,100));
base_log=2
input_output_entropy=input_output_entropy_generator(log_base=base_log);
entropy_diff=entropy_diff_generator(log_base=base_log);

orientors=list(TST=tst,ASC=ascendency, AMI=ami,EDiff=entropy_diff,Finn=finn,b=b)
plotting=c('tst','ascendency','ami','ediff_norm','finn','eigMax');
conversions=c(function(x) return(log(x)),function(x) return(log(x)),function(x) return(x),function(x) return(x),function(x) return(-log(1.00000000000001-x)),function(x) return(x))
figure=LETTERS[seq(from=1,to=length(node_ensembles))]
num_checkpoints=1

par(mfrow=c(2,2))

for (N in 1:length(node_ensembles)){
  node_folder=paste0(main_folder,'/',node_ensembles[N])
  correlations=matrix(nrow=length(orientors),ncol=length(orientors),dimnames = list(names(orientors),names(orientors)))
  deviations=matrix(nrow=length(orientors),ncol=length(orientors),dimnames = list(names(orientors),names(orientors)))
  for (i0 in 1:length(orientors)){
    print(paste0('orientors_x :',names(orientors[i0])))
    orientor_folder=paste0(node_folder,'/',names(orientors[i0]))
    setwd(orientor_folder)
    
    ens_evolve_control=readRDS('ens_evolve_control')
    results_control=ens_evolve_control$get_results()
    results_control_list=lapply(results_control,function(x) x$state$extra$history)
    
    
    variables=names(results_control_list[[1]][[1]])
    n_nets=length(results_control_list)
    
    rm(results_control,ens_evolve_control)
    gc()
    
    vis=Visualization$new(results_list=results_control_list,variables=variables)
    
    for(j in 1:length(orientors)){
      corr_net=vector(mode='numeric',length=n_nets);
      for(n in 1:n_nets){
        val=cbind(vis$data[[n]][plotting[i0]],vis$data[[n]][plotting[j]]);
        corr_net[n]=cor.test(conversions[[i0]](val[,1]),conversions[[j]](val[,2]))$estimate
      }
      deviations[i0,j]=sd(corr_net)
      #if(deviations[i0,j]<0.25){correlations[i0,j]=summary(corr_net)['Median']}
      #if(deviations[i0,j]>0.25){correlations[i0,j]=0}
      correlations[i0,j]=summary(corr_net)['Median']
      print(paste0('orientors_y :',names(orientors[j])))
    }
  }
  corrplot(correlations,diag=FALSE,title=paste0(figure[N],". ",'Control ensemble with ',node_ensembles[N],' core nodes'),cex.main=0.85, mar=c(0,0,1.5,0))
  rm(vis)
  gc()
}


compare=c()
for(n in 1:n_nets){
  compare=rbind(compare,cbind(vis$data[[n]][plotting[2]],vis$data[[n]][plotting[5]]))
}

#-------------------------------------------------------------
#DISTRIBUTION PLOT OF CORRELATIONS OF ORIENTORS-----------
source("FlowNet.R")
library(corrplot)
main_folder=getwd()

setwd(main_folder)
node_ensembles=as.integer(c(30,50,70,100));
base_log=2
input_output_entropy=input_output_entropy_generator(log_base=base_log);
entropy_diff=entropy_diff_generator(log_base=base_log);

orientors=list(TST=tst,ASC=ascendency, AMI=ami,EDiff=entropy_diff,Finn=finn,b=b)
plotting=c('tst','ascendency','ami','ediff_norm','finn','b');
conversions=c(function(x) return(log(x)),function(x) return(log(x)),function(x) return(x),function(x) return(x),function(x) return(-log(1.00000000000001-x)),function(x) return(x))
num_checkpoints=1

N=3

par(mfrow=c(length(orientors),length(orientors)),cex.axis=0.8)#, mai=c(0.1,0.1,0.2,0.1))
  node_folder=paste0(main_folder,'/',node_ensembles[N])
  correlations=matrix(nrow=length(orientors),ncol=length(orientors),dimnames = list(names(orientors),names(orientors)))
  deviations=matrix(nrow=length(orientors),ncol=length(orientors),dimnames = list(names(orientors),names(orientors)))
  for (i0 in 1:length(orientors)){
    print(paste0('orientors_x :',names(orientors[i0])))
    orientor_folder=paste0(node_folder,'/',names(orientors[i0]))
    setwd(orientor_folder)
    
    ens_evolve_control=readRDS('ens_evolve_control')
    results_control=ens_evolve_control$get_results()
    results_control_list=lapply(results_control,function(x) x$state$extra$history)
    
    variables=names(results_control_list[[1]][[1]])
    n_nets=length(results_control_list)
    
    rm(results_control,ens_evolve_control)
    gc()
    
    vis=Visualization$new(results_list=results_control_list,variables=variables)
    
    #net_statistics=vector(mode='list',length=orientors)
    for(j in 1:length(orientors)){
      if (i0==j){
        par(mar=c(1.9, 1.9, 1, 1))
        plot(x = 0:5, y = 0:5, ann = F,bty = "n",type = "n",
             xaxt = "n", yaxt = "n")
        text(x=2,y=2.5,names(orientors[j]))
      }
      else{
      corr_net=vector(mode='numeric',length=n_nets);
      for(n in 1:n_nets){
        val=cbind(vis$data[[n]][plotting[i0]],vis$data[[n]][plotting[j]]);
        corr_net[n]=cor.test(conversions[[i0]](val[,1]),conversions[[j]](val[,2]))$estimate
      }
      correlations[i0,j]=summary(corr_net)['Mean']
      deviations[i0,j]=sd(corr_net)
      print(paste0('orientors_y :',names(orientors[j])))
      h=hist(corr_net,breaks=seq(-1,1,0.1),plot=FALSE)
      h$density = h$counts/sum(h$counts)
      col_val = (correlations[i0,j]+1.0)/2.0 
      cols = rgb(1.0-col_val,0.0,col_val,1.0)
      plot(h,freq=FALSE,ylim = c(0.0,1.0),col=cols,border = rgb(1.0,1.0,1.0,alpha=1.0),main=NULL)
      lines(x=c(0.0,0.0),y=c(0.0,1.0),lty=2,col=rgb(0.0,0.0,0.0,0.3))
      lines(x=c(-1.0,1.0),y=c(0.0,0.0))
      lines(x=c(correlations[i0,j],correlations[i0,j]),y=c(0.0,1.0),lty=1,col='green')
      rect(max(-1.0,correlations[i0,j]-deviations[i0,j]),1,min(1.0,correlations[i0,j]+deviations[i0,j]),0, col= rgb(0.5,0.5,0.5,alpha=0.2),border = NA)      }
    }
  }







#------------------------------------------------------------
#DISTRIBUTION PLOT OF CORRELATIONS OF ORIENTORS ALL ENSEMBLES----
  source("FlowNet.R")
  library(corrplot)
  main_folder=getwd()
  
  node_ensembles=as.integer(c(30,50,70,100));
  base_log=2
  input_output_entropy=input_output_entropy_generator(log_base=base_log);
  entropy_diff=entropy_diff_generator(log_base=base_log);
  
  orientors=list(TST=tst,ASC=ascendency, AMI=ami,EDiff=entropy_diff,Finn=finn,b=b)
  plotting=c('tst','ascendency','ami','ediff_norm','finn','eigMax');
  conversions=c(function(x) return(log(x)),function(x) return(log(x)),function(x) return(x),function(x) return(x),function(x) return(-log(1.00000000000001-x)),function(x) return(x))
  num_checkpoints=1
  
  combine_ensembles=function(orientor){
    all_data=vector(mode="list",length=length(node_ensembles))
    for (N in 1:length(node_ensembles)){
      setwd(main_folder)
      node_folder=paste0(main_folder,'/',node_ensembles[N])
      orientor_folder=paste0(node_folder,'/',orientor)
      print(orientor_folder)
      setwd(orientor_folder)
      
      #TO CHECK: WETHER CONTROL OR NATURAL HISTORIES ARE PUT IN HISTORIES!
      ens_evolve_control=readRDS('ens_evolve_control')
      results_control=ens_evolve_control$get_results()
      results_control_list=lapply(results_control,function(x) x$state$extra$history)
      variables=names(results_control_list[[1]][[1]])
      vis=Visualization$new(results_list=results_control_list,variables=variables)
      all_data[[N]]=as.matrix(vis$data);
      gc()
    }
    out=Reduce(rbind,all_data)
    setwd(main_folder)
    return(out)
  }
  
  #all_ensembles=lapply(names(orientors),combine_ensembles)
  names(all_ensembles)=names(orientors)
  #saveRDS(all_ensembles,file='all_control_ensembles.RDS')
  #all_ensembles=readRDS('all_control_ensembles.RDS')
  n_nets=length(all_ensembles[[1]])
  
  par(mfrow=c(length(orientors),length(orientors)),cex.axis=0.8)#, mai=c(0.1,0.1,0.2,0.1))
  correlations=matrix(nrow=length(orientors),ncol=length(orientors),dimnames = list(names(orientors),names(orientors)))
  deviations=matrix(nrow=length(orientors),ncol=length(orientors),dimnames = list(names(orientors),names(orientors)))
  for(i0 in 1:length(orientors)){ 
    print(paste0('orientors_x :',names(orientors[i0])))
    for(j in 1:length(orientors)){
      if (i0==j){
        par(mar=c(1.9, 1.9, 1, 1))
        plot(x = 0:5, y = 0:5, ann = F,bty = "n",type = "n",
             xaxt = "n", yaxt = "n")
        text(x=2,y=2.5,names(orientors[j]))
      }
      else{
        corr_net=vector(mode='numeric',length=length(n_nets));
        for(n in 1:n_nets){
          val=cbind(all_ensembles[[names(orientors[i0])]][[n]][plotting[i0]],all_ensembles[[names(orientors[i0])]][[n]][plotting[j]]);
          corr_net[n]=cor.test(conversions[[i0]](val[,1]),conversions[[j]](val[,2]))$estimate
        }
        correlations[i0,j]=summary(corr_net)['Mean']
        deviations[i0,j]=sd(corr_net)
        print(paste0('orientors_y :',names(orientors[j])))
        #if(correlations[i0,j]>0){cols='blue'}
        #if(correlations[i0,j]<0){cols='red'}
        col_val = (correlations[i0,j]+1.0)/2.0 
        h=hist(corr_net,breaks=seq(-1,1,0.1),plot=FALSE)
        h$density = h$counts/sum(h$counts)
        cols = rgb(1.0-col_val,0.0,col_val,1.0)
        plot(h,freq=FALSE,ylim = c(0.0,1.0),col=cols,border = rgb(1.0,1.0,1.0,alpha=1.0),main=NULL)
        lines(x=c(0.0,0.0),y=c(0.0,1.0),lty=2,col=rgb(0.0,0.0,0.0,0.3))
        lines(x=c(-1.0,1.0),y=c(0.0,0.0))
        lines(x=c(correlations[i0,j],correlations[i0,j]),y=c(0.0,1.0),lty=1,col='green')
        rect(max(-1.0,correlations[i0,j]-deviations[i0,j]),1,min(1.0,correlations[i0,j]+deviations[i0,j]),0, col= rgb(0.5,0.5,0.5,alpha=0.2),border = NA)
      }
    }
  }
  
  
  
  
  
  
  
  