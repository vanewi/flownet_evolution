source("FlowNet.R")

Manipulation1Cp <- setRefClass("Manipulation1Cp",
                               fields = list(ensemble_size ="numeric",
                                             orientor = "character",
                                             greed_orientor='character',
                                             greed_perc='numeric',
                                             results_control='list',
                                             results_control_list='list',
                                             results_greedy='list',
                                             results_greedy_list='list',
                                             variables='character',
                                             # vis='Visualization',
                                             # vis_greedy='Visualization',
                                             main_folder='character'
                               ),
                               methods = list(
                                 initialize = function(...,ensemble_size,orientor,greed_orientor='',greed_perc=0){
                                   main_folder<<-getwd();
                                   node_folder=paste0(main_folder,'/',ensemble_size);
                                   orientor_folder=paste0(node_folder,'/',orientor);
                                   print(orientor_folder)
                                   ensemble_size<<-ensemble_size;
                                   orientor<<-orientor;
                                   greed_orientor<<-greed_orientor;
                                   greed_perc<<-greed_perc; 
                                 },
                                 open_control=function(){
                                   setwd(main_folder)
                                   node_folder=paste0(main_folder,'/',ensemble_size);
                                   orientor_folder=paste0(node_folder,'/',orientor);
                                   setwd(orientor_folder);
                                   ens_evolve_control=readRDS('ens_evolve_control')
                                   results_control<<-ens_evolve_control$get_results()
                                   results_control_list<<-lapply(results_control,function(x) x$state$extra$history)
                                   
                                   #getting rid of numeric errors
                                   eigens=lapply(results_control_list,function(x) return(unlist(map(x,'eigMax'))))
                                   
                                   bad=c()
                                   for (i in 1:length(results_control_list)){
                                     if(max(eigens[[i]])>=1){bad=c(bad,i)}
                                   }
                                   print(paste0('bad:',bad))
                                   if(length(bad)>0){
                                     for(i in 1:length(bad)){
                                       pos=min(which(eigens[[bad[i]]]>=1))
                                       
                                       results_control_list[[bad[i]]]<<-results_control_list[[bad[i]]][1:(pos-1)];
                                     }
                                   }
                                   
                                   rm(ens_evolve_control)
                                   gc()
                                   variables<<-names(results_control_list[[1]][[1]])
                                   setwd(main_folder)
                                   
                                   # vis<<-Visualization$new(results_list=results_control_list,variables=variables)
                                   
                                   return(list(results_control=results_control,results_control_list=results_control_list,variables=variables))
                                 },
                                 open_greedy=function(){
                                   setwd(main_folder)
                                   node_folder=paste0(main_folder,'/',ensemble_size);
                                   orientor_folder=paste0(node_folder,'/',orientor);
                                   print(orientor_folder)
                                   setwd(orientor_folder);
                                   ens_evolve_greedy=readRDS(paste0('ens_evolve_greedy_',greed_perc,'_',greed_orientor));
                                   results_greedy<<-ens_evolve_greedy$get_results();
                                   results_greedy_list<<-lapply(results_greedy,function(x) x$state$extra$history);
                                   
                                   #getting rid of numeric errors
                                   eigens=lapply(results_greedy_list,function(x) return(unlist(map(x,'eigMax'))))
                                   
                                   bad=c()
                                   for (i in 1:length(results_greedy_list)){
                                     if(max(eigens[[i]])>=1){bad=c(bad,i)}
                                   }
                                   print(paste0('bad:',bad))
                                   if(length(bad)>0){
                                     for(i in 1:length(bad)){
                                       pos=min(which(eigens[[bad[i]]]>=1))
                                       results_greedy_list[[bad[i]]]<<-results_greedy_list[[bad[i]]][1:(pos-1)];
                                     }
                                   }
                                   
                                   rm(ens_evolve_greedy);
                                   gc();
                                   variables<<-names(results_greedy_list[[1]][[1]]);
                                   setwd(main_folder);
                                   # vis_greedy<<-Visualization$new(results_list=results_greedy_list,variables=variables)
                                   return(list(results_greedy=results_greedy,results_greedy_list=results_greedy_list,variables=variables))
                                 },
                                 create_all_scalars=function(){},
                                 plot_both=function(net,y,iter=iterations){
                                   
                                   if(length(results_control_list)==0){open_control()}
                                   if(length(results_greedy_list)==0) {open_greedy()}
                                   
                                   vis=Visualization$new(results_list=results_control_list,variables=variables)
                                   vis_greedy=Visualization$new(results_list=results_greedy_list,variables=variables)
                                   
                                   cont = as.numeric(unlist(vis$data[[net]][y]))#[-1]#[-c(1:2)]
                                   gred = as.numeric(unlist(vis_greedy$data[[net]][y]))#[-1]#[-c(1:2)]
                                   mVal = min(c(cont,gred))
                                   MVal = max(c(cont,gred))
                                   plot(cont,col='blue',ylim=c(mVal,0),xlim=c(0,iter),cex=0.7)
                                   points(gred,col='red',cex=0.7)
                                   control = nls.control(maxiter = 5000,tol = 1e-20,minFactor = 1/9192,printEval = TRUE,warnOnly = TRUE)
                                   
                                   x_c = 1:length(cont)-1
                                   #nlmod = nls(cont ~ A-(A-cont[[1]])*exp(-x * B),start = list(A=-1,B=1),control = control)#,algorithm="port",upper=c(-1e-30,10000))
                                   nlmod = nls(cont ~ cont[[1]]*exp(-x_c * B),start = list(B=0.001),control = control)#,algorithm="port",upper=c(-1e-30,10000))
                                   
                                   x_g = 1:length(gred)-1
                                   #nlmod_g = nls(gred ~ A-(A-gred[[1]])*exp(-x * B),start = list(A=-1,B=1),control = control)#,algorithm="port",upper=c(-1e-30,10000))
                                   nlmod_g = nls(gred ~ gred[[1]]*exp(-x_g * B),start = list(B=0.001))#,control = control)#,algorithm="port",upper=c(-1e-30,10000))
                                   
                                   lines(1:6000,unlist(lapply(1:6000,function(x){cont[[1]]*exp(-x * summary(nlmod)$coefficients[[1]])})),col='green',lwd=4)
                                   lines(1:6000,unlist(lapply(1:6000,function(x){gred[[1]]*exp(-x * summary(nlmod_g)$coefficients[[1]])})),col='purple',lwd=4)
                                   
                                   # print(nlmod)
                                   # print(nlmod_g)
                                   print(summary(nlmod)$coefficients[[1]])
                                   print(summary(nlmod_g)$coefficients[[1]])
                                   print(cont[[1]]*exp(-num_fixs * summary(nlmod)$coefficients[[1]]))
                                   print(gred[[1]]*exp(-num_fixs * summary(nlmod_g)$coefficients[[1]]))
                                 }
                               )
)