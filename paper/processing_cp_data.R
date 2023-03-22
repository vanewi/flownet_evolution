source("paper/stat_util.R")

main_folder=getwd()

base_log=2
input_output_entropy=input_output_entropy_generator(log_base=base_log);
entropy_diff=entropy_diff_generator(log_base=base_log);

node_ensemble=c(30)#,50,70,100)
orientors=list(TST=tst)#, AMI=ami,EDiff=entropy_diff)
greed_perc=c(0.1,0.3,0.5)
greed_orientor=c('part_stock')#,'proportion_stock','part_finn','proportion_finn')


#----------------
#----to create resilience table and y_s:----

for (N in 1:length(node_ensemble)) {
  cp_list_scalars_all = list()
  initial_final_bases=list()
  print(paste0('ensemble: ',node_ensemble[N]))
  for (or in 1:length(orientors)) {
    print(paste0('orientor: ',names(orientors[or])))
    
    con=Manipulation1Cp$new(ensemble_size=node_ensemble[N],orientor=names(orientors[or]))
    control=con$open_control()
    #cp_list_control=control$results_control_list
    n_nets=length(control$results_control_list)
    variables=names(control$results_control_list[[1]][[1]])
    
    cp_list_scalars_all[[names(orientors[or])]]=list(control = Visualization$new(results_list=control$results_control_list,variables=variables)$data)
    
    initials = lapply(control$results_control,function(x) x$initial$adjacency_matrix)
    finals= lapply(control$results_control,function(x) x$final$adjacency_matrix)
    initial_final_bases[[names(orientors[or])]]=list(control=list(initialA=initials,finalA=finals))

    
    for (greed_or in 1:length(greed_orientor)){
      print(greed_orientor[greed_or])
      for (greed_p in 1:length(greed_perc)){
        print(greed_perc[greed_p])
        
        greedy_path=paste0('cp_greedy_',greed_perc[greed_p],'_',greed_orientor[greed_or])
        greedy_path
        
        gred=Manipulation1Cp$new(ensemble_size=node_ensemble[N],orientor=names(orientors[or]),greed_orientor=greed_orientor[greed_or],greed_perc=greed_perc[greed_p])
        greedy=gred$open_greedy()
       # cp_list_greedy=greedy$results_greedy_list
        
        cp_list_scalars_all[[names(orientors[or])]][[greedy_path]]=Visualization$new(results_list=greedy$results_greedy_list,variables=variables)$data
        initials = lapply(greedy$results_greedy,function(x) x$initial$adjacency_matrix)
        finals= lapply(greedy$results_greedy,function(x) x$final$adjacency_matrix)
        picked=lapply(greedy$results_greedy,function(x) x$initial$net_data$picked_nodes)
        initial_final_bases[[names(orientors[or])]][[greedy_path]]=list(initialA=initials,finalA=finals,picked=picked)
      }
    }
    gc()
  }
  saveRDS(cp_list_scalars_all,paste0('all_scalars_',node_ensemble[N]))
  saveRDS(initial_final_bases,paste0('initial_final_bases',node_ensemble[N]))
}



