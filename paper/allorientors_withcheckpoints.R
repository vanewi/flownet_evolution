source("paper/stat_util.R")

base_log=2
input_output_entropy=input_output_entropy_generator(log_base=base_log);
entropy_diff=entropy_diff_generator(log_base=base_log);

event_save_function = function(base,mat,iteration_num,characteristic_vals){
  event = list();
  event$tst = tst(base,mat);
  event$ami = ami(base,mat);
  eigens = base$get_eigen_mod(mat=mat);
  event$eigSD=sd(eigens);
  event$eigMax = max(eigens);
  event$b = log(event$eigMax);
  event$finn_node = base$get_finn_per_node(mat=mat);
  event$finn = sum(event$finn_node);
  event$time = iteration_num;
  event$ascendency=event$tst*event$ami;
  event$num_nodes=base$core_nodes;
  event$num_links=length(which(base$get_core(mat=mat)>0));
  event$in_links=length(which(mat[base$core_nodes+1,]>0));
  event$out_links=length(which(mat[,base$core_nodes+2]>0));
  entropy_analysis=input_output_entropy(base,mat);
  event$in_entropy=entropy_analysis$in_entropy;
  event$out_entropy=entropy_analysis$out_entropy;
  event$entropy_diff=event$out_entropy-event$in_entropy;
  event$ediff_norm=event$entropy_diff/base$get_general_entropy(rep(1,event$out_links),base_log);
  #event$mat = mat; #INCREASES MEMORY TOO MUCH
  event$finn_node_sd=sd(event$finn_node)
  event$var=names(event);
  return(event);
}

connector_f = connection_generator_with_renyi_opt(do_renyi = TRUE,renyi_connectance = 0.2);

evolution_f = evolution_function_generator(iteration_max=500,
                                           characteristics_boolean='and',
                                           event_saving_function=event_save_function
);

evolution_exp = evolution_function_generator(iteration_max=5000,
                                             characteristics_boolean='and',
                                             event_saving_function=event_save_function)

mutation_f = mutation_generator_base(loss_threshold = 0.2,gain_threshold = 0.2,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = FALSE,change_min_or_perturbation = 0.0,change_max = 0.9999);

mutation_perturbative = mutation_generator_base(loss_threshold = 0.0,gain_threshold = 0.0,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = TRUE,change_min_or_perturbation = 0.1,change_max = 0.9999);

node_ensembles=as.integer(c(100))#,70))#,100));
orientors=list(EDiff=entropy_diff)#Finn=finn,b=b,ASC=ascendency) #AMI=ami, TST=tst, 
greedy_orientors=list(part_stock=grow_part_of_stock,proportion_stock=grow_proportion_of_stock,part_finn=grow_part_finn_node,proportion_finn=grow_proportion_finn_node);
picked_greedy_perc=c(0.1,0.3,0.5)#,0.3,0.5)

num_checkpoints = 1;

PRINT_FLOWNET_DEBUG = FALSE

main_folder=getwd()

for (N in 1:length(node_ensembles)){ 
  gc();
  node_folder=paste0(main_folder,'/',node_ensembles[N])
  ifelse(!dir.exists(file.path(node_folder)), dir.create(file.path(node_folder)), FALSE)
  setwd(node_folder)
  
  ensemble=readRDS('ensemble')
# ensemble = NetEnsemble$new(ensemble_sizeMin = node_ensembles[N],
#                            ensemble_sizeMax = node_ensembles[N],
#                            size_replicate = 100L,
#                            nodes_autoloops_allowed = FALSE,
#                            size_interval= 1L);
# ensemble$generate(nodes_connect_function = connector_f,minimal_flowing = TRUE)
# saveRDS(object=ensemble,file = 'ensemble')
  
  for (or in 1:length(orientors)){ # PER ORIENTOR
    gc();
    orientor_folder=paste0(node_folder,'/',names(orientors[or]))
    ifelse(!dir.exists(file.path(orientor_folder)), dir.create(file.path(orientor_folder)), FALSE)
    setwd(orientor_folder)
    
    #NATURAL HISTORY EVOLUTION
    ens_evolve = EnsembleEvolve$new(ensemble=ensemble,
                                    mutation_func=mutation_f,
                                    evolution_func = evolution_f,
                                    characteristic_funcs = list(orientors[[or]]))
    
    system.time(ens_evolve$evolve_ensemble_par(num_cores = 8))
    
    saveRDS(object=ens_evolve,file='ens_evolve')
    
   #ens_evolve=readRDS('ens_evolve')
   # for (i in 1:length(ens_evolve_greedy$ensemble$generated_networks)){
   #   natural_ensemble$generated_networks[[i]]$net_data$picked_nodes= ens_evolve_greedy$ensemble$generated_networks[[i]]$net_data$picked_nodes
   # }
   # rm(ens_evolve_greedy)
   
    
    #CONTROL EXPERIMENT EVOLUTION
    natural_ensemble = ens_evolve$get_final_ensemble();
    rm(ens_evolve);
    gc();
    new_ensemble_control = natural_ensemble;
    dir.create(path = 'cp_control');
    for(cp in 1:num_checkpoints){
      ens_evolve_control = EnsembleEvolve$new(ensemble=new_ensemble_control,
                                              mutation_func=mutation_perturbative,
                                              evolution_func = evolution_exp,
                                              characteristic_funcs = list(orientors[[or]]))
      
      system.time(ens_evolve_control$evolve_ensemble_par(num_cores = 8))
      
      saveRDS(ens_evolve_control$get_results(),file=paste0('cp_control/',cp));
      new_ensemble_control = ens_evolve_control$get_final_ensemble();
      print(cat('CHECKPOINT: ',cp));
      gc();
    }
    
    saveRDS(ens_evolve_control,file='ens_evolve_control');
    rm(ens_evolve_control);
    gc();

 #run if only control_evolution desired:    
 # }
#} 
    # FORCED GREEDY EVOLUTION
    
    for (greed_perc in 1:length(picked_greedy_perc)){
      
      # FORCED SECTOR PICK PER NET
      pick_percent = picked_greedy_perc[greed_perc];
      for (i in 1:length(natural_ensemble$generated_networks)){
        num_nodes = natural_ensemble$generated_networks[[i]]$core_nodes
        to_pick = pick_percent * num_nodes
        if(to_pick<1)to_pick=1
        picked = sample(seq(1,num_nodes,1),to_pick,replace = FALSE)
        natural_ensemble$generated_networks[[i]]$net_data$picked_nodes=picked
      };
      gc();
      
      for(greed in 1:length(greedy_orientors)){
        new_ensemble_greedy=natural_ensemble;
        gc();
        name_cp_greedy=paste0('cp_greedy_',pick_percent,'_',names(greedy_orientors[greed]))
        dir.create(path = name_cp_greedy);
        for(cp in 1:num_checkpoints){
          gc();
          ens_evolve_greedy = EnsembleEvolve$new(ensemble=new_ensemble_greedy,
                                                 mutation_func=mutation_perturbative,
                                                 evolution_func = evolution_exp,
                                                 characteristic_funcs = list(greedy_orientors[[greed]]))
          
          system.time(ens_evolve_greedy$evolve_ensemble_par(num_cores = 8))
          
          saveRDS(ens_evolve_greedy$get_results(),file=paste0(name_cp_greedy,'/',cp));
          
          new_ensemble_greedy_tmp = ens_evolve_greedy$get_final_ensemble();
          
          for (i in 1:length(new_ensemble_greedy$generated_networks)){
            new_ensemble_greedy_tmp$generated_networks[[i]]$net_data$picked_nodes=new_ensemble_greedy$generated_networks[[i]]$net_data$picked_nodes
          }
          new_ensemble_greedy=new_ensemble_greedy_tmp
          rm(new_ensemble_greedy_tmp);
          print(cat('CHECKPOINT:',cp));
        }
        
        name_ens_evolve_greedy=paste0('ens_evolve_greedy_',pick_percent,'_',names(greedy_orientors[greed]))
        saveRDS(ens_evolve_greedy,file=name_ens_evolve_greedy);
        rm(ens_evolve_greedy);
      }
    }
  }
}

