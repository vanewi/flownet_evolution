source("FlowNet.R")

input_output_entropy=input_output_entropy_generator(log_base=2);
entropy_diff=entropy_diff_generator(log_base=2);

event_save_function = function(base,mat,iteration_num,characteristic_vals){
  event = list();
  event$tst = tst(base,mat);
  event$ami = ami(base,mat);
  event$b = base$get_exp_factor(mat=mat);
  event$eigMax = base$get_eigen_max(mat=mat);
  event$finn = finn(base,mat);
  event$time = iteration_num;
  event$ascendency=event$tst*event$ami;
  event$eigSD = sd(base$get_eigen_mod(mat=mat));
  event$num_nodes=base$core_nodes;
  event$num_links=length(which(base$get_core(mat=mat)>0));
  event$storage_pn=storage_per_node(base,mat);
  event$storage=sum(event$storage_pn);
  event$residence_time=event$storage/event$tst;
  event$finn_node = base$get_finn_per_node(mat=mat);
  entropy_analysis=input_output_entropy(base,mat);
  event$in_entropy=entropy_analysis$input_entropy;
  event$out_entropy=entropy_analysis$output_entropy;
  event$entropy_diff=entropy_diff(base,mat);
  event$c_flow_perc=(event$c_flow)/(event$tst)*100;
  event$c_flow_stock_norm=(event$c_flow_tot)/event$storage;
  event$mat = mat;
  event$var=names(event);
  return(event);
}

connector_f = connection_generator_with_renyi_opt(do_renyi = TRUE,renyi_connectance = 0.2);

evolution_f = evolution_function_generator(iteration_max=200,
                                           characteristics_boolean='and',
                                           event_saving_function=event_save_function
);

evolution_exp = evolution_function_generator(iteration_max=20,
                                             characteristics_boolean='and',
                                             event_saving_function=event_save_function)

mutation_f = mutation_generator_base(loss_threshold = 0.2,gain_threshold = 0.2,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = FALSE,change_min_or_perturbation = 0.0,change_max = 0.9999);

mutation_perturbative = mutation_generator_base(loss_threshold = 0.0,gain_threshold = 0.0,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = TRUE,change_min_or_perturbation = 0.1,change_max = 0.9999);

node_ensembles=as.integer(c(15))#,30,60,100));
orientors=list(TST=tst,AMI=ami,Finn=finn,EDiff=entropy_diff,b=b);
greedy_orientors=list(part_stock=grow_part_of_stock,proportion_stock=grow_proportion_of_stock);
picked_greedy_perc=c(0.3,0.1,0.5);

num_checkpoints = 10;

#NATURAL HISTORY EVOLUTION per orientor

#saveRDS(ensemble,file='random_ensemble')

for (N in 1:length(node_ensembles)){
  gc()
  ensemble = NetEnsemble$new(ensemble_sizeMin = node_ensembles[N],
                             ensemble_sizeMax = node_ensembles[N],
                             size_replicate = 1000L,
                             nodes_autoloops_allowed = FALSE,
                             size_interval= 1L);
  
  ensemble$generate(nodes_connect_function = connector_f,minimal_flowing = TRUE)
  
  for (or in 1:length(orientors)){
    gc()
    ens_evolve = EnsembleEvolve$new(ensemble=ensemble,
                                  mutation_func=mutation_f,
                                  evolution_func = evolution_f,
                                  characteristic_funcs = list(orientors[[or]]))
  
    system.time(ens_evolve$evolve_ensemble_par(num_cores = 8))
    results = ens_evolve$get_results()
    
    # FORCED GREEDY EVOLUTION
  
    for(greed in 1:length(greedy_orientors)){
      gc()
      for (greed_perc in 1:length(picked_greedy_perc)){
        gc()
      new_ensemble_control = ens_evolve$get_final_ensemble();
      
      new_ensemble_greedy = ens_evolve$get_final_ensemble();
  
      # FORCED SECTOR PICK PER NET
      pick_percent = picked_greedy_perc[greed_perc];
      for (i in 1:length(new_ensemble_greedy$generated_networks)){
        num_nodes = new_ensemble_greedy$generated_networks[[i]]$core_nodes
        to_pick = pick_percent * num_nodes
        if(to_pick<1)to_pick=1
        picked = sample(seq(1,num_nodes,1),to_pick,replace = FALSE)
        new_ensemble_greedy$generated_networks[[i]]$net_data$picked_nodes=picked
      }
  
    results_control=vector(mode='list',length=num_checkpoints);
    results_greedy=vector(mode='list',length=num_checkpoints);

    gc()
    for(cp in 1:num_checkpoints){
      ens_evolve_greedy = EnsembleEvolve$new(ensemble=new_ensemble_greedy,
                                           mutation_func=mutation_perturbative,
                                           evolution_func = evolution_exp,
                                           characteristic_funcs = list(greedy_orientors[[greed]]))
    
      ens_evolve_control = EnsembleEvolve$new(ensemble=new_ensemble_control,
                                            mutation_func=mutation_perturbative,
                                            evolution_func = evolution_exp,
                                            characteristic_funcs = list(orientors[[or]]))
    
      system.time(ens_evolve_greedy$evolve_ensemble_par())
    
      gc()
    
      system.time(ens_evolve_control$evolve_ensemble_par())
    
      results_greedy[[cp]] = ens_evolve_greedy$get_results()
      results_control[[cp]] = ens_evolve_control$get_results()
  
      #OJO=PUEDEN NO EVOLUCIONAR UN CHECKPOINT Y QUEDA COMO NULL
          
      new_ensemble_control = ens_evolve_control$get_final_ensemble();
      new_ensemble_greedy_tmp = ens_evolve_greedy$get_final_ensemble();
    
      for (i in 1:length(new_ensemble_greedy$generated_networks)){
        new_ensemble_greedy_tmp$generated_networks[[i]]$net_data$picked_nodes=new_ensemble_greedy$generated_networks[[i]]$net_data$picked_nodes
      }
      gc()
      new_ensemble_greedy=new_ensemble_greedy_tmp
      gc()
      print(cat('CHECKPOINT:',cp));
    }
    
    save_name=paste0('Ensemble_',ens_evolve$ensemble$name,
                 '_Orientor_',names(orientors[or]), 
                 '_GreedExp_',names(greedy_orientors[greed]),
                 '_PickedPerc_',picked_greedy_perc[greed_perc])
    
    main_dir=getwd();
    sub_dir=paste0('/',ens_evolve$ensemble$name);
    ifelse(!dir.exists(file.path(main_dir, sub_dir)), dir.create(file.path(main_dir, sub_dir)), FALSE)
    file_dir=paste0('.',sub_dir,'/',save_name,'.RDS')  
    print(file_dir)
    saveRDS(list(ensemble=ensemble,ens_evolve=ens_evolve,new_ensemble_control=new_ensemble_control,
                 new_ensemble_greedy=new_ensemble_greedy,results=results,
                 results_greedy=results_greedy,results_control=results_control),
            file = file_dir)
    }
  }
}
}

#RESULTS
#WORKING PER EXPERIMENT

open_name=opening(15,100,1,'EDiff','part_stock',0.3)
experiment <- readRDS(file=open_name)
attach(experiment)

results_list=lapply(results,function(x){x$state$extra$history})
variables=results_list[[1]][[1]]$var;
n_nets=length(ens_evolve$evolved_nets); #length(results_control[[1]])

vis = Visualization$new(results_list=results_list,variables=variables)

vis$plot_me2(x_var = 'time',
             y_var = 'finn',
             col_var='entropy_diff',
             #color_func=function(x){return(hsv(h = log((x-min_max_col$min)/(min_max_col$max-min_max_col$min)*0.8), s = 1, v = 1, 1))},
             x_func = function(x) return(exp(x)),
             y_func = function(x) return(x)
             #y_lim = list(min=4,max=15)
)


#for greedy ensembles:
cp_control=vector(mode='list',length=(length(results_control)+1))
cp_greedy=vector(mode='list',length=(length(results_control)+1))
all_control=vector(mode='list',length=n_nets)
all_greedy=vector(mode='list',length=n_nets)
for (cp in 1:(length(results_control)+1)){
  for (net in 1:n_nets){
    if(cp==1){
      end=results[[net]]$state$extra$fixations
      state_control=results[[net]]$state$extra$history[[end]]
      state_greedy=results[[net]]$state$extra$history[[end]]
    }
    else{
      end=results_control[[cp-1]][[net]]$state$extra$fixations
      state_control=results_control[[cp-1]][[net]]$state$extra$history[[end]]
      end=results_greedy[[cp-1]][[net]]$state$extra$fixations
      state_greedy=results_greedy[[cp-1]][[net]]$state$extra$history[[end]]
      all_control[[net]]=c(all_control[[net]],results_control[[cp-1]][[net]]$state$extra$history)
      all_greedy[[net]]=c(all_greedy[[net]],results_greedy[[cp-1]][[net]]$state$extra$history)
    }
    if(length(state_control)==0){
      end=results_control[[cp-2]][[net]]$state$extra$fixations
      state_control=results_control[[cp-2]][[net]]$state$extra$history[[end]]
    }
    if(length(state_greedy)==0){
      end=results_greedy[[cp-2]][[net]]$state$extra$fixations
      state_greedy=results_greedy[[cp-2]][[net]]$state$extra$history[[end]]
    }
    cp_control[[cp]][[net]]=state_control;
    cp_greedy[[cp]][[net]]=state_greedy;
  }
}

cp_list_control=vector(mode='list',length=n_nets)
cp_list_greedy=vector(mode='list',length=n_nets)
cp_list_control_scalars=vector(mode='list',length=n_nets)
cp_list_greedy_scalars=vector(mode='list',length=n_nets)
for (net in 1:n_nets){
  cp_list_control[[net]]=t(sapply(cp_control,'[[',net))
  cp_list_greedy[[net]]=t(sapply(cp_greedy,'[[',net))
  nonscalars=c()
  for (i in 1:length(colnames(cp_list_control[[net]]))){
    if(length(cp_list_control[[net]][[1,colnames(cp_list_control[[net]])[i]]])>1){
      nonscalars=c(nonscalars,i)
    }
  }
  cp_control_reduced = cp_list_control[[net]][,-nonscalars]
  cp_greedy_reduced = cp_list_greedy[[net]][,-nonscalars]
  cp_control_mat = matrix(data=unlist(cp_control_reduced),nrow=nrow(cp_control_reduced),ncol=ncol(cp_control_reduced),dimnames = dimnames(cp_control_reduced))
  cp_greedy_mat = matrix(data=unlist(cp_greedy_reduced),nrow=nrow(cp_greedy_reduced),ncol=ncol(cp_greedy_reduced),dimnames = dimnames(cp_greedy_reduced))
  cp_list_control_scalars[[net]]=as.data.frame(cp_control_mat)
  cp_list_greedy_scalars[[net]]=as.data.frame(cp_greedy_mat)
}

find_slope=function(dat){
  variables=names(dat[[1]]);
  slopes=matrix(nrow=length(dat),ncol=length(variables))
  for (var in 1:length(variables)){
    obj=vector(mode='list',length=length(data));
    for (i in 1:length(dat)){
      cp=length(dat[[i]][,1]);
      dumX=1:cp;
      colum=variables[var]
      v=dat[[i]][[colum]]
      #slopes[i,var]=lm(as.matrix(eval)~dumX)$coefficients[[2]]
      #fixed intercept (starting point):
      slopes[i,var]=lm(I(v-v[1])~0+dumX)$coefficients[[1]]
      slopes=as.data.frame(slopes)
      names(slopes)=variables
    }
  }
  return(as.data.frame(slopes))
}

slopes_control=find_slope(dat=cp_list_control_scalars)
slopes_greedy=find_slope(dat=cp_list_greedy_scalars)

slopes_results=find_slope(dat=vis$data)

x=(slopes_greedy$b-slopes_control$b)#/slopes_control$b#*100
y=(slopes_greedy$finn-slopes_control$finn)#/slopes_control$finn*100

par(mfrow=c(1,1))
plot(x,y)
abline(v=0,h=0)

