source("FlowNet.R")

# "ECOSYSTEM EVOLUTION"

#tst
#ascendency
storage_per_node=storage_per_node_generator(initial_input=1000,frequency=1L);
storage=storage_generator(initial_input=1000,frequency=1L);
#mean_time = mean_time_generator(input_percentage = 0.05);
residence_time=residence_time_generator(initial_input=1000,frequency=1L);

event_save_function = function(base,mat,iteration_num,characteristic_vals){
  event = list();
  #cycles_analyses=base$get_fast_cycles_and_flow(mat=mat,save_cycle_flow_detail = FALSE,initial_input=1000L);
  #event$cycles_num=cycles_analyses$tot_cycles;
  #event$cflow_tot=cycles_analyses$total_cycle_flow;
  #event$cycle_flow_per_link=cycles_analyses$cycle_flow_per_link;
  event$storage=storage(base,mat);
  event$tst = tst(base,mat);
  event$ln_tst=log(tst(base,mat));
  #event$mean_time = mean_time(base,mat);
  event$residence_time=residence_time(base,mat);
  event$ami = ami(base,mat);
  event$b = base$get_exp_factor(mat=mat);
  event$eigMax = base$get_eigen_max(mat=mat);
  event$finn = finn(base,mat);
  event$time = iteration_num;
  event$ascendency=event$tst*event$ami;
  event$eigSD = sd(base$get_eigen_mod(mat=mat));
  event$num_nodes=base$core_nodes;
  event$num_links=length(which(base$get_core(base$adjacency_matrix)>0));
  event$storage_pn=storage_per_node(base,mat);
  event$storage=sum(event$storage_pn);
  event$finn_node = base$get_finn_per_node(mat=mat);
  event$c_flow_perc=(event$c_flow)/(event$tst)*100;
  event$mat = mat;
  event$var=names(event);
  return(event);
}

evolution_f = evolution_function_generator(iteration_max=50,
                                           characteristics_boolean='or',
                                           event_saving_function=event_save_function
);

mutation_f = mutation_generator_base(loss_threshold = 0.2,gain_threshold = 0.2,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = FALSE,change_min_or_perturbation = 0.0,change_max = 0.9999);

connector_f = connection_generator_with_renyi_opt(do_renyi = TRUE,renyi_connectance = 0.2);

ensemble = NetEnsemble$new(ensemble_sizeMin = 8L,
                           ensemble_sizeMax = 8L,
                           size_replicate = 25L,
                           nodes_autoloops_allowed = FALSE,
                           size_interval= 1L);

ensemble$generate(nodes_connect_function = connector_f,minimal_flowing = TRUE)

ens_evolve = EnsembleEvolve$new(ensemble=ensemble,
                                mutation_func=mutation_f,
                                evolution_func = evolution_f,
                                characteristic_funcs = list(tst))

system.time(ens_evolve$evolve_ensemble_par(num_cores = 8))

# TO SAVE AND LOAD RESULTS
#saveRDS(ens_evolve, file = "15nodes_500replicates_02connect_02gainAndLoss_randomToStructure_noCycles_allorientors")
#results <- readRDS("")

#PLOTTING
n_nets=length(ens_evolve$evolved_nets)
results = ens_evolve$get_results()
results_list=lapply(results,function(x){x$state$extra$history})
variables=names(results[[1]]$state$extra$history[[1]]);

vis=Visualization$new(results_list=results_list,variables=variables)
# FOR SUBSETS
#vis = Visualization3$new(in_data=por_mat[lapply(por_mat,function(x) if(length(x)==0) 0 else x[[1]]$num_nodes)==10])

vis$plot_me2(x_var = 'finn',
            y_var = 'c_flow_perc',
            y_func= function(y) return(y*100),
            col_var = 'tst',
            size_var = NULL,
            x_lim = NULL, #list(min=0,max=7000), 
            y_lim = NULL#list(min=0.0,max=1000000)
)

# FORCED GREEDY EVOLUTION

new_ensemble_control = ens_evolve$get_final_ensemble();
new_ensemble_greedy = ens_evolve$get_final_ensemble();

evolution_exp = evolution_function_generator(iteration_max=20,
                                                characteristics_boolean='and',
                                                event_saving_function=event_save_function
);

mutation_perturbative = mutation_generator_base(loss_threshold = 0.0,gain_threshold = 0.0,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = TRUE,change_min_or_perturbation = 0.1,change_max = 0.9999);

# FORCED SECTOR PICK PER NET
#random percentage pick of nodes
pick_percent = 1.0
for (i in 1:length(new_ensemble_greedy$generated_networks)){
  num_nodes = new_ensemble_greedy$generated_networks[[i]]$core_nodes
  to_pick = floor(runif(1,max=pick_percent) * num_nodes)
  if(to_pick<1)to_pick=1
  picked = sample(seq(1,num_nodes,1),to_pick,replace = FALSE)
  new_ensemble_greedy$generated_networks[[i]]$net_data$picked_nodes=picked
}

#fixed percentage pick of nodes
pick_percent = 0.3 #0.5 $0.1
for (i in 1:length(new_ensemble_greedy$generated_networks)){
  num_nodes = new_ensemble_greedy$generated_networks[[i]]$core_nodes
  to_pick = pick_percent * num_nodes
  if(to_pick<1)to_pick=1
  picked = sample(seq(1,num_nodes,1),to_pick,replace = FALSE)
  new_ensemble_greedy$generated_networks[[i]]$net_data$picked_nodes=picked
}

num_checkpoints = 2

results_control=vector(mode='list',length=num_checkpoints);
results_greedy=vector(mode='list',length=num_checkpoints);

gc()
for(cp in 1:num_checkpoints){
  ens_evolve_greedy = EnsembleEvolve$new(ensemble=new_ensemble_greedy,
                                         mutation_func=mutation_perturbative,
                                         evolution_func = evolution_exp,
                                         characteristic_funcs = list(grow_part_of_stock))
  
  ens_evolve_control = EnsembleEvolve$new(ensemble=new_ensemble_control,
                                          mutation_func=mutation_perturbative,
                                          evolution_func = evolution_exp,
                                          characteristic_funcs = list(tst))
  
  system.time(ens_evolve_greedy$evolve_ensemble_par())
  
  gc()
  
  system.time(ens_evolve_control$evolve_ensemble_par())
  
  results_greedy[[cp]] = ens_evolve_greedy$get_results()
  results_control[[cp]] = ens_evolve_control$get_results()
  
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

cp_control=vector(mode='list',length=(length(results_control)+1))
cp_greedy=vector(mode='list',length=(length(results_control)+1))
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
  cp=length(dat[[1]][,1]);
  dumX=1:cp;
  slopes=matrix(nrow=length(dat),ncol=length(variables))
  for (var in 1:length(variables)){
    obj=vector(mode='list',length=length(data));
    for (i in 1:length(dat)){
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

      
#-------------------------------------
#de nico antes:

# SKETCHY RESULT ANALYSIS

dExpF = vector(mode='list',length=length(num_checkpoints))
dPStock = vector(mode='list',length=length(num_checkpoints))
dPFinn = vector(mode='list',length=length(num_checkpoints))
pStockI = vector(mode='list',length=length(num_checkpoints))
pStockF = vector(mode='list',length=length(num_checkpoints))
pickedNum = vector(mode='list',length=length(num_checkpoints))
finnProp = vector(mode='list',length=length(num_checkpoints))
finnPropAv = vector(mode='list',length=length(num_checkpoints))
finnPropF = vector(mode='list',length=length(num_checkpoints))
finnPropFAv = vector(mode='list',length=length(num_checkpoints))

for(cp in 1:num_checkpoints){
  dExpF[[cp]] = vector(mode='list',length=length(results_control[[cp]]))
  dPStock[[cp]] = vector(mode='list',length=length(results_control[[cp]]))
  dPFinn[[cp]] = vector(mode='list',length=length(results_control[[cp]]))
  pStockI[[cp]] = vector(mode='list',length=length(results_control[[cp]]))
  pStockF[[cp]] = vector(mode='list',length=length(results_control[[cp]]))
  pickedNum[[cp]] = vector(mode='list',length=length(results_control[[cp]]))
  finnProp[[cp]] = vector(mode='list',length=length(results_control[[cp]]))
  finnPropAv[[cp]] = vector(mode='list',length=length(results_control[[cp]]))
  finnPropF[[cp]] = vector(mode='list',length=length(results_control[[cp]]))
  finnPropFAv[[cp]] = vector(mode='list',length=length(results_control[[cp]]))
  for(n in 1:length(results_control[[cp]])){
    nodes = results_greedy[[cp]][[n]]$initial$core_nodes
    picked = results_greedy[[cp]][[n]]$initial$net_data$picked_nodes
    pickedNum[[n]]=length(picked)
    b_g = results_greedy[[cp]][[n]]$final$get_exp_factor()
    b_c = results_control[[cp]][[n]]$final$get_exp_factor()
    dExpF[[cp]][[n]] = (b_c-b_g)/abs(b_c)*100
    
    stockI = results_greedy[[cp]][[n]]$initial$get_initial_stock_after_stabilization(initial_input = 1000.0,frequency = 1L)[1:nodes]
    stockF = results_greedy[[cp]][[n]]$final$get_initial_stock_after_stabilization(initial_input = 1000.0,frequency = 1L)[1:nodes]
    
    finnI = results_greedy[[cp]][[n]]$initial$get_finn_per_node(input_multiplier = 1000.0);
    finnF = results_greedy[[cp]][[n]]$final$get_finn_per_node(input_multiplier = 1000.0);
    
    finnProp[[cp]][[n]] = sum(finnI[picked])/sum(finnI)
    
    finnPropF[[cp]][[n]] = sum(finnF[picked])/sum(finnF)
    
    finnPropAv[[cp]][[n]] = finnProp[[cp]][[n]]/length(picked)
    
    finnPropFAv[[n]] = finnPropF[[cp]][[n]]/length(picked)
    
    dPFinn[[cp]][[n]] = sum(finnF[picked])/sum(finnF)*100-sum(finnI[picked])/sum(finnI)*100
    
    pStockI[[cp]][[n]] = sum(stockI[picked])/sum(stockI)*100
    pStockF[[cp]][[n]] = sum(stockF[picked])/sum(stockF)*100
    dPStock[[cp]][[n]] = pStockF[[cp]][[n]]-pStockI[[cp]][[n]]
  }
}




h_dExpF = vector(mode='list',length=length(results_control[[1]]))

for(n in 1:length(results_control[[1]])){
  h_dExpF[[n]]=unlist(map(dExpF,n))
}

minY = min(unlist(h_dExpF))
maxY = max(unlist(h_dExpF))

plot(h_dExpF[[1]],type = 'line',ylim = c(minY,maxY))
desde = 1
cuantas = 3#length(results_control[[1]]) 
for(n in desde:(desde+cuantas)){
  lines(h_dExpF[[n]])
}

