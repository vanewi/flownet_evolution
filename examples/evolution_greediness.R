source("FlowNet.R")
library(purrr)

# "ECOSYSTEM EVOLUTION"

mean_time = mean_time_generator(input_percentage = 0.05);

event_save_function = function(base,mat,iteration_num,characteristic_vals){
  event = list();
  event$mat = mat;
  event$tst = tst(base,mat);
  event$mean_time = mean_time(base,mat);
  event$ami = ami(base,mat);
  event$b=base$get_exp_factor(mat=mat);
  event$eigMax=base$get_eigen_max(mat=mat);
  event$finn=finn(base,mat);
  event$time = iteration_num;
  event$ascendency=event$tst*event$ami;
  event$eigSD = sd(base$get_eigen_mod(mat=mat));
  return(event);
}

evolution_f = evolution_function_generator(iteration_max=30,
                                           characteristics_boolean='and',
                                           event_saving_function=event_save_function
);

mutation_f = mutation_generator_base(loss_threshold = 0.2,gain_threshold = 0.2,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = FALSE,change_min_or_perturbation = 0.0,change_max = 0.9999);

connector_f = connection_generator_with_renyi_opt(do_renyi = TRUE,renyi_connectance = 0.2);

ensemble = NetEnsemble$new(ensemble_sizeMin = 75L,
                           ensemble_sizeMax = 85L,
                           size_replicate = 4L,
                           nodes_autoloops_allowed = FALSE);

ensemble$generate(nodes_connect_function = connector_f,minimal_flowing = TRUE)

ens_evolve = EnsembleEvolve$new(ensemble=ensemble,
                                mutation_func=mutation_f,
                                evolution_func = evolution_f,
                                characteristic_funcs = list(ascendency))

system.time(ens_evolve$evolve_ensemble_par())

results = ens_evolve$get_results()

# FORCED GREEDY EVOLUTION

new_ensemble = ens_evolve$get_final_ensemble();

evolution_greedy = evolution_function_generator(iteration_max=100,
                                           characteristics_boolean='and',
                                           event_saving_function=event_save_function
);

mutation_perturbative = mutation_generator_base(loss_threshold = 0.0,gain_threshold = 0.0,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = TRUE,change_min_or_perturbation = 0.1,change_max = 0.9999);

# FORCED SECTOR PICK PER NET
pick_percent = 1.0
for (i in 1:length(new_ensemble$generated_networks)){
  num_nodes = new_ensemble$generated_networks[[i]]$core_nodes
  to_pick = floor(runif(1,max=pick_percent) * num_nodes)
  if(to_pick<1)to_pick=1
  picked = sample(seq(1,num_nodes,1),to_pick,replace = FALSE)
  new_ensemble$generated_networks[[i]]$net_data$picked_nodes=picked
}

grow_part_of_stock=function(base,mat){
  stock=base$get_initial_stock_after_stabilization(mat=mat,initial_input=1000.0,frequency=1L)[1:base$core_nodes];
  return(sum(stock[base$net_data$picked_nodes]))
}

grow_proportion_of_stock=function(base,mat){
  stock=base$get_initial_stock_after_stabilization(mat=mat,initial_input=1000.0,frequency=1L)[1:base$core_nodes];
  return(sum(stock[base$net_data$picked_nodes])/sum(stock[1:base$core_nodes]))
}

ens_evolve_greedy = EnsembleEvolve$new(ensemble=new_ensemble,
                                 mutation_func=mutation_perturbative,
                                 evolution_func = evolution_greedy,
                                 characteristic_funcs = list(grow_part_of_stock,grow_proportion_of_stock))

ens_evolve_control = EnsembleEvolve$new(ensemble=new_ensemble,
                                       mutation_func=mutation_perturbative,
                                       evolution_func = evolution_greedy,
                                       characteristic_funcs = list(ascendency))

system.time(ens_evolve_greedy$evolve_ensemble_par())

system.time(ens_evolve_control$evolve_ensemble_par())

results_greedy = ens_evolve_greedy$get_results()
results_control = ens_evolve_control$get_results()

