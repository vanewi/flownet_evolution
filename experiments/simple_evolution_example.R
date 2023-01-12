source("FlowNet.R")

#mean_time=mean_time_generator(input_percentge=0.5)
input_output_entropy=input_output_entropy_generator(log_base=2);
entropy_diff=entropy_diff_generator(log_base=2);

# AS THE EVENT SAVING CALCULATES CYCLES, DO NOT USE LARGER THAN 30 NODE NETS
event_saving_function=function(base,mat,iteration_num,characteristic_vals){
  event = list();
  #ent = c_ent(base,mat);
  #event$trajectory_num = length(ent$netdata$trajectories);
  #event$ent_tot = ent$link_type_entropy$total_link_ent;
  #event$ent_av = ent$link_type_entropy$average_ent_per_link;
  #event$mean_time = mean_time(base,mat);
  #event$residence_time=residence_time(base,mat);
  #event$c_analysis=base$get_fast_cycles_and_flow(mat=mat,initial_input=1000,save_cycle_flow_detail = FALSE);
  #event$c_flow=event$c_analysis$total_cycle_flow;
  #event$c_flow_tot=event$c_analysis$total_cycle_flow;
  #event$cycles_detail=event$c_analysis$cycles;
  #event$c_flow_per_cycle=event$c_analysis$cycle_flow_per_cycle;
  #event$c_flow_per_link=event$c_analysis$cycle_flow_per_link;
  #event$c_flow_per_node=colSums(event$c_flow_per_link);
  entropy_analysis=input_output_entropy(base,mat);
  event$in_entropy=entropy_analysis$input_entropy;
  event$out_entropy=entropy_analysis$output_entropy;
  event$entropy_diff=entropy_diff(base,mat);
  event$storage=storage(base,mat)
  event$mat = mat;
  event$tst = tst(base,mat);
  event$ami = ami(base,mat);
  event$finn_node = base$get_finn_per_node(mat=mat);
  event$b=base$get_exp_factor(mat=mat);
  event$eigMax=base$get_eigen_max(mat=mat);
  event$cycles_num = event$c_analysis$tot_cycles;
  event$finn=finn(base,mat);
  event$time = iteration_num;
  event$storage_pn=storage_per_node(base,mat);
  event$storage=sum(event$storage_pn);
  event$ascendency=event$tst*event$ami;
  event$eigSD = sd(base$get_eigen_mod(mat=mat));
  event$c_flow_perc=(event$c_flow)/(event$tst)*100;
  event$c_flow_stock_norm=(event$c_flow_tot)/event$storage;
  event$num_nodes=base$core_nodes;
  event$num_links=length(which(base$get_core(base$adjacency_matrix)>0));
  event$var=names(event);
  return(event);
}

evolution_f = evolution_function_generator(iteration_max=100,
                                           characteristics_boolean='and',
                                           event_saving_function = event_saving_function
  );

mutation_f = mutation_generator_base(loss_threshold = 0.2,gain_threshold = 0.2,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = FALSE,change_min_or_perturbation = 0.0,change_max = 0.9999);

connector_f = connection_generator_with_renyi_opt(do_renyi = TRUE,renyi_connectance = 0.2);

ensemble = NetEnsemble$new(ensemble_sizeMin = 15L,
                           ensemble_sizeMax = 15L,
                           size_replicate = 50L,
                           size_interval=1,
                           nodes_autoloops_allowed = FALSE);

ensemble$generate(nodes_connect_function = connector_f,minimal_flowing = TRUE)

ens_evolve = EnsembleEvolve$new(ensemble=ensemble,
                                mutation_func=mutation_f,
                                evolution_func = evolution_f,
                                characteristic_funcs = list(tst))

system.time(ens_evolve$evolve_ensemble())

results = ens_evolve$get_results()

# TO SAVE AND LOAD RESULTS
#saveRDS(list(results=results,ens_evolve=ens_evolve), file = "10_20_100_5_Orientor_Finn_Ensemble_NH_Evolution.RDS")
#results <- readRDS("resPrueba.RDS")

results_list=lapply(results,function(x){x$state$extra$history})
n_nets=length(ens_evolve$evolved_nets)
variables=results_list[[1]][[1]]$var

vis = Visualization$new(results_list=results_list,variables=variables)

vis$plot_me2(x_var = 'entropy_diff',
             y_var = 'finn',
             col_var='time',
             #color_func=function(x){return(hsv(h = log((x-min_max_col$min)/(min_max_col$max-min_max_col$min)*0.8), s = 1, v = 1, 1))},
             x_func = function(x) return(x),
             y_func = function(x) return(x)
             #y_lim = list(min=4,max=15)
)



