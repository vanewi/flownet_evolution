source("FlowNet.R")

base_log=2
input_output_entropy=input_output_entropy_generator(log_base=base_log);
entropy_diff=entropy_diff_generator(log_base=base_log);

# WHEN CALCULATING CYCLES, ITS RECOMMENDED NOT TO USE LARGER THAN 30 NODE NETS
event_saving_function=function(base,mat,iteration_num,characteristic_vals){
  event = list();
  ###################################################################
  # Different cycle stats. Depends on searching for cycles. 
  # Not recommended for larger Nets or too many links.
  ###################################################################
  #event$c_analysis=base$get_fast_cycles_and_flow(mat=mat,initial_input=1000,save_cycle_flow_detail = FALSE);
  #event$cycles_num = event$c_analysis$tot_cycles;
  #event$c_flow=event$c_analysis$total_cycle_flow;
  #event$c_flow_tot=event$c_analysis$total_cycle_flow;
  #event$cycles_detail=event$c_analysis$cycles; # Depends on save_cycle_details parameter
  #event$c_flow_per_cycle=event$c_analysis$cycle_flow_per_cycle;
  #event$c_flow_per_link=event$c_analysis$cycle_flow_per_link;
  #event$c_flow_per_node=colSums(event$c_flow_per_link);
  ####################################################################
  entropy_analysis=input_output_entropy(base,mat);
  event$storage=storage(base,mat)
  event$mat = mat;
  event$tst = tst(base,mat);
  event$ami = ami(base,mat);
  event$finn_node = base$get_finn_per_node(mat=mat);
  event$b=base$get_exp_factor(mat=mat);
  event$eigMax=base$get_eigen_max(mat=mat);
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
  event$in_links=length(which(mat[base$core_nodes+1,]>0));
  event$out_links=length(which(mat[,base$core_nodes+2]>0));
  entropy_analysis=input_output_entropy(base,mat);
  event$in_entropy=entropy_analysis$in_entropy;
  event$out_entropy=entropy_analysis$out_entropy;
  event$entropy_diff=event$out_entropy-event$in_entropy;
  event$ediff_norm=event$entropy_diff/base$get_general_entropy(rep(1,event$out_links),base_log);
  event$var=names(event);
  return(event);
}

evolution_f = evolution_function_generator(iteration_max=500,
                                           characteristics_boolean='and',
                                           event_saving_function = event_saving_function
  );

# General mutation function with chances of adding or removing links
mutation_f = mutation_generator_base(loss_threshold = 0.2,gain_threshold = 0.2,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = FALSE,change_min_or_perturbation = 0.0,change_max = 0.9999);

# Perturbative mutation function. Does not allow adding/removing links and the new value is close to the previous one.
#mutation_perturbative = mutation_generator_base(loss_threshold = 0.0,gain_threshold = 0.0,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = TRUE,change_min_or_perturbation = 0.1,change_max = 0.9999);

connector_f = connection_generator_with_renyi_opt(do_renyi = TRUE,renyi_connectance = 0.2);

ensemble = NetEnsemble$new(ensemble_sizeMin = 10L,
                           ensemble_sizeMax = 20L,
                           size_replicate = 10L,
                           size_interval=5,
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

vis$plot_me2(x_var = 'time',
             y_var = 'tst',
             col_var='num_nodes',
             #color_func=function(x){return(hsv(h = log((x-min_max_col$min)/(min_max_col$max-min_max_col$min)*0.8), s = 1, v = 1, 1))},
             x_func = function(x) return(log(x)),
             y_func=function(y) return(y),
             with_raster = TRUE
             #y_lim = list(min=4,max=15)
)



