source("FlowNet.R")
library(purrr)

# "ECOSYSTEM EVOLUTION"

mean_time = mean_time_generator(input_percentage = 0.05);

event_save_function = function(base,mat,iteration_num,characteristic_vals){
  event = list();
  #event$mat = mat;
  event$tst = tst(base,mat);
  event$mean_time = mean_time(base,mat);
  event$ami = ami(base,mat);
  event$b = base$get_exp_factor(mat=mat);
  event$eigMax = base$get_eigen_max(mat=mat);
  event$finn = finn(base,mat);
  event$time = iteration_num;
  event$ascendency=event$tst*event$ami;
  event$eigSD = sd(base$get_eigen_mod(mat=mat));
  return(event);
}

evolution_f = evolution_function_generator(iteration_max=100,
                                           characteristics_boolean='and',
                                           event_saving_function=event_save_function
);

mutation_f = mutation_generator_base(loss_threshold = 0.2,gain_threshold = 0.2,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = FALSE,change_min_or_perturbation = 0.0,change_max = 0.9999);

connector_f = connection_generator_with_renyi_opt(do_renyi = TRUE,renyi_connectance = 0.2);

ensamble = NetEnsamble$new(ensemble_sizeMin = 100L,
                           ensemble_sizeMax = 100L,
                           size_replicate = 16L,
                           nodes_autoloops_allowed = FALSE);

ensamble$generate(nodes_connect_function = connector_f,minimal_flowing = TRUE)

ens_evolve = EnsambleEvolve$new(ensamble=ensamble,
                                mutation_func=mutation_f,
                                evolution_func = evolution_f,
                                characteristic_funcs = list(ascendency))

system.time(ens_evolve$evolve_ensamble_par(num_cores = 8))

results = ens_evolve$get_results()

cada_cuantos = 0
offset = 1
inicial = offset
final = inicial+cada_cuantos#length(results)
for(k in inicial:final){
  bn = unlist(map(results[[k]]$state$extra$history,'b'))
  #bn = (bn-min(bn))/(max(bn)-min(bn))-0.1*k
  if(k==inicial)
    plot(bn,type='l')
  else
    lines(bn,type='l')
}
results[[1]]$final$get_series_of_node_stock_after_stabilization(initial_input = 1000.0,frequency = 30,do_plot = TRUE)
map(results[[1]]$state$extra$history,'b')
# FORCED GREEDY EVOLUTION
results[[1]]$final$get_exp_factor()

new_ensamble_control = ens_evolve$get_final_ensamble();
new_ensamble_greedy = ens_evolve$get_final_ensamble();

evolution_exp = evolution_function_generator(iteration_max=20,
                                                characteristics_boolean='and',
                                                event_saving_function=event_save_function
);

mutation_perturbative = mutation_generator_base(loss_threshold = 0.0,gain_threshold = 0.0,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = TRUE,change_min_or_perturbation = 0.1,change_max = 0.9999);

# FORCED SECTOR PICK PER NET
pick_percent = 1.0
for (i in 1:length(new_ensamble_greedy$generated_networks)){
  num_nodes = new_ensamble_greedy$generated_networks[[i]]$core_nodes
  to_pick = floor(runif(1,max=pick_percent) * num_nodes)
  if(to_pick<1)to_pick=1
  picked = sample(seq(1,num_nodes,1),to_pick,replace = FALSE)
  new_ensamble_greedy$generated_networks[[i]]$net_data$picked_nodes=picked
}

grow_part_of_stock=function(base,mat){
  stock=base$get_initial_stock_after_stabilization(mat=mat,initial_input=1000.0,frequency=1L)[1:base$core_nodes];
  return(sum(stock[base$net_data$picked_nodes]))
}

grow_proportion_of_stock=function(base,mat){
  stock=base$get_initial_stock_after_stabilization(mat=mat,initial_input=1000.0,frequency=1L)[1:base$core_nodes];
  return(sum(stock[base$net_data$picked_nodes])/sum(stock[1:base$core_nodes]))
}

num_checkpoints = 10

results_control=vector(mode='list',length=num_checkpoints);
results_greedy=vector(mode='list',length=num_checkpoints);
gc()
for(cp in 1:num_checkpoints){
  ens_evolve_greedy = EnsambleEvolve$new(ensamble=new_ensamble_greedy,
                                         mutation_func=mutation_perturbative,
                                         evolution_func = evolution_exp,
                                         characteristic_funcs = list(grow_part_of_stock,grow_proportion_of_stock))
  
  ens_evolve_control = EnsambleEvolve$new(ensamble=new_ensamble_control,
                                          mutation_func=mutation_perturbative,
                                          evolution_func = evolution_exp,
                                          characteristic_funcs = list(ascendency))
  
  system.time(ens_evolve_greedy$evolve_ensamble_par())
  
  gc()
  
  system.time(ens_evolve_control$evolve_ensamble_par())
  
  results_greedy[[cp]] = ens_evolve_greedy$get_results()
  results_control[[cp]] = ens_evolve_control$get_results()
  
  new_ensamble_control = ens_evolve_control$get_final_ensamble();
  new_ensamble_greedy_tmp = ens_evolve_greedy$get_final_ensamble();
  
  for (i in 1:length(new_ensamble_greedy$generated_networks)){
    new_ensamble_greedy_tmp$generated_networks[[i]]$net_data$picked_nodes=new_ensamble_greedy$generated_networks[[i]]$net_data$picked_nodes
  }
  gc()
  new_ensamble_greedy=new_ensamble_greedy_tmp
  gc()
  print(cat('CHECKPOINT:',cp));
}


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
desde = 60
cuantas = 30#length(results_control[[1]]) 
for(n in desde:(desde+cuantas)){
  lines(h_dExpF[[n]])
}

