source("FlowNet.R")

# evolution algorithm functions

evolution_simple=function(iteration_num,base,current_state,current_net,mutated_net,characteristic_vals){
  cont = TRUE
  if(iteration_num>30) 
    cont = FALSE
  new_state = current_state;
  new_net = current_net;
  saveit = FALSE
  for(i in 1:length(characteristic_vals)){
    #natural selection rule: condition for saving
    if(characteristic_vals[[i]]>current_state$vals[[i]]){
      saveit = TRUE;
      break;
    }
  }
  if(saveit){
    new_state$vals = characteristic_vals;
    if(is.null(new_state$extra$fixations)){
      new_state$extra$fixations = 1
    }
    else{
      new_state$extra$fixations = new_state$extra$fixations + 1;
    }
    #saves the history of selected networks
    if(is.null(new_state$extra$history)){
      new_state$extra$history = list()
    }
    new_state$extra$history = append(new_state$extra$history,list(mutated_net-current_net));
    new_net = mutated_net;
  }
  return(list(new_state = new_state,
              new_net = new_net,
              continue = cont));
};

mutation_simple = function(base,mat){
  over_0 = which((mat[1:base$core_nodes,1:base$core_nodes])>0);
  if(length(over_0)==0)
    return(mat);
  link_to_change = sample(over_0,1);
  from = link_to_change%%(base$core_nodes);
  if(from==0)from = base$core_nodes; 
  to = ceiling(link_to_change/(base$core_nodes));
  mat_out = mat;
  mat_out[[from,to]]=runif(1,min=0.00001,max=0.99999);
  mat_out=base$normalize_rows(mat=mat_out,set_adjacency=FALSE);
  return(mat_out);
};


only_mean_time = function(base,mat){
  return(base$mean_time_discrete(input_percentage=0.05,mat=mat));
};

ensamble = NetEnsamble$new(ensemble_sizeMin = 4L,
                           ensemble_sizeMax = 7L,
                           size_replicate = 10L,
                           nodes_autoloops_allowed = FALSE);
ensamble$generate()

ens_evolve = EnsambleEvolve$new(ensamble=ensamble,
                                mutation_func=mutation_simple,
                                evolution_func = evolution_simple,
                                characteristic_funcs = list(only_mean_time))

ens_evolve$evolve_ensamble()


results = ens_evolve$get_results()

results[[4]]$final$adjacency_matrix

for(i in 1:length(results)){
  tst_0 = results[[i]]$initial$get_TST()
  tst_1 = results[[i]]$final$get_TST()
  print((tst_1-tst_0)/tst_0)
}
results[[6]]$initial$get_TST()
results[[1]]$initial$net_plot()
