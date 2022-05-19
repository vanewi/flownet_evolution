source("FlowNet.R")

# evolution algorithm functions

evolution_simple=function(iteration_num,base,current_state,current_net,mutated_net,characteristic_vals){
  cont = TRUE
  if(iteration_num>30) 
    cont = FALSE
  new_state = current_state;
  new_net = current_net;
  saveit = TRUE
  # "AND" OF ALL CHARACTERISTICS
  for(i in 1:length(characteristic_vals)){
    #natural selection rule: condition for saving
    if(characteristic_vals[[i]]<current_state$vals[[i]]){
      saveit = FALSE;
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
    event = list();
    event$mat = mutated_net;
    event$time = iteration_num;
    event$mean_time = characteristic_vals[[1]]
    new_state$extra$history = append(new_state$extra$history,list(event));
    new_net = mutated_net;
  }
  return(list(new_state = new_state,
              new_net = new_net,
              continue = cont));
};

mutation_fixed_links = function(base,mat){
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
                                mutation_func=mutation_fixed_links,
                                evolution_func = evolution_simple,
                                characteristic_funcs = list(only_mean_time))

system.time(ens_evolve$evolve_ensamble_par())

results = ens_evolve$get_results()

results_list=lapply(results,function(x){x$state$extra$history})

por_mat=lapply(results_list,function(x){lapply(x,function(h){list(mean_time=h$mean_time,
                                                                  time=h$time,
                                                                  num_nodes=ncol(h$mat)-2)})})

vis = Visualization$new(in_data=por_mat)

vis$plot_me(x_var = 'time',y_var = 'mean_time',col_var = 'num_nodes')
