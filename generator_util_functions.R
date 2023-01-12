# BASE CONNECT FUNCTION GENERATOR WITH RENYI OPTION

connection_generator_with_renyi_opt = function(do_renyi=TRUE,renyi_connectance=0.2){
  output = function(N,autoloop,respiration){
    n_input=sample(1:N,1);
    n_output=sample(1:N,1);
    non_core_nodes = if(respiration) 3 else 2;
    mat = matrix(0, nrow = N+non_core_nodes, ncol = N+non_core_nodes);
    mat[N+1,1:n_input]=1;
    mat[(N-n_output):N,N+2]=1;
    if(respiration){
      mat[1:N,N+3]=1;
    }
    core = mat[1:N,1:N];
    if(do_renyi){
      #THIS IS EXTRA APPART FROM THE DEFAULT FUNCTION 
      core = as_adjacency_matrix(erdos.renyi.game(N, renyi_connectance, directed=TRUE),sparse=FALSE);
      #############################################
    }
    mat[1:N,1:N]=core;
    return(mat);
  }
}


# ANY LINK MUTATION FUNCTION WITH LOSS CHANCE
# CONSTANTS : loss_threshold (Probability of loosing a link), 
#             gain_threshold (probability of gaining a link). 
#             These should add up to a value less than 1.0. 1.0-loss_threshold-gain_threshold is the probability of modifying an existing link

# DESCRIPTION : Mutates link values randomly, selects links from any link in core depending on the random value obtained using the probabilities described before
# OBSERVATIONS : 

mutation_generator_base = function(loss_threshold=0.2,gain_threshold=0.2,gain_min_val=0.0,gain_max_val=0.9999,perturbative_change=FALSE,change_min_or_perturbation=0.0,change_max=0.99999){
  output = function(base,mat){
    diag_pos=which(row(mat[1:base$core_nodes,1:base$core_nodes])==col(mat[1:base$core_nodes,1:base$core_nodes]))
    non_zero_links = which((mat[1:base$core_nodes,1:base$core_nodes])>0);
    zero_links = which((mat[1:base$core_nodes,1:base$core_nodes])==0);
    non_zero=non_zero_links[-diag_pos]
    zero=zero_links[-diag_pos]
    type_of_mutation = runif(1,min=0.0,max=1.0)
    new_value = 0.0;
    if(type_of_mutation<loss_threshold){
      link_to_change = sample(non_zero,1);
    }else if(type_of_mutation<(gain_threshold+loss_threshold)){
      link_to_change = sample(zero,1);
      new_value=runif(1,min=gain_min_val,max=gain_max_val);
    }else{
      link_to_change = sample(non_zero,1);
      if(perturbative_change){
        from = link_to_change%%(base$core_nodes);
        if(from==0)from = base$core_nodes; 
        to = ceiling(link_to_change/(base$core_nodes));
        min_val = mat[[from,to]]-change_min_or_perturbation;
        max_val = mat[[from,to]]+change_min_or_perturbation;
        if(min_val<0)min_val=0.0;
        if(max_val>1)max_val=0.9999;
        new_value = runif(1,min=min_val,max=max_val);
      }else{
        new_value = runif(1,min=change_min_or_perturbation,max=change_max);
      }
    }
    from = link_to_change%%(base$core_nodes);
    if(from==0)from = base$core_nodes; 
    to = ceiling(link_to_change/(base$core_nodes));
    
    mat_out = mat;
    mat_out[[from,to]]=new_value;
    mat_out = base$turn_to_flowing_network(mat=mat_out,allow_autoloops=FALSE,set_adjacency=FALSE,new_val=0.001);
    mat_out = base$normalize_rows(mat=mat_out,set_adjacency=FALSE);
    return(mat_out);
  }
  return(output)
}

tst = function(base,mat){
  return(base$get_TST(mat=mat,input_multiplier=1000));
};

storage_per_node=function(base,mat){
    return(base$get_initial_stock_after_stabilization(mat=mat,initial_input = 1000,frequency = 1)[1:base$core_nodes])
  };

storage=function(base,mat){
    return(sum(base$get_initial_stock_after_stabilization(mat=mat,initial_input = 1000,frequency = 1)[1:base$core_nodes]))
};

c_ent = function(base,mat){
  return(base$get_link_type_entropy(mat=mat));
};

ami=function(base,mat){
  return(base$get_ami(mat=mat))
};

b=function(base,mat){
  return(base$get_exp_factor(mat=mat))
}

input_output_entropy_generator=function(log_base=2){
  return(function(base,mat){
    return(base$get_input_output_entropy(mat=mat,input_multiplier=1000,logbase=log_base))
  });
}

entropy_diff_generator=function(log_base=2){
  return (function(base,mat){
    comparison=base$get_input_output_entropy(mat=mat,input_multiplier=1000L,logbase=log_base);
    return(comparison$out_entropy-comparison$in_entropy)
  })
}

ascendency=function(base,mat){
  return(base$get_ascendency(mat=mat))
};

residence_time=function(base,mat){
  return(sum(base$get_initial_stock_after_stabilization(mat=mat,initial_input = 1000,frequency = 1)[1:base$core_nodes])/base$get_TST(input_multiplier=1000))
};

mean_time_generator = function(input_percentage=0.5){
  output=function(base,mat){
    return(base$mean_time_discrete(mat=mat,input_percentage=input_percentage));
  }
  return(output)
};

finn=function(base,mat){
  return(base$get_finn(mat=mat));
}

grow_part_of_stock=function(base,mat){
  stock=base$get_initial_stock_after_stabilization(mat=mat,initial_input=1000.0,frequency=1L)[1:base$core_nodes];
  return(sum(stock[base$net_data$picked_nodes]))
}

grow_proportion_of_stock=function(base,mat){
  stock=base$get_initial_stock_after_stabilization(mat=mat,initial_input=1000.0,frequency=1L)[1:base$core_nodes];
  return(sum(stock[base$net_data$picked_nodes])/sum(stock[1:base$core_nodes]))
}

grow_part_finn_node=function(base,mat){
  finn_node=base$get_finn_per_node(mat=mat,input_multiplier=1000)
  return(sum(finn_node[base$net_data$picked_nodes]))
}

grow_proportion_finn_node=function(base,mat){
  finn_node=base$get_finn_per_node(mat=mat,input_multiplier=1000)
  return(sum(finn_node[base$net_data$picked_nodes]/sum(finn_node)))
}

evolution_function_generator = function(iteration_max=50,characteristics_boolean='and',fixed_iters=10,event_saving_function=function(base,mat,iteration_num,characteristic_vals){return(list())}){
  output=function(iteration_num,base,current_state,current_net,mutated_net,characteristic_vals){
    cont = FALSE
    if(iteration_num<iteration_max) 
      cont = TRUE
    new_state = current_state;
    new_net = current_net;
    saveit = FALSE
    if(characteristics_boolean=='and'){
      saveit = TRUE;
      # "AND" OF ALL CHARACTERISTICS
      for(i in 1:length(characteristic_vals)){
        #natural selection rule: condition for saving
        if(characteristic_vals[[i]]<current_state$vals[[i]]){
          saveit = FALSE;
          break;
        }
      }
    }else if(characteristics_boolean=='or'){
      # "OR" OF ALL CHARACTERISTICS
      for(i in 1:length(characteristic_vals)){
        #natural selection rule: condition for saving
        if(characteristic_vals[[i]]>current_state$vals[[i]]){
          saveit = TRUE;
          break;
        }
      }
    }else if(characteristics_boolean=='fixed_iters'){
      saveit = (iteration_num%%fixed_iters) == 0;
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
      event=event_saving_function(base,mutated_net,iteration_num);
      new_state$extra$history = append(new_state$extra$history,list(event));
      new_net = mutated_net;
    }
    
    return(list(new_state = new_state,
                new_net = new_net,
                continue = cont));
  };
  return(output);
}

opening=function(nodes,replicates,interval,orientor,greed_experiment,picked_percent){
  if(length(nodes)==1){
    open_name=paste0(nodes,'_',nodes,'_',replicates,'_',interval,'/',
                     'Ensemble_',nodes,'_',nodes,'_',replicates,'_',interval,
                     '_Orientor_',orientor,
                     '_GreedExp_',greed_experiment,
                     '_PickedPerc_',picked_percent,
                     '.RDS')
    
  }
  else{
    open_name=paste0(head(nodes,1),'_',tail(nodes,1),'_',replicates,'_',interval,
                     '_Ensemble_',nodes,'_',nodes,'_',replicates,'_',interval,
                     '_Orientor_',orientor,
                     '_GreedExp_',greed_experiment,
                     'PickedPerc_',picked_percent,
                     '.RDS')
  }
  return(open_name)
}