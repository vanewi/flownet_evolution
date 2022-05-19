source("FlowNet.R")

evolution_c_ent=function(iteration_num,base,current_state,current_net,mutated_net,characteristic_vals){
  cont = TRUE
  if(iteration_num>50) 
    cont = FALSE
  new_state = current_state;
  new_net = current_net;

  current_ent = characteristic_vals[[1]]$link_type_entropy$average_ent_per_link;
  pre_ent = current_state$vals[[1]]$link_type_entropy$average_ent_per_link;
  
  if(current_ent>pre_ent){
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
    event$tst = characteristic_vals[[2]];
    event$mean_time = characteristic_vals[[3]];
    event$ami = characteristic_vals[[4]];
    event$cycle_num = length(characteristic_vals[[1]]$netdata$cycles);
    event$trajectory_num = length(characteristic_vals[[1]]$netdata$trajectories);
    event$ent_tot = characteristic_vals[[1]]$link_type_entropy$total_link_ent;
    event$ent_av = characteristic_vals[[1]]$link_type_entropy$average_ent_per_link;
    event$time = iteration_num;
    new_state$extra$history = append(new_state$extra$history,list(event));
    new_net = mutated_net;
  }

  return(list(new_state = new_state,
              new_net = new_net,
              continue = cont));
};

mutation_any_link = function(base,mat){
  diag_pos=which(row(mat[1:base$core_nodes,1:base$core_nodes])==col(mat[1:base$core_nodes,1:base$core_nodes]))
  
  mod_links = which((mat[1:base$core_nodes,1:base$core_nodes])>=0);
  modlinks=mod_links[-diag_pos]
  if(length(mod_links)==0)
    return(mat);
  link_to_change = sample(mod_links,1);
  from = link_to_change%%(base$core_nodes);
  if(from==0)from = base$core_nodes; 
  to = ceiling(link_to_change/(base$core_nodes));
  mat_out = mat;
  loss_chance = runif(1,min=0.0,max=1.0);
  if(mat[[from,to]]>0 && loss_chance>0.8)
    mat_out[[from,to]]=0.0
  else
    mat_out[[from,to]]=runif(1,min=0.00001,max=0.99999); # CONSIDERAR EN UN RANGO EN TORNO AL VALOR ANTERIOR TAMBIEN, NO SOLO FULL RANDOM
  mat_out=base$normalize_rows(mat=mat_out,set_adjacency=FALSE);
  return(mat_out);
};

tst = function(base,mat){
  return(base$get_TST(mat=mat));
};

c_ent = function(base,mat){
  return(base$get_link_type_entropy(mat=mat));
};

ami=function(base,mat){
  return(base$get_ami(mat=mat))
}

mean_time = function(base,mat){
  return(base$mean_time_discrete(input_percentage=0.05,mat=mat));
};

ensamble = NetEnsamble$new(ensemble_sizeMin = 5L,
                           ensemble_sizeMax = 13L,
                           size_replicate = 2L,
                           nodes_autoloops_allowed = FALSE);
ensamble$generate()

ens_evolve = EnsambleEvolve$new(ensamble=ensamble,
                                mutation_func=mutation_any_link,
                                evolution_func = evolution_c_ent,
                                characteristic_funcs = list(c_ent,tst,mean_time,ami))

system.time(ens_evolve$evolve_ensamble_par())

results = ens_evolve$get_results()

results_list=lapply(results,function(x){x$state$extra$history})

por_mat=lapply(results_list,function(x){lapply(x,function(h){list(tst=h$tst,
                                                                  mean_time=h$mean_time,
                                                                  ami=h$ami,
                                                                  asc = h$ami*h$tst,
                                                                  time=h$time,
                                                                  trajectory_num=h$trajectory_num,
                                                                  cycles_num=h$cycle_num,
                                                                  ent_tot=h$ent_tot ,
                                                                  ent_av=h$ent_av,
                                                                  num_nodes=ncol(h$mat)-2)})})

vis = Visualization$new(in_data=por_mat)

vis$plot_me(x_var = 'time',
            y_var = 'ent_av',
            col_var = 'num_nodes',
            size_var = 'cycles_num',
            y_lim = NULL,#list(min=0,max=25000),
            x_lim = NULL#list(min=0,max=50)
            )
