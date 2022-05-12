source("FlowNet.R")

evolution_simple=function(iteration_num,base,current_state,current_net,mutated_net,characteristic_vals){
  cont = TRUE
  if(iteration_num>50) 
    cont = FALSE
  new_state = current_state;
  new_net = current_net;
  if(iteration_num==1 && is.null(new_state$ent)){
    new_state$ent = base$get_link_type_entropy()$link_type_entropy$average_ent_per_link;
  }
  saveit = FALSE
  cycle_num=0;
  trajectory_num=0;
  ent_av=0;
  ent_tot=0;
  trajectory_num=0;
  if(characteristic_vals[[1]]>current_state$vals[[1]]
  ){
      ent = base$get_link_type_entropy(mat=mutated_net);
      saveit = TRUE;
      ent_av = ent$link_type_entropy$average_ent_per_link;
      ent_tot = ent$link_type_entropy$total_link_ent;
      cycle_num = length(ent$netdata$cycles);
      trajectory_num = length(ent$netdata$trajectories);
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
    event$vals = characteristic_vals;
    event$cycle_num = cycle_num;
    event$trajectory_num = trajectory_num;
    event$ent_tot = ent_tot;
    event$ent_av = ent_av;
    event$time = iteration_num;
    new_state$extra$history = append(new_state$extra$history,list(event));
    new_net = mutated_net;
  }
  return(list(new_state = new_state,
              new_net = new_net,
              continue = cont));
};


evolution_double=function(iteration_num,base,current_state,current_net,mutated_net,characteristic_vals){
  cont = TRUE
  if(iteration_num>50) 
    cont = FALSE
  new_state = current_state;
  new_net = current_net;
  if(iteration_num==1 && is.null(new_state$ent)){
    new_state$ent = base$get_link_type_entropy()$link_type_entropy$average_ent_per_link;
  }
  saveit = FALSE
  cycle_num=0;
  trajectory_num=0;
  ent_av=0;
  ent_tot=0;
  trajectory_num=0;
    if(characteristic_vals[[1]]>current_state$vals[[1]]
       ){
      ent = base$get_link_type_entropy(mat=mutated_net);
      if(ent$link_type_entropy$average_ent_per_link>new_state$ent){
        saveit = TRUE;
        ent_av = ent$link_type_entropy$average_ent_per_link;
        ent_tot = ent$link_type_entropy$total_link_ent;
        cycle_num = length(ent$netdata$cycles);
        trajectory_num = length(ent$netdata$trajectories);
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
    event$vals = characteristic_vals;
    event$cycle_num = cycle_num;
    event$trajectory_num = trajectory_num;
    event$ent_tot = ent_tot;
    event$ent_av = ent_av;
    event$time = iteration_num;
    new_state$extra$history = append(new_state$extra$history,list(event));
    new_net = mutated_net;
  }
  return(list(new_state = new_state,
              new_net = new_net,
              continue = cont));
};

mutation_simple = function(base,mat){
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
  return(base$get_link_type_entropy(mat=mat)$link_type_entropy$average_ent_per_link);
};

ami=function(base,mat){
  return(base$get_ami(mat=mat))
}

ensamble = NetEnsamble$new(ensemble_sizeMin = 5L,
                           ensemble_sizeMax = 7L,
                           size_replicate = 10L,
                           nodes_autoloops_allowed = FALSE);
ensamble$generate()

ens_evolve = EnsambleEvolve$new(ensamble=ensamble,
                                mutation_func=mutation_simple,
                                evolution_func = evolution_simple,
                                characteristic_funcs = list(ami))

system.time(ens_evolve$evolve_ensamble())
system.time(ens_evolve$evolve_ensamble_par())

results = ens_evolve$get_results()

results_list=lapply(results,function(x){x$state$extra$history})
por_mat=lapply(results_list,function(x){lapply(x,function(h){list(ami=h$vals[[1]],
                                                                  #c_ent=h$vals[[2]],
                                                                  time=h$time,
                                                                  trajectory_num=h$trajectory_num,
                                                                  cycles_num=h$cycle_num,
                                                                  ent_tot=h$ent_tot ,
                                                                  ent_av=h$ent_av,
                                                                  num_nodes=ncol(h$mat)-2)})})
                                                                
dfs=lapply(por_mat,function(x){return(data.frame(t(matrix(unlist(x),nrow = 7, dimnames=list(c("ami","time","trajectory_num","cycles_num","entropy_tot","entropy_av","num_nodes"),vector(mode="expression"))))))})

get_min_max=function(dfs,variable){
  maximum=vector(mode="numeric",length=length(dfs))
  minimum=vector(mode="numeric",length=length(dfs))
  for(i in 1:length(dfs)){
    maximum[i]=max(dfs[[i]][[variable]])
    minimum[i]=min(dfs[[i]][[variable]])
  }
  y_max=max(maximum)
  y_min=min(minimum)
  output=list();
  output$max=y_max;
  output$min=y_min;
  return(output);
}

plot_x = 'ami'
plot_color = 'num_nodes'
plot_variable='cycles_num'
plot_size='trajectory_num'

min_max_col = get_min_max(dfs,plot_color);
min_max_size = get_min_max(dfs,plot_size);
min_max_x = get_min_max(dfs,plot_x);
color_func=function(x){return(hsv(h = (x-min_max_col$min)/(min_max_col$max-min_max_col$min)*0.8, s = 1, v = 1, 1))}
size_func = function(x){return(0.5+(x-min_max_size$min)/(min_max_size$max-min_max_size$min))}
min_max_var = get_min_max(dfs,plot_variable);

col = unlist(lapply(dfs[[1]][[plot_color]],color_func));
size = unlist(lapply(dfs[[1]][[plot_size]],size_func));
plot(dfs[[1]][[plot_x]],dfs[[1]][[plot_variable]], type = "b",col = col,cex=size, xlab = "x", ylab = "y",xlim=c(min_max_x$min,min_max_x$max),ylim=c(min_max_var$min,min_max_var$max))
for(i in 2:length(dfs)){
  col = unlist(lapply(dfs[[i]][[plot_color]],color_func));
  size = unlist(lapply(dfs[[i]][[plot_size]],size_func));
  lines(dfs[[i]][[plot_x]],dfs[[i]][[plot_variable]], type = "b",col = col,cex=size, xlab = "x", ylab = "y")
}

length(results[[1]]$final$navigate_net(save_trajectories = TRUE)$cycles)
