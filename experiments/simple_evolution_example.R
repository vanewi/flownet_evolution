source("FlowNet.R")
library(purrr)
# AS THE EVENT SAVING CALCULATES CYCLES, DO NOT USE LARGER THAN 20 NODE NETS

event_saving_function=function(base,mat,iteration_num,characteristic_vals){
  event = list();
  ent = c_ent(base,mat);
  c_analysis=base$get_cycles_flow(mat=mat,initial_input=1000)$cycles_flow;
  event$c_flow=c_analysis;
  event$mat = mat;
  event$tst = tst(base,mat);
  event$mean_time = mean_time(base,mat);
  event$ami = ami(base,mat);
  event$b=base$get_exp_factor(mat=mat);
  event$eigMax=base$get_eigen_max(mat=mat);
  event$cycle_num = length(ent$netdata$cycles);
  event$trajectory_num = length(ent$netdata$trajectories);
  event$ent_tot = ent$link_type_entropy$total_link_ent;
  event$ent_av = ent$link_type_entropy$average_ent_per_link;
  event$finn=finn(base,mat);
  event$time = iteration_num;
  event$ascendency=event$tst*event$ami;
  event$eigSD = sd(base$get_eigen_mod(mat=mat));
  return(event);
}

mean_time = mean_time_generator(input_percentage = 0.05);

evolution_f = evolution_function_generator(iteration_max=30,
                                           characteristics_boolean='and',
                                           event_saving_function = event_saving_function
  );

mutation_f = mutation_generator_base(loss_threshold = 0.2,gain_threshold = 0.2,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = FALSE,change_min_or_perturbation = 0.0,change_max = 0.9999);

connector_f = connection_generator_with_renyi_opt(do_renyi = TRUE,renyi_connectance = 0.2);

ensamble = NetEnsamble$new(ensemble_sizeMin = 5L,
                           ensemble_sizeMax = 8L,
                           size_replicate = 4L,
                           nodes_autoloops_allowed = FALSE);

ensamble$generate(nodes_connect_function = connector_f,minimal_flowing = TRUE)

ens_evolve = EnsambleEvolve$new(ensamble=ensamble,
                                mutation_func=mutation_f,
                                evolution_func = evolution_f,
                                characteristic_funcs = list(tst))

system.time(ens_evolve$evolve_ensamble_par())

results = ens_evolve$get_results()

# TO SAVE AND LOAD RESULTS
#saveRDS(results, file = "resPrueba.RDS")
#results <- readRDS("resPrueba.RDS")

results_list=lapply(results,function(x){x$state$extra$history})

por_mat=lapply(results_list,function(x){lapply(x,function(h){list(ltst=log(h$tst),
                                                                 mean_time=h$mean_time,
                                                                  ami=h$ami,
                                                                  finn=h$finn,
                                                                  eigMax=h$eigMax,
                                                                  b=h$b,
                                                                  eEM=exp(h$eigMax),
                                                                 lnFinn=log(h$finn),
                                                                  #asc = h$ascendency
                                                                  c_flowln=log(sum(unlist(map(h$c_flow,"cycle_flow")))),
                                                                  c_flow=sum(unlist(map(h$c_flow,"cycle_flow"))),
                                                                  c_flow_max=max(unlist(map(h$c_flow,"cycle_flow"))),
                                                                  c_flow_max_n=max(unlist(map(h$c_flow,"cycle_flow")))/sum(unlist(map(h$c_flow,"cycle_flow"))),
                                                                  c_flow_sd=sd(unlist(map(h$c_flow,"cycle_flow"))),
                                                                  c_inhom=mean(unlist(map(h$c_flow,"cycle_inhom"))),
                                                                  c_length=mean(unlist(map(h$c_flow,"cycle_length"))),
                                                                  c_length_sd=sd(unlist(map(h$c_flow,"cycle_length"))),
                                                                 asc=h$tst*h$ami,
                                                                 eAsc=exp(h$tst*h$ami),
                                                                  time=h$time,
                                                                 tst=h$tst,
                                                                 eigSD=h$eigSD,
                                                                  trajectory_num=log10(h$trajectory_num),
                                                                  cycles_num=log10(h$cycle_num),
                                                                 prop_vals=h$cycle_num/h$trajectory_num,
                                                                 prop_vals_i=h$trajectory_num/h$cycle_num,
                                                                  ent_tot=h$ent_tot ,
                                                                  ent_av=h$ent_av,
                                                                  num_nodes=ncol(h$mat)-2,
                                                                  links_prop=length(which(h$mat[1:(ncol(h$mat)-2),1:(ncol(h$mat)-2)]>0))/((ncol(h$mat)-2)*(ncol(h$mat)-3)),
                                                                  finn=h$finn,
                                                                  num_links=length(which(h$mat[1:(ncol(h$mat)-2),1:(ncol(h$mat)-2)]>0)))})})

for(i in 1:length(por_mat))for(j in 1:length(por_mat[[i]]))if(is.na(por_mat[[i]][[j]]$c_flow_sd))por_mat[[i]][[j]]$c_flow_sd=0
for(i in 1:length(por_mat))for(j in 1:length(por_mat[[i]]))if(is.na(por_mat[[i]][[j]]$c_length))por_mat[[i]][[j]]$c_length=0
for(i in 1:length(por_mat))for(j in 1:length(por_mat[[i]]))if(is.na(por_mat[[i]][[j]]$c_length))por_mat[[i]][[j]]$c_length=0

# FOR SUBSETS
#vis = Visualization$new(in_data=por_mat[lapply(por_mat,function(x) if(length(x)==0) 0 else x[[1]]$num_nodes)==10])

vis = Visualization$new(in_data=por_mat)


vis$plot_me(x_var = 'eigSD',
            y_var = 'ent_av',
            col_var = 'time',
            size_var = 'cycles_num',
            y_lim = NULL,#list(min=1,max=10), 
            x_lim = NULL#list(min=0,max=2000)
)

vis$plot_me(x_var = 'b',
            y_var = 'finn',
            col_var = 'time',
            size_var = 'num_nodes',
            y_lim = NULL,#list(min=0,max=25000),
            x_lim = NULL#list(min=0,max=50)
)

  
  
