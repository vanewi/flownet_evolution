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

ensamble = NetEnsamble$new(ensemble_sizeMin = 30L,
                           ensemble_sizeMax = 30L,
                           size_replicate = 300L,
                           nodes_autoloops_allowed = FALSE);

ensamble$generate(nodes_connect_function = connector_f,minimal_flowing = TRUE)

ens_evolve = EnsambleEvolve$new(ensamble=ensamble,
                                mutation_func=mutation_f,
                                evolution_func = evolution_f,
                                characteristic_funcs = list(ascendency))

system.time(ens_evolve$evolve_ensamble_par())

results = ens_evolve$get_results()

# FORCED GREEDY EVOLUTION

new_ensamble = ens_evolve$get_final_ensamble();

evolution_greedy = evolution_function_generator(iteration_max=100,
                                           characteristics_boolean='and',
                                           event_saving_function=event_save_function
);

mutation_perturbative = mutation_generator_base(loss_threshold = 0.0,gain_threshold = 0.0,gain_min_val = 0.0,gain_max_val = 0.9999,perturbative_change = TRUE,change_min_or_perturbation = 0.1,change_max = 0.9999);

# FORCED SECTOR PICK PER NET
pick_percent = 1.0
for (i in 1:length(new_ensamble$generated_networks)){
  num_nodes = new_ensamble$generated_networks[[i]]$core_nodes
  to_pick = floor(runif(1,max=pick_percent) * num_nodes)
  if(to_pick<1)to_pick=1
  picked = sample(seq(1,num_nodes,1),to_pick,replace = FALSE)
  new_ensamble$generated_networks[[i]]$net_data$picked_nodes=picked
}

grow_part_of_stock=function(base,mat){
  stock=base$get_initial_stock_after_stabilization(mat=mat,initial_input=1000.0,frequency=1L)[1:base$core_nodes];
  return(sum(stock[base$net_data$picked_nodes]))
}

ens_evolve_greedy = EnsambleEvolve$new(ensamble=new_ensamble,
                                 mutation_func=mutation_perturbative,
                                 evolution_func = evolution_greedy,
                                 characteristic_funcs = list(grow_part_of_stock))

ens_evolve_control = EnsambleEvolve$new(ensamble=new_ensamble,
                                       mutation_func=mutation_perturbative,
                                       evolution_func = evolution_greedy,
                                       characteristic_funcs = list(ascendency))

system.time(ens_evolve_greedy$evolve_ensamble_par())

system.time(ens_evolve_control$evolve_ensamble_par())

results_greedy = ens_evolve_greedy$get_results()
results_control = ens_evolve_control$get_results()


# SKETCHY RESULT ANALYSIS

dExpF = vector(mode='list',length=length(results_control))
dPStock = vector(mode='list',length=length(results_control))
dPFinn = vector(mode='list',length=length(results_control))
pStockI = vector(mode='list',length=length(results_control))
pStockF = vector(mode='list',length=length(results_control))
pickedNum = vector(mode='list',length=length(results_control))
for(n in 1:length(results_control)){
  nodes = results_greedy[[n]]$initial$core_nodes
  picked = results_greedy[[n]]$initial$net_data$picked_nodes
  pickedNum[[n]]=length(picked)
  b_g = results_greedy[[n]]$final$get_exp_factor()
  b_c = results_control[[n]]$final$get_exp_factor()
  dExpF[[n]] = (b_c-b_g)/abs(b_c)*100
  
  stockI = results_greedy[[n]]$initial$get_initial_stock_after_stabilization(initial_input = 1000.0,frequency = 1L)[1:nodes]
  stockF = results_greedy[[n]]$final$get_initial_stock_after_stabilization(initial_input = 1000.0,frequency = 1L)[1:nodes]
  
  finnI = results_greedy[[n]]$initial$get_finn_per_node(input_multiplier = 1000.0);
  finnF = results_greedy[[n]]$final$get_finn_per_node(input_multiplier = 1000.0);
  
  dPFinn[[n]] = sum(finnF[picked])/sum(finnF)*100-sum(finnI[picked])/sum(finnI)*100
  
  pStockI[[n]] = sum(stockI[picked])/sum(stockI)*100
  pStockF[[n]] = sum(stockF[picked])/sum(stockF)*100
  dPStock[[n]]=pStockF[[n]]-pStockI[[n]]
}
plot(unlist(dExpF))
plot(unlist(pickedNum))
abline(a=0,b=0,col="red")

plot(unlist(pickedNum),unlist(dExpF))
abline(a=0,b=0,col="red")

df=data.frame(cbind(unlist(dPStock),unlist(dExpF),unlist(dPFinn),unlist(pickedNum)))
names(df)=c('dPStock','dExpF','dPFinn','pickedNum')

p <- ggplot(df, aes(pickedNum, dExpF, colour = dPStock)) 
p <- p + geom_point()
p


abline(a=0,b=0,col="red")

GoodGReedN=which(dExpF<0)
BadGreedN=which(dExpF>=0)
length(GoodGReedN)
length(BadGreedN)

# BETTER FOR GREEDY
cual_bien = 10
elected=results_greedy[[GoodGReedN[[cual_bien]]]]$initial$net_data$picked_nodes
results_greedy[[GoodGReedN[[cual_bien]]]]$initial$get_finn_per_node()[elected]
summary(results_greedy[[GoodGReedN[[cual_bien]]]]$initial$get_finn_per_node())
results_greedy[[GoodGReedN[[cual_bien]]]]$initial$get_finn()
results_greedy[[GoodGReedN[[cual_bien]]]]$initial$get_finn_per_node()[elected]
summary(results_greedy[[GoodGReedN[[cual_bien]]]]$initial$get_finn_per_node())
plot(unlist(map(results_greedy[[GoodGReedN[[cual_bien]]]]$state$extra$history,'b')))
results_greedy[[GoodGReedN[[cual_bien]]]]$initial$get_initial_stock_after_stabilization(initial_input = 1000,frequency = 1)[elected]
summary(as.vector(results_greedy[[GoodGReedN[[cual_bien]]]]$initial$get_initial_stock_after_stabilization(initial_input = 1000,frequency = 1)))


# WORSE FOR GREEDY
cual_mal = 2
results_greedy[[BadGreedN[[cual_mal]]]]$initial$net_data$picked_nodes
results_greedy[[BadGreedN[[cual_mal]]]]$initial$get_finn_per_node()
results_greedy[[BadGreedN[[cual_mal]]]]$initial$get_finn()
summary(results_greedy[[BadGreedN[[cual_mal]]]]$initial$get_finn_per_node())
plot(unlist(map(results_greedy[[BadGreedN[[cual_mal]]]]$state$extra$history,'b')))
results_greedy[[BadGreedN[[cual_mal]]]]$initial$get_initial_stock_after_stabilization(initial_input = 1000,frequency = 1)

length(GoodGReedN)

plot(unlist(map(results_greedy[[6]]$state$extra$history,'b')))
plot(unlist(map(results_control[[6]]$state$extra$history,'b')))


