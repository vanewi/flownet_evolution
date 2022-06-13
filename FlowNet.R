library(igraph)
library(sets)
library(Rcpp)
library(ggplot2)
library(reshape2)
library(parallel)
sourceCpp("flownet.cpp")

NetBase <- setRefClass("NetBase",
											 fields = list(core_nodes='numeric',
											               #nodes_connect_function='function',
											               adjacency_matrix='matrix',
											               allow_autoloop='logical',
											               respiration='logical',
											               total_nodes='numeric',
											               net_data='list'),
											 methods = list(
											 	initialize=function(core_nodes=5L,
											 											nodes_connect_function=function(N,autoloop,respiration){
												 												n_input=sample(1:N,1);
												 												n_output=sample(1:N,1);
												 												non_core_nodes = if(respiration) 3 else 2;
												 												mat = matrix(0, nrow = N+non_core_nodes, ncol = N+non_core_nodes);
												 												input_ready = vector(mode="list",length = N);
												 												output_ready = vector(mode="list",length = N);
												 												mat[N+1,1:n_input]=1;
												 												input_ready[1:n_input]=TRUE;
												 												mat[(N-n_output):N,N+2]=1;
												 												output_ready[(N-n_output):N]=TRUE;
												 												if(respiration){
												 													mat[1:N,N+3]=1;
												 												}
												 												return(mat);
											 												},
											 											allow_autoloop=FALSE,
											 											respiration = FALSE,
											 											make_flowing = TRUE,
											 											random_weights = TRUE,
											 											normalized_rows = TRUE,
											 											fully_connected = FALSE,
											 											...){
											 		core_nodes<<-core_nodes;
											 		#nodes_connect_function <<- nodes_connect_function;
											 		allow_autoloop <<- allow_autoloop;
											 		respiration <<- respiration;
											 		adjacency_matrix <<- nodes_connect_function(core_nodes,autoloop,respiration);
											 		total_nodes <<- nrow(adjacency_matrix);
											 		net_data<<- list();
											 		if(fully_connected){
											 		  make_fully_connected();
											 		}
											 		if(make_flowing){
											 		  turn_to_flowing_network();
											 		}
											 		if(random_weights){
											 		  set_random_weights();
											 		}
											 		if(normalized_rows){
											 		  normalize_rows();
											 		}
											 		callSuper(...);
											 	},
											 	create_from_base_and_adjacency=function(base,adjacency){
											 	  core_nodes<<- base$core_nodes;
											 	  #nodes_connect_function <<- base$nodes_connect_function;
											 	  allow_autoloop <<- base$allow_autoloop;
											 	  respiration <<- base$respiration;
											 	  adjacency_matrix <<- adjacency
											 	  total_nodes <<- nrow(adjacency_matrix);
											 	  net_data<<- list();
											 	},
											 	make_fully_connected=function(){
											 	  mat = matrix(1,nrow=core_nodes+2,ncol=core_nodes+2);
											 	  diag(mat)=0;
											 	  mat[core_nodes+2,]=0;
											 	  mat[,core_nodes+1]=0;
											 	  mat[core_nodes+1,core_nodes+2]=0;
											 	  adjacency_matrix <<- mat;
											 	},
											 	get_max_number_cycles_core=function(min_nodes=2,max_nodes=core_nodes){
											 	  return(sum(unlist(lapply(min_nodes:max_nodes,function(x)choose(core_nodes,x)*factorial(x-1)))));
											 	  },
											 	get_max_number_flowing_trajectories=function(min_nodes=0,max_nodes=core_nodes){
											 	  return(sum(unlist(lapply(min_nodes:max_nodes,function(x)choose(core_nodes,x)*factorial(x)))));
											 	},
											 	get_ami=function(mat=adjacency_matrix,input_multiplier=1000,consider_input=TRUE,consider_output=TRUE,log_base=2){
											 	  final_node_in = core_nodes;
											 	  final_node_out = core_nodes;
											 	  if(consider_input)
											 	    final_node_in = core_nodes + 1
											 	  if(consider_output)
											 	    final_node_out = core_nodes + 2
											 	  m = get_quantity_matrix(mat=mat,multiplier=input_multiplier)[1:final_node_in,1:final_node_out];
											 	  tst=get_TST(mat = mat,input_multiplier=input_multiplier,consider_input = consider_input,consider_output = consider_output)
											 	  ami=0
											 	  for (i in 1:final_node_in){
											 	    for (j in 1:final_node_out){
											 	      if(m[i,j]==0) next
											 	      val = (m[i,j]*tst)/(sum(m[i,])*sum(m[,j]));
											 	      ami=ami+(m[i,j]/tst)*(log2(val)/log2(log_base));
									            }
											 	  }
											 	  return(ami)
											 	},
											 	mean_time_discrete=function(input_percentage=0.5,mat=adjacency_matrix){
											 		X=rep(0,ncol(mat));
											 		X[[core_nodes+1]]=1.0;
											 		innet = 1.0
											 		iters = 0
											 		while(innet > input_percentage){
											 				X=(X)%*%mat
											 				iters = iters+1
											 				innet = sum(X[1:core_nodes])
											 		}
											 		return(iters) 
											 	},
											 	net_plot=function(mat=adjacency_matrix,margin=0.0){
											 		#plot.network(as.network(adjacency_matrix),interactive=TRUE,label = seq(1:total_nodes));
											 		g = graph_from_adjacency_matrix(mat,weighted=TRUE)
											 	
											 		names=as.character(c(seq(1,core_nodes,1),"IN","OUT"))
											 		color_pallete=c("orange","green",'red')
											 		colors=as.factor(c(rep(1,core_nodes),2,3))
											 		
											 		set_vertex_attr(g,"names",value=names)
											 		set_vertex_attr(g,"colors",value=colors)
				
											 		plot(g,edge.arrow.size=0.1,edge.label=round(E(g)$weight, 3),edge.label.cex=0.7,vertex.size=20,vertex.color=color_pallete[colors],vertex.label=names,margin=margin)#, edge.color=edge.col);
											 		},
											 	iterate_through_mat_inout=function(core,nodes_ready,node,in_or_out){
											 		if(!is.null(nodes_ready[[node]])){
											 			return(nodes_ready);
											 		}
											 		nodes_ready_temp = nodes_ready;
											 		nodes_ready_temp[[node]]=TRUE;
											 		connected = list();
											 		if(in_or_out==TRUE){
											 			connected=which(core[node,]==1);
											 		}else{
											 			connected=which(core[,node]==1);
											 		}
											 		if(length(connected)==0)
											 			return(nodes_ready_temp);
											 		for(i in 1:length(connected)){
											 			nodes_ready_temp = iterate_through_mat_inout(core,nodes_ready_temp,connected[[i]],in_or_out);
											 		}
											 		return(nodes_ready_temp);
											 	},
											 	get_connected_from=function(mat=adjacency_matrix,node){
											 		connected = vector(mode="list",length = nrow(mat));
											 		return(iterate_through_mat_inout(mat,connected,node,TRUE));
											 	},
											 	get_connected_to=function(mat=adjacency_matrix,node){
											 		connected = vector(mode="list",length = nrow(mat));
											 		return(iterate_through_mat_inout(mat,connected,node,FALSE));
											 	},
											 	turn_to_flowing_network=function(mat=adjacency_matrix,allow_autoloops=allow_autoloop,set_adjacency=TRUE,new_val=1.0){
											 		input_ready = (get_connected_from(mat,core_nodes+1))[1:core_nodes];
											 		output_ready = (get_connected_to(mat,core_nodes+2))[1:core_nodes];
											 		core = mat[1:core_nodes,1:core_nodes];
											 		while(length(which(sapply(input_ready,is.null)))!=0 || length(which(sapply(output_ready,is.null)))!=0){
											 			zeros = which(core==0);
											 			if(length(zeros)==0){
											 			  print("WE SHOULD NOT BE HERE!!!");
											 			  break;
											 			}
											 		  non_connected=sample(zeros,1);
											 			from = non_connected%%core_nodes;
											 			if(from==0){from = core_nodes;}
											 			to = ceiling(non_connected/core_nodes)
											 			core[from,to]=new_val;
											 			if(!is.null(input_ready[[from]]) && is.null(input_ready[[to]])){
											 				input_ready=iterate_through_mat_inout(core,input_ready,to,TRUE);
											 			}
											 			if(!is.null(output_ready[[to]]) && is.null(output_ready[[from]])){
											 				output_ready=iterate_through_mat_inout(core,output_ready,from,FALSE);
											 			}
											 		}
											 		if(!allow_autoloops){
											 			diag(core)=0;
											 		}
											 		mat_aux = mat;
											 		mat_aux[1:core_nodes,1:core_nodes] = core;
											 		if(set_adjacency){
											 			adjacency_matrix <<- mat_aux;
											 		}
											 		return(mat_aux);
											 	},
											 	set_random_weights=function(min_val=0.0,max_val=1.0){
											 		mat = matrix(0,nrow=nrow(adjacency_matrix),ncol=ncol(adjacency_matrix));
											 		non_zero_components = which(adjacency_matrix>0);
											 		mat[non_zero_components]=runif(length(non_zero_components), min=min_val, max=max_val);
											 		adjacency_matrix <<- mat;
											 	},
											 	normalize_rows=function(mat=adjacency_matrix,set_adjacency=TRUE){
											 		sums = apply(mat,1,sum);
											 		sums[which(sums==0)]=1;
											 		res = sweep(mat,1,sums,"/");
											 		if(set_adjacency){
											 			adjacency_matrix<<-res;
											 		}else{
											 			return(res);
											 		}
											 	},
											 	normalize_cols=function(mat=adjacency_matrix,set_adjacency=FALSE){
											 	  sums = apply(mat,2,sum);
											 	  sums[which(sums==0)]=1;
											 	  res = sweep(mat,2,sums,"/");
											 	  if(set_adjacency){
											 	    adjacency_matrix<<-res;
											 	  }else{
											 	    return(res);
											 	  }
											 	},
											 	get_binary_adjacency=function(mat=adjacency_matrix){
											 		mat_aux = matrix(0,nrow=nrow(mat),ncol=ncol(mat));
											 		mat_aux[which(mat>0)]=1;
											 		return(mat_aux);
											 	},
											 	operate_adjacency_matrix=function(mat=adjacency_matrix,initial_input=1.0,frequency=1L,iterations=100L,do_plot=FALSE){
											 		X=rep(0,ncol(mat));
											 		if(do_plot){
											 		  out=vector(mode='list',length=iterations);
											 		}
											 		input=X;
											 		input[[core_nodes+1]]=initial_input;
											 		for(i in 1:iterations){
											 			if((i%%frequency)==0){
											 				X=(X+input)%*%mat
											 			}else{
											 				X=(X)%*%mat
											 			}
											 		  if(do_plot){
											 		    out[[i]]=X
											 		  }
											 		}
											 		if(do_plot){
											 		  serie_mat = t(matrix(unlist(out), ncol=length(out), byrow=FALSE))
											 		  serie_mat=serie_mat[,-c(length(out[[1]]),length(out[[1]])-1)]
											 		  df <- data.frame(serie_mat)
											 		  df['index']=seq(1,iterations,1)
											 		  df <- melt(df ,  id.vars = 'index', variable.name = 'series')
											 		  print(ggplot(df, aes(index, value)) + geom_line(aes(colour = series)))
											 		}
											 		return(X) 
											 	},
											 	check_is_mat_same_as_adjacency=function(mat){
											 	  return(length(which(abs(mat-adjacency_matrix)>0))==0);
											 	},
											 	iteration_series_of_adjacency_operations=function(mat=adjacency_matrix,initial_input=1.0,frequency=2L,after_pulses=30L,do_plot=FALSE){
											 		base_iterations = frequency*after_pulses;
											 		out=vector(mode='list',length=frequency);
											 		for(i in 1:frequency){
											 			out[[i]]=operate_adjacency_matrix(mat=mat,initial_input=initial_input,frequency=frequency,iterations=base_iterations+(i-1));
											 		}
											 		if(do_plot){
											 		  serie_mat = t(matrix(unlist(out), ncol=length(out), byrow=FALSE))
											 		  serie_mat=serie_mat[,-c(length(out[[1]]),length(out[[1]])-1)]
											 		  df <- data.frame(serie_mat)
											 		  df['index']=seq(1,frequency,1)
											 		  df <- melt(df ,  id.vars = 'index', variable.name = 'series')
											 		  print(ggplot(df, aes(index, value)) + geom_line(aes(colour = series)))
											 		}
											 		return(out);
											 	},
											 	get_quantity_matrix=function(mat=adjacency_matrix,multiplier=1000,iterations=100L){
											 		X = operate_adjacency_matrix(mat=mat,initial_input=multiplier,frequency=1L,iterations=iterations);
											 		X[[core_nodes+1]]=multiplier;
											 		return(sweep(mat,1,X,"*"));
											 	},
											 	get_as_igraph=function(mat=adjacency_matrix){
											 		return(graph_from_adjacency_matrix(mat,weighted=TRUE));
											 	},
											 	get_core_edge_set=function(mat=adjacency_matrix){
											 		core = mat[1:core_nodes,1:core_nodes];
											 		edges=which(core>0);
											 		out=set();
											 		for(i in 1:length(edges)){
											 			from = edges[[i]]%%core_nodes;
											 			if(from==0){from = core_nodes;}
											 			to = as.integer(ceiling(edges[[i]]/core_nodes));
											 			out=set_union(out,set(tuple(from,to)));
											 		}
											 		return(out);
											 	},
											 	get_core=function(mat=adjacency_matrix){
											 		return(mat[1:core_nodes,1:core_nodes]);
											 	},
											 	get_edges_that_connect_with=function(edges,ed){
											 		out=set();
											 		for(i in 1:length(edges)){
											 			if(edges[[i]][[1]]==ed[[2]]){
											 				out=set_union(out,set(edges[[i]]));
											 			}
											 		}
											 		return(out);
											 	},
											 	get_paths_from=function(edges,paths,started_node){
											 		out = list();
											 		cycles = set();
											 		for(p in paths){
											 			if(length(p[[2]])==0)
											 				next;
											 			connect = get_edges_that_connect_with(edges,p[[2]][[length(p[[2]])]]);
											 			for(ed in connect){
											 				if(ed[[2]]==started_node){
											 					p2_new=append(p[[2]],list(ed));
											 					cycles=set_union(cycles,set(as.set(p2_new)));
											 				}else if(set_contains_element(p[[1]], ed[[2]]) || ed[[2]]<started_node){
											 					# NO GOOD PATH
											 				}else{
											 					p2_new=append(p[[2]],set(ed));
											 					p1_new=set_union(p[[1]],set(ed[[2]]));
											 					out=append(out,list(list(p1_new,p2_new)));
											 				}
											 			}
											 		}
											 		return(list(out,cycles));
											 	},
											 	get_cycles_R=function(){
											 		out=set();
											 		edges = as.list(get_core_edge_set());
											 		for(i in 1:(length(edges)-1)){
											 			if(edges[[i]][[2]]<edges[[i]][[1]]){
											 				next;
											 			}
											 			should_we_search = FALSE;
											 			for(j in (i+1):length(edges)){
											 				if(edges[[i]][[1]] == edges[[j]][[2]]){
											 					if(edges[[i]][[2]] == edges[[j]][[1]]){
											 						out=set_union(out,set(set(edges[[i]],edges[[j]])));
											 						next;
											 					}
											 					should_we_search = TRUE;
											 				}
											 			}
											 			if(should_we_search){
											 				current = get_edges_that_connect_with(edges[(i+1):length(edges)],edges[[i]]);
											 				paths = list();
											 				for(ed in current){
											 					if(ed[[2]]!=edges[[i]][[1]]){
												 					paths=append(paths,
												 								 list(list(
												 								 	set(edges[[i]][[2]],ed[[2]]), #1 its the set of nodes
												 								 	list(edges[[i]],ed)))          #2 is the set of edges or tuples
												 								 );
											 					}
											 				}
											 				res=get_paths_from(edges[(i+1):length(edges)],paths,edges[[i]][[1]]);
											 				if(length(res[[2]])>0){
											 					out=set_union(out,res[[2]]);
											 				}
											 				while(length(res[[1]])>0){
											 					res=get_paths_from(edges[(i+1):length(edges)],res[[1]],edges[[i]][[1]]);
											 					if(length(res[[2]])>0){
											 						out=set_union(out,res[[2]]);
											 					}
											 				}
											 			}
											 		}
											 		return(out);
											 	},
											 	get_xdegree_by_nodes=function(mat=adjacency_matrix,in_or_out='in'){
											 	  nc = ncol(mat);
											 	  X = rep(0,nc);
											 	  for(i in 1:nc){
											 	    if(in_or_out=='in'){
											 	      X[i] = length(which(mat[,i]>0));
											 	    }else{
											 	      X[i] = length(which(mat[i,]>0));
											 	    }
											 	  }
											 	  return(X);
											 	},
											 	get_TST=function(mat=adjacency_matrix,input_multiplier=1000,consider_input=FALSE,consider_output=FALSE,proportional_to_input=FALSE){
											 		final_node_in = core_nodes;
											 		final_node_out = core_nodes;
											 		if(consider_input)
											 			final_node_in = core_nodes + 1
											 		if(consider_output)
											 			final_node_out = core_nodes + 2
											 		m = get_quantity_matrix(mat=mat,multiplier=input_multiplier)[1:final_node_in,1:final_node_out];
											 		total = sum(m)#sum(m[which(m>0)])
											 		if(proportional_to_input)
											 			return(total/input_multiplier)
											 		return(total)
											 	},
											 	navigate_net=function(mat=adjacency_matrix,save_trajectories=FALSE){
											 	  system.time(res<-paths(mat,core_nodes+1,core_nodes+2,save_trajectories))
											 	  if(check_is_mat_same_as_adjacency(mat)){
  											    net_data$cycle_participation <<- t(matrix(unlist(res[[4]]), ncol = core_nodes+2, nrow = core_nodes+2));
  											    net_data$trajectory_participation <<- t(matrix(unlist(res[[2]]), ncol = core_nodes+2, nrow = core_nodes+2));
  											    if (save_trajectories){
  											      net_data$cycles <<- res[[3]];
  											      net_data$trajectories <<- res[[1]];
  											    }
  											    return(net_data)
											 	  }else{
											 	    output=list()
											 	    output$cycle_participation = t(matrix(unlist(res[[4]]), ncol = core_nodes+2, nrow = core_nodes+2));
											 	    output$trajectory_participation = t(matrix(unlist(res[[2]]), ncol = core_nodes+2, nrow = core_nodes+2));
											 	    if (save_trajectories){
											 	      output$cycles = res[[3]];
											 	      output$trajectories = res[[1]];
											 	    }
											 	    return(output)
											 	  }
											 	},
											 	get_cycles_flow=function(mat=adjacency_matrix,x=vector(mode='numeric'),cycles=list()){
											 	  cycles_to_use = cycles;
											 	  X = x;
											 	  if(length(x)==0){
											 	    X = as.vector(operate_adjacency_matrix(mat = mat));
											 	  }
											 	  if(length(cycles)==0){
											 	    if(is.null(net_data$cycles)){
											 	      navigate_net(save_trajectories = TRUE);
											 	    }
											 	    cycles_to_use = net_data$cycles;
											 	  }
											 	  
											 	  output = vector(mode='list',length=length(cycles_to_use));
											 	  
											 	  for(k in 1:length(cycles_to_use)){
											 	    current_element = list();
											 	    current_element$cycle_length = length(cycles_to_use[[k]]); 
											 	    mult = 1.0;
											 	    for(i in 1:(length(cycles_to_use[[k]])-1)){
											 	      mult = mult * mat[cycles_to_use[[k]][[i]],cycles_to_use[[k]][[i+1]]];
											 	    }
											 	    mult = mult * mat[cycles_to_use[[k]][[length(cycles_to_use[[k]])]],cycles_to_use[[k]][[1]]];
											 	    cycle_val = sum(X[cycles_to_use[[k]]])*mult;
											 	    current_element$cycle_flow=cycle_val;
											 	    output[[k]] = current_element;
											 	  }
											 	  return(output);
											 	},
											 	get_trajectories_from_link=function(from,to){
											 	if(is.null(net_data$trajectories)) return(list())
											 	out=list()
											 	for (i in 1:length(net_data$trajectories)){
											 	  positions=match(c(from,to),net_data$trajectories[[i]])
											 	  if(!is.na(positions[1]) && !is.na(positions[2])){
											 	    dif = positions[2]-positions[1];
											 	    if(dif==1){
											 	      out=append(out,list(net_data$trajectories[[i]]))
											 	    }
											 	  }
											 	}
											 	return(out)
											 	},
											 	calculate_link_type_entropy_from_participation=function(cycle_part,traj_part){
											 	  sum_of_parts = (cycle_part+traj_part);
											 	  cycle_prop = cycle_part / sum_of_parts;
											 	  nlink = length(which(!is.na(cycle_prop[1:core_nodes,1:core_nodes])));
											 	  cycle_prop[which(is.na(cycle_prop))]=0; 
											 	  traj_prop = traj_part / sum_of_parts;
											 	  traj_prop[which(is.na(traj_prop))]=0;
											 	  s = -cycle_prop*log2(cycle_prop)-traj_prop*log2(traj_prop);
											 	  s[which(is.na(s))] = 0; # NUEVA LINEA OJO!!! 
											 	  sum_of_ent = sum(s);#sum(s[which(!is.na(s))]);
											 	  if(nlink==0){
											 	    nlink=1;
											 	  }
											 	  output=list();
											 	  output$average_ent_per_link = sum_of_ent/nlink;
											 	  output$total_link_ent = sum_of_ent;
											 	  output$mat = s;
											 	  return(output);
											 	},
											 	get_link_type_entropy=function(mat=adjacency_matrix,recalculate=TRUE,adjoin_netdata=TRUE){
											 	  netdata = net_data;
											 	  is_adjacency = check_is_mat_same_as_adjacency(mat);
											 	  if(is_adjacency && !recalculate && !is.null(net_data$cycle_participation) && !is.null(net_data$trajectory_participation)){
											 	    
											 	  }else{
											 	    netdata = navigate_net(mat = mat,save_trajectories = TRUE);
											 	  }
											 	  
											 	  mean_s = calculate_link_type_entropy_from_participation(netdata$cycle_participation,netdata$trajectory_participation);
											 	  
											 	  if(is_adjacency){
											 	    net_data$link_type_entropy <<- mean_s;
											 	  }
											 	  output=list();
											 	  output$link_type_entropy = mean_s;
											 	  if(adjoin_netdata)
											 	    output$netdata = netdata;
											 	  return(output);
											 	}
											 	
											 	)
)

NetEnsamble <- setRefClass("NetEnsamble",
													 fields = list(ensemble_sizeMin = "numeric",
													 							ensemble_sizeMax = "numeric", 
													 							#size_interval = "numeric", 
													 							size_replicate = "numeric",
													 							nodes_autoloops_allowed = "logical",
													 							random_seed = "numeric",
													 							nodes_connect_function = "function",
													 							node_seq = "vector",
													 							generated_networks = "list"
													 							),
													 methods = list(
													 	initialize = function(...,
													 																			ensemble_sizeMin = 2L,
													 																			ensemble_sizeMax = 10L,
													 																			#size_interval = 1L,
													 																			size_replicate = 10L,
													 																			nodes_autoloops_allowed = FALSE,
													 																			random_seed = 55L){
													 	ensemble_sizeMin <<- ensemble_sizeMin;
													 	ensemble_sizeMax <<- ensemble_sizeMax;
													 	size_replicate <<- size_replicate;
													 	nodes_autoloops_allowed <<- nodes_autoloops_allowed;
													 	node_seq <<- as.integer(seq(ensemble_sizeMin,ensemble_sizeMax,1));
													 	callSuper(...)
														},
														generate = function(){
															num_total = (ensemble_sizeMax-ensemble_sizeMin)*size_replicate;
															#print(num_total);
															output = vector(mode="list",length = num_total);
															for(i in ensemble_sizeMin:ensemble_sizeMax){
																for(j in 1:size_replicate){
																	current = NetBase$new(core_nodes = i,allow_autoloop=nodes_autoloops_allowed,respiration=FALSE)
																	#current$turn_to_flowing_network()
																	#current$set_random_weights()
																	#current$normalize_rows()
																	output[size_replicate*(i-ensemble_sizeMin)+j] = current
																}
															}
															generated_networks <<-output;
														},
														navigate_all_paths=function(save_trajectories=FALSE){
														  for(i in 1:length(generated_networks)){
														    generated_networks[[i]]$navigate(save_trajectories=save_trajectories)
														    cat('Navigated',i,'of',length(generated_networks),'\n')
														  }
														}
														)
)
						
NetEvolve <- setRefClass("NetEvolve",
												fields = list(initial_net = "NetBase",
																			current_net = "matrix",
																			current_state = "list",
																			mutation_func = "function",
																			characteristic_funcs = "list",
																			evolution_func = "function"
												),
												methods = list(
													initialize = function(...,
																								initial_net = NetBase$new(),
																								mutation_func = function(base,mat){
																									
																								},
																								characteristic_funcs = list(function(base,mat){return(1)}),
																								evolution_func = function(iteration_num,base,current_state,mutated_net,characteristic_vals){
																									return(list(new_state = current_state,
																															new_net = current_net,
																															continue = FALSE));
																								}){
														# TODO : CHECK VARIABLES ARE COHERENT
														initial_net <<- initial_net;
														mutation_func <<- mutation_func;
														current_state <<- list(vals=list(),extra=list());
														current_net <<- initial_net$adjacency_matrix;
														evolution_func <<- evolution_func;
														characteristic_funcs <<- characteristic_funcs;
														callSuper(...);
													},
													evolve = function(){
														should_continue = TRUE;
														iteration = 1;
														characteristic_vals = vector(mode="list",length=length(characteristic_funcs));
														for(i in 1:length(characteristic_funcs)){
															characteristic_vals[[i]] = characteristic_funcs[[i]](initial_net,initial_net$adjacency_matrix);
														}
														current_state$vals <<- characteristic_vals;
														current_state$initial_vals <<- characteristic_vals;
														while(should_continue){
														  system(sprintf('echo "\n%s\n"', paste0('Iteration num for net with ',initial_net$core_nodes,' : ',iteration, collapse="")));
															mutated_net = mutation_func(initial_net,current_net);
															for(i in 1:length(characteristic_funcs)){
																characteristic_vals[[i]] = characteristic_funcs[[i]](initial_net,mutated_net);
															}
															result = evolution_func(iteration_num = iteration,
																											base = initial_net,
																											current_state = current_state,
																											current_net = current_net,
																											mutated_net = mutated_net,
																											characteristic_vals = characteristic_vals
																											);
															current_net <<- result$new_net;
															current_state <<- result$new_state;
															should_continue = result$continue;
															iteration = iteration + 1;
															#if(!(result$continue))
															#	break;
														}
													}
												)
)

EnsambleEvolve <- setRefClass("EnsambleEvolve",
                              fields = list(ensamble = "NetEnsamble",
                                            mutation_func = "function",
                                            characteristic_funcs = "list",
                                            evolution_func = "function",
                                            evolved_nets = "list"
                              ),
                              methods = list(
                                initialize = function(...,
                                                      ensamble = NetEnsamble$new(),
                                                      mutation_func = function(base,mat){
                                                        
                                                      },
                                                      characteristic_funcs = list(function(base,mat){return(1)}),
                                                      evolution_func = function(iteration_num,base,current_state,mutated_net,characteristic_vals){
                                                        return(list(new_state = current_state,
                                                                    new_net = current_net,
                                                                    continue = FALSE));
                                                      }){
                                  ensamble <<- ensamble;
                                  mutation_func <<- mutation_func;
                                  evolution_func <<- evolution_func;
                                  characteristic_funcs <<- characteristic_funcs;
                                  evolved_nets <<- vector(mode="list",length=length(ensamble$generated_networks));
                                  callSuper(...);
                                  if(length(ensamble$generated_networks)>0){
                                    for(i in 1:length(ensamble$generated_networks)){
                                      evolved_nets[[i]] <<- NetEvolve$new(initial_net=ensamble$generated_networks[[i]],
                                                                          characteristic_funcs=characteristic_funcs,
                                                                          mutation_func=mutation_func,
                                                                          evolution_func=evolution_func);
                                    }
                                  }
                                },
                                evolve_ensamble = function(){
                                  for(i in 1:length(evolved_nets)){
                                    evolved_nets[[i]]$evolve();
                                    cat('Ensamble evolve',i,'of',length(evolved_nets),'\n')
                                  }
                                },
                                evolve_ensamble_par = function(){
                                    num_cores=detectCores();
                                    #lapply(evolved_nets,function(x){x$evolve()});
                                    
                                    res<-mclapply(evolved_nets,function(x){x$evolve();return(x);},mc.cores=num_cores);
                                    evolved_nets<<-res;
                                },
                                get_results = function(){
                                  output = vector(mode="list",length=length(evolved_nets));
                                  for(i in 1:length(evolved_nets)){
                                    output[[i]]$initial = evolved_nets[[i]]$initial_net
                                    final = NetBase$new()
                                    final$create_from_base_and_adjacency(base = evolved_nets[[i]]$initial_net,
                                                                         adjacency = evolved_nets[[i]]$current_net)
                                    output[[i]]$final = final
                                    output[[i]]$state = evolved_nets[[i]]$current_state
                                  }
                                  return(output)
                                }
)
)

Visualization <- setRefClass("Visualization",
                     fields = list(data = "list"),
                     methods = list(
                       initialize=function(...,in_data){
                         data<<-lapply(in_data,function(x){
                           if(length(x)==0) return(NULL) 
                           else return(
                             data.frame(t(matrix(unlist(x),nrow = length(names(in_data[[1]][[1]])),
                                                 dimnames=list(names(in_data[[1]][[1]]),
                                                               vector(mode="expression"))))))
                           })
                         data<<-data[lengths(data) != 0]
                       },
                       plot_me=function(x_var,
                                        y_var,
                                        col_var = NULL,
                                        size_var = NULL,
                                        color_func=function(x){return(hsv(h = (x-min_max_col$min)/(min_max_col$max-min_max_col$min)*0.8, s = 1, v = 1, 1))},
                                        size_func=function(x){return(0.5+(x-min_max_size$min)/(min_max_size$max-min_max_size$min))},
                                        x_lim=NULL,
                                        y_lim=NULL){
                         
                         if(!is.null(col_var)){
                           min_max_col = get_min_max(data,col_var);
                         }
                         if(!is.null(size_var)){
                           min_max_size = get_min_max(data,size_var);
                         }
                         
                         min_max_x = x_lim;
                         min_max_y = y_lim;
                         if(is.null(x_lim)){
                           min_max_x = get_min_max(data,x_var);
                         }
                         if(is.null(y_lim)){
                           min_max_y = get_min_max(data,y_var);
                         }
                         layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
                         par(mar = c(5.1, 4.1, 4.1, 2.1));
                         for(i in 1:length(data)){
                           col = 'black';
                           size = 0.5;
                           if(!is.null(col_var)){
                            col = unlist(lapply(data[[i]][[col_var]],color_func));
                           }
                           if(!is.null(size_var)){
                            size = unlist(lapply(data[[i]][[size_var]],size_func));
                           }
                           if(i==1){
                             plot(data[[i]][[x_var]],data[[i]][[y_var]], 
                                  col = col,
                                  cex = size, 
                                  xlab = x_var, 
                                  ylab = y_var,
                                  xlim=c(min_max_x$min,min_max_x$max),
                                  ylim=c(min_max_y$min,min_max_y$max))
                           }
                           else{
                             points(data[[i]][[x_var]],
                                    data[[i]][[y_var]], 
                                    col = col,
                                    cex = size)
                           }
                           line_col = col;
                           if(!is.null(col_var)){
                             line_col = col[-1L];
                           }
                           segments(data[[i]][[x_var]][-length(data[[i]][[x_var]])],
                                    data[[i]][[y_var]][-length(data[[i]][[y_var]])],
                                    data[[i]][[x_var]][-1L],data[[i]][[y_var]][-1L],
                                    col=line_col);
                         }
                         col_seq_range = seq(min_max_col$min,min_max_col$max,l=5);
                         legend_image <- as.raster(matrix(unlist(lapply(col_seq_range,color_func)), ncol=1))
                         par(mar = c(2, 0, 2, 0));
                         plot(c(0,10),c(min_max_col$min,min_max_col$max),type = 'n', axes = F,xlab = '', ylab = '', main = col_var)
                         text(x=5, y = col_seq_range, labels =  format(col_seq_range, scientific = TRUE, digits = 3));
                         rasterImage(legend_image, 0, min_max_col$max, 2,min_max_col$min)
                         par(mar = c(5.1, 4.1, 4.1, 2.1));
                       },
                       get_min_max=function(dfs,variable){
                         maximum=vector(mode="numeric",length=length(dfs))
                         minimum=vector(mode="numeric",length=length(dfs))
                         for(i in 1:length(dfs)){
                           maximum[i]=max(dfs[[i]][variable])
                           minimum[i]=min(dfs[[i]][variable])
                           if(is.infinite(maximum[i])) maximum[i] = -Inf;
                           if(is.infinite(minimum[i])) minimum[i] = Inf;
                         }
                         y_max=max(maximum)
                         y_min=min(minimum)
                         output=list();
                         output$max=y_max;
                         output$min=y_min;
                         return(output);
                       }
                     )
)
