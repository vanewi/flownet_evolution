library(igraph)
library(sets)


NetBase <- setRefClass("NetBase",
											 fields = list(core_nodes='numeric',nodes_connect_function='function',adjacency_matrix='matrix',allow_autoloop='logical',respiration='logical',total_nodes='numeric'),
											 methods = list(
											 	initialize=function(core_nodes=5L,
											 											nodes_connect_function=function(N,autoloop,respiration){
												 												n_input=sample(1:N,1);
												 												n_output=sample(1:N,1);
												 												non_core_nodes = if(respiration) 3 else 2;
												 												mat = matrix(0, nrow = N+non_core_nodes, ncol = N+non_core_nodes);
												 												input_ready = vector(mode="list",length = N);
												 												output_ready = vector(mode="list",length = N);
												 												mat[N+1,i:n_input]=1;
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
											 											...){
											 		core_nodes<<-core_nodes;
											 		nodes_connect_function <<- nodes_connect_function;
											 		allow_autoloop <<- allow_autoloop;
											 		respiration <<- respiration;
											 		adjacency_matrix <<- nodes_connect_function(core_nodes,autoloop,respiration);
											 		total_nodes <<- nrow(adjacency_matrix);
											 		callSuper(...);
											 	},
											 	ami=function(){},
											 	mean_time_discrete=function(input_percentage=0.5,mat=adjacency_matrix){
											 		X=rep(0,ncol(mat));
											 		X[[core_nodes+1]]=1.0;
											 		innet = 1.0
											 		iters = 0
											 		while(innet > input_percentage){
											 				X=(X)%*%mat
											 				iters = iters+1
											 				innet = sum(X[1:N])
											 		}
											 		return(iters) 
											 	},
											 	net_plot=function(mat=adjacency_matrix){
											 		#plot.network(as.network(adjacency_matrix),interactive=TRUE,label = seq(1:total_nodes));
											 		g = graph_from_adjacency_matrix(mat,weighted=TRUE)
				
											 		plot(g,edge.arrow.size=0.1,edge.label=round(E(g)$weight, 3),edge.label.cex=0.7);
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
											 	get_connected_from=function(mat,node){
											 		connected = vector(mode="list",length = nrow(mat));
											 		return(iterate_through_mat_inout(mat,connected,node,TRUE));
											 	},
											 	get_connected_to=function(mat,node){
											 		connected = vector(mode="list",length = nrow(mat));
											 		return(iterate_through_mat_inout(mat,connected,node,FALSE));
											 	},
											 	turn_to_flowing_network=function(){
											 		input_ready = (get_connected_from(adjacency_matrix,core_nodes+1))[1:core_nodes];
											 		output_ready = (get_connected_to(adjacency_matrix,core_nodes+2))[1:core_nodes];
											 		core = adjacency_matrix[1:core_nodes,1:core_nodes];
											 		while(length(which(sapply(input_ready,is.null)))!=0 || length(which(sapply(output_ready,is.null)))!=0){
											 			non_connected=sample(which(core==0),1);
											 			from = non_connected%%core_nodes;
											 			if(from==0){from = core_nodes;}
											 			to = ceiling(non_connected/core_nodes)
											 			core[from,to]=1
											 			if(!is.null(input_ready[[from]]) && is.null(input_ready[[to]])){
											 				input_ready=iterate_through_mat_inout(core,input_ready,to,TRUE);
											 			}
											 			if(!is.null(output_ready[[to]]) && is.null(output_ready[[from]])){
											 				output_ready=iterate_through_mat_inout(core,output_ready,from,FALSE);
											 			}
											 		}
											 		if(!allow_autoloop){
											 			diag(core)=0;
											 		}
											 		adjacency_matrix[1:core_nodes,1:core_nodes] <<- core;
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
											 	get_binary_adjacency=function(){
											 		mat = matrix(0,nrow=nrow(adjacency_matrix),ncol=ncol(adjacency_matrix));
											 		mat[which(adjacency_matrix>0)]=1
											 		return(mat)
											 	},
											 	operate_adjacency_matrix=function(initial_input=1.0,frequency=1L,iterations=100L){
											 		X=rep(0,ncol(adjacency_matrix));
											 		input=X;
											 		input[[core_nodes+1]]=initial_input;
											 		for(i in 1:iterations){
											 			if((i%%frequency)==0){
											 				X=(X+input)%*%adjacency_matrix
											 			}else{
											 				X=(X)%*%adjacency_matrix
											 			}
											 		}
											 		return(X) 
											 	},
											 	iteration_series_of_adjacency_operations=function(initial_input=1.0,frequency=2L,after_pulses=30L){
											 		base_iterations = frequency*after_pulses;
											 		out=vector(mode='list',length=frequency);
											 		for(i in 1:frequency){
											 			out[[i]]=operate_adjacency_matrix(initial_input,frequency,base_iterations+(i-1));
											 		}
											 		return(out);
											 	},
											 	get_quantity_matrix=function(multiplier=1000){
											 		X = operate_adjacency_matrix(initial_input=multiplier,frequency=1L,iterations=100L);
											 		X[[core_nodes+1]]=multiplier;
											 		return(sweep(adjacency_matrix,1,X,"*"));
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
											 	get_core=function(){
											 		return(adjacency_matrix[1:core_nodes,1:core_nodes]);
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
											 	get_cycles=function(){
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
											 	})
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
																	current = NetBase$new(core_nodes =i,allow_autoloop=nodes_autoloops_allowed,respiration=TRUE)
																	current$turn_to_flowing_network()
																	current$set_random_weights()
																	current$normalize_rows()
																	output[size_replicate*(i-ensemble_sizeMin)+j] = current
																}
															}
															generated_networks <<-output;
														}
														)
)
						

