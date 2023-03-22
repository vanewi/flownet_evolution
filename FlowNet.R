library(igraph)
library(sets)
library(Rcpp)
library(ggplot2)
library(reshape2)
library(bettermc)
library(purrr)
sourceCpp("flownet.cpp")
source("generator_util_functions.R")

PRINT_FLOWNET_DEBUG = TRUE

NetBasicNodesConnectFunction = connection_generator_with_renyi_opt(do_renyi = FALSE);

NetBase <- setRefClass("NetBase",
                       fields = list(core_nodes='numeric',
                                     adjacency_matrix='matrix',
                                     allow_autoloop='logical',
                                     respiration='logical',
                                     total_nodes='numeric',
                                     net_data='list',
                                     name='character'),
                       methods = list(
                           initialize=function(core_nodes=5L,
                                               nodes_connect_function=NetBasicNodesConnectFunction,
                                               allow_autoloop=FALSE,
                                               respiration = FALSE,
                                               make_flowing = TRUE,
                                               minimal_flowing = FALSE,
                                               random_weights = TRUE,
                                               normalized_rows = TRUE,
                                               fully_connected = FALSE,
                                               name='GenericName',
                                               ...){
                               name <<- name;
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
                                   turn_to_flowing_network(minimal=minimal_flowing);
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
                           get_ascendency=function(mat=adjacency_matrix,input_multiplier=1000,consider_input=TRUE,consider_output=TRUE,log_base=2){
                               final_node_in = core_nodes;
                               final_node_out = core_nodes;
                               if(consider_input)
                                   final_node_in = core_nodes + 1
                               if(consider_output)
                                   final_node_out = core_nodes + 2
                               m = get_quantity_matrix(mat=mat,multiplier=input_multiplier)[1:final_node_in,1:final_node_out];
                               tst=get_TST(mat = mat,input_multiplier=input_multiplier,consider_input = consider_input,consider_output = consider_output)
                               asc=0
                               for (i in 1:final_node_in){
                                   for (j in 1:final_node_out){
                                       if(m[i,j]==0) next
                                       val = (m[i,j]*tst)/(sum(m[i,])*sum(m[,j]));
                                       asc=asc+(m[i,j])*(log2(val)/log2(log_base));
                                   }
                               }
                               return(asc)
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
                           net_plot=function(mat=adjacency_matrix,margin=0.0,interactive=TRUE){
                               g = graph_from_adjacency_matrix(mat,weighted=TRUE)
                               
                               names=as.character(c(seq(1,core_nodes,1),"IN","OUT"))
                               color_pallete=c("orange","green",'red')
                               colors=as.factor(c(rep(1,core_nodes),2,3))
                               
                               set_vertex_attr(g,"names",value=names)
                               set_vertex_attr(g,"colors",value=colors)

                               if(interactive){
                                   tkplot(g,
                                          edge.label=round(E(g)$weight, 3),
                                          vertex.color=color_pallete[colors],vertex.label=names,
                                          margin=margin,layout=layout_with_lgl(g))#, edge.color=edge.col);
                               }else{
                                   plot(g,
                                        edge.arrow.size=0.1,edge.label=round(E(g)$weight, 3),edge.label.cex=0.7,
                                        vertex.size=20,vertex.color=color_pallete[colors],vertex.label=names,
                                        margin=margin,layout=layout_with_lgl(g))#, edge.color=edge.col);
                               }
                           },
                           iterate_through_mat_inout=function(core,nodes_ready,node,in_or_out){
                               if(!is.null(nodes_ready[[node]])){
                                   return(nodes_ready);
                               }
                               nodes_ready_temp = nodes_ready;
                               nodes_ready_temp[[node]]=TRUE;
                               connected = list();
                               if(in_or_out=='in'){
                                   connected=which(core[node,]>0);
                               }else{
                                   connected=which(core[,node]>0);
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
                               return(iterate_through_mat_inout(mat,connected,node,'in'));
                           },
                           get_connected_to=function(mat=adjacency_matrix,node){
                               connected = vector(mode="list",length = nrow(mat));
                               return(iterate_through_mat_inout(mat,connected,node,'out'));
                           },
                           turn_to_flowing_network=function(mat=adjacency_matrix,allow_autoloops=allow_autoloop,set_adjacency=TRUE,new_val=1.0,minimal=FALSE){
                                        # TODO : MAKE ALGORITHM DEPENDANT OPTIONS FOR THIS (MINIMAL VERSUS RANDOM SEARCH FOR EXAMPLE)
                               input_ready = (get_connected_from(mat,core_nodes+1))[1:core_nodes];
                               output_ready = (get_connected_to(mat,core_nodes+2))[1:core_nodes];
                               missing_in = which(sapply(input_ready,is.null));
                               missing_out = which(sapply(output_ready,is.null));
                                        # TODO : FIX IN CASE READY IN OR OUT IS EMPTY!!
                               core = mat[1:core_nodes,1:core_nodes];
                               full = mat;
                               if(minimal){ # SELECT ALL FAILING IN AND OUT AND CONNECT RANDOMLY TO READY IN EACH SET
                                   while(length(missing_in)>0 || length(missing_out)>0){
                                       ready_in = which(sapply(input_ready,isTRUE));
                                       ready_out = which(sapply(output_ready,isTRUE));
                                       ready_in = append(ready_in,core_nodes+1);
                                       ready_out = append(ready_out,core_nodes+2);
                                       if(length(missing_in)>0){
                                           to = missing_in[1];
                                           from = sample(ready_in,1);
                                           full[from,to]=new_val;
                                           input_ready=iterate_through_mat_inout(full[1:core_nodes,1:core_nodes],input_ready,to,'in');
                                           missing_in = which(sapply(input_ready,is.null));
                                       }
                                       if(length(missing_out)>0){
                                           from = missing_out[1];
                                           to = sample(ready_out,1);
                                           full[from,to]=new_val;
                                           output_ready=iterate_through_mat_inout(full[1:core_nodes,1:core_nodes],output_ready,from,'out');
                                           missing_out = which(sapply(output_ready,is.null));
                                       }
                                   }
                                   core = full[1:core_nodes,1:core_nodes];
                               }
                               else{ #FULL ALGO, RANDOM TILL EVERY NODE IS FLOWING
                                   while(length(which(sapply(input_ready,is.null)))!=0 || length(which(sapply(output_ready,is.null)))!=0){
                                       zeros = which(core==0);
                                       if(length(zeros)==0){
                                           print("WE SHOULD NOT BE HERE!!!");
                                           print(adjacency_matrix);
                                        #break;
                                       }
                                       non_connected=sample(zeros,1);
                                       from = non_connected%%core_nodes;
                                       if(from==0){from = core_nodes;}
                                       to = ceiling(non_connected/core_nodes)
                                       core[from,to]=new_val;
                                       if(!is.null(input_ready[[from]]) && is.null(input_ready[[to]])){
                                           input_ready=iterate_through_mat_inout(core,input_ready,to,'in');
                                       }
                                       if(!is.null(output_ready[[to]]) && is.null(output_ready[[from]])){
                                           output_ready=iterate_through_mat_inout(core,output_ready,from,'out');
                                       }
                                   }
                               }
                               if(!allow_autoloops){
                                   diag(core)=0;
                               }
                               mat_aux = full;
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
                           check_is_mat_bin_same_as_adjacency_bin=function(mat){
                               return(length(which(abs(get_binary_adjacency(mat=mat)-get_binary_adjacency())>0))==0);
                           },
                           check_is_mat_bin_same_as_mat_bin=function(mat,mat2){
                               return(length(which(abs(get_binary_adjacency(mat=mat)-get_binary_adjacency(mat=mat2))>0))==0);
                           },
                           plot_series_of_vectors=function(vector_serie,xmax=0){
                               serie_mat = t(matrix(unlist(vector_serie), ncol=length(vector_serie), byrow=FALSE))
                               serie_mat=serie_mat[,-c(length(vector_serie[[1]]),length(vector_serie[[1]])-1)]
                               df <- data.frame(serie_mat)
                               limit = length(vector_serie)
                               if(xmax>0 && xmax<limit){
                                   limit=xmax;
                               }
                               df['index']=seq(1,length(vector_serie),1)
                               df <- melt(df[1:limit,] ,  id.vars = 'index', variable.name = 'series')
                               print(ggplot(df, aes(index, value)) + geom_line(aes(colour = series)))
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
                           get_initial_stock_after_stabilization=function(mat=adjacency_matrix,initial_input=1.0,frequency=2L){
                               E=operate_adjacency_matrix(mat=mat,initial_input = initial_input,frequency = 1,iterations = 1)[1:core_nodes];
                               output = cbind(E%*%get_leontief(mat = get_mat_nth_power(mat = mat,n=frequency),in_or_out_normalized = 'out',renormalize_rows = FALSE,input_multiplier = initial_input),0,0);
                               output[[core_nodes+2]]=(output%*%get_mat_nth_power(mat = mat,n=frequency,keep_in_out = FALSE))[[core_nodes+2]];
                               return(output)
                           },
                           get_series_of_node_stock_after_stabilization=function(mat=adjacency_matrix,initial_input=1.0,frequency=2L,do_plot=FALSE){
                               out=vector(mode='list',length=frequency);
                               E=operate_adjacency_matrix(mat=mat,initial_input = initial_input,frequency = 1,iterations = 1)[1:core_nodes];
                               out[[1]]=cbind(E%*%get_leontief(mat = get_mat_nth_power(mat = mat,n=frequency),in_or_out_normalized = 'out',renormalize_rows = FALSE,input_multiplier = initial_input),0,0);
                               for(n in 2:(frequency)){
                                   out[[n]]=(out[[n-1]]%*%mat);
                               }
                               out[[1]][[core_nodes+2]]=(out[[frequency]]%*%mat)[[core_nodes+2]];
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
                           get_quantity_matrix=function(mat=adjacency_matrix,multiplier=1000){
                               X = get_initial_stock_after_stabilization(mat=mat,initial_input=multiplier,frequency=1L);#operate_adjacency_matrix(mat=mat,initial_input=multiplier,frequency=1L,iterations=iterations);
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
                           get_TST=function(mat=adjacency_matrix,input_multiplier=1000,consider_input=FALSE,consider_output=FALSE,proportional_to_input=FALSE,return_quatity_matrix=FALSE){
                               final_node_in = core_nodes;
                               final_node_out = core_nodes;
                               if(consider_input)
                                   final_node_in = core_nodes + 1
                               if(consider_output)
                                   final_node_out = core_nodes + 2
                               m = get_quantity_matrix(mat=mat,multiplier=input_multiplier);
                               total = sum(m[1:final_node_in,1:final_node_out])#sum(m[which(m>0)])
                               if(proportional_to_input){
                                   if(return_quatity_matrix){
                                       return(list(TST=total/input_multiplier,QM=m))
                                   }
                                   return(total/input_multiplier)
                               }
                               if(return_quatity_matrix)
                                   return(list(TST=total,QM=m))
                               return(total)
                           },
                           navigate_net=function(mat=adjacency_matrix,save_trajectories=FALSE,save_cycles=FALSE,return_data=TRUE,check_pre_data=TRUE,pre_net_data=net_data,do_parallel=FALSE,limit_for_parallel=1000){
                               if(check_pre_data
                                  && !is.null(pre_net_data$cycle_participation)
                                  && !is.null(pre_net_data$trajectory_participation)
                                  && !is.null(pre_net_data$binary_mat)
                                  && check_is_mat_bin_same_as_mat_bin(mat,pre_net_data$binary_mat)){
                                   if(return_data){
                                       return(pre_net_data);
                                   }
                                   else{
                                       return();
                                   }
                               }
                               system.time(res<-paths(mat,core_nodes+1,core_nodes+2,save_trajectories,use_parallel = do_parallel,limit_for_parallel=limit_for_parallel))
                               if(check_is_mat_same_as_adjacency(mat)){
                                   net_data$cycle_participation <<- t(matrix(unlist(res[[4]]), ncol = core_nodes+2, nrow = core_nodes+2));
                                   net_data$trajectory_participation <<- t(matrix(unlist(res[[2]]), ncol = core_nodes+2, nrow = core_nodes+2));
                                   net_data$binary_mat <<- get_binary_adjacency();
                                   if (save_trajectories){
                                       net_data$trajectories <<- res[[1]];
                                   }
                                   if(save_cycles){
                                       net_data$cycles <<- res[[3]];
                                   }
                                   if(return_data)
                                       return(net_data)
                               }else{
                                   output=list()
                                   output$cycle_participation = t(matrix(unlist(res[[4]]), ncol = core_nodes+2, nrow = core_nodes+2));
                                   output$trajectory_participation = t(matrix(unlist(res[[2]]), ncol = core_nodes+2, nrow = core_nodes+2));
                                   if(save_trajectories || save_cycles){
                                       output$binary_mat = get_binary_adjacency(mat=mat);
                                   }
                                   if (save_trajectories){
                                       output$trajectories = res[[1]];
                                   }
                                   if(save_cycles){
                                       output$cycles = res[[3]];
                                   }
                                   if(return_data)
                                       return(output)
                               }
                           },
                           get_fast_cycles=function(mat=adjacency_matrix,save_cycles=TRUE,return_data=TRUE,check_pre_data=TRUE,pre_net_data=net_data,do_parallel=FALSE){
                               if(check_pre_data
                                        #&& !is.null(pre_net_data$cycle_participation)
                                  && !is.null(pre_net_data$binary_mat)
                                  && check_is_mat_bin_same_as_mat_bin(mat,pre_net_data$binary_mat)){
                                   if(return_data){
                                       return(pre_net_data);
                                   }
                                   else{
                                       return();
                                   }
                               }
                                        #system.time(res<-fast_cycles(get_core(mat=get_binary_adjacency(mat=mat))))#,use_parallel = do_parallel))
                               res<-fast_cycles(m = get_core(mat=mat),save_cycles = save_cycles);
                               if(check_is_mat_same_as_adjacency(mat)){
                                        #net_data$cycle_participation <<- t(matrix(unlist(res[[2]]), ncol = core_nodes+2, nrow = core_nodes+2));
                                   net_data$binary_mat <<- get_binary_adjacency();
                                   net_data$cycles <<- res;
                                   if(return_data)
                                       return(net_data)
                               }else{
                                   output=list()
                                        #output$cycle_participation = t(matrix(unlist(res[[2]]), ncol = core_nodes+2, nrow = core_nodes+2));
                                   output$binary_mat = get_binary_adjacency(mat=mat);
                                   output$cycles = res;
                                   if(return_data)
                                       return(output)
                               }
                           },
                           get_cycles_flow=function(mat=adjacency_matrix,x=vector(mode='numeric'),cycles=list(),initial_input=1000.0){
                               cycles_to_use = cycles;
                               X = x;
                               if(length(x)==0){
                                   X = as.vector(operate_adjacency_matrix(mat = mat,initial_input = initial_input));
                               }
                               if(length(cycles)==0){
                                   is_adjacency = check_is_mat_same_as_adjacency(mat);
                                   if(is_adjacency){
                                       if(is.null(net_data$cycles)){
                                           netdata=navigate_net(save_trajectories = TRUE);
                                       }else{
                                           netdata=net_data;
                                       }
                                   }else{
                                       netdata = navigate_net(mat = mat,save_trajectories = TRUE);
                                   }
                                   
                                   cycles_to_use = netdata$cycles;
                               }
                               
                               if(length(cycles_to_use)==0) {
                                   return(list())
                               }
                               link_flow = matrix(data=0,nrow=nrow(mat),ncol=ncol(mat));
                               output = vector(mode='list',length=length(cycles_to_use));
                               
                               for(k in 1:length(cycles_to_use)){
                                   current_element = list();
                                   n_nodes = length(cycles_to_use[[k]]); 
                                   current_element$cycle_length=n_nodes;
                                   mult = 1.0;
                                   mult_dummy=vector(mode="numeric",length=n_nodes);
                                   mult_sd=vector(mode="numeric",length=n_nodes);
                                   for(i in 1:(n_nodes-1)){
                                       link=mat[cycles_to_use[[k]][[i]],cycles_to_use[[k]][[i+1]]]
                                       mult = mult * link;
                                       mult_dummy[i]=link;
                                   }
                                   link = mat[cycles_to_use[[k]][[n_nodes]],cycles_to_use[[k]][[1]]];
                                   mult = mult * link;
                                   mult_norm = mult^(1/n_nodes);
                                   mult_dummy[n_nodes]=link;
                                   mult_sd=sd(mult_dummy)
                                   
                                   d2_mult = 0;
                                   for(i in 1:(length(cycles_to_use[[k]])-1)){
                                       d2_mult = d2_mult + (mult_dummy[[i]]-mult_norm)^2;
                                   }
                                   d2_mult = d2_mult + (mult_dummy[[n_nodes]]-mult_norm)^2;
                                   d2_mult = d2_mult/n_nodes;
                                   
                                   stock = sum(X[cycles_to_use[[k]]]);
                                   stock_sd=sd(X[cycles_to_use[[k]]]);
                                   cycle_val = stock*mult;
                                   cycle_val_per_link = cycle_val/n_nodes;
                                   for(i in 1:(n_nodes-1)){
                                       link_flow[cycles_to_use[[k]][[i]],cycles_to_use[[k]][[i+1]]] = cycle_val_per_link + link_flow[cycles_to_use[[k]][[i]],cycles_to_use[[k]][[i+1]]]
                                   }
                                   link_flow[cycles_to_use[[k]][[n_nodes]],cycles_to_use[[k]][[1]]] = cycle_val_per_link + link_flow[cycles_to_use[[k]][[n_nodes]],cycles_to_use[[k]][[1]]];
                                   
                                   current_element$cycle_flow=cycle_val;
                                   current_element$cycle_stock_sd=stock_sd;
                                   current_element$cycle_mult_sd=mult_sd;
                                   current_element$cycle_mult_av=mean(mult_dummy)
                                   current_element$cycle_stock_av=stock/n_nodes;
                                   current_element$cycle_adjacency=mult_dummy;
                                   current_element$cycle_stock=stock;
                                   current_element$cycle_mult=mult;
                                   current_element$cycle_inhom = d2_mult;
                                   current_element$cycle_mult_norm = mult_norm;
                                   current_element$cycle_flownorm=cycle_val_per_link;
                                   output[[k]] = current_element;
                               }
                               
                               return(list(cycles_flow=output,links_flow=link_flow));
                           },
                           get_fast_cycles_flow=function(mat=adjacency_matrix,x=vector(mode='numeric'),cycles=list(),initial_input=1000.0,save_cycle_flow_detail = TRUE){
                               cycles_to_use = cycles;
                               X = x;
                               if(length(x)==0){
                                   X = get_initial_stock_after_stabilization(mat = mat,initial_input = initial_input,frequency = 1L)
                               }
                               if(length(cycles)==0){
                                   is_adjacency = check_is_mat_same_as_adjacency(mat);
                                   if(is_adjacency){
                                       if(is.null(net_data$cycles)){
                                           netdata=get_fast_cycles(save_cycles = TRUE);
                                       }else{
                                           netdata=net_data;
                                       }
                                   }else{
                                       netdata = get_fast_cycles(mat = mat,save_cycles = TRUE);
                                   }
                                   
                                   cycles_to_use = netdata$cycles$cycles;
                               }
                               
                               if(length(cycles_to_use)==0) {
                                   return(list())
                               }
                               res = fast_cycle_flow(m = mat,x = X,cycles = cycles_to_use,save_cycle_flows = save_cycle_flow_detail);
                               return(res);
                           },
                           get_fast_cycles_and_flow=function(mat=adjacency_matrix,x=vector(mode='numeric'),initial_input=1000.0,save_cycle_flow_detail = TRUE){
                               X = x;
                               if(length(x)==0){
                                   X = get_initial_stock_after_stabilization(mat = mat,initial_input = initial_input,frequency = 1L)
                               }
                               res = fast_cycles(m = get_core(mat),save_cycles=save_cycle_flow_detail,calculate_flow = TRUE,x = X);
                               return(res);
                           },
                           get_link_cycle_flow_entropy=function(links_flow,mat=adjacency_matrix,multiplier=1000.0){
                               qm_full = get_quantity_matrix(mat = mat,multiplier = multiplier);
                               elements = which(links_flow>0);
                               lf = links_flow[elements];
                               qm_nz = qm_full[elements];
                               lf_p = lf/qm_nz;
                               prob_per_link = links_flow;
                               prob_per_link[elements]=lf_p;
                               entr_per_link = links_flow;
                               entr = -(lf_p*log2(lf_p)+(1.0-lf_p)*log2(1.0-lf_p));
                               entr_per_link[elements]=entr;
                               ent_tot = sum(entr);
                               ent_av = ent_tot/length(which(qm_full[1:core_nodes,1:core_nodes]>0));
                               return(list(probs_per_link=prob_per_link,entropy_per_link=entr_per_link,entropy_tot = ent_tot,entropy_mean = ent_av));
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
                           },
                           get_leontief=function(mat=adjacency_matrix,in_or_out_normalized='in',renormalize_cols = TRUE, renormalize_rows = TRUE,input_multiplier=1000){
                               if(in_or_out_normalized=='in'){
                                   nCols = mat;
                                   if(renormalize_cols)
                                       nCols = normalize_cols(mat = get_quantity_matrix(mat = mat, multiplier = input_multiplier));
                                   return(solve(diag(core_nodes)-nCols[1:core_nodes,1:core_nodes],tol=0));
                               }else{
                                   if(renormalize_rows){
                                       return(solve(diag(core_nodes)-normalize_rows(mat=mat,set_adjacency = FALSE)[1:core_nodes,1:core_nodes],tol=0)); 
                                   }else{
                                       return(solve(diag(core_nodes)-mat[1:core_nodes,1:core_nodes],tol=0)); 
                                   }
                               }
                           },
                           get_finn=function(mat=adjacency_matrix,input_multiplier=1000,renormalize_cols = TRUE,consider_input=FALSE,consider_output=FALSE){
                               L = get_leontief(mat = mat,in_or_out_normalized = 'in',renormalize_cols = renormalize_cols, input_multiplier = input_multiplier);
                               TSTQM = get_TST(mat = mat, input_multiplier = input_multiplier,consider_input = consider_input, consider_output = consider_output,return_quatity_matrix = TRUE);
                               TST = TSTQM$TST;
                               QM = TSTQM$QM;
                               Si = apply(QM,2,sum);
                               out = 0;
                               for(i in 1:core_nodes){
                                   out = out + Si[[i]]/TST*(L[[i,i]]-1.0)/L[[i,i]];
                               }
                               return(out);
                           },
                           get_finn_per_node=function(mat=adjacency_matrix,input_multiplier=1000,renormalize_cols = TRUE,consider_input=FALSE,consider_output=FALSE){
                               L = get_leontief(mat = mat,in_or_out_normalized = 'in',renormalize_cols = renormalize_cols, input_multiplier = input_multiplier);
                               TSTQM = get_TST(mat = mat, input_multiplier = input_multiplier,consider_input = consider_input, consider_output = consider_output,return_quatity_matrix = TRUE);
                               TST = TSTQM$TST;
                               QM = TSTQM$QM;
                               Si = apply(QM,2,sum);
                               finn = 0;
                               out = vector(mode='numeric',length=core_nodes);
                               for(i in 1:core_nodes){
                                   out[i] = Si[[i]]/TST*(L[[i,i]]-1.0)/L[[i,i]];
                                        #finn = finn + out[i];
                               }
                                        #if(finn==0)
                                        # return(rep(0, core_nodes))
                               return(out);
                           },
                           get_finn2=function(mat=adjacency_matrix,input_multiplier=1000,renormalize_cols = TRUE,renormalize_rows=FALSE,consider_input=FALSE,consider_output=FALSE){
                               L = get_leontief(mat = mat,in_or_out_normalized = 'out',renormalize_cols = renormalize_cols,renormalize_rows = renormalize_rows, input_multiplier = input_multiplier);
                               TSTQM = get_TST(mat = mat, input_multiplier = input_multiplier,consider_input = consider_input, consider_output = consider_output,return_quatity_matrix = TRUE);
                               TST = TSTQM$TST;
                               QM = TSTQM$QM;
                               Si = apply(QM,1,sum);
                               out = 0;
                               for(i in 1:core_nodes){
                                   out = out + Si[[i]]/TST*(L[[i,i]]-1.0)/L[[i,i]];
                               }
                               return(out);
                           },
                           get_mat_nth_power=function(mat=adjacency_matrix,n=2,keep_in_out=TRUE){
                               output = mat;
                               iters = n-1;
                               while(iters>0){
                                   output = output %*% mat;
                                   iters=iters-1;
                               }
                               if(keep_in_out){
                                   output[,(core_nodes+1):(core_nodes+2)]=mat[,(core_nodes+1):(core_nodes+2)];
                                   output[(core_nodes+1):(core_nodes+2),]=mat[(core_nodes+1):(core_nodes+2),];
                               }
                               return(output)
                           },
                           get_eigen=function(mat=adjacency_matrix){
                               return(eigen(get_core(mat=mat)));
                           },
                           get_eigen_mod=function(mat=adjacency_matrix){
                               ev = get_eigen(mat=mat)$values;
                               return(sqrt(Re(Conj(ev)*ev)));
                           },
                           get_eigen_max=function(mat=adjacency_matrix){
                               ev = get_eigen(mat=mat)$values;
                               return(sqrt(max(Re(Conj(ev)*ev))));
                           },
                           get_exp_factor=function(mat=adjacency_matrix){
                               return(log(get_eigen_max(mat=mat)));
                           },
                           get_general_entropy=function(data=c(1,1,1),logbase=2){
                             if(sum(data)!=1.0){
                               data=data/sum(data)
                             }
                             n=length(data)
                             entr=vector(mode='numeric',length=n)
                             for (i in 1:n){
                               entr[i]=data[i]*log(data[i],base=logbase)
                             }
                             out=sum(entr)*-1
                             return(out)
                           },
                           get_input_output_entropy=function(mat=adjacency_matrix,input_multiplier=1000L,logbase=2){
                             A=get_quantity_matrix(mat = mat,multiplier = input_multiplier)
                             L=list()
                             in_links=which(A[core_nodes+1,]>0)
                             if(length(in_links)==1) {in_storages=sum(A[,in_links])}
                             else {in_storages=colSums(A[,in_links])}
                             L$in_entropy=get_general_entropy(in_storages,logbase=logbase)
                             out_links=A[which(A[,core_nodes+2]>0),core_nodes+2]
                             L$out_entropy=get_general_entropy(out_links,logbase=logbase)
                             return(L)
                           }
                       )
                       )

NetEnsemble <- setRefClass("NetEnsemble",
                           fields = list(ensemble_sizeMin = "numeric",
                                         ensemble_sizeMax = "numeric", 
                                         size_interval = "numeric", 
                                         size_replicate = "numeric",
                                         nodes_autoloops_allowed = "logical",
                                         random_seed = "numeric",
                                         nodes_connect_function = "function",
                                         respiration = "logical",
                                         make_flowing = "logical",
                                         minimal_flowing = "logical",
                                         random_weights = "logical",
                                         normalized_rows = "logical",
                                         fully_connected = "logical",
                                         node_seq = "vector",
                                         generated_networks = "list",
                                         name = 'character'
                                         ),
                           methods = list(
                               initialize = function(...,
                                                     ensemble_sizeMin = 2L,
                                                     ensemble_sizeMax = 10L,
                                                     size_interval = 1L,
                                                     size_replicate = 10L,
                                                     random_seed = 55L,
                                                     name = paste0(ensemble_sizeMin,'_',ensemble_sizeMax,'_',
                                                                   size_replicate,'_',size_interval)
                                                     ){
                                   ensemble_sizeMin <<- ensemble_sizeMin;
                                   ensemble_sizeMax <<- ensemble_sizeMax;
                                   size_replicate <<- size_replicate;
                                   size_interval <<- size_interval;
                                   name <<- name;
                                   node_seq <<- as.integer(seq(ensemble_sizeMin,ensemble_sizeMax,size_interval));
                                   callSuper(...)
                               },
                               generate = function(nodes_connect_function = NetBasicNodesConnectFunction,
                                                   nodes_autoloops_allowed = FALSE,
                                                   respiration=FALSE,
                                                   make_flowing = TRUE,
                                                   minimal_flowing = FALSE,
                                                   random_weights = TRUE,
                                                   normalized_rows = TRUE,
                                                   fully_connected = FALSE){
                                   num_total = length(node_seq) * size_replicate;
                                   output = vector(mode="list",length = num_total);
                                   nodes_connect_function <<- nodes_connect_function;
                                   respiration     <<- respiration;
                                   make_flowing    <<- make_flowing; 
                                   minimal_flowing <<- minimal_flowing;
                                   random_weights  <<- random_weights;
                                   normalized_rows <<- normalized_rows;
                                   fully_connected <<- fully_connected;
                                   nodes_autoloops_allowed <<- nodes_autoloops_allowed;
                                   num_generated=1;
                                   for(i in seq(ensemble_sizeMin,ensemble_sizeMax,size_interval)){
                                       for(j in 1:size_replicate){
                                           current = NetBase$new(
                                                                 core_nodes = i,
                                                                 allow_autoloop = nodes_autoloops_allowed,
                                                                 respiration = respiration,
                                                                 nodes_connect_function = nodes_connect_function,
                                                                 make_flowing = make_flowing,
                                                                 minimal_flowing = minimal_flowing,
                                                                 random_weights = random_weights,
                                                                 normalized_rows = normalized_rows,
                                                                 fully_connected = fully_connected,
                                                                 name = paste0(name,"-",num_generated)
                                                             );
                                           output[num_generated] = current;
                                           num_generated = num_generated + 1;
                                       }
                                   }
                                   generated_networks <<-output;
                               },
                               generate_from_ensemble=function(pre_ensemble=NetEnsemble$new(),adjacencies=list()){
                                   ensemble_sizeMin <<- pre_ensemble$ensemble_sizeMin;
                                   ensemble_sizeMax <<- pre_ensemble$ensemble_sizeMax;
                                   size_replicate <<- pre_ensemble$size_replicate;
                                   node_seq <<- pre_ensemble$node_seq;
                                   nodes_connect_function <<- pre_ensemble$nodes_connect_function;
                                   respiration     <<- pre_ensemble$respiration;
                                   make_flowing    <<- pre_ensemble$make_flowing; 
                                   minimal_flowing <<- pre_ensemble$minimal_flowing;
                                   random_weights  <<- pre_ensemble$random_weights;
                                   normalized_rows <<- pre_ensemble$normalized_rows;
                                   fully_connected <<- pre_ensemble$fully_connected;
                                   nodes_autoloops_allowed <<- pre_ensemble$nodes_autoloops_allowed;
                                   
                                   if(length(adjacencies)==0 || (length(adjacencies)!=length(pre_ensemble$generated_networks))){
                                       generated_networks <<-pre_ensemble$generated_networks;
                                   }
                                   else{
                                       num_total = length(node_seq);
                                       output = vector(mode="list",length = num_total);
                                       
                                       for(i in 1:length(pre_ensemble$generated_networks)){
                                           current = NetBase$new()
                                           current$create_from_base_and_adjacency(pre_ensemble$generated_networks[[i]],adjacencies[[i]]);
                                           output[[i]] = current;
                                       }
                                       generated_networks <<- output;
                                   }
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
                                     if(PRINT_FLOWNET_DEBUG){
                                      system(sprintf('echo "\n%s\n"', paste0('Iteration num for net ',initial_net$name,' with ',initial_net$core_nodes,' : ',iteration, collapse="")));
                                     }
                                     mutated_net = mutation_func(initial_net,current_net);
                                     tryCatch(
                                       expr = {
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
                                       },
                                       finally = {
                                         iteration = iteration + 1;
                                       })
                                        #if(!(result$continue))
                                        #	break;
                                 }
                             }
                         )
                         )

EnsembleEvolve <- setRefClass("EnsembleEvolve",
                              fields = list(ensemble = "NetEnsemble",
                                            mutation_func = "function",
                                            characteristic_funcs = "list",
                                            evolution_func = "function",
                                            evolved_nets = "list"
                                            ),
                              methods = list(
                                  initialize = function(...,
                                                        ensemble = NetEnsemble$new(),
                                                        mutation_func = function(base,mat){
                                                            
                                                        },
                                                        characteristic_funcs = list(function(base,mat){return(1)}),
                                                        evolution_func = function(iteration_num,base,current_state,mutated_net,characteristic_vals){
                                                            return(list(new_state = current_state,
                                                                        new_net = current_net,
                                                                        continue = FALSE));
                                                        }){
                                      ensemble <<- ensemble;
                                      mutation_func <<- mutation_func;
                                      evolution_func <<- evolution_func;
                                      characteristic_funcs <<- characteristic_funcs;
                                      evolved_nets <<- vector(mode="list",length=length(ensemble$generated_networks));
                                      callSuper(...);
                                      if(length(ensemble$generated_networks)>0){
                                          for(i in 1:length(ensemble$generated_networks)){
                                              evolved_nets[[i]] <<- NetEvolve$new(initial_net=ensemble$generated_networks[[i]],
                                                                                  characteristic_funcs=characteristic_funcs,
                                                                                  mutation_func=mutation_func,
                                                                                  evolution_func=evolution_func);
                                          }
                                      }
                                  },
                                  evolve_ensemble = function(){
                                      for(i in 1:length(evolved_nets)){
                                          print(i)
                                          evolved_nets[[i]]$evolve();
                                          cat('Ensemble evolve',i,'of',length(evolved_nets),'\n')
                                      }
                                  },
                                  evolve_ensemble_par = function(num_cores=4){
                                        #lapply(evolved_nets,function(x){x$evolve()});
                                      
                                        #counter=1;
                                      res<-mclapply(evolved_nets,function(x){
                                        # system(sprintf('echo "\n%s\n"', paste0("generated network",":", counter,collapse = "")));
                                        #  counter <<- counter + 1;
                                          x$evolve();
                                          return(x);
                                      }
                                     ,mc.cores=num_cores);
                                      evolved_nets<<-res;
                                  },
                                  get_final_ensemble = function(){
                                      output = vector(mode="list",length=length(evolved_nets));
                                      for(i in 1:length(evolved_nets)){
                                          output[[i]]=evolved_nets[[i]]$current_net;
                                      }
                                      output_ens = NetEnsemble$new()
                                      output_ens$generate_from_ensemble(pre_ensemble = ensemble,adjacencies = output)
                                      return(output_ens);
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
                                 initialize = function(...,results_list,variables){
                                     
                                     nonscalar=c();
                                     for (i in 1:length(results_list[[1]][[1]])){
                                         if(!length(results_list[[1]][[1]][[i]])==1) nonscalar=c(nonscalar,i)
                                     };
                                    
                                     if(is.null(nonscalar)) {variables=variables}
                                     else {variables=variables[-nonscalar]}
                                     
                                     to_plot=vector(mode='list',length(results_list));
                                     for (k in 1:length(results_list)){
                                         columns=list();
                                         for (i in 1:length(variables)){
                                             columns[[i]]=unlist(map(results_list[[k]],variables[[i]]));
                                         };
                                         columns=set_names(columns,variables)
                                         to_plot[[k]]=columns;
                                     };
                                     
                                     data<<-lapply(to_plot,function(x){
                                         if(length(x)==0) return(NULL)
                                         else return(
                                                  as.data.frame(cbind(matrix(unlist(x),
                                                                             nrow=length(x[[1]]),
                                                                             dimnames=list(NULL,as.list(names(to_plot[[1]])))))
                                                                ))
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
                                                  cex.lab=1.5,
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
                                 plot_me2=function(x_var,
                                                  y_var,
                                                  with_raster=TRUE,
                                                  x_func=function(x) return(x),
                                                  y_func=function(y) return(y),
                                                  col_var = NULL,
                                                  size_var = NULL,
                                                  color_func=function(x){return(hsv(h = (x-min_max_col$min)/(min_max_col$max-min_max_col$min)*0.8, s = 1, v = 1, 1))},
                                                  size_func=function(x){return(0.5+(x-min_max_size$min)/(min_max_size$max-min_max_size$min))},
                                                  x_lim=NULL,
                                                  y_lim=NULL){
                                   
                                   x_apply= function(dat,x) {for (i in 1:length(dat)) dat[[i]][x]=x_func(dat[[i]][x]);return(dat)}
                                   y_apply= function(dat,y) {for (i in 1:length(dat)) dat[[i]][y]=y_func(dat[[i]][y]);return(dat)}
                                   
                                   dat=x_apply(data,x_var)
                                   dat=y_apply(dat,y_var)
                                   
                                   if(!is.null(col_var)){
                                     min_max_col = get_min_max(dat,col_var);
                                   }
                                   if(!is.null(size_var)){
                                     min_max_size = get_min_max(dat,size_var);
                                   }
                                   
                                  
                                   min_max_x = x_lim;
                                   min_max_y = y_lim;
                                   if(is.null(x_lim)){
                                     min_max_x = get_min_max(dat,x_var);
                                   }
                                   if(is.null(y_lim)){
                                     min_max_y = get_min_max(dat,y_var);
                                   }
                                   
                                   for(i in 1:length(dat)){
                                     col = 'black';
                                     size = 0.5;
                                     if(!is.null(col_var)){
                                       col = unlist(lapply(dat[[i]][[col_var]],color_func));
                                     }
                                     if(!is.null(size_var)){
                                       size = unlist(lapply(dat[[i]][[size_var]],size_func));
                                     }
                                     if(i==1){
                                       if(!is.null(col_var)){
                                           if(with_raster){
                                           layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
                                           par(mar = c(5.1, 4.1, 4.1, 2.1));
                                           }
                                        };
                                       plot(dat[[i]][[x_var]],dat[[i]][[y_var]], 
                                            col = col,
                                            cex = size, 
                                            xlab = x_var, 
                                            ylab = y_var,
                                            cex.lab=1.3,
                                            xlim=c(min_max_x$min,min_max_x$max),
                                            ylim=c(min_max_y$min,min_max_y$max))
                                     }
                                     else{
                                       points(dat[[i]][[x_var]],
                                              dat[[i]][[y_var]], 
                                              col = col,
                                              cex = size)
                                     }
                                     line_col = col;
                                     if(!is.null(col_var)){
                                       line_col = col[-1L];
                                     }
                                     segments(dat[[i]][[x_var]][-length(dat[[i]][[x_var]])],
                                              dat[[i]][[y_var]][-length(dat[[i]][[y_var]])],
                                              dat[[i]][[x_var]][-1L],dat[[i]][[y_var]][-1L],
                                              col=line_col);
                                   }
                                   if(!is.null(col_var)){
                                     if(with_raster){
                                     col_seq_range = seq(min_max_col$min,min_max_col$max,l=5);
                                     legend_image <- as.raster(matrix(unlist(lapply(col_seq_range,color_func)), ncol=1))
                                     par(mar = c(2, 0, 2, 0));
                                     plot(c(0,10),c(min_max_col$min,min_max_col$max),type = 'n', axes = F,xlab = '', ylab = '', main = col_var)
                                     text(x=5, y = col_seq_range, labels =  format(col_seq_range, scientific = TRUE, digits = 3));
                                     rasterImage(legend_image, 0, min_max_col$max, 2,min_max_col$min)
                                     par(mar = c(5.1, 4.1, 4.1, 2.1))
                                     }
                                   }
                                 },
                                 get_min_max=function(dfs,variable){
                                     maximum=vector(mode="numeric",length=length(dfs))
                                     minimum=vector(mode="numeric",length=length(dfs))
                                        for(i in 1:length(dfs)){
                                            infinite=which(is.infinite(as.matrix(dfs[[i]][variable])))
                                            if(length(infinite)>0) dfs[[i]][variable][infinite]=NA
                                            maximum[i]=max(dfs[[i]][variable],na.rm=TRUE)
                                            minimum[i]=min(dfs[[i]][variable],na.rm = TRUE)
                                        }
                                     y_max=max(maximum,na.rm = TRUE)
                                     y_min=min(minimum,na.rm=TRUE)
                                     output=list();
                                     output$max=y_max;
                                     output$min=y_min;
                                     return(output);
                                 }
                             )
                    )

