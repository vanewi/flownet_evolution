source("FlowNet.R")

netest = NetEnsamble$new(ensemble_sizeMin = 5L,ensemble_sizeMax = 14L,size_replicate = 30L,nodes_autoloops_allowed = FALSE);
netest$generate();
total_nets = length(netest$generated_networks);

evolved = vector(mode='list',length=total_nets);
evolved_info = vector(mode='list',length=total_nets);

temp = NetBase$new()

use_only_existing_links = FALSE
cumulative_evolve = TRUE

#TO DO
#tiempo continuo ?
#ami, fin index. 
#ciclo

#entropia: por ciclo, 
#distribucion de energia por ciclos vs energia por flujo
#proporcion nodos en ciclos vs trayectoria
#trayectorias en la redes vs. ciclos ()
#funcion: componentes desconectados de la red. deberiamos arreglar esas patologias con otra funcion. esta funcion esta asociada a la conexion de modulos entre redes tb

#pulsos: tiempo de independencia de nodos de los pulsos y su dinamica de ritmos de entrada (incluyendo random) e interacciones de ciclos
#cuanto tiempo puede vivir un nodo sin nada? 

#Leontief - structure mat

# grandes numeros para la teoria. 

#LISTOS
#indegree y outdegree, falta ver que tipo de analisis hacerle a esos valores

for(i in 1:total_nets){
	base = (netest$generated_networks)[[i]];
	evolved_info[[i]] = list(-1);
	start_mean_time = base$mean_time_discrete(input_percentage = 0.01);
	current_mean_time = 0;
	reps = 0;
	while(start_mean_time>current_mean_time && reps<500){
		if((cumulative_evolve && reps==0) || !cumulative_evolve)
			evolved[[i]]=base$adjacency_matrix;
		link_to_change = 1;
		if(use_only_existing_links){
			over_0 = which((evolved[[i]][1:base$core_nodes,1:base$core_nodes])>0);
			if(length(over_0)==0)
				break;
			link_to_change = sample(over_0,1);
		}else{
			link_to_change = sample(base$core_nodes * base$core_nodes,1);
		}
		
		from = link_to_change%%(base$core_nodes);
		if(from==0){from = base$core_nodes;}
		to = ceiling(link_to_change/(base$core_nodes));
		
		evolved[[i]][[from,to]]=runif(1,min=0.00001,max=0.99999);
		evolved[[i]]=base$normalize_rows(mat=evolved[[i]],set_adjacency=FALSE);
		current_mean_time = base$mean_time_discrete(input_percentage = 0.01,mat=evolved[[i]]);
		reps = reps + 1;
		if(current_mean_time>start_mean_time){
			evolved_info[[i]] = list(reps,from,to,current_mean_time,start_mean_time,base$core_nodes);
			print("EN :")
			print(i);
			print(reps);
			next;
		}
	}
}
