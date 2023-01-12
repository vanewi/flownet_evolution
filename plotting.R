source("FlowNet.R")
#-------------------
#fig. comparing orientors between evolution scenarios.

input_output_entropy=input_output_entropy_generator(log_base=2);
entropy_diff=entropy_diff_generator(log_base=2);

orientors=list(TST=tst,TSS=storage,ASC=ascendency,Tau=residence_time,AMI=ami,Finn=finn,EDiff=entropy_diff,b=b);
plotting=c('tst','storage','residence_time','ascendency','ami','finn','entropy_diff','b');
ylims_max=c(NULL,NULL,NULL,NULL,NULL,NULL,NULL)
ylims_min=c(NULL,NULL,NULL,NULL,NULL,NULL,NULL)
transforming=list((function(x) return(log(x))),
              (function(x) return(log(x))),
              (function(x) return(log(x))),
              (function(x) return(log(x))),
              (function(x) return(x)),
              (function(x) return(x)),
              (function(x) return(x)),
              (function(x) return(x))
              )
node_ensembles=as.integer(c(15))#,30,60,100));

par(mfrow=c(length(orientors),length(orientors)),cex.axis=0.8)#, mai=c(0.1,0.1,0.2,0.1))
for (i0 in 1:length(orientors)){
  print(paste0('orientors_x :',names(orientors[i0])))
  open_name=opening(nodes = node_ensembles[1],
                    replicates = 100,
                    interval = 1,
                    orientor =names(orientors[i0]),
                    greed_experiment = 'part_stock',
                    picked_percent =0.3 
                    )
  experiment <- readRDS(file=open_name)

  results_control=experiment$results_control;
  n_nets=length(results_control[[1]])
  variables=names(results_control[[1]][[1]]$state$extra$history[[1]])
  rm(experiment);
  gc()

  all_control=vector(mode='list',length=n_nets)
  for (cp in 1:length(results_control)){
    for (net in 1:n_nets){
        all_control[[net]]=c(all_control[[net]],results_control[[cp]][[net]]$state$extra$history)
    }
  }
  
  for(j in 1:length(orientors)){
  vis=Visualization$new(results_list=all_control,variables=variables)
   if (i0==j){
     plot(x = 0:5, y = 0:5, ann = F,bty = "n",type = "n",
          xaxt = "n", yaxt = "n")
     text(x=2,y=2.5,plotting[j])
   }
   else{
     
     par(mar=c(1.9, 1.9, 1, 1))
     vis$plot_me2(x_var = plotting[i0],
                  y_var = plotting[j],
                  col_var=NULL,
                  x_func = transforming[[i0]],
                  y_func = transforming[[j]]
                  #y_lim=list(min=ylims_min[j],max=ylims_max[j]), 
                  )
     
   }
    print(paste0('orientors_y :',names(orientors[j])))
  }
}

##--------------------------------------
out=vector(mode='list',length=length(results_list))
for (i in 1:length(results_list)){
  end=length(results_list[[i]])
  out[[i]]$cflown=results_list[[i]][[end]]$c_flow_per_node
  out[[i]]$cflow_tot=results_list[[i]][[end]]$c_analysis$total_cycle_flow
  out[[i]]$cflown_perc=(out[[i]]$cflown)/(out[[i]]$cflow_tot)*100
  out[[i]]$num_nodes=results_list[[i]][[end]]$num_nodes
  out[[i]]$finnn=results_list[[i]][[end]]$finn_node
  out[[i]]$finn=results_list[[i]][[end]]$finn
  out[[i]]$finnn_perc=(out[[i]]$finnn)/(out[[i]]$finn)*100
  out[[i]]$storage=results_list[[i]][[end]]$storage
  out[[i]]$cflown_stock_norm=out[[i]]$cflown/out[[i]]$storage
  out[[i]]$cflow_times_finn=log(results_list[[i]][[end]]$finn/(1.0-results_list[[i]][[end]]$finn+0.000000000001))
}


x=unlist(map(out,'finnn'))
y=unlist(map(out,'cflown_stock_norm'))
num_nodes=unlist(map(out,'num_nodes'))

color_func=function(x){return(hsv(h = (x-min(num_nodes))/(max(num_nodes)-min(num_nodes))*0.8, s = 1, v = 1, 1))}
col=color_func(num_nodes)

par(mfrow=c(1,1))
plot(x,y,col=col,ylab='Normalized Cycles Flow (per Total Storage) per node',xlab='Finn Index (per node)',cex.lab=1.3)
legend('topleft',inset=0.05,legend=c("10","15",'20'), fill=c(col[1],col[105],col[205]),title = 'core nodes',bty='n')

#----------------------------------------
#fig. orientors relations for natural history evolution (changes in gaining/loosing links)

open_name=opening(nodes = node_ensembles[1],
                  replicates = 100,
                  interval = 1,
                  orientor = 'TST', #names(orientors[i0]),
                  greed_experiment = 'part_stock',
                  picked_percent =0.3 
)
experiment <- readRDS(file=open_name)

results=experiment$results;
n_nets=length(results)
variables=names(results[[1]]$state$extra$history[[1]])
rm(experiment);
gc()

all_control=vector(mode='list',length=n_nets)
  for (net in 1:n_nets){
    all_control[[net]]=c(all_control[[net]],results[[net]]$state$extra$history)
  }

vis=Visualization$new(results_list=all_control[-63],variables=variables)

vis$plot_me2(x_var = 'tst',
             y_var = 'storage',
             col_var = 'num_links',
             #color_func = ,
             size_var = 'cycles_num',
        #    x_func = function(x) return(log(x)),
      #       y_func = function(y) return(log(y)) ,
             x_lim = NULL, #list(min=0,max=7000), 
             y_lim = NULL#list(min=0.0,max=1000000)
)


#-----------------
#Fig. decay analysis of resilience
test=NetBase$new(core_nodes=12L);
input=1000;
freq=40;
serieL=test$get_series_of_node_stock_after_stabilization2(initial_input = input,frequency = freq,do_plot = TRUE)

buscaValorSerie = function(s){
  esta = TRUE;
  comienza = 1;
  for(n in 2:length(s)){
    dif = (s[[n]]-s[[n-1]])
    if(!esta && dif<=0){
      esta = TRUE
      comienza = n
    }else if(dif>0){
      esta = FALSE
    }
  }
  if(esta)return(comienza)
  return(length(s))
}

limite = vector(mode='list',length=test$core_nodes);
for(n in 1:test$core_nodes){
  limite[[n]]=buscaValorSerie(unlist(map(serieL,n)))
}
thres = max(unlist(limite))

coeffs=vector(mode='list',length=length(serieL))
coeffsA=vector(mode='list',length=length(serieL))

yMax=max(t(matrix(unlist(serieL),nrow=test$core_nodes+2))[thres:freq,1:test$core_nodes])

A=(serieL[[1]]%*%test$get_mat_nth_power(n=thres-1))[1:test$core_nodes]
b=test$get_exp_factor()

vectTime = function(k,freq,A,b){
  output=vector(mode='list',length=freq)
  for(x in 1:freq){
    output[[x]]=A[[k]]*exp(b*(x-1))
  }
  return(output)
}

first=TRUE
for(k in 1:test$core_nodes){
  data0=unlist(map(serieL,k))
  data=data0[thres:length(data0)]
  s=seq(thres,thres+(length(data)-1),1)
  df=data.frame(s,data)
  names(df)=c('x','y')
  exp.eq <- function(x, a, b) {
    a*exp(b*x)
  }
  if(data[[1]]>0){
    m=nls(y ~ exp.eq(x, a, b), data = df, start = list(a = data[[1]], b = -0.1),trace = T)  
    if(first){
      first=FALSE
      plot(df,ylim=c(0,yMax),cex=0.5,ylab="energy accumulation", xlab="decay time")
    }
    else{
      points(df,cex=0.5)
    }
    lines(s, predict(m, list(x = s)), col = "blue")
    lines(s, vectTime(k,freq = length(data),A,b), col = "red")
    coeffs[[k]]=coefficients(m)[[2]]
    coeffsA[[k]]=coefficients(m)[[1]]
  }
}

media=mean(unlist(coeffs),na.rm = TRUE)
desv=sd(unlist(coeffs),na.rm = TRUE)
ev=test$get_exp_factor();

legend('topright',legend=c("nls fitted decay", "eigen simulated decay"), 
       fill = c("blue","red"),bty = "n")

text(36,85,paste0("mean= ", round(media, 4),'  sd=',round(desv,4)),cex=0.8,col='blue')
text(36,75,paste0('log(max eigen)= ',round(ev,4)),cex=0.8,col='red')
