source("FlowNet.R")

table=list()

greed_experiment=c('part_stock','proportion_stock')
picked_percent=c(0.1,0.3,0.5)

for (greed in 1:length(greed_experiment)){
  for (pickp in 1:length(picked_percent)){
open_name=opening(30,100,1,'TST',greed_experiment[greed],picked_percent[pickp]) 
experiment <- readRDS(file=open_name)
#attach(experiment)
results=experiment$results;
results_control=experiment$results_control;
results_greedy=experiment$results_greedy;
ens_evolve=experiment$ens_evolve;

results_list=lapply(results,function(x){x$state$extra$history})
variables=results_list[[1]][[1]]$var;
n_nets=length(ens_evolve$evolved_nets); #length(results_control[[1]])

#for greedy ensembles:
cp_control=vector(mode='list',length=(length(results_control)+1))
cp_greedy=vector(mode='list',length=(length(results_control)+1))
all_control=vector(mode='list',length=n_nets)
all_greedy=vector(mode='list',length=n_nets)
for (cp in 1:(length(results_control)+1)){
  for (net in 1:n_nets){
    if(cp==1){
      end=results[[net]]$state$extra$fixations
      state_control=results[[net]]$state$extra$history[[end]]
      state_greedy=results[[net]]$state$extra$history[[end]]
    }
    else{
      end=results_control[[cp-1]][[net]]$state$extra$fixations
      state_control=results_control[[cp-1]][[net]]$state$extra$history[[end]]
      end=results_greedy[[cp-1]][[net]]$state$extra$fixations
      state_greedy=results_greedy[[cp-1]][[net]]$state$extra$history[[end]]
      all_control[[net]]=c(all_control[[net]],results_control[[cp-1]][[net]]$state$extra$history)
      all_greedy[[net]]=c(all_greedy[[net]],results_greedy[[cp-1]][[net]]$state$extra$history)
    }
    if(length(state_control)==0){
      end=results_control[[cp-2]][[net]]$state$extra$fixations
      state_control=results_control[[cp-2]][[net]]$state$extra$history[[end]]
    }
    if(length(state_greedy)==0){
      end=results_greedy[[cp-2]][[net]]$state$extra$fixations
      state_greedy=results_greedy[[cp-2]][[net]]$state$extra$history[[end]]
    }
    cp_control[[cp]][[net]]=state_control;
    cp_greedy[[cp]][[net]]=state_greedy;
  }
}

cp_list_control=vector(mode='list',length=n_nets)
cp_list_greedy=vector(mode='list',length=n_nets)
cp_list_control_scalars=vector(mode='list',length=n_nets)
cp_list_greedy_scalars=vector(mode='list',length=n_nets)
for (net in 1:n_nets){
  cp_list_control[[net]]=t(sapply(cp_control,'[[',net))
  cp_list_greedy[[net]]=t(sapply(cp_greedy,'[[',net))
  nonscalars=c()
  for (i in 1:length(colnames(cp_list_control[[net]]))){
    if(length(cp_list_control[[net]][[1,colnames(cp_list_control[[net]])[i]]])>1){
      nonscalars=c(nonscalars,i)
    }
  }
  cp_control_reduced = cp_list_control[[net]][,-nonscalars]
  cp_greedy_reduced = cp_list_greedy[[net]][,-nonscalars]
  cp_control_mat = matrix(data=unlist(cp_control_reduced),nrow=nrow(cp_control_reduced),ncol=ncol(cp_control_reduced),dimnames = dimnames(cp_control_reduced))
  cp_greedy_mat = matrix(data=unlist(cp_greedy_reduced),nrow=nrow(cp_greedy_reduced),ncol=ncol(cp_greedy_reduced),dimnames = dimnames(cp_greedy_reduced))
  cp_list_control_scalars[[net]]=as.data.frame(cp_control_mat)
  cp_list_greedy_scalars[[net]]=as.data.frame(cp_greedy_mat)
}

find_slope=function(dat){
  variables=names(dat[[1]]);
  slopes=matrix(nrow=length(dat),ncol=length(variables))
  for (var in 1:length(variables)){
    obj=vector(mode='list',length=length(data));
    for (i in 1:length(dat)){
      cp=length(dat[[i]][,1]);
      dumX=1:cp;
      colum=variables[var]
      v=dat[[i]][[colum]]
      #slopes[i,var]=lm(as.matrix(eval)~dumX)$coefficients[[2]]
      #fixed intercept (starting point):
      slopes[i,var]=lm(I(v-v[1])~0+dumX)$coefficients[[1]]
      slopes=as.data.frame(slopes)
      names(slopes)=variables
    }
  }
  return(as.data.frame(slopes))
}

slopes_control=find_slope(dat=cp_list_control_scalars)
slopes_greedy=find_slope(dat=cp_list_greedy_scalars)


#vis = Visualization$new(results_list=results_list,variables=variables)
#slopes_results=find_slope(dat=vis$data)

x=(slopes_greedy$b-slopes_control$b)/slopes_control$b*100
y=(slopes_greedy$tst-slopes_control$tst)/slopes_control$tst*100

Gperc=length(which(x>0))/length(x)*100
better=which(x>0)
if(!length(better)==0){
XG=summary(x[better])
xrangeG=c(format(round(XG[['Min.']],2)),format(round(XG[['Median']],2)),format(round(XG[['Max.']],2)))
XC=summary(-1*x[-better])
xrangeC=c(format(round(XC[['Min.']],2)),format(round(XC[['Median']],2)),format(round(XC[['Max.']],2)))
}
else{
  XG=0;
  xrangeG=c(0,0,0)
  XC=summary(-1*x)
  xrangeC=c(format(round(XC[['Min.']],2)),format(round(XC[['Median']],2)),format(round(XC[['Max.']],2)))
}

table[[open_name]]=rbind(Gperc,xrangeG,xrangeC)
  }
}


df=data.frame(x,y)

proportional_finn=function(base){
  finn_node=base$get_finn_per_node()/base$get_finn()*100
  pick=base$net_data$picked_nodes
  return(sum(finn_node[pick]))
}

new_ensemble_greedy=experiment$new_ensemble_greedy

df$greed=unlist(lapply(new_ensemble_greedy$generated_networks,proportional_finn))
df$order = findInterval(df$greed, sort(df$greed))

library(RColorBrewer)
cols = brewer.pal(4, "Blues")
pal = colorRampPalette(cols)

par(mfrow=c(1,1))
plot(x ~ y, df, pch=19, col=pal(nrow(df))[df$order])
legend("topleft", col=pal(2), pch=19,
       legend=c(round(range(df$greed), 1)))
abline(v=0,h=0)

