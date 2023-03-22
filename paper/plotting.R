source("paper/stat_util.R")

#-----Fig 5. Finn per node vs. cycle flow----
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
------

#----Fig 3. decay analysis of resilience----
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



-----
  
#----Fig 11. resilience/performance relation ----

source("FlowNet.R")
all_scalars30_100_all=readRDS('all_scalars_30')
i_f_bases30_100_all=readRDS('initial_final_bases30')
all_scalars30_100=list()
i_f_bases30_100=list()
all_scalars30_100$TST=all_scalars30_100_all$TST
i_f_bases30_100$TST=i_f_bases30_100_all$TST
rm(all_scalars30_100_all,i_f_bases30_100_all)
gc()

all_scalars30_300_1=readRDS('all_scalars_30_primeros300')
i_f_bases30_300_1=readRDS('initial_final_bases30_primeros300')

all_scalars30_300_2=readRDS('all_scalars_30_segundos300')
i_f_bases30_300_2=readRDS('initial_final_bases30_segundos300')

all_scalars30_300_3=readRDS('all_scalars_30_terceros300')
i_f_bases30_300_3=readRDS('initial_final_bases30_terceros300')

getB_GOAL_mats = function(container,goal,goal_min){
  bExp=vector(mode='list',length=length(container[[goal]]$control))
  bM=vector(mode='list',length=length(container[[goal]]$control))
  goalM=vector(mode='list',length=length(container[[goal]]$control))
  for(i in 1:length(container[[goal]]$control)){
    print(i)
    k=1
    row_b = vector(mode='numeric',length=length(names(container[[goal]])))
    row_bExp = vector(mode='numeric',length=length(names(container[[goal]])))
    row_goal = vector(mode='numeric',length=length(names(container[[goal]])))
    row_names = vector(mode='character',length=length(names(container[[goal]])))
    for(a in names(container[[goal]])){
      #print(a)
      if(grepl('finn',a) || grepl('proportion',a)){
        row_bExp[k] = 0
      }else{
        row_bExp[k] = minimalFit(as.numeric(unlist(container[[goal]][[a]][[i]]$b)),do_plot = FALSE,print_eval = FALSE,num_fixs = 4000)
      }
      row_b[k] = container[[goal]][[a]][[i]]$b[length(container[[goal]][[a]][[i]]$b)]
      row_goal[k] = container[[goal]][[a]][[i]]$tst[length(container[[goal]][[a]][[i]][[goal_min]])]
      row_names[k] = a
      k=k+1
    }
    bExp[[i]]=row_bExp
    bM[[i]]=row_b
    goalM[[i]]=row_goal
  }
  M_b = as.matrix(Reduce(rbind,bM))
  M_bExp = as.matrix(Reduce(rbind,bExp))
  M_goal = as.matrix(Reduce(rbind,goalM))
  colnames(M_b)=row_names
  colnames(M_goal)=row_names
  colnames(M_bExp)=row_names
  #M_b
  #M_goal
  diff_b=matrix(0,nrow=nrow(M_b),ncol=ncol(M_b))
  diff_goal=matrix(0,nrow=nrow(M_goal),ncol=ncol(M_goal))
  for(i in 1:ncol(M_b)){
    diff_b[,i]=(M_b[,1]-M_b[,i])/abs((M_b[,1]+M_b[,i]))
  }
  for(i in 1:ncol(M_goal)){
    diff_goal[,i]=(M_goal[,1]-M_goal[,i])/abs((M_goal[,1]+M_goal[,i]))
  }
  colnames(diff_b)=row_names
  colnames(diff_goal)=row_names
  return(list(b=M_b,g=M_goal,d_b=diff_b,d_g=diff_goal,bExp=M_bExp))
}

get_FinnAndStock=function(container,goal,goal_min,nodes){
  nflow = NetBase$new(core_nodes = nodes)
  initialFinnPN=vector(mode='list',length=length(container[[goal]]$control$initialA))
  finalFinnPN=vector(mode='list',length=length(container[[goal]]$control$initialA))
  initialFinnTOT=vector(mode='list',length=length(container[[goal]]$control$initialA))
  finalFinnTOT=vector(mode='list',length=length(container[[goal]]$control$initialA))
  initialFinnPNSD=vector(mode='list',length=length(container[[goal]]$control$initialA))
  initialFinnTOTSD=vector(mode='list',length=length(container[[goal]]$control$initialA))
  
  initialStockPN=vector(mode='list',length=length(container[[goal]]$control$initialA))
  finalStockPN=vector(mode='list',length=length(container[[goal]]$control$initialA))
  initialStockTOT=vector(mode='list',length=length(container[[goal]]$control$initialA))
  finalStockTOT=vector(mode='list',length=length(container[[goal]]$control$initialA))
  initialStockPNSD=vector(mode='list',length=length(container[[goal]]$control$initialA))
  initialStockTOTSD=vector(mode='list',length=length(container[[goal]]$control$initialA))
  
  for(i in 1:length(container[[goal]]$control$initialA)){
    k=1
    row_initialPN = vector(mode='numeric',length=length(names(container[[goal]])))
    row_finalPN = vector(mode='numeric',length=length(names(container[[goal]])))
    row_initialTOT = vector(mode='numeric',length=length(names(container[[goal]])))
    row_finalTOT = vector(mode='numeric',length=length(names(container[[goal]])))
    
    row_initialStockPN = vector(mode='numeric',length=length(names(container[[goal]])))
    row_finalStockPN = vector(mode='numeric',length=length(names(container[[goal]])))
    row_initialStockTOT = vector(mode='numeric',length=length(names(container[[goal]])))
    row_finalStockTOT = vector(mode='numeric',length=length(names(container[[goal]])))
    
    row_initialStockPNSD = vector(mode='numeric',length=length(names(container[[goal]])))
    row_initialStockTOTSD = vector(mode='numeric',length=length(names(container[[goal]])))
    row_initialFinnTOTSD = vector(mode='numeric',length=length(names(container[[goal]])))
    row_initialFinnPNSD = vector(mode='numeric',length=length(names(container[[goal]])))
    
    row_names = vector(mode='character',length=length(names(container[[goal]])))
    for(a in names(container[[goal]])){
      iFinnPN = nflow$get_finn_per_node(mat=container[[goal]][[a]]$initialA[[i]],input_multiplier = 1000.0)
      fFinnPN = nflow$get_finn_per_node(mat=container[[goal]][[a]]$finalA[[i]],input_multiplier = 1000.0)
      
      iStockPN = nflow$get_initial_stock_after_stabilization(mat=container[[goal]][[a]]$initialA[[i]],initial_input = 1000.0,frequency = 1L)
      fStockPN = nflow$get_initial_stock_after_stabilization(mat=container[[goal]][[a]]$finalA[[i]],initial_input = 1000.0,frequency = 1L)
      
      iTot = sum(iFinnPN)
      fTot = sum(fFinnPN)
      iPN = iTot
      fPN = fTot
      
      iTOTSD = sd(iFinnPN)
      iPNSD = iTOTSD
      
      iStockTot = sum(iStockPN)
      fStockTot = sum(fStockPN)
      iSPN = iStockTot
      fSPN = fStockTot
      
      iSTOTSD = sd(iStockPN)
      iSPNSD = iSTOTSD
      
      if(k>1){
        picked = container[[goal]][[a]]$picked[[i]]
        iPN = sum(iFinnPN[picked])
        fPN = sum(fFinnPN[picked])
        
        iSPN = sum(iStockPN[picked])
        fSPN = sum(fStockPN[picked])
        
        iPNSD = sd(iFinnPN[picked])
        
        iSPNSD = sd(iStockPN[picked])
      }
      row_initialPN[k] = iPN
      row_finalPN[k] = fPN
      row_initialTOT[k] = iTot
      row_finalTOT[k] = fTot
      
      row_initialStockPN[k] = iSPN
      row_finalStockPN[k] = fSPN
      row_initialStockTOT[k] = iStockTot
      row_finalStockTOT[k] = fStockTot
      
      row_initialStockPNSD[k] = iSPNSD
      row_initialStockTOTSD[k] = iSTOTSD
      row_initialFinnTOTSD[k] = iTOTSD
      row_initialFinnPNSD[k] = iPNSD
      
      row_names[k] = a
      
      k=k+1
    }
    initialFinnPN[[i]] = row_initialPN
    finalFinnPN[[i]] = row_finalPN
    initialFinnTOT[[i]] = row_initialTOT
    finalFinnTOT[[i]] = row_finalTOT
    
    initialStockPN[[i]] = row_initialStockPN
    finalStockPN[[i]] = row_finalStockPN
    initialStockTOT[[i]] = row_initialStockTOT
    finalStockTOT[[i]] = row_finalStockTOT
    
    initialStockPNSD[[i]] = row_initialStockPNSD
    initialFinnPNSD[[i]] = row_initialFinnPNSD
    initialStockTOTSD[[i]] = row_initialStockTOTSD
    initialFinnTOTSD[[i]] = row_initialFinnTOTSD
  }
  M_i_PN = as.matrix(Reduce(rbind,initialFinnPN))
  M_f_PN = as.matrix(Reduce(rbind,finalFinnPN))
  M_i_TOT = as.matrix(Reduce(rbind,initialFinnTOT))
  M_f_TOT = as.matrix(Reduce(rbind,finalFinnTOT))
  
  M_i_s_PN = as.matrix(Reduce(rbind,initialStockPN))
  M_f_s_PN = as.matrix(Reduce(rbind,finalStockPN))
  M_i_s_TOT = as.matrix(Reduce(rbind,initialStockTOT))
  M_f_s_TOT = as.matrix(Reduce(rbind,finalStockTOT))
  
  M_i_PNSD = as.matrix(Reduce(rbind,initialFinnPNSD))
  M_i_s_PNSD = as.matrix(Reduce(rbind,initialStockPNSD))
  M_i_TOTSD = as.matrix(Reduce(rbind,initialFinnTOTSD))
  M_i_s_TOTSD = as.matrix(Reduce(rbind,initialStockTOTSD))
  
  colnames(M_i_PN)=row_names
  colnames(M_f_PN)=row_names
  colnames(M_i_TOT)=row_names
  colnames(M_f_TOT)=row_names
  
  colnames(M_i_s_PN)=row_names
  colnames(M_f_s_PN)=row_names
  colnames(M_i_s_TOT)=row_names
  colnames(M_f_s_TOT)=row_names
  
  colnames(M_i_PNSD)=row_names
  colnames(M_i_s_PNSD)=row_names
  colnames(M_i_TOTSD)=row_names
  colnames(M_i_s_TOTSD)=row_names
  
  return(list(i_pn=M_i_PN,f_pn=M_f_PN,i_tot=M_i_TOT,f_tot=M_f_TOT,
              i_s_pn=M_i_s_PN,f_s_pn=M_f_s_PN,i_s_tot=M_i_s_TOT,f_s_tot=M_f_s_TOT,
              i_pn_sd=M_i_PNSD,i_s_pn_sd=M_i_s_PNSD,i_tot_sd=M_i_TOTSD,i_s_tot_sd=M_i_s_TOTSD))
}

get_InOuts=function(container,goal,goal_min,nodes){
  nflow = NetBase$new(core_nodes = nodes)
  ins=vector(mode='list',length=length(container[[goal]]$control$initialA))
  outs=vector(mode='list',length=length(container[[goal]]$control$initialA))
  m_val_outs=vector(mode='list',length=length(container[[goal]]$control$initialA))
  non_outs=vector(mode='list',length=length(container[[goal]]$control$initialA))
  weigh_outs=vector(mode='list',length=length(container[[goal]]$control$initialA))
  for(i in 1:length(container[[goal]]$control$initialA)){
    k=1
    #JUST SUM OF 1 or 0 FOR EACH NODE IF IT HAS OUT
    row_ins = vector(mode='numeric',length=length(names(container[[goal]])))
    row_outs = vector(mode='numeric',length=length(names(container[[goal]])))
    
    #SUM OF OUT VALUES FOR ALL NODES
    val_outs = vector(mode='numeric',length=length(names(container[[goal]])))
    
    #SUM OF NON OUTS
    num_non_outs_of_outs = vector(mode='numeric',length=length(names(container[[goal]])))
    
    #SUM OF OUT VALUE MULTIPLIED BY NODE STOCK
    wheighted_outs = vector(mode='numeric',length=length(names(container[[goal]])))
    
    row_names = vector(mode='character',length=length(names(container[[goal]])))
    
    for(a in names(container[[goal]])){
      cmat =container[[goal]][[a]]$initialA[[i]]
      row_ins[k] = length(which(cmat[nodes+1,]>0))
      row_outs[k] = length(which(cmat[,nodes+2]>0))
      iStockPN = nflow$get_initial_stock_after_stabilization(mat=cmat,initial_input = 1000.0,frequency = 1L)
      val_outs[k] = sum(cmat[which(cmat[,nodes+2]>0),nodes+2])
      num_non_outs_of_outs[k] = length(which(cmat[which(cmat[,nodes+2]>0),1:nodes]>0))
      has_out = which(cmat[,nodes+2]>0)
      wheighted_outs[k] = sum(as.vector(cmat[has_out,nodes+2])*iStockPN[has_out])
      if(k>1){
        picked = container[[goal]][[a]]$picked[[i]]
        has_out = picked[which(cmat[picked,nodes+2]>0)]
        row_ins[k] = length(which(container[[goal]][[a]]$initialA[[i]][nodes+1,picked]>0))
        row_outs[k] = length(which(container[[goal]][[a]]$initialA[[i]][picked,nodes+2]>0))
        val_outs[k] = sum(cmat[picked,nodes+2][which(cmat[picked,nodes+2]>0)])
        num_non_outs_of_outs[k] = length(which(cmat[which(cmat[picked,nodes+2]>0),1:nodes]>0))
        wheighted_outs[k] = sum(as.vector(cmat[has_out,nodes+2])*iStockPN[has_out])
      }
      
      row_names[k] = a
      
      k=k+1
    }
    ins[[i]] = row_ins
    outs[[i]] = row_outs
    m_val_outs[[i]]= val_outs
    non_outs[[i]]=num_non_outs_of_outs
    weigh_outs[[i]]=wheighted_outs
  }
  M_ins = as.matrix(Reduce(rbind,ins))
  M_outs = as.matrix(Reduce(rbind,outs))
  M_val_outs = as.matrix(Reduce(rbind,m_val_outs))
  M_non_outs = as.matrix(Reduce(rbind,non_outs))
  M_weigh_outs = as.matrix(Reduce(rbind,weigh_outs))
  
  
  colnames(M_ins)=row_names
  colnames(M_outs)=row_names
  colnames(M_val_outs)=row_names
  colnames(M_non_outs)=row_names
  colnames(M_weigh_outs)=row_names
  
  
  return(list(ins=M_ins,outs=M_outs,val_outs=M_val_outs,non_outs=M_non_outs,w_outs=M_weigh_outs))
}

test100 = getB_GOAL_mats(all_scalars30_100,'TST','tst')
testFS100 = get_FinnAndStock(i_f_bases30_100,'TST','tst',30)
testIO100 = get_InOuts(i_f_bases30_100,'TST','tst',30)

test1 = getB_GOAL_mats(all_scalars30_300_1,'TST','tst')
testFS1 = get_FinnAndStock(i_f_bases30_300_1,'TST','tst',30)
testIO1 = get_InOuts(i_f_bases30_300_1,'TST','tst',30)

test2 = getB_GOAL_mats(all_scalars30_300_1,'TST','tst')
testFS2 = get_FinnAndStock(i_f_bases30_300_1,'TST','tst',30)
testIO2 = get_InOuts(i_f_bases30_300_1,'TST','tst',30)

test3 = getB_GOAL_mats(all_scalars30_300_1,'TST','tst')
testFS3 = get_FinnAndStock(i_f_bases30_300_1,'TST','tst',30)
testIO3 = get_InOuts(i_f_bases30_300_1,'TST','tst',30)

test1000 = test100
testFS1000 = testFS100
testIO1000 = testIO100
for(i in 1:5){
  test1000[[i]] = rbind(test1000[[i]][,1:7],rbind(rbind(test1[[i]],test2[[i]]),test3[[i]]))
}
for(i in 1:12){
  testFS1000[[i]] = rbind(testFS1000[[i]][,1:7],rbind(rbind(testFS1[[i]],testFS2[[i]]),testFS3[[i]]))
}
for(i in 1:5){
  testIO1000[[i]] = rbind(testIO1000[[i]][,1:7],rbind(rbind(testIO1[[i]],testIO2[[i]]),testIO3[[i]]))
}

length(which(test1000$d_b[,4]<0))
length(which(test1000$d_g[,4]<0))

test = test1000
testFS = testFS1000
testIO = testIO1000
#to_test=(testFS$i_pn[,4]/testFS$i_tot[,4])
#rbPal <- colorRampPalette(c('red','yellow','blue'))
#cols <- rbPal(100)[as.numeric(cut(to_test,breaks = 100))]

exper = 4
plot(testIO$non_outs[,exper],test$d_b[,exper])
#abline(v=0,h=0)
binary_vals = unlist(lapply(test$d_b[,exper],FUN=function(x){if(x>0)return(1)else return(-1)}))
resglm = glm(test$d_b[,exper] ~ testIO$non_outs[,exper] + testIO$val_outs[,exper] + testFS$i_pn[,exper])
summary(resglm)

#FINN PERCENTAGE
dat=data.frame(Finn_greed_perc=testFS$i_pn[,exper]/testFS$i_tot[,exper],b_diff=test$d_b[,exper])

tst_diff=(test$d_g[,exper])#testIO$outs[,exper]#

#principal plot
p1=qplot(Finn_greed_perc, b_diff, data=dat, colour=tst_diff) + scale_colour_gradient2(low=hsv(0,0.8,1.0,0.8),mid = hsv(0.333,0.8,1.0,0.5), high=hsv(0.666,0.8,1.0,0.8))+theme_classic()+
  geom_vline(xintercept = 0,linetype="dotted")+geom_hline(yintercept = 0,linetype="dotted")+
  theme(legend.position='left')+
  #labs(x=str2lang("(sum(Finn(Ngreed))/Finn(Ntot))"),y=str2lang('(b_control - b_greedy)/abs((b_control+b_greedy))'),colour='TST Diff')+
  labs(x='Finn percentage (Ngreed)',y='Resilience Difference (b)',colour='TST Diff')+
  scale_x_continuous(expand = c(0,0),lim=c(0,1.01))+
  scale_y_continuous(expand = c(0,0),lim=c(-1,1),breaks = scales::pretty_breaks(n = 10))

gred=which(test$d_b[,exper]<0)
condition=rep('control',1000)#length(all_scalars30_300_1$TST$control))
condition[gred]='greedy'
dat$condition=condition

#histogram x
p2=ggplot(dat,aes(x=Finn_greed_perc))+
  geom_histogram(data=subset(dat,condition == 'control'),aes(y = after_stat(count / sum(count))),bins=21,boundary=0,fill = "blue", alpha = 0.2,color='white') +
  geom_histogram(data=subset(dat,condition == 'greedy'),aes(y = after_stat(count / sum(count))),bins=21,boundary=0,fill = "red", alpha = 0.25,color='white') +
  theme_classic()+
  labs(y='Frequency',x='')+theme(legend.position='right',axis.text=element_text(size=5),axis.title=element_text(size=8))+
  scale_y_continuous(expand = c(0,0),lim=c(0,0.3))+
  scale_x_continuous(expand = c(0,0),lim=c(0,1),breaks = scales::pretty_breaks(n = 11))+
  annotate(geom='text',x=0.9,y=0.165,label='better control',color='blue',size=3)+
  annotate(geom='text',x=0.9,y=0.125,label='better greedy',color='red',size=3)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

#histogram y
dat$positive=as.factor(dat$b_diff>0)

p3=ggplot(dat,aes(y=b_diff,fill=positive))+
  geom_histogram(aes(x = after_stat(count / sum(count)),fill=condition),bins=21,boundary=-0.000001, alpha = 0.3,color="white") +
  #geom_histogram(data=subset(dat,condition == 'greedy'),aes(x = after_stat(count / sum(count))),bins=21,boundary=0,fill = "blue", alpha = 0.2,color='white') +
  theme_classic()+
  scale_fill_manual(values=c('red','blue'))+
  labs(x='Frequency',y='')+theme(legend.position='top',axis.title=element_text(size=8),axis.text=element_text(size=5))+
  scale_x_continuous(expand = c(0,0),lim=c(0,0.2))+
  scale_y_continuous(expand = c(0,0),lim=c(-1,1),breaks = scales::pretty_breaks(n = 11))+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = 'none')

#multiplot:
library(aplot)

p1 %>%
  insert_top(p2,height=0.2) %>%
  insert_right(p3,width=0.2)

#OUTPUT PERCENTAGE

#exper=2
bins=15
freq_max = 0.2 #exp 2 0.4, exp 3 0.3, exp 4 0.2

dat=data.frame(Out_greed_perc=(testIO$outs[,exper]),b_diff=test$d_b[,exper])

tst_diff=(test$d_g[,exper])

#principal plot
p1=qplot(Out_greed_perc, b_diff, data=dat, colour=tst_diff) + scale_colour_gradient2(low=hsv(0,0.8,1.0,0.8),mid = hsv(0.333,0.8,1.0,0.5), high=hsv(0.666,0.8,1.0,0.8))+theme_classic()+
  geom_hline(yintercept = 0,linetype="dotted")+
  theme(legend.position='left')+
  #labs(x=str2lang("(sum(Finn(Ngreed))/Finn(Ntot))"),y=str2lang('(b_control - b_greedy)/abs((b_control+b_greedy))'),colour='TST Diff')+
  labs(x='Output Percentage (Ngreed)',y='Resilience Difference (b)',colour='TST Diff')+
  scale_x_continuous(expand = c(0,0),lim=c(-0.5,bins+0.5))+
  scale_y_continuous(expand = c(0,0),lim=c(-1,1),breaks = scales::pretty_breaks(n = bins))

gred=which(test$d_b[,exper]<0)
condition=rep('control',1000)#length(all_scalars30_300_1$TST$control))
condition[gred]='greedy'
dat$condition=condition

#histogram x
p2=ggplot(dat,aes(x=Out_greed_perc))+
  geom_histogram(center = 0,binwidth=1,data=subset(dat,condition == 'control'),aes(y = after_stat(count / sum(count))),bins=bins,fill = "blue", alpha = 0.2,color='white') +
  geom_histogram(center = 0,binwidth=1,data=subset(dat,condition == 'greedy'),aes(y = after_stat(count / sum(count))),bins=bins,fill = "red", alpha = 0.25,color='white') +
  theme_classic()+
  labs(y='Frequency',x='')+theme(legend.position='right',axis.text=element_text(size=5),axis.title=element_text(size=8))+
  scale_y_continuous(expand = c(0,0),lim=c(0,freq_max))+
  scale_x_continuous(expand = c(0,0),lim=c(-0.5,bins+0.5),breaks = scales::pretty_breaks(n = bins))+
  annotate(geom='text',x=bins-2.2,y=0.17,label='better control',color='blue',size=3)+
  annotate(geom='text',x=bins-2.2,y=0.13,label='better greedy',color='red',size=3)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

#histogram y
dat$positive=as.factor(dat$b_diff>0)

p3=ggplot(dat,aes(y=b_diff,fill=positive))+
  geom_histogram(aes(x = after_stat(count / sum(count)),fill=condition),bins=21,boundary=-0.000001, alpha = 0.3,color="white") +
  #geom_histogram(data=subset(dat,condition == 'greedy'),aes(x = after_stat(count / sum(count))),bins=21,boundary=0,fill = "blue", alpha = 0.2,color='white') +
  theme_classic()+
  scale_fill_manual(values=c('red','blue'))+
  labs(x='Frequency',y='')+theme(legend.position='top',axis.title=element_text(size=8),axis.text=element_text(size=5))+
  scale_x_continuous(expand = c(0,0),lim=c(0,0.2))+
  scale_y_continuous(expand = c(0,0),lim=c(-1,1),breaks = scales::pretty_breaks(n = 11))+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = 'none')

#multiplot:
library(aplot)

p1 %>%
  insert_top(p2,height=0.2) %>%
  insert_right(p3,width=0.2)

#percStock_i=testFS$i_s_pn[,4]/testFS$i_s_tot[,4]
#percFinn_i=testFS$i_pn[,4]/testFS$i_tot[,4]
#plot((testIO$outs[,4])/(testIO$outs[,1]),test$d_b[,4])
#plot(percFinn_i,test$d_b[,4])
  

#-----Fig SI. orientor scaling with network size----
#using control evolution 
source('FlowNet.R')
input_output_entropy=input_output_entropy_generator(log_base=2);
entropy_diff=entropy_diff_generator(log_base=2);

node_ensemble=c(30,50,70,100)
orientors=list(TST=tst, ASC=ascendency, AMI=ami,Finn=finn,EDiff=entropy_diff,b=b)
plotting=c('tst','ascendency','ami','finn','entropy_diff','b');

or=6

all_controls=vector(mode='list', length=length(node_ensemble))
for (N in 1:length(node_ensemble)){
  setwd(main_folder)
  node_folder=paste0(main_folder,'/',node_ensemble[N])
  orientor_folder=paste0(node_folder,'/',names(orientors[or]))
  orientor_folder
  setwd(orientor_folder)
  
  ens_evolve_control=readRDS('ens_evolve_control')
  results_control=ens_evolve_control$get_results()
  results_control_list=lapply(results_control,function(x) x$state$extra$history)
  
  all_controls[[N]]=results_control_list
}

controls=Reduce(rbind,all_controls)
n_nets=length(controls)
variables=all_controls[[1]][[1]][[1]]$var

vis = Visualization$new(results_list=controls,variables=variables)

plot=vis$plot_me2(x_var = 'time',
             y_var = 'ediff_norm',
             col_var='num_nodes',
             #color_func=function(x){return(hsv(h = log((x-min_max_col$min)/(min_max_col$max-min_max_col$min)*0.8), s = 1, v = 1, 1))},
             x_func = function(x) return(x),
             y_func=function(y) return(y),
             with_raster = TRUE
             #y_lim = list(min=4,max=15)
)



-----
  
#-----Fig SI. stability evolution----

source("FlowNet.R")
main_folder=getwd()
input_output_entropy=input_output_entropy_generator(log_base=2);
entropy_diff=entropy_diff_generator(log_base=2);

node_ensemble=c(30,50,70,100)
orientors=list(TST=tst, ASC=ascendency, AMI=ami,Finn=finn,EDiff=entropy_diff,b=b)
plotting=c('tst','ascendency','ami','finn','entropy_diff','b');
greed_perc=c(0.1,0.3,0.5)
greed_orientor=c('part_stock','proportion_stock')#,'part_finn','proportion_finn')

N=3
or=1
greed_or=1
greed_p=1

con=Manipulation1Cp$new(ensemble_size=node_ensemble[N],orientor=names(orientors[or]))
control=con$open_control()

greedy_path=paste0('cp_greedy_',greed_perc[greed_p],'_',greed_orientor[greed_or])
print(greedy_path)

gred=Manipulation1Cp$new(ensemble_size=node_ensemble[N],orientor=names(orientors[or]),greed_orientor=greed_orientor[greed_or],greed_perc=greed_perc[greed_p])
greedy=gred$open_greedy()

#gred$plot_both(net = 7,y = 'b',iter=7000)

vis=Visualization$new(results_list=con$results_control_list,variables=con$variables)
vis_greedy=Visualization$new(results_list=gred$results_greedy_list,variables=gred$variables)

plot_both_history = function(cual,var,num_fixs = 5000){
  cont = as.numeric(unlist(vis$data[[cual]][var]))#[-1]#[-c(1:2)]
  gred = as.numeric(unlist(vis_greedy$data[[cual]][var]))#[-1]#[-c(1:2)]
  mVal = min(c(cont,gred))
  MVal = max(c(cont,gred))
  plot(cont,col='blue',ylim=c(mVal,0),xlim=c(0,num_fixs),xlab='evolutionary fixations (~time)',ylab='flux decay rate (b)')
  points(gred,col='red')
  control = nls.control(maxiter = 5000,tol = 1e-20,minFactor = 1/9192,printEval = TRUE,warnOnly = TRUE)
  
  x_c = 1:length(cont)-1
  #nlmod = nls(cont ~ (A)*exp(-x_c * B),start = list(A=cont[[1]],B=0.001),control = control)#,algorithm="port",upper=c(-1e-30,10000))
  nlmod = nls(cont ~ cont[[1]]*exp(-x_c * B),start = list(B=0.0001),control = control)#,algorithm="port",upper=c(-1e-30,10000))
  
  x_g = 1:length(gred)-1
  #nlmod_g = nls(gred ~ (A)*exp(-x_g * B),start = list(A=gred[[1]],B=0.001))#,control = control)#,algorithm="port",upper=c(-1e-30,10000))
  nlmod_g = nls(gred ~ gred[[1]]*exp(-x_g * B),start = list(B=0.0001))#,control = control)#,algorithm="port",upper=c(-1e-30,10000))
  
  lines(1:num_fixs,unlist(lapply(1:num_fixs,function(x){cont[[1]]*exp(-x * summary(nlmod)$coefficients[[1]])})),col='green',lwd=2)
  
  lines(1:num_fixs,unlist(lapply(1:num_fixs,function(x){gred[[1]]*exp(-x * summary(nlmod_g)$coefficients[[1]])})),col='purple',lwd=2)
  #print(nlmod)
  #print(nlmod_g)
  print(summary(nlmod)$coefficients[[1]])
  print(summary(nlmod_g)$coefficients[[1]])
  print(cont[[1]]*exp(-num_fixs * summary(nlmod)$coefficients[[1]]))
  print(gred[[1]]*exp(-num_fixs * summary(nlmod_g)$coefficients[[1]]))
}

plot_both_history(22,'b',num_fixs = 3000)


minimalFit = function(values,do_plot=TRUE,print_eval=TRUE,num_fixs=4000){
  control = nls.control(maxiter = 5000,tol = 1e-20,minFactor = 1/9192,printEval = print_eval,warnOnly = TRUE)
  x_v = 1:length(values)-1
  nlmod = nls(values ~ values[[1]]*exp(-x_v * B),start = list(B=0.001),control = control)#,algorithm="port",upper=c(-1e-30,10000))
  if(do_plot){
    mVal = min(values)
    MVal = max(values)
    plot(values,col='blue',ylim=c(mVal,0),xlim=c(0,num_fixs))
    lines(1:num_fixs,unlist(lapply(1:num_fixs,function(x){values[[1]]*exp(-x * summary(nlmod)$coefficients[[1]])})),col='green',lwd=4)
  }
  return(summary(nlmod)$coefficients[[1]])
}
minimalFit(as.numeric(unlist(vis$data[[32]]['b'])),do_plot = FALSE,print_eval = FALSE,num_fixs = 4000)

  
