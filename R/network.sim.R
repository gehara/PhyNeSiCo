## simualte network
sim.sp.tree<-function(tree,
                      Ne.prior,
                      biffurcating=T,
                      migration=F,
                      admixture=F,
                      mig,
                      hib.clade,
                      hib.priors,
                      major.sister,
                      minor.sister,
                      bp,
                      mi,
                      nsims,
                      ntrees,
                      segsites,
                      gen.time,
                      time.modif,
                      coaltrees,
                      spectral,
                      distance,
                      time.scalar=1)
{
  
  #### exclude branch lengths and change taxon to numbers
  tree.b<-getrid(tree)
  tree.b<-taxons.to.numbers(tree.b)
  
  ### get all joints of the tree
  input<-tree.b[[2]]
  join<-NULL
  while(length(grep("(",input,fixed=T))>=1){
    xx<-get.node(input)
    join<-rbind(join,xx[[2]])
    input<-xx[[1]]
  }
  
  ### get all joint times of the tree
  input<-tree
  tree2<-tree
  join.t<-NULL
  while(length(grep("(",input,fixed=T))>=1){
    xx<-get.time(input,tree2)
    join.t<-rbind(join.t,xx[[3]])
    input<-xx[[1]]
  }
  
  ####### generate ej (nodes) flag string. Nodes ages are rescaled.
  ej<-cbind(join.t,join)
  ej<-cbind(rep("-ej",nrow(ej)),ej)
  ej[,2]<-as.numeric(ej[,2])*time.scalar
  
  ###### generate en (ancestral Ne) flag. Sample random Nes. Ne values will change in the simulation cicle
  en<-cbind(rep("-en",nrow(ej)),ej[,2],ej[,3],runif(nrow(ej),Ne.prior[1],Ne.prior[2]))
  for(i in 1:nrow(en)){
    en[i,3]<-strsplit(ej[i,3]," ")[[1]][2]
  }
  
  #######
  ms.string<-list()
  simulated<-NULL
  ## master Ne
  master.Ne<-mean(Ne.prior)
  ## get nodes
  nodes<-as.numeric(ej[,2])
  
  #######################
  #######################
  # start simulations ###
  for(j in 1:nsims){
    sim.t<-NULL
    trees<-list()
    
    # put the original dates back on time strings
    ej[,2]<-nodes
    en[,2]<-nodes
    # rescale to coalescent (Ne proportion)
    ej[,2]<-as.numeric(ej[,2])/(4*master.Ne)/gen.time
    en[,2]<-as.numeric(ej[,2])#*1.01
    # sample time modifier
    time.mod<-runif(1,time.modif[1],time.modif[2])
    # modify time
    ej[,2]<-as.numeric(ej[,2])*time.mod
    en[,2]<-as.numeric(en[,2])*time.mod
    
    # sample an Ne mean
    Ne.mean<-runif(1,Ne.prior[1],Ne.prior[2])
    # sample Ne Standart deviation
    Ne.SD<-Ne.mean*runif(1,0.1,1)
    # sample contemporary Nes (truncated in 100 individuals)
    Nes<-rtnorm((nrow(ej)+1),Ne.mean,Ne.SD,lower=100)
    # Sample ancestral Nes (truncatedd 100 individuals)
    anc.Nes<-rtnorm(nrow(ej),Ne.mean,Ne.SD,lower=100)
    if(migration==T){
      Mig.prop<-runif(1,mig[1],mig[2])
      ntrees.mig<-floor(ntrees*Mig.prop)
      }else{
        Mig.prop<-0
      }
    if(admixture==T){
      Admix.prob.minor<-runif(1,hib.priors[4],hib.priors[5])
    }else{
      Admix.prob.minor<-0
    }
    ### simulate n-trees
    count<-0
    for(i in 1:ntrees){
      master.theta<-0
      while(master.theta<0.000001){
      # sample mutation rate
      mi.mean<-runif(1,mi[1],mi[2])
      mi.SD<-mi.mean*runif(1,0.1,1)
      #mi.SD<-mi.mean*mi.SD
      rate<-rnorm(1,mi.mean,mi.SD)
      while(rate<=0){
        rate<-rnorm(1,mi.mean,mi.SD)
      }
      # sample sequence length
      seq.length<-rnorm(1,bp[1],bp[2])
      # generate theta
      master.theta<-master.Ne*4*seq.length*(rate*gen.time)
      }
      # theta string
      ms.string[[1]]<-paste("-t",master.theta)
      ### pop structure
      ms.string[[2]]<-paste(c("-I",(nrow(ej)+1),(rep(1,nrow(ej)+1))),collapse=" ")
      ### current popsize strings
      n<-cbind(rep("-n",(nrow(ej)+1)),c(1:(nrow(ej)+1)),rep(0,(nrow(ej)+1)),Nes)
      n[,3]<-as.numeric(n[,4])/master.Ne
      ms.string[[3]]<-paste(apply(n[,1:3],1,paste,collapse=" "),collapse=" ")
      ### ancestral pop sizes string
      en[,4]<-anc.Nes
      en[,4]<-as.numeric(en[,4])/master.Ne
      ms.string[[4]]<-paste(apply(en,1,paste,collapse=" "),collapse=" ")
      ###  node string
      ms.string[[5]]<-paste(apply(ej,1,paste,collapse=" "),collapse=" ")
      
      
      #### network string
      
      if(biffurcating==F){
        if(admixture==T){
        split.time<-runif(1,hib.priors[1],hib.priors[2])*time.scalar
        ej.hib<-runif(1,split.time,(hib.priors[2]*time.scalar))
        split.time<-((split.time/gen.time)*time.mod)/(4*master.Ne)
        ej.hib<-((ej.hib/gen.time)*time.mod)/(4*master.Ne)
        es<-paste("-es",split.time,max(hib.clade),(1-Admix.prob.minor))
        ejh<-paste("-ej",ej.hib,(length(tree.b[[1]])+1),max(minor.sister))
        ms.string[[6]]<-paste(es,ejh)
      }
      if(migration==T){
        count<-count+1
        if(count>ntrees.mig){
          Mig<-0
        }else{
          Mig=1
        }
        #Mig<-rtnorm(1,Mig.mean,Mig.SD,0)
        mig.time<-((hib.priors[1]/gen.time)*time.mod)/(4*master.Ne)
        em<-paste("-em",mig.time,max(minor.sister),max(hib.clade),Mig)
        ms.string[[6]]<-em
      }
      }
      
      ######
      # combining all string peaces in one ms string
      ms.string.final<-paste(unlist(ms.string),collapse=" ")
      #print(ms.string.final)
      ########################################
      ########################################
      ########################################
      #### get the coalescent trees
      if(coaltrees==T){
        t<-ms(nreps = 1, nsam=(nrow(ej)+1),opts=paste("-T",ms.string.final))
        t<-strsplit(t,"//")[[3]]
        t<-read.tree(text=t)
        trees[[i]]<-t
      }
      
      #### build trees from simulated segregating sites
      if(segsites==T){
        fas<-ms.to.DNAbin(ms(nreps = 1, nsam=(nrow(ej)+1),opts=ms.string.final),bp.length = 0)
        
        while(length(fas)==0){
          fas<-ms.to.DNAbin(ms(nreps = 1, nsam=(nrow(ej)+1),opts=ms.string.final),bp.length = 0)  
        }
        
        d<-dist.dna(fas, model="N")/seq.length
        
        if(distance==T){
          sim.t<-rbind(sim.t,as.vector(d))
        }
       
        if(spectral==T){
          upgma.tree<-upgma(d)
          upgma.tree$edge.length<-upgma.tree$edge.length*1000
          trees[[i]]<-upgma.tree
        }
        
      }
      rm(fas)
    }
    
    if(spectral==T){
      eigen<-get.panda(trees,method="standard")
      sim.t<-cbind(eigen$principal_eigenvalue,eigen$asymmetry,eigen$eigenvalues)
      if(nrow(sim.t)>1){
        sim<-c(apply(sim.t,2,mean),apply(sim.t,2,var),apply(sim.t,2,kurtosis),apply(sim.t,2,skewness))
      } else {sim<-t(sim.t)} 
    }
    if(distance==T){
      if(segsites==F){
        stop("you need to have segsites=TRUE to get the distance")} else {
          if(nrow(sim.t)>1){
            nam<-t(combn(attr(d,"Labels"),2))
            nam<-apply(nam,1,paste,collapse="_")
            colnames(sim.t)<-nam
            sim<-c(apply(sim.t,2,mean))#,apply(sim.t,2,var),apply(sim.t,2,kurtosis),apply(sim.t,2,skewness))
            mean<-paste("mean",names(sim[1:length(nam)]),sep="_")
            #var<-paste("var",names(sim[(length(nam)+1):(2*length(nam))]),sep="_")
            #kur<-paste("kur",names(sim[((2*length(nam))+1):(3*length(nam))]),sep="_")
            #skew<-paste("skew",names(sim[((3*length(nam))+1):(4*length(nam))]),sep="_")
            names(sim)<-c(mean)#,var,kur,skew)
          } else {
            nam<-t(combn(attr(d,"Labels"),2))
            nam<-apply(nam,1,paste,collapse="_")
            colnames(sim.t)<-nam
            sim<-t(sim.t)} 
        }
    }
    #print(i)
    sim<-cbind(time.mod,Ne.mean,Ne.SD,mi.mean,mi.SD,Mig.prop,Admix.prob.minor,t(sim))
    simulated<-rbind(simulated,sim)
    rm(time.mod,Ne.mean,Ne.SD,mi.mean,mi.SD,Mig.prop,Admix.prob.minor,sim)
    print(j)
  }
  return(simulated)
}

