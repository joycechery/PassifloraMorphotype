#function to computer posterior probabilities at the nodes (L88-93) by liam revell from http://blog.phytools.org/2012/12/plotting-node-piecharts-on-top-of.html
#function to plot simmaps and node posterior probabilities (117-229--"Plot_simmap") by Dr. Michael May (UC Berkeley)

#this code plots stochastic character maps along the branches of the phylogeny and piecharts at the nodes
library(phytools)
library(grDevices)
library(geiger)

setwd("~/Desktop/Passiflora/")
tree <-read.nexus("Passiflora_maximum_clade_credibility_ultrametric.tree")
plotTree(tree)

#tree$tip.label
tree <- reroot(tree, 1, interactive = TRUE)

# upload data 
data  <- read.csv("Passiflora_morphotypes.csv", row.names = 1) 
datum <- as.data.frame(cbind(rownames(data), data[,1]), row.names = FALSE) #select the column with the trait of interest
datum <- datum[complete.cases(datum), ]
datum

# species
species <- datum$V1

# reformat data
datas <- datum$V2
names(datas) <- species

# drop tips

KeepMe <-c("Passiflora_actinia",
           "Passiflora_alata",
           "Passiflora_amethystina",
           "Passiflora_biflora",
           "Passiflora_caerulea",
           "Passiflora_capsularis",
           "Passiflora_cincinnata",
           "Passiflora_coccinea",
           "Passiflora_coriacea",
           "Passiflora_edmundoi",
           "Passiflora_edulis",
           "Passiflora_foetida",
           "Passiflora_galbana",
           "Passiflora_gibertii",
           "Passiflora_hatschbachii",
           "Passiflora_kermesina",
           "Passiflora_ligularis",
           "Passiflora_malacophylla",
           "Passiflora_maliformis",
           "Passiflora_micropetala",
           "Passiflora_miersii",
           "Passiflora_misera",
           "Passiflora_mollissima",
           "Passiflora_morifolia",
           "Passiflora_mucronata",
           "Passiflora_nitida",
           "Passiflora_organensis",
           "Passiflora_pohlii",
           "Passiflora_racemosa",
           "Passiflora_rubra",
           "Passiflora_setacea",
           "Passiflora_sidifolia",
           "Passiflora_tenuifila",
           "Passiflora_tricuspis",
           "Passiflora_villosa")

pruned_tree <- keep.tip(tree, KeepMe)
pruned_tree=ladderize(pruned_tree)
plot(pruned_tree)

####character transition model
fitER<-fitDiscrete(pruned_tree,datas,model= "ER")
fitSYM<-fitDiscrete(pruned_tree,datas,model= "SYM")
fitARD<-fitDiscrete(pruned_tree,datas,model= "ARD")

#Compare the AIC scores from the model outputs below, picking the model with the lowest AIC
fitER$opt$aic
fitSYM$opt$aic
fitARD$opt$aic

# stochastic mapping, use the model with the lowest AIC score
simmap.trees <- make.simmap(pruned_tree, datas, model = "ER", pi = "estimated", nsim=1000)

# set colors & character states
col_vec <- c("steelblue3", "darkorange","forestgreen", "red2", "mediumpurple3", "sienna", "plum1")
states<- c("Morphotype_A", "Morphotype_B", "Morphotype_C",  "Morphotype_D", "Morphotype_E", "Morphotype_F", "Morphotype_G")

# function to compute the node states
foo<-function(x){
  y<-sapply(x$maps,function(x) names(x)[1])
  names(y)<-x$edge[,1]
  y<-y[as.character(length(x$tip)+1:x$Nnode)]
  return(y)
}

XX<-sapply(simmap.trees,foo)
pies<-t(apply(XX,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=states,Nsim=1000))

#generate summary of stochastic maps with pies of posterior at nodes..
#see below (line 119) for plot_simmap function

plot_simmap(time_tree = simmap.trees[[1]], 
            tree = simmap.trees[[1]], 
            simmaps = simmap.trees, 
            states = states,
            show.tip.label = T,
            lwd = 8,
            label.cex = 1.2,
            label.offset=.8,
            colors = col_vec, edge.width=0.1, nt=10001)

add.simmap.legend(colors=col_vec, prompt=FALSE,x=105,y=180, fsize =.9)
nodelabels(pie=pies,cex=0.8,piecol=col_vec, lwd=1)
legend(x=0,y=1, legend = states, col = col_vec, pch = 20, yjust = 0, bty = 'n', cex =.7, pt.cex = 4)


#####this is the function being use to plot density plot + pie chart  
###written by Dr. Mike May 2021

plot_simmap = function(time_tree, tree, simmaps, states, colors, nt=1001, show.tip.label=FALSE, edge.width=0, plot_pie=TRUE, pie_size=1.5, label.offset=0, label.cex=1, lwd=1, ...) {
  
  # compute dt
  dt = sum(time_tree$edge.length) / nt
  
  # for each branch, compute the probability at each time slice
  num_branches = nrow(tree$edge)
  num_sims     = length(simmaps)
  num_states   = length(states)
  
  if ( missing(colors) ) {
    colors = 1:num_states
  }
  col_rgb = col2rgb(colors)
  
  branch_x0 = vector("list", num_branches)
  branch_x1 = vector("list", num_branches)
  branch_sample_colors = vector("list", num_branches)
  for(i in 1:num_branches) {
    # compute the scale factor
    factor = time_tree$edge.length[i] / tree$edge.length[i]
    
    # compute the time points
    sampled_times = seq(0, time_tree$edge.length[i], by=dt)
    #if (length(sampled_times) == 1) {
    #    sampled_times = c(0, time_tree$edge.length[i])
    #}
    
    # compute the states per time slice
    sampled_states = vector("list", num_sims)
    for(j in 1:num_sims) {
      this_map = simmaps[[j]]$maps[[i]] * factor
      this_map_cum = cumsum(this_map)
      this_map_states = names(this_map_cum)
      sampled_states[[j]] = this_map_states[findInterval(sampled_times, this_map_cum, left.open=TRUE)+1]
    }
    sampled_states = do.call(rbind, sampled_states)
    # compute the probabilities per time slice
    state_probs = matrix(NA, nrow=num_states, ncol=length(sampled_times))
    for(j in 1:num_states) {
      state_probs[j,] = colMeans(sampled_states == states[j])
    }
    
    # compute the colors per time slice
    segment_colors = numeric(length(sampled_times))
    for(j in 1:length(sampled_times)) {
      these_cols = col_rgb %*% state_probs[,j]
      segment_colors[j] = rgb(these_cols[1,1], these_cols[2,1], these_cols[3,1], maxColorValue = 255)
    }
    
    # store computed values
    branch_x0[[i]] = sampled_times
    branch_x1[[i]] = c(sampled_times[-1], time_tree$edge.length[i])
    branch_sample_colors[[i]] = segment_colors
    
  }
  
  # get the coordinates for the tree
  yy = node.height(time_tree)
  xx = node.depth.edgelength(time_tree)
  
  # plot the tree
  pp = plot.phylo(time_tree, type = "phylogram",
                  use.edge.length = TRUE,
                  show.tip.label  = show.tip.label,
                  edge.width      = edge.width,
                  direction       = "rightwards",
                  cex             = label.cex,
                  label.offset    = label.offset)
  
  # plot the segments
  for(i in 1:num_branches) {
    
    # get the horiztonal segments
    this_x0  = branch_x0[[i]] + xx[ time_tree$edge[i,1] ]
    this_x1  = branch_x1[[i]] + xx[ time_tree$edge[i,1] ]
    this_y0  = yy[ time_tree$edge[i,2] ]
    this_col = branch_sample_colors[[i]]
    segments(x0=this_x0, x1=this_x1, y0=this_y0, col=this_col, lwd=lwd, ...)
    
    # get the vertical segment
    this_x0  = xx[ time_tree$edge[i,1] ]
    this_y0  = yy[ time_tree$edge[i,1] ]
    this_y1  = yy[ time_tree$edge[i,2] ]
    
    this_col = branch_sample_colors[[i]][1]
    segments(x0=this_x0, y0=this_y0, y1=this_y1, col=this_col, lwd=lwd, ...)
    
  }
  
  # plot data pie charts
  if ( plot_pie == FALSE ) {
    return()
  }
  
  for(i in 1:nrow(data)) {
    
    this_taxon = rownames(data)[i]
    this_data  = data[i,]
    this_index = which(tree$tip.label == this_taxon)
    
    this_x = xx[ this_index ]
    this_y = yy[ this_index ]
    
    floating.pie(xpos = this_x, ypos = this_y, x=this_data+1e-10, col=colors, radius=pie_size, lwd=ifelse(edge.width > lwd, edge.width - lwd, lwd)  )
    
  }
  
  
}


