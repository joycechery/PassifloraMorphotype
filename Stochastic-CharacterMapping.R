#function to computer posterior probabilities at the nodes (L89-95) by liam revell from http://blog.phytools.org/2012/12/plotting-node-piecharts-on-top-of.html
#function to plot simmaps and node posterior probabilities in the code file named "Plot_simmap" by Dr. Michael May

#this code plots stochastic character maps along the branches of the phylogeny and piecharts at the nodes
library(phytools)
library(grDevices)
library(geiger)

setwd("~/Desktop/Passiflora/")
tree <-read.nexus("Passiflora_maximum_clade_credibility.tree")
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

#generate summary of stochastic maps with pies of posterior at nodes.
plot_simmap(time_tree = simmap.trees[[1]], 
            tree = simmap.trees[[1]], 
            simmaps = simmap.trees, 
            states = states,
            show.tip.label = T,
            lwd = 9,
            label.cex = 1.2,
            label.offset=.007,
            colors = col_vec, edge.width=0.1, nt=10001)

add.simmap.legend(colors=col_vec, prompt=FALSE,x=105,y=180, fsize =.9)
nodelabels(pie=pies,cex=0.7,piecol=col_vec, lwd=1)
legend(x=0,y=1, legend = states, col = col_vec, pch = 20, yjust = 0, bty = 'n', cex =.7, pt.cex = 4)

