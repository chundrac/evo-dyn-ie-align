require(rstan)
require(phytools)
require(phangorn)
require(expm)
require(tikzDevice)
require(RColorBrewer)

#rate graphs

load('all_rates.Rdata')

for (f in 1:9) {
rates <- data.frame(extract(total.list[[f]])$R)
tree.ind <- rep(1:50,each=nrow(rates)/50)

F <- length(levels(diacl.data[,f])) #careful!!

charnames <- c('Present-past',
               'Noun, present progressive',
               'Noun, simple past',
               'Pronoun, present progressive',
               'Pronoun, simple past',
               'Reflexive with agent',
               'Reflexive with object',
               'Verb, present progressive',
               'Verb, simple past')
              

likelihood <- to.matrix(diacl.data[,f],seq=levels(diacl.data[,f]))
rownames(likelihood) <- rownames(diacl.data)
likelihood[rowSums(likelihood)==0,] <- likelihood[rowSums(likelihood)==0,]+1

fill.matrix <- function(R) {
    Q <- matrix(nrow=F,ncol=F)
    k <- 1
    for (i in 1:F) {
        for (j in 1:F) {
            if (i != j) {
                Q[i,j] <- R[[k]]
                k <- k + 1
            }
        }
    }
    for (i in 1:F) {
        z <- 0
        for (j in 1:F) {
            if (i != j) {
                z <- z - Q[i,j]
            }
        }
        Q[i,i] <- z
    }
    return(Q)
}



anc.states <- function(t,f) {
    rate.slice <- as.vector(rates[t,])
    tree <- IE.trees[[tree.ind[t]]]
    tree <- reorder.phylo(tree,'pruningwise')
    parent <- tree$edge[,1]
    child <- tree$edge[,2]
    Q <- fill.matrix(rate.slice)
    brlens <- tree$edge.length/1000
    B <- length(brlens)
    P <- list() #each branch's transition probability matrix
    #get likelihoods	
    curr.states <- likelihood[tree$tip.label,]
    lambda <- matrix(data=0,nrow=N,ncol=F)
    for (n in 1:T) {
        for (f in 1:F) {
            lambda[n,f] <- log(curr.states[n,f])
        }
    }
    for (b in 1:B) {
        P[[b]] <- expm(Q*brlens[b])
        for (f in 1:F) {
            lambda[parent[b],f] <- lambda[parent[b],f] + log(P[[b]][f,]%*%exp(lambda[child[b],]))
        }
    }
    pi = expm(10000*Q)[1,]
    lambda[parent[B],] <- lambda[parent[B],] + log(pi)
    #normalize likelihoods
    phi <- prop.table(exp(lambda),margin=1)
    #draw state and root and move down the tree
    states.sim <- matrix(data=0,nrow=N,ncol=F)
    states.sim[parent[B],] <- t(rmultinom(1,1,phi[parent[B],]))
    for (b in B:1) {
        states.sim[child[b],] <- t(rmultinom(1,1,(states.sim[parent[b],]%*%P[[b]]) * phi[child[b],]))
    }
    return(states.sim)
    #    return(states.sim[(T+1):N,])
}


set.seed(1)
states.list <- list()
for (i in 1:1000) {
    j = sample(1:nrow(rates),1)
    states.list[[i]] <- anc.states(j,f)
}




#states.agg <- data.frame(prop.table(Reduce('+',states.list),margin=1))
#states.agg <- data.frame(Reduce('+',states.list))
states.agg <- prop.table(Reduce('+',states.list),margin=1)
#rownames(states.agg) <- tree$node.label
colnames(states.agg) <- levels(diacl.data[,f])
#states.agg$node <- c(1:N)

saveRDS(states.agg,file=paste(f,' ',charnames[f],'.rds',sep=''))

node.colors <- brewer.pal(n = 8, name = "Set2")
#tikz('state_tree')
cairo_pdf(paste(f,' ',charnames[f],'.pdf',sep=''),height=10,width=15)
mcc.tree <- maxCladeCred(IE.trees)


boxlabel<-function(x,y,text,cex=1,bg="transparent",offset=0){
    w<-strwidth(text)*cex*1.1
    h<-strheight(text)*cex*1.4
    os<-offset*strwidth("W")*cex
    rect(x+os,y-0.5*h,x+w+os,y+0.5*h,col=bg,border=0)
    text(x,y,text,cex=cex,pos=4,offset=offset,font=1)
}

#par(fg="transparent")
plotTree(mcc.tree,title='charnames[f]',direction="downwards")

#add.scale.bar()
#tiplabels(pie=states.agg[1:T,],cex=.5,piecol=node.colors)
#nodelabels(pie=states.agg[(T+1):N,],cex=.5,piecol=node.colors)
non.advent<-c((T+1):N)[sapply(c((T+1):N),function(x) {return(!any(mcc.tree$edge.length[mcc.tree$edge[,1]==x]==1))})]
nodelabels(pie=states.agg[non.advent,],node=non.advent,cex=.5,piecol=node.colors)
tiplik<-as.matrix(likelihood[mcc.tree$tip.label,]/rowSums(likelihood[mcc.tree$tip.label,]))
#tiplabels(pie=tiplik,node=mcc.tree$tip.label,cex=.5,piecol=node.colors)

#pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
#N<-Ntip(mcc.tree)
#par(fg="black")
#for(i in 1:Ntip(mcc.tree)) {
#    if (sum(likelihood[mcc.tree$tip.label[i],])==1) {
#    #if !(is.na(diacl.data[mcc.tree$tip.label[i],f])) {
#  boxlabel(pp$xx[i],pp$yy[i],tree$tip.label[i],bg=node.colors[as.numeric(diacl.data[mcc.tree$tip.label[i],f])])
#  }
#  else {
#      boxlabel(pp$xx[i],pp$yy[i],tree$tip.label[i],bg='transparent')
#  }
#}


#CAREFUL!!
legend('topleft',legend=levels(diacl.data[,f]),fill=node.colors,bty='n',cex=1)
#legend('topleft',legend=levels(states$V2),fill=node.colors,bty='n',cex=1)
dev.off()
}
