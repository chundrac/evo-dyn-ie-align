require(rstan)
require(phytools)
require(phangorn)
require(expm)
require(tikzDevice)
require(RColorBrewer)

#rate graphs

load('all_rates.Rdata')

charnames <- c('Present-past',
               'Noun, present progressive',
               'Noun, simple past',
               'Pronoun, present progressive',
               'Pronoun, simple past',
               'Reflexive with agent',
               'Reflexive with object',
               'Verb, present progressive',
               'Verb, simple past')

charlist <- c()
deltas <- c()

for (f in 1:9) {
  
  F <- length(levels(diacl.data[,f])) #careful!!
  
  node.probs <- readRDS(file=paste(f,' ',charnames[f],'.rds',sep=''))
  
model <- "data {
  int<lower=1> N;
  int<lower=1> F;
  matrix[N,F] pi;
}
transformed data {
  real e[N];
  real F_;
  F_ = F;
  for (i in 1:N) {
    e[i] = 0;
    for (f in 1:F) {
      if (pi[i,f] <= 1./F_) {
        e[i] += pi[i,f];
      }
      else {
        e[i] += ((1./(1.-(F_)))*pi[i,f]) - (1./(1.-(F_)));
      }
    }
    if (e[i] == 0) {
      e[i] += 1e-10;
    }
  }
  //print(e);
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> lambda_0;
}
model {
  lambda_0 ~ gamma(1,1);
  alpha ~ exponential(lambda_0);
  beta ~ exponential(lambda_0);
  for (i in 1:N) {
    e[i] ~ beta(alpha,beta);
  }
}"

data.list <- list(pi=node.probs,N=nrow(node.probs),F=ncol(node.probs))

fit <- stan(model_code=model,data=data.list)

charlist <- c(charlist,rep(charnames[f],4000))
deltas <- c(deltas,extract(fit)$beta/extract(fit)$alpha)


#deltas[[f]] <- extract(fit)$beta/extract(fit)$alpha

}

deltas <- data.frame(character=charlist,delta=deltas)
#deltas$character <- factor(deltas$character,levels=rev(levels(deltas$character)))
new_order=with(deltas, reorder(character,delta,median))

require(ggplot2)
pdf('character_deltas.pdf')
ggplot(data=deltas,aes(x=new_order,y=delta)) + geom_violin() + coord_flip() + labs(x = "character")
dev.off()