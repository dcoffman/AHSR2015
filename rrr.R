setwd("/Users/dlc30/Documents/Conferences/AHSR conference")
library(survey)
library(twang)
########################################################################################################
# Mediation
########################################################################################################

dat <- read.csv("rrr-cjdats.csv", header=T)

#Obtain weights for a continuous mediator
num.mod <- lm(RRA_30 ~ COND, data=dat)
den.mod <- lm(RRA_30 ~ COND + CBLACK + ABUSE6M + suprob + drgprb + alcprb + ownhouse + nevermarried + married + fulltime +
                parttime + LIVSP + JAIL30D + SATXIMP + ARST30D + MEDINS + HS + DGTXJLC + ALCDYS + FINSUP1 + condom_intake,
              data=dat)
sigma.n <- summary(num.mod)$sigma
sigma.d <- summary(den.mod)$sigma
num.p <- dnorm((dat$RRA_30-num.mod$fitted)/sigma.n, 0,1)
den.p <- dnorm((dat$RRA_30-den.mod$fitted)/sigma.d, 0,1)
dat$wt.cont <- num.p/den.p

#Outcome analysis for continuous mediator
design.ps <- svydesign(ids= ~1, weights= ~wt.cont, data=dat)
mod.y.cont <- svyglm(condomyn ~ RRA_30 + COND, design=design.ps, family=quasibinomial(link="logit"))
summary(mod.y.cont)

mod.y.mod <- svyglm(condomyn ~ RRA_30 + COND + RRA_30:COND, design=design.ps, family=quasibinomial(link="logit"))
summary(mod.y.mod)

mod.m.cont <- lm(RRA_30 ~ COND, data=dat)
summary(mod.m.cont)

#Check balance
design.b <- svydesign(ids= ~1, weights=~wt.cont, data=dat)
covmat1 <- svyvar(RRA_30 ~ CBLACK + ABUSE6M + suprob + drgprb + alcprb + ownhouse + nevermarried + married + fulltime +
                 parttime + LIVSP + JAIL30D + SATXIMP + ARST30D + MEDINS + HS + DGTXJLC + ALCDYS + FINSUP1 + condom_intake, design.b)
var1 <- diag(covmat1)
corr1 <- covmat1[1,-1]/sqrt(var1[1]*var1[-1])
AAC1 <- mean(abs(corr1))

##################################################################################
#Using boosting
##################################################################################

library(energy)
library(gbm)
library(polycor)

F.aac.iter=function(i,data,ps.model,ps.num,rep,criterion) { 
  # i: number of iterations (trees)
  # data: dataset containing the treatment and the covariates
  # ps.model: the boosting model to estimate p(T_i|X_i)
  # ps.num: the estimated p(T_i)
  # rep: number of replications in bootstrap
  # criterion: the correlation metric used as the stopping criterion
  GBM.fitted=predict(ps.model,newdata=data,n.trees=floor(i),
                     type="response") 
  ps.den=dnorm((data$T-GBM.fitted)/sd(data$T-GBM.fitted),0,1)
  wt=ps.num/ps.den
  aac_iter=rep(NA,rep)
  for (i in 1:rep){
    bo=sample(1:dim(data)[1],replace=TRUE,prob=wt)
    newsample=data[bo,]
    j.drop=match(c("T"),names(data))
    j.drop=j.drop[!is.na(j.drop)]
    x=newsample[,-j.drop]
    if(criterion=="spearman"|criterion=="kendall"){
      ac=apply(x, MARGIN=2, FUN=cor, y=newsample$T,
               method=criterion)
    } else if (criterion=="distance"){
      ac=apply(x, MARGIN=2, FUN=dcor, y=newsample$T) 
    } else if (criterion=="pearson"){
      ac=matrix(NA,dim(x)[2],1)
      for (j in 1:dim(x)[2]){
        ac[j]=ifelse (!is.factor(x[,j]), cor(newsample$T, x[,j], method=criterion),
                      polyserial(newsample$T, x[,j]))
      }
    } else print("The criterion is not correctly specified")
    aac_iter[i]=mean(abs(1/2*log((1+ac)/(1-ac))),na.rm=TRUE)
  }
  aac=mean(aac_iter)
  return(aac)
}

dat$T <- dat$RRA_30
num.mod <- lm(RRA_30 ~ COND, data=dat)
num.p <- dnorm((dat$RRA_30-num.mod$fitted)/(summary(num.mod))$sigma, 0,1)
den.mod <- gbm(RRA_30 ~ COND + CBLACK + ABUSE6M + suprob + drgprb + alcprb + ownhouse + nevermarried + married + fulltime +
                 parttime + LIVSP + JAIL30D + SATXIMP + ARST30D + MEDINS + HS + DGTXJLC + ALCDYS + FINSUP1 + condom_intake, 
               data=dat, shrinkage = 0.0005,
               interaction.depth=4, distribution="gaussian", n.trees=5000)
opt <- optimize(F.aac.iter, interval=c(1,5000), data=dat, ps.model= den.mod, 
                ps.num=num.p, rep=50, criterion="pearson")
best.aac.iter <- opt$minimum
best.aac <- opt$objective

# Calculate the inverse probability weights
den.mod$fitted <- predict(den.mod, newdata=dat,
                         n.trees=floor(best.aac.iter), type="response")
den.p <- dnorm((dat$RRA_30 - den.mod$fitted)/sd(dat$RRA_30 - den.mod$fitted), 0, 1)
dat$weight.gbm <- num.p/den.p

# Outcome analysis using survey package
design.b <- svydesign(ids= ~1, weights= ~weight.gbm, data = dat)
fit <- svyglm(condomyn ~ RRA_30 + COND, family=quasibinomial(), design=design.b)
summary(fit)

