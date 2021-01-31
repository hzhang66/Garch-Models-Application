library(readr)
library(DataAnalytics) 
library(pracma)
library(stats)
library(fGarch)
library(tidyr)
library(rugarch)

#  step 1  (a)
data <- read_csv("Downloads/etfglobal/week2&3/SP500_data.csv",col_names = FALSE)
colnames(data)<-c("date", "price")
data$return<-(data$price-back(data$price))/back(data$price)
data$return2<-data$return^2
data<-drop_na(data)
acf(data$return)
pacf(data$return)
#determine mean model
armaSearch<-function(data){
  armacoef<-data.frame()
  for (p in 0:7){
    for (q in c(0:7)) {
      data.arma = arima(data, order = c(p, 0, q))
      armacoef<-rbind(armacoef,c(p,q,data.arma$aic))
    
    }
  }
  colnames(armacoef)<-c("p","q","AIC")
  pos<-which(armacoef$AIC==min(armacoef$AIC))
  cat("the min AIC=",armacoef$AIC[pos], ",p=",armacoef$p[pos],",q=",armacoef$q[pos])
  return(armacoef,armacoef$p[pos],armacoef$q[pos])
}
armaSearch(data$return)


garch1=myspec=ugarchspec(
  
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1), submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
  
  mean.model = list(armaOrder = c(5, 4), include.mean = TRUE, archm = FALSE, archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE),
  
  distribution.model = "norm"
  
)
fit1=ugarchfit(garch1,data=data$return,solver="solnp")
mu=mean(data$return)
sigma=fit1@fit[["sigma"]]
residual=(data$return-mu)/sigma
corr(as.matrix(cbind(data$return2[1:length(data$return2)],back(data$return2)[1:length(data$return2)])))
residual2=residual^2
corr(as.matrix(cbind(residual2[1:length(residual2)],back(residual2)[1:length(residual2)])))
Box.test(data$return, type="Box-Pierce")

#  step 1  (b)

garch2=myspec=ugarchspec(
  
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1), submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
  
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE, archm = FALSE, archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE),
  
  distribution.model = "ged"
  
)

fit2=ugarchfit(garch2,data=data$return,solver="solnp")

#  step 1  (c)

garch3 = ugarchspec(
      
  variance.model = list(model="fGARCH",submodel="TGARCH", garchOrder=c(1,1),external.regressors = NULL, variance.targeting = FALSE),　  
                    
  mean.model = list(armaOrder=c(0,0), include.mean = TRUE, archm = FALSE, archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE),
                    
  distribution.model = "norm"
)  
fit3=ugarchfit(garch3,data=data$return,solver="solnp")

#step2
sumlike2=sum(fit2@fit[["log.likelihoods"]])
sumlike3=sum(fit3@fit[["log.likelihoods"]])

infocriteria(fit1)
infocriteria(fit3)
#model1 is the best

#step3

rf=0.0051/365
div=0.023/365
dt=1
omega=fit3@fit[["coef"]][["omega"]]
alpha=fit3@fit[["coef"]][["alpha1"]]
gamma=fit3@fit[["coef"]][["eta11"]]
beta=fit3@fit[["coef"]][["beta1"]]
S0=2126.41+3.58
forc=ugarchforecast(fit3,n.ahead=15)
sigmaT=c(fit3@fit[["sigma"]][length(fit3@fit[["sigma"]])],forc@forecast[["sigmaFor"]])
GBMTGARCH<-function (day){
  a_new=data$return[length(data$return)]-fit3@fit[["coef"]][["mu"]]
  rand=rnorm(day, mean=0, sd=1)
  sigmat=c(fit3@fit[["sigma"]][length(fit3@fit[["sigma"]])])
  St=c(S0)
  for (t in 1:day){
    s_new=St[t]*exp((rf-div)*dt-sigmaT[t]^2/2+sigmaT[t]*rand[t])
    St=c(St,s_new)
    return_new=(s_new-St[t])/St[t]
    a_new=return_new-mu
    if (a_new<0){s=1} else{s=0}
    sigma_new=sqrt(omega+alpha*a_new^2+gamma*s*a_new^2+beta*sigmat[t]^2)
    sigmat=c(sigmat,sigma_new)
  }
  return(St[day+1])
}

optionprice<- function(K,day,type,N=50000){
  
   ST=c()
   
   for (k in 1:N){ST=c(ST,GBMTGARCH(day))}
   
   if (type=="call") {payoff= pmax(ST-K,0)}
   if (type=="put") {payoff= pmax(K-ST,0)}
   
   price=mean(payoff)*exp(-rf*day*dt)
   
   return(price)
}
p1<-optionprice(2070,15,"put")
p2<-optionprice(2100,15,"put")
p3<-optionprice(2130,15,"call")
p4<-optionprice(2160,15,"call")
p5<-optionprice(2190,15,"call")
result<-data.frame(option=c("2070 puts", "2100 puts", "2130 calls", "2160 calls", "2190 calls"))
result$marketprice=c(15.2,21.9,28.3,12.9,4)
result$tgarchprice=c(p1,p2,p3,p4,p5)
for (i in 1:5){
  if (result$tgarchprice[i]<result$marketprice[i]) {result$assesment[i]="overvalue"} 
  else {result$assesment[i]="undervalue"}
}


