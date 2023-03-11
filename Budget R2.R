#--------------------------------------------------------------------------------#
## Simulation für R^2-Budget ##
#--------------------------------------------------------------------------------#

###
p=14
N=1000
R=0.49
Qse=10

#Modell:b*x1+b*l*x2+b*u*x3+b*a*X4+b*e*X1:X3+b*d*X1:X4+b*c*X2:X4+b*w*X1:X2+b*y*X3:X4+b*z*X2:X3+b*t*X1:X1+b*g*X2:X2+b*s*X3:X3+b*q*X4:X4
#Annahmen: multivariat nv, Erwartungswerte alle gleich 0, var=1


#Kovarianz/Korrelationsmatrix
Sigma <- lavaan::lav_matrix_vech_reverse(c(0.3,0.2,0.1,0.4,0.3,0.1), diag = FALSE)
diag(Sigma) <- c(1)

#Budgetfunktion
Budget<-function(QSE,Rq,N,p,
                 KOR12,KOR13,KOR14,KOR23,KOR24,KOR34,
                 l,u,a,t,g,s,q,w,e,d,z,c,y){
  QSY=QSE/(1-Rq)
  VARY=QSY/(N-1)
  VARE=QSE/(N-p-1)
  beta<-sqrt(VARY-VARE)/sqrt(1+l^2+u^2+a^2+e^2*(1+KOR13^2)+d^2*(1+KOR14^2)+c^2*(1+KOR24^2)+w^2*(KOR12^2+1)+y^2*(1+KOR34^2)+z^2*(1+KOR23^2)+2*t^2+2*g^2+2*s^2+2*q^2+
                               2*(l*KOR12+u*KOR13+a*KOR14)+
                               2*l*(u*KOR23+a*KOR24)+
                               2*u*a*KOR34+
                               2*e*(d*(KOR34+KOR14*KOR13)+c*(KOR12*KOR34+KOR14*KOR23)+w*(KOR23+KOR13*KOR12)+y*(KOR14+KOR34*KOR13)+z*(KOR12+KOR23*KOR13)+2*t*KOR13+2*g*KOR12*KOR23+2*s*KOR13+2*q*KOR14*KOR34)+
                               2*d*(c*(KOR12+KOR14*KOR24)+w*(KOR24+KOR12*KOR14)+y*(KOR13+KOR14*KOR13)+z*(KOR12*KOR34+KOR13*KOR24)+2*t*KOR14+2*g*KOR12*KOR24+2*s*KOR13*KOR14+2*q*KOR14)+
                               2*c*(w*(KOR14+KOR12*KOR24)+y*(KOR23+KOR24*KOR34)+z*(KOR34+KOR23*KOR24)+2*t*KOR12*KOR14+2*g*KOR24+2*s*KOR23*KOR34+2*q*KOR24)+
                               2*w*(y*(KOR13*KOR24+KOR23*KOR14)+z*(KOR13+KOR12*KOR23)+2*t*KOR12+2*g*KOR12+2*s*KOR13*KOR23+2*q*KOR14*KOR24)+
                               2*y*(z*(KOR24+KOR23*KOR34)+2*t*KOR13*KOR14+2*g*KOR23*KOR24+2*s*KOR34+2*q*KOR34)+
                               2*z*(2*t*KOR12*KOR13+2*g*KOR23+2*s*KOR23+2*q*KOR34*KOR24)+
                               2*t*(2*g*KOR12^2+2*s*KOR13^2+2*q*KOR14^2)+
                               2*g*(2*s*KOR23^2+2*q*KOR24^2)+
                               2*s*2*q*KOR34^2)
return(beta)
}

Bmod<-Budget(QSE = 10, Rq = 0.49, N = 1000, p = 14,
             KOR12 = 0.3, KOR13 = 0.2, KOR14 = 0.1, KOR23 = 0.4, KOR24 = 0.3, KOR34 = 0.1,
             l = 2, u = 3, a = 4, t = 5, g = 6, s = 7, q = 8, w = 9, e = 10, d = 11, z = 12, c = 13, y = 14)
betas<-c(Bmod*1:14)

#Simulation
sample<-1:15000
SimulationRQuadrat<-sapply(sample,FUN=function(isample)  {            
  
  X <- mvtnorm::rmvnorm(N , mean =c(0,0,0,0),sigma = Sigma)
  X <- model.matrix(~ X1 + X2 + X3 + X4 + I(X1^2) + I(X2^2) + I(X3^2) + I(X4^2)+ X1*X2 + X1*X3 + X1*X4 + X2*X3 + X2*X4 + X3*X4 , data = data.frame(X))
  Xoi<-X[,-1]
  y<- Xoi %*% betas+ rnorm(N ,mean = 0, sd= sqrt(Qse/(N-p-1)))
  data<-data.frame(y,X)
  
  fitR<-lm(y~X,data)
  fitRsum<-summary(fitR)
  list(fitR$coefficients,fitRsum$sigma^2*(N-p-1),fitRsum$r.squared)
  
}, simplify="array")

#Auswertung
RQuadrat<-sapply(SimulationRQuadrat[3,],sum)
sum(RQuadrat)/15000
QSe<-sapply(SimulationRQuadrat[2,],sum)
sum(QSe)/15000
#passt 

