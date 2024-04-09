install.packages("rgl")
install.packages("psych",dependencies=TRUE)
install.packages(c("quantmod", "ggplot2"))
install.packages("copula")
install.packages("VineCopula")
install.packages("openxlsx")
install.packages("VC2copula")

library(MASS)
set.seed(100)
m <- 3
n <- 2000
corr_matrix_norm <- matrix(c(1, 0.5, 0.1,
                  0.5, 1, -0.8,
                  0.1, -0.8, 1), 
                nrow=3)

z <- mvrnorm(n,mu=rep(0, m),Sigma=corr_matrix_norm,empirical=T)
cor(z,method='pearson')
head(z)


#Is mvrnorm doing a Choolesky decomposition? Is done by eigen -> deepen this
help(mvrnorm, package = "MASS")

library(psych)
cor(z,method='pearson')
pairs.panels(z,hist.col = 'steelblue')


#transformation in Uniform with CDF - ITT
u <- pnorm(z)
pairs.panels(u, hist.col = 'steelblue', )

#3D plot of the dependence structure
library(rgl)
plot3d(u[,1],u[,2],u[,3],pch=30,col='navyblue')

#Choosing marginals. I assumed marginal are already Known.
x1 <- qgamma(u[,1],shape=2,scale=1)
x2 <- qbeta(u[,2],2,2)
x3 <- qt(u[,3],df=5)

#scatter after apply the inverse CDF
df <- cbind(x1,x2,x3)
pairs.panels(df,method = "pearson", hist.col = 'steelblue')
cor(df,meth='kendall')

help("pairs.panels")

#3D plot of the final dependence structure using x1, x2, x3 as marginals
plot3d(x1,x2,x3,pch=30,col='steelblue')



# bivariate density extracted with a T copula with t marginals
persp (copula_dist,dMvdc,xlim=c(-1,1),ylim=c(0,2),main="Multivariate density (extracted from a t Copula with t marginals)")

#Countor plot, measure of dependence, difficult to interpret
contour(copula_dist,dMvdc,xlim=c(-4,4),ylim=c(0,2))

#Simulating data from our Copula
sim <- rMvdc(4000, copula_dist)

ts.plot(sim)

###########################################################################################

#Real applications - 

#APPLE AND EXXON MULTIVARIATE DISTRIBUTION WITH A SURVIVAL GUMBEL AND T-STUDENT COPULA

#LIBRARIES
library(tidyquant)
library(quantmod)
library(VineCopula)
library(copula)
library(MASS)


#Start and end dates
start_date <- as.Date("2019-01-01")
end_date <- as.Date("2023-12-18")

#Download the data from Yahoo Finance
getSymbols("AAPL", from = start_date, to = end_date,auto.assign = TRUE)
getSymbols("XOM", from = start_date, to = end_date,auto.assign = TRUE)


#Extracting the Ad. Closed Prices
apple_price <- AAPL$AAPL.Adjusted
datax<-index(apple_price)
plot(datax, apple_price, type = "l", main = "Apple Price", col = "blue", xlab = "Date", ylab = "Adjusted Price")

exxon_price <- XOM$XOM.Adjusted
plot(datax, exxon_price,main = "Exxon Price",type = "l", col = "blue", xlab = "Date", ylab = "Adjusted Price")


# Computing daily returns
exxon_returns <- dailyReturn(exxon_price)
ts.plot(exxon_returns)

apple_returns <- dailyReturn(apple_price)
ts.plot(apple_returns)

#Join returns into a Dataframe
emp_data <- data.frame(apple_returns,exxon_returns)


#Scatter Plot with regression line
plot(emp_data,col=4, pch= 20, main = "Daily Returns", xlab = "Apple", ylab = "Exxon");grid();
abline(coefficients(mod1),col=2,lwd=3 );

pearson_data <- cor(exxon_returns, apple_returns, method='pearson')
kendall_data <- cor(exxon_returns, apple_returns, method = 'kendall')
c(pearson_data,kendall_data)

####FIRST STEP

####ESTIMATING THE MARGINAL PARAMETER FROM THE DATA WITH MLE

#Due to fat tails of returns, I reckon a suitable distribution to fit the data can be a t student
hist(apple_returns, col='steelblue',breaks = 200)
hist(exxon_returns, col='steelblue',breaks = 200, xlab = 'Exxon daily returns' , main =' Empirical Exxon returns')
fitted_t1 <- fitdistr(apple_returns,densfun = "t") #can be a log normal? no bcs ret are negative
fitted_t2 <- fitdistr(exxon_returns,densfun = 't')
print(fitted_t1)
print(fitted_t2)

#Simulation from our T student
set.seed(1233)
simulated_t1 <- rt(1000, df = 3.6, ncp = 0.00037)
hist(simulated_t1,col='salmon',breaks = 200)

set.seed(126)
simulated_t2 <- rt(2000, df = 4, ncp = 0.00037)
hist(simulated_t2,col='salmon',breaks = 200, xlab = 'Exxon daily simulated returns' , main =' Simulated Exxon returns from the estimated t distribution')


#SECOND STEP
#Estimate of Copula parameter!

#A) Transform your data in into uniform distributed values

#Pseudo observation are calculated as rank(X_i)/(n+1). They are in the [0,1] intervals
#This proces transforms the original data into uniform distributed values.
u <- pobs(as.matrix(cbind(exxon_returns,apple_returns)))[,1]
v <- pobs(as.matrix(cbind(exxon_returns,apple_returns)))[,2]


# B) CHOOSING THE BEST COPULA and ESTIMATE THE PARAMETERS - with BiCopSelelect()

#At first, all copulas are fitted with MLE, and we have a parameter estimates.
# Then BiCopSelelect function pins down the most appropriate copula family according to either AIC or BIC.
#We recall that the Survival Gumbel is a type of Archimedian Copula.

selectedCopula <- BiCopSelect(u,v,familyset=NA, selectioncrit = "BIC", method = "mle")
summary(selectedCopula)


##Two-dimensional examples : Density of suvival Gumbel(1.24)
survC <- rotCopula(gumbelCopula(1.22)) 
persp(survC, dCopula,xlim=c(0,1),ylim=c(0,1), main = "Survival Gumbel (1.22) density")


#3 STEP: APPLY SKLAR THEOREM- COPULA + MARGINALS = MULT DISTRIB

#MULTIVARIATE DISTIBUTION WITH SUVIVAL GUMBELL COPULA
plot(selectedCopula,type = "surface", size = 100L, xlim = c(-3,3), ylim =c(-2,2), zlim= c(0,0.20), margins = "norm", main = 'Survival Copula with estimated t-marginals')



#Example1:Sampling from a Survival Gumbel, here we can see there is left tail dependence
dat<-BiCopSim(35000,14, 1.24) 
u2<-dat[,1] 
v2<-dat[,2] 
plot(u2,v2,pch='.',col='blue', main ="Dependence structure of Survival Gumbel Copula with theta  = 1.24")



##NON FUNZIONA LA SIMULAZIONE DELLA SURV GUMBELL
#mv.NE <- mvdc(rotCopula(gumbelCopula(1.24),c("norm","exp"), list(list(mean= 0,sd=2),list(rate=2)))
#sim_SG_cop <- rMvdc(300,mv.NE)




#2 ESTIMATE A T-COPULA PARAMETER WITH THE KENDALL TAU.
#Employing kendall Tau to estimate the parameter.
#Two parameter families (other than the t-copula) cannot be handled by method "itau". 
#They are automatically removed from the familyset.
#For the t- copula the second parameter is found by a crude profile likelihood optimisation in (2,10]
#Here we want to fit a t-student copulas
selectedCopula <- BiCopSelect(u, v,selectioncrit = "BIC", method = "itau", familyset=c(1,2,3,5,6,7,8,9) )
summary(selectedCopula)

#Check if the parameter suggested by the function above are coherent - Double check is always useful in model validation

tcopula <- tCopula(dim=2)
pseudo_obs <- pobs(as.matrix(cbind(exxon_returns,apple_returns)))
est_para <- fitCopula(tcopula, pseudo_obs, method='ml')
coef(est_para)

#It is nice to see that the parameters of the fitted copula 
#are the same as those suggested by the BiCopSelect()
#Density of the t student copula with dCopula
rho <- coef(est_para)[1]
df <- coef(est_para)[2]

#change rho and see how it changes
persp(tCopula(dim = 2, rho, df=df),dCopula, main= "Density of a t stud Copula ρ = 0.25 with df=5 ")

 
 
 ##Construct a Bivariate Distribution whose marginals 
 ##are both t di student with df1= 4 and df2=3.6 respectively,coupled 
 ##together via a t copula with df=5
 
 copula_dist <- mvdc(copula=tCopula(rho,dim=2,df=5), margins = c("t","t"),
                     paramMargins=list(list(df = 3.6),
                                   list(df = 4)))
 
 # Bivariate density extracted with a T copula with t marginals
 persp (copula_dist,dMvdc,xlim=c(-4,4),ylim=c(-2,2),main=" T-stud Copula with t student marginals ")
 

 
 ###Comparing measure of Dependence before and after employing copulas
 
 #Measures of Dependence on original data
 kendall_data
 pearson_data
 
 #Measure of dependence after employing copula - stability of rank correlation
 cor(sim[,1],sim[,2],method='kendall')
 cor(sim[,1],sim[,2],method='pearson')
 
 
 
 #Simulating SYNTHETIC RETURNS from our T Copula
 set.seed(1123)
 sim_t_cop<- rMvdc(3000, copula_dist)
 
 #Copula value in percentual terms to get percentual returns
 sim_t_cop_n <- sim_t_cop/100
 
 #compare plots of syntethic returns and empirical stock returns
 plot(emp_data, xlab ='  Apple daily returns', ylab =' Exxon daily returns', main =' Empirical daily returns ')
 plot(sim_t_cop_n, xlab =' sim Apple daily returns', ylab ='sim Exxon daily returns', main =' Simulated returns from a T Copula with t marginals')
 
 
 #Overlapped Plot - Empirical vs Simulated
 plot(emp_data[[1]], emp_data[[2]], col = "blue", xlab ='Apple daily returns', ylab ='Exxon daily returns', main = 'Empirical vs Simulated returns' )
 points(sim_t_cop_n[, 1], sim_t_cop_n[, 2], col = 'red')
 # legend
 legend('bottomright', c('Observed', 'Simulated'), col = c('blue', 'red'), pch = 21)



 
 
 
 
 
 
 
 
 
