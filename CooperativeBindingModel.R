#### Thermodynamic modeling and estimation of cooperative binding of transcription factors ###
#### Based on the Senear and Brenowitz model, 1991 ####

library("nls2")
library("ggplot2")

###### functions ##############################################################

cal_fbound <- function(abound0,abound1,abound2){
  sum_fbound = abound0 + abound1 + abound2
  fbound0 <- abound0 /sum_fbound
  fbound1 <- abound1/sum_fbound
  fbound2 <- abound2/sum_fbound
  list(Unbound = fbound0, Mono_bound = fbound1, Dimer_bound = fbound2)
}

cal_fbound_stat <- function(fbound,num_conc){
  fbound_r1 = fbound[1:num_conc] 
  if(num_reps == 1){
    fbound.mean = fbound_r1
    fbound.se = rep(0,num_conc)
  }
  else{
    fbound_r2 = fbound[(num_conc+1):(num_conc*num_reps)]
    fbound.mean = (fbound_r1 + fbound_r2) / 2
    fbound.se = NULL
    for(i in 1:length(fbound_r1)){ 
      fbound.se <- c(fbound.se, sd(c(fbound_r1[i], fbound_r2[i])/sqrt(2))) 
    }
  }
  list(fbound.mean = fbound.mean, fbound.se = fbound.se)
}

combined_fit_sumf1_f2_fhalf <- function(fbound0_fit, fbound1_2_fit, fboundhalf_fit, ki,k12,b,m,n){   
  conc = conc[1:num_conc]
  conc2 = conc*conc 
  Z = 1 + ((2*ki*conc)^n) + (ki^2*k12*conc2)^n
  bounds = num_conc*num_reps
  fbound0_fit_eq = fbound0_fit*(1/Z)
  fbound0_fit_eq = fbound0_fit_eq[1:bounds] 
  fbound1_2_fit_eq = fbound1_2_fit*(  ((2*ki*conc)^n + (ki^2*k12*conc2)^n) /Z )
  fbound1_2_fit_eq = fbound1_2_fit_eq[(bounds+1):(bounds*2)]  
  fboundhalf_fit_eq = fboundhalf_fit*(b + (m-b)/( 1+ (1/(2*ki*conc)^n)))
  fboundhalf_fit_eq = fboundhalf_fit_eq[((bounds*2)+1):(bounds*3)]
  return(c(fbound0_fit_eq,fbound1_2_fit_eq,fboundhalf_fit_eq)) 
}

#init values for nls model: 

run_nls <- function(data,fbound0_fit, fbound1_2_fit, fboundhalf_fit, ki,k12,b,m,n){
  com.model <- data ~ combined_fit_sumf1_f2_fhalf(fbound0_fit, fbound1_2_fit, fboundhalf_fit, ki,k12,b,m,n)
  fit <- nls(com.model, start=c(ki=0.009,k12=50,b=0.01,m =0.9,n=3),
             trace=TRUE,algorithm="port",
             lower=list(ki=0.0009,k12=1,b=0.01,m=0.8,n=1),
             control=nls.control(maxiter=1000)
             )
  return(fit)
}

define_fit <- function(num_conc){
  fbound0_fit <- rep(c(1,0,0), rep(num_conc*num_reps,3))
  fbound1_2_fit <- rep(c(0,1,0), rep(num_conc*num_reps,3))
  fboundhalf_fit <- rep(c(0,0,1), rep(num_conc*num_reps,3))
  list(fbound0_fit = fbound0_fit, fbound1_2_fit=fbound1_2_fit, fboundhalf_fit = fboundhalf_fit)
}

############### init values #################
num_conc = 10
conc <- c(8.2,11.09,14.97,20.22,27.29,36.85,49.75,67.16,90.67,122.4,165.25,223.09,301.17,406.5,548.88)
conc <- conc*0.9
num_reps = 2

#######   test data ############################
#ERR erre61 binding on a palindormic response element 
pal.abound0 <- c(12359.024,9405.853,11112.56,7658.439,7089.832,4066.296,1304.426,487.042,
                 rep(0,7),
                 12322.095,12232.995,10659.439,9768.196,9162.974,6318.317,1869.79,rep(0,8)
                 )
pal.abound1 <- c(1117.426,1715.74,2938.983,3069.64,4773.175,4038.276,2888.962,487.042, rep(0,7),
                  838.184,1242.548,1684.497,2483.205,3786.69,4601.447,3328.548,470.042, rep(0,7)
                  )
pal.abound2 <- c(0,153.021,389.163,495.163,1574.841,2678.497,6090.397,10138.731,11715.267,10958.116,
                  12456.936,10363.037,10405.765,9264.865,8770.087,
                  0,0,0,286.92,645.113,1985.719,5389.619,8532.983,10845.782,7359.146,9535.986,
                  12760.856,10259.271,9473.513,10074.291
                  )

##ERR binding on a half-site erre61 response element
half.abound0 <- c(15906.995,16339.945,14734.48,12946.702,7244.388,2468.69,344.506, rep(0,8),
                  8971.095,11141.731,8521.196,5886.903,4469.368,758.184,262.87, rep(0,8))
half.abound1 <- c(736.648,1316.719,2546.033,5028.882,8340.61,14309.388,14841.752,15196.702,
                  14197.187,15397.714,16212.806,15605.957,13712.291,13628.563,13860.978,
                  880.426,1885.861,2597.104,4285.347,7862.832,10138.731,12329.853,12057.752,
                  8895.803,6177.56,3880.024,5712.179,8408.383,7733.291,10947.546)

#######################################################################################

## get data in the proper format for nls fitting; this includes 2-site and half-site data

fbound_pal = cal_fbound(pal.abound0,pal.abound1,pal.abound2)
fbound_half = cal_fbound(half.abound0,half.abound1,rep(0,(num_conc*num_reps)) )

fbound0_pal.stat = cal_fbound_stat(fbound_pal$Unbound,num_conc)
fbound1_pal.stat = cal_fbound_stat(fbound_pal$Mono_bound,num_conc)
fbound2_pal.stat = cal_fbound_stat(fbound_pal$Dimer_bound,num_conc)
fbound1_half.stat = cal_fbound_stat(fbound_half$Mono_bound,num_conc)

combineddata <- c(fbound_pal$Unbound, fbound_pal$Mono_bound+fbound_pal$Dimer_bound, 
                          fbound_half$Mono_bound)


encode_fit_def <- define_fit(num_conc)

Globalfit <- run_nls(combineddata,
                     encode_fit_def$fbound0_fit,
                     encode_fit_def$fbound1_2_fit,
                     encode_fit_def$fboundhalf_fit,
                     ki,k12,b,m,n)

conc <- conc[1:num_conc]

Protein_summary <- data.frame(Protein= "AncERR",
                      RE=c(rep("EREpal.fbound0",num_conc),
                          rep("EREpal.fbound1_2",num_conc),
                         rep("EREhalf", num_conc)
                        ),
                     Conc = rep(conc,3),
                     frac_bound = c(fbound0_pal.stat$fbound.mean,
                                    fbound1_pal.stat$fbound.mean + fbound2_pal.stat$fbound.mean, 
                                    fbound1_half.stat$fbound.mean
                                    ),
                     CI = c(fbound0_pal.stat$fbound.se,
                            fbound1_pal.stat$fbound.se+fbound2_pal.stat$fbound.se,
                            fbound1_half.stat$fbound.se
                            )
                    )
conc = rep(conc,3)
p = ggplot(Protein_summary, aes(x=conc,y=frac_bound, colour=RE)) + geom_point(size=3) + scale_x_log10() + xlab(expression("AncERR log10[nM]")) + ylab(expression("Fraction bound")) 
colours <- c(EREhalf = "orange", EREpal.fbound0 = "red", 
             EREpal.fbound1_2 = "green")
p = p + aes(colour = RE) + scale_colour_manual(values = colours)

#draw out fits
r <- range(conc)
conc_grid <- seq(r[1],r[2],1)
yNew0 <- expand.grid(conc = conc_grid)
yNew1_2 <- expand.grid(conc = conc_grid)
yNewhalf <- expand.grid(conc = conc_grid)

ki = summary(Globalfit)$param[1]
k12 = summary(Globalfit)$param[2]
b = summary(Globalfit)$param[3]
m = summary(Globalfit)$param[4]
n = summary(Globalfit)$param[5]
Z = 1 + ((2*ki*conc_grid)^n) + ((ki^2*k12*conc_grid^2)^n)

yNew0$frac_bound = 1/Z
yNew0$RE = rep("EREpal.fbound0",length(conc_grid))
yNew1_2$frac_bound = (((2*ki*conc_grid)^n) + ((ki^2*k12*conc_grid^2)^n))/Z
yNew1_2$RE = rep("EREpal.fbound1_2",length(conc_grid))
yNewhalf$frac_bound = b + (m-b)/ (1+(1/(2*ki*conc_grid)^n))
yNewhalf$RE = rep("EREhalf",length(conc_grid))


p = p + geom_line(data = yNew0, colour = "red", size= 0.8) 
limits = aes(ymax = frac_bound + CI, ymin=frac_bound - CI)
p = p + geom_errorbar(limits, width = 0.08, size=0.4) 

p = p + geom_line(data = yNew1_2, colour = "green", size= 0.8) 
limits = aes(ymax = frac_bound + CI, ymin=frac_bound - CI)
p = p + geom_errorbar(limits, width = 0.08, size=0.4) 

p = p + geom_line(data = yNewhalf, colour = "orange", size= 0.8) 
limits = aes(ymax = frac_bound + CI, ymin=frac_bound - CI)
p = p + geom_errorbar(limits, width = 0.08, size=0.4) 

##generate expected binding curves for a given macroscopic binding constant ######
#[M] DBD
conc <- c(1e-12,5e-12,1e-11,5e-11,1e-10,5e-10, 1e-09, 5e-09, 1e-08, 5e-08, 1e-07, 5e-07, 1e-06, 5e-06, 1e-05, 5e-05, 1e-04, 5e-04, 1e-03,
          5e-03, 1e-02, 5e-02, 1e-01)
test = c(1:3)
#conc_grid = sapply(conc,function(x) x*test, simplify="Array")
conc_grid = NULL
for(i in 1:length(conc)){conc_grid = c(conc_grid,conc[i]*test)}

#following is Kds in M on 6,9,12 bp RE resp; ERREhalf, EREhalf, EREpal
Kd_mac_ancerr <- c(1/5e+06, 1/4.5e+07,1/7.47e+06,1/5,55e+07)

plot(log(conc,10),(1/(1+(Kd_mac_ancerr/conc))), col="magenta", xlab="log DBD [M]", ylab="E[Fraction bound]", main="AncSR1Chordate SR2 CTE Chimera")
lines(log(conc_grid,10),1/(1+(Kd_mac_ancerr[1]/conc_grid)), col="magenta")
points(log(conc,10),(1/(1+(Kd_mac_ancerr[2]/conc))), col="magenta", pch=16)
lines(log(conc_grid,10),(1/(1+(Kd_mac_ancerr[2]/conc_grid))), col="magenta")
points(log(conc,10),(1/(1+(Kd_mac_ancerr[3]/conc))), col="magenta", pch=22)
lines(log(conc_grid,10),(1/(1+(Kd_mac_ancerr[3]/conc_grid))), col="magenta",pch=22)
RE <- c("6bp","9bp","12bp")
legend(-4,0.5,pch=c(1,16,22),legend=RE, col="magenta")


#simulate frac bound ERE-pal, two-site
#Kd1 (uM) on EREpal for ancerr, ancnr3 and ancsr1
Kd1 <- c(0.134,0.0868,0.0376)
#Ka1 in [M]
Ka1 <- 1/(Kd1*10^-6)
#w on EREpal for ancerr,ancnr3 and ancsr1
w <- c(1,1.66,1.69)
monomer = 2*Ka1[1]*conc
dimer = w[1]*Ka1[1]^2*conc^2
Z = 1 + monomer + dimer
frac_bound_2site = (monomer + dimer)/Z
