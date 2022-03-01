library(peccary)
library(QSPVP)
library(RxODE)
library(progress)
library(R6)
library(crayon)
# Step 1 create or load project ---------------------------------------------------
create_VT_Project("D:/these/Second_project/QSP/modeling_work/VT_simeoni")

# Step 2: fulfil the file

shell.exec("D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config.r")
shell.exec("D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/2_verif.r")

# test
ind_param = tibble(k1 = 3, k2 = 0.4, ke = 27.5/52.8, lambda0 = 0.06, lambda1 = 6, Vd = 52.8)


celltheque


dose = 100
add_events <- tibble(cmt = c("Central")) %>%  # Administration sampling, with "concX" being replaced by concentration
  mutate(time = 0, amt = dose)

#
res <- simulations(ind_param = ind_param, add_events = add_events,returnSim = T)
simulations(ind_param = ind_param, add_events = add_events,returnSim = T) %>%
  filter(time >23)
simulations(ind_param = ind_param, add_events = add_events,returnSim = F)
#
# res %>%
#   gather("key", "value", Conc, tumVol) %>%
#   ggplot()+
#   geom_line(aes(time, value))+
#   facet_wrap(~key, scales = "free")+
#   scale_y_log10()
read.table("D:/Peccary_Annexe/Exemple_demo/Simeoni/closeIV/IndividualParameters/estimatedIndividualParameters.txt", header = T, sep = ",") %>%
  as_tibble %>%
  filter(id == 218)

# Step3: create celltheque
 #%>%
  # slice(1:(3*6000))# what you want to add in your celltheque as individuals


# targets <-  tribble(~protocols, ~time, ~cmt, ~ min, ~max,
#                     "dose0",10,"tumVol", 50, 120,
#                     "dose0",25,"tumVol", 70,220,
#                     "dose0",0,"Conc",0,1,
#                     "dose50",10,"tumVol", 40,125,
#                     "dose50",25,"tumVol", 30,70,
#                     "dose50",10,"Conc", 1e-6,1,
#                     "dose100",25,"tumVol", 3,167)


# Two compart ---------------------------------------------------------------


source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()



self$set_targets(filter = cmt == "tumVol", ntime = 6)
self$set_targets(filter = Dose ==50 ,timeforce = c(12,19, 30,45))

VP_df <- crossing(k1 = c(0.5),
         k2 = seq(0,8,0.2),
         ke = 1 ,#*  seq(0.6,1.4,0.2),
         lambda0 =seq(0,0.16,0.025),
         lambda1 = c(12),
         Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

self$add_VP(VP_df, fillatend = F, reducefilteratend = T, testnofV = F)



self$plot_VP()

self$plot_2D(x = k2, y = lambda0, plotMain = T)

self$plot_2D(x = k2, y = lambda0,toaddneg = VP_df, plotMain = T)
self$n_filter_reduc()
# self$fill_simul()
self$filters_neg_above
self$filters_neg_below
self$filters_pos_below
self$filters_pos_above


# 2cm big -----------------------------------------------------------------




source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()



self$set_targets(filter = Dose == 50 & cmt == "tumVol", ntime = 8)
self$set_targets(filter = Dose==50 & cmt == "tumVol",timeforce = c(12,19, 30,45))

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,8,0.01),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,0.16,0.0025),
                  lambda1 = c(12),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )



self$add_VP(VP_df, fillatend = F, reducefilteratend = F,  npersalve = 2000,  time_compteur = F)

self$n_filter_reduc()
self$plot_VP()

self$targets

self$poolVP %>%
  unnest(simul) %>%
  filter(time == 12 & tumVol < 21) %>%
  pull(cellid) -> idtarget

self$plot_2D(x = k2, y = lambda0, plotMain = F)

self$n_filter_reduc()

ti <- Sys.time()
self$fill_simul()

difftime( Sys.time(), ti, units = "sec")
self$filters_neg_above
self$filters_neg_below
self$filters_pos_above
self$filters_pos_below



# three compart ---------------------------------------------------------------



source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()



# self$set_targets(filter = Dose == 50 & cmt == "tumVol", ntime = 8)
self$set_targets(filter = Dose == 50 & cmt == "tumVol",timeforce = c(12,19, 30,45))



VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,8,0.5),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0.04,0.16,0.025),
                  lambda1 = c(12),
                  Vd =  c(20:40)) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

# VP_df <- crossing(k1 = c(0.5),
#                   k2 = seq(0,8,0.1),
#                   ke = 1 ,#*  seq(0.6,1.4,0.2),
#                   lambda0 =seq(0,0.16,0.01),
#                   lambda1 = c(12),
#                   Vd =  c(0:40)) %>% #c(0.8,1,1.2)) %>%
#   map_df(function(x){
#
#     if(is.character(x)) return(x)
#     round(x,3)
#
#   } )

self$add_VP(VP_df, fillatend =F, reducefilteratend = F)

self$plot_VP()
self$n_filter_reduc()
# self$plot_2D(k2, lambda0, plotMain = T)
self$plot_3D(k2, lambda0, Vd)

self$filters_neg_above
self$filters_neg_below
self$filters_pos_above
self$filters_pos_below


# 4 dimensions and more

# Infinite dimension ------------------------------------------------------


source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()
#
# saveRDS(self, "D:/these/Second_project/QSP/modeling_work/VT_simeoni/longrun.RDS")

# self$set_targets(filter = Dose == 50 & cmt == "tumVol", ntime = 8)
self$set_targets(filter = Dose == 50 & cmt == "tumVol",timeforce = c(12,19, 30,45))

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,8,0.1),
                  ke = 1 ,#*  seq(0.6,1.4,0.2),
                  lambda0 =seq(0,0.16,0.01),
                  lambda1 = c(10,12,14,33),
                  Vd =  c(0:40)) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

t0 <- Sys.time()
self$add_VP(VP_df, fillatend = F, reducefilteratend = F,use_green_filter = F, npersalve = 2000, time_compteur = F, pctActivGreen = 0.75)
Sys.time() - t0

self$plot_VP()
self$n_filter_reduc()

self$filters_neg_above
self$filters_neg_below
self$filters_pos_above
self$filters_pos_below



# PK AND PD ---------------------------------------------------------------





source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")


self <- VP_proj_creator$new()



# self$set_targets(filter = cmt == "tumVol", ntime = 6)
self$set_targets(filter = Dose ==50 ,timeforce = c(3,19, 30,45))

tar <- self$targets ; tar$min[[5]] <- 0.005; tar$max[[2]] <-100; tar$max[[3]] <- 200; tar$max[[5]] <- 300
self$targets  <- tar

VP_df <- crossing(k1 = c(0.5),
                  k2 = seq(0,8,0.2),
                  ke = seq(0,3,0.2) ,#*  seq(0.6,1.4,0.2),
                  lambda0 = c(0.025),
                  lambda1 = c(12),
                  Vd =  40) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

self$add_VP(VP_df, fillatend = F, reducefilteratend =F, use_green_filter = F)

self$plot_VP()
self$poolVP
# rxode n optimal ---------------------------------------------------------


#>     intersect, setdiff, setequal, union
ev1 <- eventTable(amount.units="mg", time.units="hours") %>%
  add.dosing(dose=10000, nbr.doses=1, dosing.to=2) %>%
  add.sampling(seq(0,48,length.out=10));

ev2 <- eventTable(amount.units="mg", time.units="hours") %>%
  add.dosing(dose=5000, nbr.doses=1, dosing.to=2) %>%
  add.sampling(seq(0,48,length.out=8));

dat <- rbind(data.frame(ID=1, ev1$get.EventTable()),
             data.frame(ID=2, ev2$get.EventTable()))


rxSolve(model_RxODE, dat )
 rxSolve(mod, theta, dat, omega=omega, sigma=sigma)




 mod <- RxODE({
   eff(0) = 1
   C2 = centr/V2;
   C3 = peri/V3;
   CL =  TCl ## This is coded as a variable in the model
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
   e = eff
   cp = centr*(1)
 })

 theta <- c(KA=2.94E-01, TCl=1.86E+01, V2=4.02E+01,  # central
            Q=1.05E+01, V3=2.97E+02,                # peripheral
            Kin=1, Kout=1, EC50=200)                # effects
 omega <- lotri(eta.Cl ~ 0.4^2)
 sigma <- lotri(eff.err ~ 0.1, cp.err ~ 0.1)
 tmp <- matrix(rnorm(8^2), 8, 8)
 tMat <- tcrossprod(tmp, tmp) / (8 ^ 2)
 dimnames(tMat) <- list(NULL, names(theta))
 ev <- et(amount.units="mg", time.units="hours") %>%
   et(amt=10000, cmt="centr")

 a <- Sys.time()
 sim  <- rxSolve(mod, theta, ev, omega=omega, nSub=100, sigma=sigma,
                 thetaMat=tMat, nStud=10,
                 simVariability=FALSE)
 difftime(Sys.time(), a)



 # obs <- tibble


 theta %>% as.data.frame()

 ### 1
 aa_1 <- tibble( KA = 2,  TCl  = 2,    V2  = 2,     Q  = 2,    V3  = 2,   Kin = 2,   Kout  = 2,  EC50= 2) %>%
   crossing(id = 1:1)
 ev3_1 <- ev2 %>% select(-id) %>% crossing (id = 1:1) %>% arrange(id) %>%  select(id, everything())

 test1 <- function(){


 # a <- Sys.time()
 b <- mod$solve(aa_1, ev3_1)
 # difftime(Sys.time(), a)
 }

 ### 10
 aa_10 <- tibble( KA = rnorm(10, 2,0.1),  TCl  =  rnorm(10, 2,0.1),    V2  =  rnorm(10, 2,0.1),     Q  =  rnorm(10, 2,0.1),
                  V3  =  rnorm(10, 2,0.1),   Kin =  rnorm(10, 2,0.1),   Kout  =  rnorm(10, 2,0.1),  EC50=  rnorm(10, 2,0.1)) %>%
   rowid_to_column("id")
 ev3_10 <- ev2 %>% select(-id) %>% crossing (id = 1:10) %>% arrange(id) %>%  select(id, everything())

 test10 <- function(){


   # a <- Sys.time()
   b <- mod$solve(aa_10, ev3_10)
   # difftime(Sys.time(), a)
 }

 b %>%
   as_tibble() %>%
   ggplot()+
   geom_line(aes(time, C2, group= factor(id)))+
   scale_y_log10()
 ### 100
 aa_100 <- tibble( KA = 2,  TCl  = 2,    V2  = 2,     Q  = 2,    V3  = 2,   Kin = 2,   Kout  = 2,  EC50= 2) %>%
   crossing(id = 1:100)
 ev3_100 <- ev2 %>% select(-id) %>% crossing (id = 1:100) %>% arrange(id) %>%  select(id, everything())

 test100 <- function(){


   # a <- Sys.time()
   b <- mod$solve(aa_100, ev3_100)
   # difftime(Sys.time(), a)
 }


 ### 1000
 aa_1000 <- tibble( KA = 2,  TCl  = 2,    V2  = 2,     Q  = 2,    V3  = 2,   Kin = 2,   Kout  = 2,  EC50= 2) %>%
   crossing(id = 1:1000)
 ev3_1000 <- ev2 %>% select(-id) %>% crossing (id = 1:1000) %>% arrange(id) %>%  select(id, everything())

 test1000 <- function(){


   # a <- Sys.time()
   b <- mod$solve(aa_1000, ev3_1000)
   # difftime(Sys.time(), a)
 }


 ### 1000
 aa_10000 <- tibble( KA = 2,  TCl  = 2,    V2  = 2,     Q  = 2,    V3  = 2,   Kin = 2,   Kout  = 2,  EC50= 2) %>%
   crossing(id = 1:10000)
 ev3_10000 <- ev2 %>% select(-id) %>% crossing (id = 1:10000) %>% arrange(id) %>%  select(id, everything())

 test10000 <- function(){


   # a <- Sys.time()
   b <- mod$solve(aa_10000, ev3_10000)
   # difftime(Sys.time(), a)
 }


 test(3)

 library(microbenchmark)


 mb <- microbenchmark(test1(),test10(), test100(), test1000(), times = 30L, unit = "s"); mb
 plot(mb)

 microbenchmark(test10000(), times = 30L, unit = "s")


 tribble(~nsimn, ~mean,
         1,0.0027,
         10,0.0030,
         100, 0.0056,
         1000,0.035,
         2000, 0.539) %>%
   mutate(per_sim = mean / nsimn) %>%
   ggplot()+
   geom_line(aes(nsimn, per_sim))+
   scale_y_log10()+
   scale_x_log10()

 tribble(~nsimn, ~mean,
         1,0.0027,
         10,0.0030,
         100, 0.0056,
         1000,0.035,
         2000, 0.539) %>%
   mutate(per_sim = mean / nsimn) %>%
   ggplot()+
   geom_line(aes(nsimn, mean))+
   scale_x_log10()
   scale_y_log10()


# test two compartment whole time -----------------------------------------



   model <- model_RxODE <-  RxODE({
     d/dt(Central) <- -ke * Central
     Conc <- Central/Vd

     tumVol <- X1 + X2 + X3 + X4
     growth <- lambda0 * X1/((1 + (lambda0 * tumVol/lambda1)^20)^(1/20))
     X1(0) <- 50
     d/dt(X1) <- growth - X1 * Conc * k2
     d/dt(X2) <- X1 * Conc * k2 - k1 * X2
     d/dt(X3) <- k1 * (X2 - X3)
     d/dt(X4) <- k1 * (X3 - X4)
   })


   VP_df <- crossing(k1 = c(0.5),
                     k2 = seq(0,8,0.2),
                     ke = 1 ,#*  seq(0.6,1.4,0.2),
                     lambda0 =seq(0,0.16,0.025),
                     lambda1 = c(12),
                     Vd =  40) %>%
     rowid_to_column("id")

   VP_df <-  crossing(k1 = c(0.5),
            k2 = seq(0,8,0.01),
            ke = 1 ,#*  seq(0.6,1.4,0.2),
            lambda0 =seq(0,0.16,0.0025),
            lambda1 = c(12),
            Vd =  40)%>%
     rowid_to_column("id")

   VP_df <-  crossing(k1 = c(0.5),
                      k2 = 4,
                      ke = 1 ,#*  seq(0.6,1.4,0.2),
                      lambda0 =0.16,
                      lambda1 = c(12),
                      Vd =  1:1000)%>%
     rowid_to_column("id")

  eve <-  tibble(cmt = "Central", time = 0, amt = 50, evid = 1) %>%
     bind_rows(tibble(cmt = NA, time = seq(0,54,1), amt = NA, evid = 0)) %>%
    crossing(id = VP_df$id)


  b <- model$solve(VP_df, eve)


  mb2 <- microbenchmark( model$solve(VP_df, eve), unit = "s", times = 30);mb2

  ####♠ Il prend 0.15 sec pour 1000 -> 30 sec
 30/6

  mb2 <- microbenchmark( model$solve(VP_df, eve), unit = "s")
     crossing(id = 1:287)

     bind_rows(mb, mb2)

found_opti_n_pack_simul <- function(model_RxODE){

params <-   model_RxODE$params

matrix

}



# Lindner !  --------------------------------------------------------------



source("D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner.r")


# Defining targets


source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner.r")

self$targets <- tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                        "unique",740,"Pore", 20, 23
                  )


VP_df <- crossing(Bcl20 = seq(100,1000,100),
                  Bclxl0 = seq(100,1000,100),
                  Mcl10 = 50 ,#*  seq(0.6,1.4,0.2),
                  BIM0 = seq(100,1000,100),
                  PUMA0 =seq(100,1000,100),
                  NOXA0 =  seq(100,1000,100),
                  BAXc0 = 1000,
                  BAK0 = 1000) %>% #c(0.8,1,1.2)) %>%
  map_df(function(x){

    if(is.character(x)) return(x)
    round(x,3)

  } )

self$add_VP(VP_df, fillatend = F, reducefilteratend = F)



self$plot_VP()

self$plot_2D(x = k2, y = lambda0, plotMain = T)

self$plot_2D(x = k2, y = lambda0,toaddneg = VP_df, plotMain = T)


# Massive screening -------------------------------------------------------



source("D:/these/Second_project/QSP/QSPVP/R/R6object.R")
self <- VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner.r")

self$targets <- tribble(~protocol, ~time, ~cmt, ~ min, ~max,
                        "unique",740,"Pore", 20, 20.0001
)


domain <- tribble(~param, ~from, ~to, ~step,
                  "Bcl20", 0, 1000, 1 ,
                  "Bclxl0", 0, 1000, 1,
                  "Mcl10", 0,1000,1,
                  "BIM0", 0, 1000, 1 ,
                  "PUMA0", 0, 1000, 1,
                  "NOXA0", 0,1000,1,
                  "BAXc0", 1000,1000,0,
                  "BAK0" , 1000,1000,0

)

demo <- self$big_screening(domain)


demo2 <- self$big_screening(domain = demo[[1]])

demo3 <- self$big_screening(domain = demo2[[1]])

demo3 %>%
  unnest() %>%
  filter(time == 740) %>%
  pull(Pore)
