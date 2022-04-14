
# Administration ----------------------------------------------------------

model_RxODE <-RxODE({
  d/dt(Admin) <- -ka * Admin
  d/dt(Central) <- ka * Admin - ke * Central

})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,100, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "Admin", time = 0, amt = 50))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "ka", 0.01,0.1,1,
                  "ke", 0.01,0.1,1
                  )

find_relative(Central, protocol = "unique", model = model_RxODE, values = domain)


# homostastie ----------------------------------------------------------

model_RxODE <-RxODE({

  IRM(0) <- kin/kout

  d/dt(IRM) <- - kout * IRM + kin
  d/dt(IRM2) <- kout * IRM - kout2 * IRM2

})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,100, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "IRM", time = 0, amt = 0))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "kout", 0.01,0.1,1,
                  "kout2", 0.01,0.1,1,
                  "kin", 0.01,0.1,1
)

find_relative(IRM, IRM2, protocol = "unique", model = model_RxODE, values = domain,time_simul = times)
find_relative(IRM, IRM2, protocol = "unique", model = model_RxODE, values = domain,time_simul = times, deepAnalysis = T)
# find_relative(Central, protocol = "unique", model = model_RxODE, values = domain,time_simul = c(0,seq(5000,10000,100)))


# inhib kin------------------------------------------------------

model_RxODE <-RxODE({

IRM(0) <- kin/kout

d/dt(Drug) <- - ke * Drug
d/dt(IRM) <- - kout * IRM + kin * ( 1 - Emax * Drug / (Drug + D50))
d/dt(IRM2) <- kout * IRM - kout2 * IRM2

})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,200, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "Drug", time = 50, amt = 50))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "kout", 0.01,0.1,1,
                  "kout2", 0.01,0.1,1,
                  "kin", 0.01,0.1,1,
                  "ke",  0.01,0.1,1,
                  "Emax", 0.1,0.5,1,
                   "D50", 1,10,30
)

find_relative(IRM, IRM2, protocol = "unique", model = model_RxODE, values = domain,time_simul = times)


# increase kout------------------------------------------------------

model_RxODE <-RxODE({

  IRM(0) <- kin/kout

  d/dt(Drug) <- - ke * Drug
  d/dt(IRM) <- - kout * ( 1 + Emax * Drug / (Drug + D50)) * IRM + kin
  d/dt(IRM2) <- kout * IRM - kout2 * IRM2

})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,200, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "Drug", time = 50, amt = 50))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "kout", 0.01,0.1,1,
                  "kout2", 0.01,0.1,1,
                  "kin", 0.01,0.1,1,
                  "ke",  0.01,0.1,1,
                  "Emax", 0.1,0.5,1,
                  "D50", 1,10,30
)

find_relative(IRM, IRM2, protocol = "unique", model = model_RxODE, values = domain,time_simul = times)


# increase kin------------------------------------------------------

model_RxODE <-RxODE({

  IRM(0) <- kin/kout

  d/dt(Drug) <- - ke * Drug
  d/dt(IRM) <- - kout * IRM + kin * ( 1 + Emax * Drug / (Drug + D50))
  d/dt(IRM2) <- kout * IRM - kout2 * IRM2

})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,200, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "Drug", time = 50, amt = 50))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "kout", 0.01,0.1,1,
                  "kout2", 0.01,0.1,1,
                  "kin", 0.01,0.1,1,
                  "ke",  0.01,0.1,1,
                  "Emax", 0.1,0.5,1,
                  "D50", 1,10,30
)

find_relative(IRM, IRM2, protocol = "unique", model = model_RxODE, values = domain,time_simul = times)

# distribution ------------------------------------------------------------

##Pk.2comp.ke.k12.k21.difEq
#k12 <- Q/V1
#k21 <- Q/V2
#ke <- Cl/V1
model_RxODE <-RxODE({
d/dt(Central) <- k21 * Periph - k12 * Central - ke * Central
d/dt(Periph) <- - k21 * Periph + k12 * Central
Conc <- Central / V1

})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,100, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "Central", time = 0, amt = 50))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "k21", 0.01,0.1,1,
                  "k12", 0.01,0.1,1,
                  "ke",  0.01,0.1,1,
                  "V1", 1,5, 10)

find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain)


# Distribution 2 ----------------------------------------------------------



##Pk.2comp.ke.k12.k21.difEq

#ke <- Cl/V1
model_RxODE <-RxODE({
  k12 <- Q/V1
  k21 <- Q/V2
  d/dt(Central) <- k21 * Periph - k12 * Central - ke * Central
  d/dt(Periph) <- - k21 * Periph + k12 * Central
  Conc <- Central / V1

})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,100, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "Central", time = 0, amt = 50))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "V2", 1,5, 10,
                  "Q", 0.01,0.1,1,
                  "ke",  0.01,0.1,1,
                  "V1", 1,5, 10)

find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain)


# Distribution 3 ----------------------------------------------------------



##Pk.2comp.ke.k12.k21.difEq

#ke <- Cl/V1
model_RxODE <-RxODE({



  k21 <- ratio / k12
  d/dt(Central) <- k21 * Periph - k12 * Central - ke * Central
  d/dt(Periph) <- - k21 * Periph + k12 * Central
  Conc <- Central / V1

})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,100, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "Central", time = 0, amt = 50))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "ratio", 0.01,1,10,
                  "k12", 0.01,0.1,1,
                  "ke",  0.01,0.1,1,
                  "V1", 1,5, 10)

find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain)


# Distribution 4 ----------------------------------------------------------





##Pk.2comp.ke.k12.k21.difEq

#ke <- Cl/V1
model_RxODE <-RxODE({



  k21 <- ratio / k12
  d/dt(Central) <- k21 * Periph - k12 * Central - ke * Central
  d/dt(Periph) <- - k21 * Periph + k12 * Central
  Conc <- Central / V1

})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,100, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "Central", time = 0, amt = 50))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "ratio", 0.01,1,10,
                  "k12", 0.01,0.1,1,
                  "ke",  0.01,0.1,1,
                  "V1", 1,5, 10)

find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain)


# dos elim ------------------------------------------------------------

##Pk.2comp.ke.k12.k21.difEq
#k12 <- Q/V1
#k21 <- Q/V2
#ke <- Cl/V1
model_RxODE <-RxODE({
  d/dt(Central) <- k21 * Periph - k12 * Central - ke * Central
  d/dt(Periph) <- - k21 * Periph + k12 * Central - ke2 * Periph
  Conc <- Central / V1

})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,100, 2) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "Central", time = 0, amt = 50))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "k21", 0.01,0.1,1,
                  "k12", 0.01,0.1,1,
                  "ke",  0.01,0.11,1,
                  "ke2",  0.01,0.1,1,
                  "V1", 1,5, 10)

find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain)
find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain, deepAnalysis = T)


#ke > ke2 -> k12 / k21 non utilisable
domain <- tribble(~param, ~min, ~ref,  ~max,
                  "k21", 0.01,0.1,1,
                  "k12", 0.01,0.1,1,
                  "ke",  0.01,0.1,1,
                  "ke2",  0,0.001,0.1,
                  "V1", 1,5, 10)
find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain)

#ke2 > ke -> k12 / k21 utilisable
domain <- tribble(~param, ~min, ~ref,  ~max,
                  "k21", 0.01,0.1,1,
                  "k12", 0.01,0.1,1,
                  "ke",  0.01,0.1,1,
                  "ke2",  0.01,1,10,
                  "V1", 1,5, 10)
find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain)

#ke2 = ke -> k12 / k21 utilisable
domain <- tribble(~param, ~min, ~ref,  ~max,
                  "k21", 0.01,0.1,1,
                  "k12", 0.01,0.1,1,
                  "ke",  0.01,0.1,1,
                  "ke2", 0.01,0.1,1,
                  "V1", 1,5, 10)
find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain,time_simul = seq(0,500,1))



# Michaelis elimination ---------------------------------------------------



model_RxODE <-RxODE({


    d/dt(Central) <-  - Central * Vmax / (Central + km)
 })


parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Central = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,52, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "Central", time = 0, amt = 50))
# tribblecreator(model_RxODE)

domain <- tribble(~param, ~min, ~ref,  ~max,
                  "Vmax", 0.1,1,10 ,
                  "km", 1,5,25)

find_relative(Central, protocol = "unique", model = model_RxODE, values = domain)


# TMDD --------------------------------------------------------------------

model_RxODE <- RxODE({

  R(0) <- ksyn/kdeg

  d/dt(L) <- -ke * L - kon * L * R + koff * P - k12 * L + k21 * A
  d/dt(R) <- ksyn - kdeg * R - kon * L * R + koff * P
  d/dt(P) <- kon * L * R - koff * P - kint * P
  d/dt(A) <- k12 * L - k21 * A
  Conc <- L/Vd
  LP <- (L + P) / Vd
})

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(R = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,200, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "L", time = 0, amt = 50))
# tribblecreator(model_RxODE)

#kint > ke
domain <- domain <- tribble(~param, ~min, ~ref,  ~max,
                            "ke", 0.01,0.1,1,
                            "kon",  0.01,0.1,1,
                            "koff", 0.01,0.1,1 ,
                            "k12", 0.01,0.1,1 ,
                            "k21",  0.01,0.1,1 ,
                            "ksyn",  0.01,0.1,1,
                            "kdeg", 0.01,0.1,1 ,
                            "kint", 0,0.2,10 ,
                            "Vd", 1,5,10)


find_relative(LP,Conc, protocol = "unique", model = model_RxODE, values = domain)

# ke > kint
domain <- domain <- tribble(~param, ~min, ~ref,  ~max,
                            "ke", 0.01,0.6,1,
                            "kon",  0.01,0.1,1,
                            "koff", 0.01,0.1,1 ,
                            "k12", 0.01,0.1,1 ,
                            "k21",  0.01,0.1,1 ,
                            "ksyn",  0.01,0.1,1,
                            "kdeg", 0.01,0.1,1 ,
                            "kint", 0,0.1,10 ,
                            "Vd", 1,5,10)


find_relative(LP,Conc, protocol = "unique", model = model_RxODE, values = domain)

# ke > kint
domain <- domain <- tribble(~param, ~min, ~ref,  ~max,
                            "ke", 0.01,0.1,1,
                            "kon",  0.01,0.1,1,
                            "koff", 0.01,0.1,1 ,
                            "k12", 0.01,0.1,1 ,
                            "k21",  0.01,0.1,1 ,
                            "ksyn",  0.01,0.1,1,
                            "kdeg", 0.01,0.1,1 ,
                            "kint", 0,0.01,10 ,
                            "Vd", 1,5,10)


find_relative(LP,Conc, protocol = "unique", model = model_RxODE, values = domain)


# ke == kint
domain <- domain <- tribble(~param, ~min, ~ref,  ~max,
                            "ke", 0.01,0.1,1,
                            "kon",  0.01,0.1,1,
                            "koff", 0.01,0.1,1 ,
                            "k12", 0.01,0.1,1 ,
                            "k21",  0.01,0.1,1 ,
                            "ksyn",  0.01,0.1,1,
                            "kdeg", 0.01,0.1,1 ,
                            "kint", 0.01,0.1,1 ,
                            "Vd", 1,5,10)


find_relative(LP,Conc, protocol = "unique", model = model_RxODE, values = domain)



# QE ----------------------------------------------------------------------


dLtot <- - kin * Ltot - (kel - kint) * L
dRtot <- ksyn - kdeg * Rtot - (kint - kdeg) * ( Ltot - L )
L <- 1 / 2 * ( ( Ltot - Rtot - Kd ) + sqrt(( Ltot - Rtot - kd ) ** 2 + 4 * kd * Ltot ))
Conc_plot <- L / Vd


# QSS ---------------------------------------------------------------------

dLtot <- - kin * Ltot - (kel - kint) * L
dRtot <- ksyn - kdeg * Rtot - ( kint - kdeg ) * ( Ltot - L )
L <- 1 / 2 * ( ( Ltot - Rtot - Kss ) + sqrt(( Ltot - Rtot - kss ) ** 2 + 4 * kss * Ltot ))
Conc_plot <- L / Vd


# Rtot --------------------------------------------------------------------
model_RxODE <- RxODE({

dL <- - ( kel + kon * R0 ) * L + ( koff + kon * L ) * P
dP <- kon * R0 * L - ( koff + kint + kon * L ) * P
Conc_plot <- L / Vd

})

# Ireversible binding -----------------------------------------------------
model_RxODE <- RxODE({

  R(0) <- ksyn/kdeg

d/dt(L) <- - kel * L - kon * L * R
d/dt(R) <- ksyn - kdeg * R - kon * L * R
d/dt(P) <- kon * L * R - kint * P
Conc<- L / Vd
})


parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(L = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,200, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "L", time = 0, amt = 50))
# tribblecreator(model_RxODE)

domain <- domain <- tribble(~param, ~min, ~ref,  ~max,
                            "kel", 0.01,0.1,1,

                            "kon", 0.01,0.1,1 ,
                            "ksyn",  0.01,0.1,1 ,
                            "kdeg",  0.01,0.1,1 ,
                            "kint",  0.01,0.1,1,
                            "Vd", 1,5,10)
find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain, sensitivity = 0.001)



# TMDD MM--------------------------------------------------------------------
model_RxODE <- RxODE({

d/dt(L) <- - ke * L - Vm * L / ( kd + L ) - k12 * L + k21 * A
d/dt(A) <- k12 * L - k21 * A
Conc <- L / Vd

})


parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(L = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,200, 1) # times you want to see your observations
protocols <- list( unique = tibble(cmt = "L", time = 0, amt = 50))
# tribblecreator(model_RxODE)

domain <- domain <- tribble(~param, ~min, ~ref,  ~max,
                            "ke", 0.01,0.1,1,

                            "k12", 0.01,0.1,1 ,
                            "k21",  0.01,0.1,1 ,
                            "Vm", 0.1,1,10,
                            "kd",  1,5,25,
                            "Vd", 1,5,10)
find_relative(Conc, protocol = "unique", model = model_RxODE, values = domain)




