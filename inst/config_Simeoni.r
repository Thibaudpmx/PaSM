
model_RxODE <-  RxODE({
  d/dt(Central) <- -ke *  max(Central,0)
  Conc <- Central/Vd

  tumVol <- X1 + X2 + X3 + X4
  growth <- lambda0 *  max(X1,0)/((1 + (lambda0 * max(tumVol,0)/lambda1)^psi)^(1/psi))

  X1(0) <- w0

  d/dt(X1) <- growth - X1 * max(Conc,0) * k2
  d/dt(X2) <- X1 * max(Conc,0) * k2 - k1 * X2
  d/dt(X3) <- k1 * (X2 - X3)
  d/dt(X4) <- k1 * (X3 - X4)

})



# paremeters value to be used by default (NA if you want to force providing value)
parameters_default_values <- c(psi = 20)

# initial compartment values. At least one, and every missing cmt name would be set to 0
initial_cmt_values <- c(X1 = 50)

# times you want to see your observations
times <- seq(0,52, 1)


# Protocols ----------------------------------

protocols <- list( dose0 = tibble(cmt = "Central", time = 0, amt = 0),
                   dose50 = tibble(cmt = "Central", time = 0, amt = 50),
                   dose100 = tibble(cmt = "Central", time = 0, amt = 100)

)

##### Fill the parameter #####

param_reduce<- list(tumVol = c("k2"), Conc = c("Vd", "ke")) #conc1 not here because....

param_increase <- list(tumVol = c("lambda0", "lambda1", "Vd", "ke", "w0" ), Conc = character())

param_no_impact <- list(tumVol = character(), Conc = c("lambda0", "lambda1", "k2", "k1" ))


# Data used -------------------------------------

# data shoud have at least ID, Value and concX columns, X being replace by drug number (one col by drug concentration)
# Avoid any column starting with "conc" if it is not a drug concentration / dose column

data_VT <- read.table("D:/Peccary_Annexe/Exemple_demo/DATA/Simeoni.txt", header = T, sep = ";", na.strings = ".") %>%
  mutate(protocol = paste0("dose", Dose), cmt = if_else(YTYPE == 2, "tumVol", "Conc")) %>%
  as_tibble %>%
  filter(!is.na(cmt))



