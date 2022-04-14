
# Model ---------------------------------------

#### Write the equations of your model, followng RxODE format


# model <- model_RxODE <-  RxODE({
#   d/dt(Central) <- -ke * Central
#   Conc <- Central/Vd
#   
#   tumVol <- X1 + X2 + X3 + X4
#   growth <- lambda0 * X1/((1 + (lambda0 * tumVol/lambda1)^psi)^(1/psi))
#   
#   X1(0) <- w0
#   
#   d/dt(X1) <- growth - X1 * Conc * k2
#   d/dt(X2) <- X1 * Conc * k2 - k1 * X2
#   d/dt(X3) <- k1 * (X2 - X3)
#   d/dt(X4) <- k1 * (X3 - X4)
# 
# })



model <- model_RxODE <-  RxODE({
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

# Various ---------------------------------------

parameters_default_values <- c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(X1 = 50) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,52, 1) # times you want to see your observations


# Protocols ----------------------------------

protocols <- list( dose0 = tibble(cmt = "Central", time = 0, amt = 0),
                   dose50 = tibble(cmt = "Central", time = 0, amt = 50),
                   dose100 = tibble(cmt = "Central", time = 0, amt = 100)
  
)


# Parameter influence -----------------------------------------------------


#### HELPER - TO COMMENT AFTER USE ############


# Step1: launch the followind command (tribblecreator) 
#               tribblecreator(model_RxODE)

# Step2 copy paste code printed in the console, then 
# 1) fill the domain table, 
# 2) replace your_output_without_quotes by your outputes (eg. tumVol, Conc, ...), and give a protocol
# name provided in protocols

 ### replace this line by the code procuded in step1 
    
# domain <- tribble(~param, ~min, ~ref,  ~max,
#                   "ke", 0.1,0.5,2,
#                   "Vd", 10,30,50 ,
#                   "lambda0", 0.02,0.05,1 ,
#                   "lambda1", 20,50,100,
#                   "psi", 10,20,30 ,
#                   "k2", 0.05, 1,2,
#                   "k1", 0.1,1,2 )
# 
# find_relative(tumVol, Conc, protocol = "dose50", model = model_RxODE, values = domain)
 
# Step3: verify the plot, the table summarising the influence, and if you find it
# plausible copy paste the three line (param_reduce, param_increase, param_no_impact)
# in the next bloc. Feel free to modify manually if needed
#### END HELPER - Please comment the full bloc above ############

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


