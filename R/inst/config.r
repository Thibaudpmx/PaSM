
# Part 1: write the model equations ---------------------------------------

#### Write the equations of your model, followng RxODE format


model_RxODE <- RxODE({

  # write your equations here !
  # Example (to replace):
  d/dt(Venetoclax) <- -ke_Venetoclax * Venetoclax
  #

})


## Verification: perform verification number1


# Part 2: parameters, initial states and time measurement ---------------------------------------

## You need to fill the following section.
## Easiest way is to first eval "model_extract()" and copy paste the output
## To directly have pre-fille the right parameter defaults values and initial states
## You can of course do modifications if needed

parameters_default_values <- c(ke_Venetoclax = 0.3) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(Venetoclax = 0) # initial compartment values. At least one, and every missing cmt name would be set to 0
times <- seq(0,48, 1) # times you want to see your observations

## Verification: perform verification number2


# Part 3: determine the fate of the cell ----------------------------------

# res will be the results of the simulations (do "res <- simulations()" to help fill the file)
# criteria should be an expression working with res with final output being a TRUE of FALSE
# TRUE being cell death, FALSE being cell survival

criteria <- expr ({

  min(res$Venetoclax[res$time <6]) > 1E-5

})

## Verification: perform verification number3

# parametre / conc qui s'ils sont augmenté et mort alors morts, si diminué et survie alors survie
param_death <- c("conc1", "conc2", "conc3", "conc4")

# parametre / conc qui s'ils sont augmenté et mort alors morts, si diminué et survie alors survie
param_survive <- c("ke_Venetoclax")



# Part 4: data and concentration to test -------------------------------------

# data shoud have at least ID, Value and concX columns, X being replace by drug number (one col by drug concentration)
# Avoid any column starting with "conc" if it is not a drug concentration / dose column

data_VT <- read.table("D:/these/Second_project/QSP/VirtualTumor/datademo_rework.csv", sep = ";", header = T) %>%
  as_tibble


# Normally the following code will automatically detect numbers of drug and
# extract the different concentration levels

ndrug <- sum(grepl("^conc", names(data_VT)))

for(a in 1:ndrug){

drug <- paste0("conc", a)

expr(!!parse_expr(drug) <- unique(data_VT[[!!drug]])) %>%
  eval

}

# You can always modify manually with such code
# ndrug <- 4
# conc1 <- c(0,0.08,0.16,0.32,0.64,1.3,2.60,5,10,20)
# conc2 <- c(0,0.08,0.16,0.32,0.64,1.3,2.60,5,10,20)
# conc3 <- c(0,0.08,0.16,0.32,0.64,1.3,2.60,5,10,20)
# conc4 <- c(0,5,10,15)





