## ---------------------------
## Script name: VP_proj_creator
##
## Purpose of script: create the R6 class of a project
##
## Author: Thibaud Derippe
##
## Date Created: 2022-04-04 (day of code repartition from R6object.R)
##
## Under GPL-3 License
## Email: thibaud.derippe@gmail.com
## GitHub: https://github.com/Thibaudpmx/QSPVP
## ---------------------------
##
## Notes:
##
## see https://adv-r.hadley.nz/r6.html for more details on R6 objects
##
## ---------------------------


#' Initiate a R6 object project
#' @author Thibaud Derippe
#' @import RxODE tibble dplyr forcats crayon ggplot2 progress purrr tidyr stringr magrittr
#' @export


VP_proj_creator <- R6Class("VT",

                           public = list( model = NULL,

                                          param = NULL,
                                          filters_neg_above = data.frame(),
                                          filters_neg_below  = data.frame(),
                                          filters_pos_above = NULL,
                                          filters_pos_below = NULL,
                                          data = NULL,
                                          parameters_default_values = NULL,# c(psi = 20) # paremeters value to be used by default (NA if you want to force providing value)
                                          initial_cmt_values = NULL, #c(X1 = 50) # initial compartment values. At least one, and every missing cmt name would be set to 0
                                          times = NULL, #seq(0,52, 1)
                                          poolVP = data.frame(),
                                          errorComputation = data.frame(),
                                          protocols = NULL,
                                          param_reduce = NULL,
                                          param_increase = NULL,
                                          param_no_impact = NULL,
                                          targets = NULL,
                                          timesaver = NULL,
                                          zone_maybe = NULL,
                                          zone_sure = NULL,
                                          timeTrack = NULL,
                                          algo2list = list(),
                                          action_programmed = list(),

                                          initialize = function(sourcefile= "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config.r"){
                                            myEnv <- new.env()
                                            source(sourcefile, local=myEnv)
                                            # attach(myEnv, name="sourced_scripts")

                                            self$model <- myEnv$model_RxODE

                                            self$data <- myEnv$data_VT
                                            self$parameters_default_values <- myEnv$parameters_default_values
                                            self$initial_cmt_values <- myEnv$initial_cmt_values
                                            self$times <- myEnv$times
                                            self$protocols <- myEnv$protocols
                                            self$param_reduce <- myEnv$param_reduce
                                            self$param_increase <- myEnv$param_increase
                                            self$param_no_impact <- myEnv$param_no_impact
                                            self$filters_neg_above <- tibble()
                                            self$filters_neg_below <- tibble()
                                            self$poolVP <- tibble()
                                            self$zone_maybe <- tibble()
                                            self$zone_sure <- tibble()
                                            # get param

                                            param <- myEnv$model_RxODE$params
                                            param  <-  param[!param%in%  names(myEnv$parameters_default_values)]

                                            lines_model <- str_split(deparse(myEnv$model_RxODE$model), pattern = "\\\\n")[[1]] %>%
                                              gsub(pattern = "==", replacement = "nop") %>%
                                              gsub(pattern = "^structure\\(\"",replacement =  "")



                                            already_computed <- lines_model[grepl("=", lines_model)] %>%
                                              gsub(pattern = "=.+", replacement = "")

                                            param <- param[!param %in% already_computed]
                                            # param <- param[!param %in% gsub("^d", "", already_computed)]

                                            self$param <- param


                                          },
                                          print = function(...) {
                                            cat(green(paste0("Number of VP found: ", nrow(self$poolVP))))
                                            cat(red(paste0("\nNumber of red filter above: ", nrow(self$filters_neg_above))))
                                            cat(red(paste0("\nNumber of red filter below: ", nrow(self$filters_neg_below))))
                                            invisible(self)
                                          }

                           )
)
