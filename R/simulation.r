# model_equations <- function(){
#
#   test <-  exists("model_RxODE")
#
#   if(! test){
#
#
#
#
#   }
#
#   return(model_RxODE)
#
# }


# source("D:/these/Second_project/QSP/VirtualTumor/inst/model.r")


# simulations(add_events = add_events)

#' Title
#'
#' @param ind_param
#' @param add_events
#'
#' @return
#' @export
#'
#' @examples
#'
# simulations(ind_param = celltheque %>% slice(1), add_events = add_events) %>%
#   ggplot()+
#   geom_line(aes(time, Bcl2_I))
simulations <- function(ind_param = double(), add_events = tibble(), returnSim = T){


  events <-  tibble(cmt = names(initial_cmt_values)[[1]], time = 0, amt = 0) %>%
    bind_rows(add_events) %>%
    mutate(evid = 1) %>%
    bind_rows(tibble(time = times, evid = 0, cmt = c(NA)))



  parameter <- parameters_default_values[!parameters_default_values %in% names(ind_param)]

  parameterinput <- as.double(unlist(ind_param[1, ]))
  names(parameterinput) <- names(ind_param)
  parameter <- c(parameter, parameterinput)


  # initial_cmt_values
  states <- initial_cmt_values

  # see if some states need to be evaluated
  line_to_eval <-  which(is.na(as.double(states)))
  # evaluate the character
  if(length(line_to_eval) > 0){

    for(a in line_to_eval){

      states[a] <- with(data = as.tibble(parameter) %>%
                          mutate(rowname = names(parameter)) %>%
                          spread(key = "rowname", value = "value"), expr = eval(parse_expr(as.character(states[a])) ))

    }
  }

  states <- as.double(states)
  names(states) <- names(initial_cmt_values)


  res <- as_tibble(model_RxODE$solve(parameter, events, states))



if(returnSim == T ) return(res)

eval(criteria)

}

#
# TV <- res %>% filter(time == 25) %>% pull(tumVol)
# targets <- tribble(~Dose, ~min, ~max, 0, 141, 208, 50, 22,
#                    216, 100, 3, 167)
# mint <- targets[targets$Dose == dose, "min"]
# maxt <- targets[targets$Dose == dose, "min"]
# case_when(TV > maxt ~ "above", TV < mint ~ "below",
#           T ~ "good")
