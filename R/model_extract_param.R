

#' model_RxODE
#'
#' @param model
#'
#' @return
#' @export
#'
#' @examples
model_extract <- function(model = model_RxODE){

  param <- model$params
  equa <- model$model
  states <- model$state

  param <- param[!map_lgl(param, ~ grepl(paste0(.x, " *="), equa))]

  output <- paste0("parameters_default_values <- c(\n\n", paste0(param, " = 0") %>% paste0(collapse = ",\n"),"\n\n)")

  state <-   states[paste0(states,"0") %in% param]

  output2 <-  paste0("initial_cmt_values <- c(\n\n", paste0(state, " = ", paste0("\"",state,"0\"") ) %>% paste0(collapse = ",\n"), "\n\n)")

  paste0(output, "\n\n# Now the initial values, with can be parameter quoted\n\n", output2) %>% cat

}
