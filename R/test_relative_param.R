#' Title
#'
#'
#' @return
#' @export
#'
#' @examples
#'
find_relative <- function(..., model = model_RxODE, time_simul = times, states = initial_cmt_values,   values = domain, protocol = "dose50", sensitivity = 1E-6,
                          printPlot = T, deepAnalysis = F, printCode = T){

  # output  = exprs(tumVol, Conc)
  # output = exprs(Pore)
  output <- enexprs(...)
  model_RxODE <- model


  if(deepAnalysis == T){

    print("here")

  base <-  find_relative( !!!output, model = model, time_simul = time_simul, states = states,   values = values, protocol = protocol, sensitivity = sensitivity,
                  printPlot = F, deepAnalysis = F, printCode = F ) %>%
    mutate(which = "base")

  for(a in values$param){
    # min
    valuestemp <- values
    valuestemp$ref[valuestemp$param == a] <-  valuestemp$min[valuestemp$param == a]

   temp <- find_relative( !!!output, model = model, time_simul = time_simul, states = states,   values = valuestemp, protocol = protocol, sensitivity = sensitivity,
                   printPlot = F, deepAnalysis = F , printCode = F) %>%
      mutate(which = paste0(a, "_min"))

   base <- bind_rows(base, temp)

   # max
   valuestemp <- values
   valuestemp$ref[valuestemp$param == a] <-  valuestemp$max[valuestemp$param == a]

   temp <- find_relative( !!!output, model = model, time_simul = time_simul, states = states,   values = valuestemp, protocol = protocol, sensitivity = sensitivity,
                          printPlot = F, deepAnalysis = F, printCode = F ) %>%
     mutate(which = paste0(a, "_max"))

   base <- bind_rows(base, temp)
  }

# base %>%
#   spread(which, Pore) %>%
#   View

base %>%
  group_by(param, !!!output) %>%
  slice(1) %>%
  gather("ytype", "res", !!!output)-> analyse



  }else{

  # get the name of parameters
  param <- values$param



  # get the reference profile (set all parameter value as their "ref" columns)
  ref <- values %>%
    distinct(param, ref) %>%
    spread(param, ref) %>%
    mutate(what = "ref", param = "ref")


  # Create all the id
  # for each param, create two rows equal to the ref
  # then replace the param in the map by respectively its min and max value
  # at the end add the ref raw
  map(param, function(x){

    min <- max <- ref
    min[[x]] <- values$min[values$param == x]
    max[[x]] <- values$max[values$param == x]

    bind_rows(

      min %>% mutate(what = "min"),
      max %>% mutate(what = "max")
    ) %>%
      mutate(param = x)

  }) %>%
    bind_rows() %>%
    bind_rows(ref) %>%
    rowid_to_column("id") -> ind_params

  # add the default values if not added

  for(a in names(parameters_default_values)[! names(parameters_default_values) %in% names(ind_params)]){

    ind_params[[a]] <- parameters_default_values[[a]]

  }

  # Select the protocol to use
  protocol2 <- protocols[[protocol]]

  # Create the protocol dataset
  # first generate the observation raws
  # combine with   the administraiton one
  # then copy this df for each id
  protocol2 <-   protocol2 %>%
    select(-time) %>%
    crossing(time = time_simul) %>%
    mutate(amt = 0, evid = 0) %>%
    bind_rows(protocol2 %>% mutate(evid = 1)) %>%
    crossing(id = unique(ind_params$id)) %>%
    arrange(id, time, desc(amt))


  # Perform the simulation, transform into tibble and join the parameter for use after
  simulations <- model_RxODE$solve(ind_params %>% select(-what, - param),protocol2 , states) %>%
    as_tibble %>%
    left_join(ind_params)


  # Time to analyse the influence of each parameter
  # first select th esimulation of reference
  refsim <- simulations %>%
    filter(what == "ref")

  # for each parameter
  tibble(param) %>%
    mutate(res = map(param, function(x){

      # select the two id with low and high current param value
      simulations %>%
        filter(param == x) %>%
        select(time, !!!output, x) %>% # select the output/ytype desired
        gather("ytype", "value", !!!output) %>% # one row per ytype observation
        arrange(x) %>%
        spread( x, value) -> temp # but one column for min and max (easier for comparison)

      names(temp) <- c("time","ytype", "min", "max") # modify the name of the column for easier manipulation

      temp %>%
        # join the reference for each time and ytype
        left_join(refsim %>% distinct(time, !!!output) %>% gather("ytype", "int", !!!output),by = c("time", "ytype")) %>%
        # compare 1) the intermerdiate value minus the lowest, and 2) the highest versus the intermediate
        mutate(test1 = int-min, test2 = max-int) %>%
        # if the difference are smaller than the sensitivity, replace by 0
        mutate(test1 = if_else(abs(test1) < sensitivity, 0, test1),
               test2 = if_else(abs(test2) < sensitivity, 0, test2)) %>%
        # gather to have each test in a separate column
        gather("test", "value", test1, test2) %>%
        # compute the result of the test (increased, decreased, or unchanged)
        mutate(res= case_when( value == 0  ~ "same",
                               value <= 0  ~ "dec",
                               value >= 0 ~ "inc")) %>%
        distinct(ytype, res) -> temp

      #  For each ytype we can have several res at different time
      #  We need to extract only one value per ytype and param
      #  so let's do a loop for each ytype
      for (a in unique(temp$ytype)){

        # see the results of the tests
        line <- temp %>% filter(ytype == a) %>% pull(res)

        # if "same" + another one, remove the "same" line
        if(length(line) > 1 & "same" %in% line){
          temp <- temp %>%
            filter(!(ytype == a & res == "same"))
        }

        # if both decreased and increased value, then put to "alt" ("alternating")
        if("dec" %in% line & "inc" %in% line){

          temp <- temp %>%
            filter(!ytype %in% a ) %>%
            bind_rows(tibble(ytype = a, res = "alt") )

        }
      }

      # finally return the temp
      temp


    })) %>%
    unnest(res) -> analyse

  }




  analyse_temp <- analyse

  # in case of deepAnalysis, need to remove the inconsistent one
  if(deepAnalysis == T){



    # rem NON-USABLE
    analyse_temp %>%
      filter(res == "NON-USABLE") %>% pull(param) -> torem

    analyse_temp <- analyse_temp %>%
      filter(! param %in% torem)

    # once dec, once inc
    analyse_temp %>%
      filter(res %in% c("dec", "inc")) %>%
      group_by(param) %>%
      tally %>%
      filter(n == 2) %>%
      pull(param) -> torem

    analyse_temp <- analyse_temp %>%
      filter(! param %in% torem)

    # dec / in + ind

    analyse_temp %>%
      # filter(res %in% c("dec", "inc")) %>%
      group_by(param) %>%
      tally %>%
      filter(n == 2) %>%
      pull(param) -> torem_ind

    analyse_temp <- analyse_temp %>%
      filter( !( param %in% torem_ind & res == "ind"))

  }



  # Here the goal is to produce the code for the user
  # Because we need to produce three code (for reducer, increaser and non-modifier)
  # Let's create a generic function
  # here "sens" represent the possible variation analysed above ("inc", "dec" or "alt)
  # "name" is the name of the vector we want to create


  func_temp <- function(sens, name){





    # for each ytype
    temp <- map_chr(unique(analyse_temp$ytype),function(x){

      # look at the parameters that increase/decrease/unchange (use of sens!) that YTYPE
      para <- analyse_temp$param[analyse_temp$ytype == x & analyse_temp$res == sens]


      # if there is none, simply output "character"
      if(length(para) == 0) return(  paste0(x, " =character()" ) )

      # else paste all the paramer to have a code like " "param1", "param2",..."
      temp <- paste0("\"", para, "\"") %>% paste0(collapse = ", ")

      # then create the final vector code for that YTYPE
      paste0(x, " = c(",temp, ")" )

    }) %>%
      paste0(collapse = ", ") # separate all YTYPE vector by a comma

    paste0(name, " <- list(", temp, ")\n\n") # finally generate the final vector

  }

 if( printCode ){

  cat(red("Code generator - use carefully"))

  paste0("\n\n", func_temp("inc", "param_increase" ),
         func_temp("dec", "param_reduce" ),
         func_temp("same", "param_no_impact" ), "\n"
  ) %>%
    cat

  cat(red("End of code generation"))
 }

  # if(deepAnalysis == T) return()
# Then Print the plot to help the users verify/understand the dynamic
# One plot per ytype/output in the end gather with plot_grid
  if(printPlot == T & deepAnalysis == F){
  print(  map(output, ~   simulations %>%
                filter(what != "ref") %>%
                # filter(tumVol > 1e-5) %>%
                ggplot()+
                geom_line(aes(time, !!.x, col = what))+
                geom_line(data = refsim %>% select(-param), aes(time, !!.x, col = "ref"))+
                facet_wrap(~param, scales = "free")+
                scale_y_log10()+
                theme_bw()+
                ggtitle(deparse(.x))
  ) %>%
    invoke(.fn = plot_grid))
  }
# As a output, a table summarising the influence of each paramter for each ytype

  analyse %>%
    mutate(res = if_else(res == "same", "ind", res)) %>%
    mutate(res = if_else(res == "alt", "NON-USABLE", res)) %>%
    spread(ytype, res)
}


# find_relative(Pore, protocol = "unique", model = model_RxODE, values = domain, sensitivity = 0.1, deepAnalysis = T)


# helper to create --------------------------------------------------------

#' Title
#'
#'
#' @return
#' @export
#'
#' @examples
#'
tribblecreator <- function(model_RxODE, rm_default_param = T){

  # take all the paramters of the model as detected by RxODE

  params <- model_RxODE$params


  # if the user wan's to remove the default values as put into the parameters_default_values vector
  if(rm_default_param == T) params <- params[!params %in% names(parameters_default_values)]


  # then we need to remove what RxODE detected as parameter even if
  # the value is actually computed inside the code


  str_split(deparse(model_RxODE$model), ";")[[1]] %>%
    gsub(pattern = "\\\\n", replacement = "") %>%
    gsub(pattern = "structure\\(\"", replacement = "") %>%
    gsub(pattern = "=.+", replacement =  "") ->torem

  params <- params[! params %in% torem]

  # Finally print the code to help the users

  paste0("domain <- tribble(~param, ~min, ~ref,  ~max,\n",
         paste0("\"", params, "\",") %>% paste0(collapse = " ,\n"), " )\n\n",
         "find_relative(your_output_without_quotes, protocol = \"protocol_to_use\", model = model_RxODE, values = domain)") %>%
    cat

}

