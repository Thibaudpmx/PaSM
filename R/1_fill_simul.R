## ---------------------------
## Script name: fill_simul
##
## Purpose of script: Perform ODE solving for VP extrapolated as accepted
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
##
## ---------------------------

VP_proj_creator$set("public", "fill_simul", function(nsalve = 500){

  self$poolVP %>%
    mutate(type = map_chr(simul, ~ class(.x)[[1]]))

  self$poolVP %>%
    mutate(type = map_chr(simul, ~ class(.x)[[1]])) %>%
    filter(type == "NULL") %>%
    select(-type, -simul) -> todo

  if("id" %in% names(todo)) todo <- todo %>% select(-id)

  ntodo <- nrow(todo)

  pb <- progress_bar$new(
    format = "  VP creation [:bar] :current/:total (:percent) in :elapsed",
    total = ntodo, clear = FALSE, width= 60)



  while(nrow(todo) > 0){

    pb$update(ratio = (ntodo-nrow(todo))/ntodo)

    line <- todo %>%
      slice(1:(min(nsalve, nrow(todo)))) %>%
      rowid_to_column("id")


    protocol <-  line %>%
      mutate(protocol2 = map(protocol, ~ self$protocols[[.x]])) %>%
      select(id, protocol2) %>%
      unnest() %>%
      mutate(evid = 1) %>%
      bind_rows( crossing(id = unique(line$id), time = self$times, amt = 0, evid = 0,  cmt = self$protocols[[1]]$cmt[[1]] ))%>%
      arrange(id, time)

    # add_events_line$amt[is.na(add_events_line$amt )] <- 0

    # And now we can make the simulation and extract the result !

    b <- Sys.time()
    res <- simulations2(ind_param = line, add_events = protocol, returnSim = T,
                        icmt = self$initial_cmt_values, time_vec =self$times,
                        pardf = self$parameters_default_values, model = self$model) %>%
      as.data.frame()

    res %>%

      left_join(line %>% distinct(id,rowid), by = "id") %>%
      select(-id) %>%
      group_by(rowid) %>%
      nest() %>%
      rename(simul = data)-> newsimuls



    bind_rows(

      self$poolVP %>% filter(! rowid %in% newsimuls$rowid),

      self$poolVP %>%
        filter(rowid %in% newsimuls$rowid) %>%
        select(-simul) %>%
        left_join(newsimuls, by = "rowid")

    ) %>%
      arrange(rowid) -> output

    self$poolVP <- output

    todo <-  todo %>%
      slice(-(1:(min(nsalve, nrow(todo)))))

  }

})

