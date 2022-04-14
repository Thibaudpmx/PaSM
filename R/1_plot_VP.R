## ---------------------------
## Script name: plot_VP
##
## Purpose of script: plot VP that have been accepted and already solved
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

VP_proj_creator$set("public", "plot_VP", function(nmax = Inf , rowids = NULL){




  self$poolVP %>%
    mutate(type = map_chr(simul, ~ typeof(.x)[[1]])) %>%
    filter(type != "NULL")

  if(is.null(rowids)){

    temp <-  self$poolVP %>%
      mutate(type = map_chr(simul, ~ typeof(.x)[[1]])) %>%
      filter(type != "NULL")

    temp <- temp %>%
      sample_n( min(nrow(temp), nmax ))


  }else{

    temp <-  self$poolVP %>%
      filter(rowid %in%rowids)


  }





  # slice() %>%
  temp %>%
    unnest(simul) %>%
    gather("cmt", "value", !!!parse_exprs(unique(self$targets$cmt))) %>%
    filter( ! (cmt == "Conc" & value == 0)) %>%
    # filter( ! (cmt == "Conc" & time > 15)) %>%
    ggplot()+
    geom_line(aes(time, value, group = id))+
    facet_wrap(protocol~cmt, scales = "free")+
    scale_y_log10()+
    geom_segment(data = self$targets ,
                 aes(x = time, xend = time, y = min, yend = max), col ="red", size = 2)


})
