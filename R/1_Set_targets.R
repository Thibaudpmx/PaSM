## ---------------------------
## Script name: set_targets
##
## Purpose of script: Attribute targets to the main R6 object, either by
## analyzing a dataset or by manual attribution
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



VP_proj_creator$set("public", "set_targets", function(..., filter = NULL, ntime = 3, manual = NULL, timeforce = NULL){

  filter <- enexpr(filter)

  if(!is.null(self$self$poolVP))  return("Virtual Patients have already been generated - It is thus no possible to modify the targets.")
  if(!is.null(self$filters_neg_below)) if(map_dbl(self$filters_neg_below, nrow) %>% sum > 0)  return("Filters have already been generated - It is thus no possible to modify the targets.")
  if(!is.null(self$filters_neg_above)) if(map_dbl(self$filters_neg_above, nrow) %>% sum > 0)  return("Filters have already been generated - It is thus no possible to modify the targets.")
  if(!is.null(self$filters_pos_above)) if(map_dbl(self$filters_neg_below, nrow) %>% sum > 0)  return("Filters have already been generated - It is thus no possible to modify the targets.")
  if(!is.null(self$filters_pos_below)) if(map_dbl(self$filters_pos_below, nrow) %>% sum > 0)  return("Filters have already been generated - It is thus no possible to modify the targets.")



  if(!is.null(manual)){

    targets <- manual

    print(

      ggplot()+
        facet_wrap(cmt~protocol, scales = "free")+
        geom_segment(data = targets,
                     aes(x = time, xend = time, y = min, yend = max), col ="red", size = 2)+
        scale_y_log10()+
        geom_point(data = targets %>%
                     gather("key", "value", min, max),
                   aes(x = time, y = value), col = "red")+
        theme_bw()
    )



  }else{

    targets <- data_segment(data = self$data, protocol, cmt, filter = !!filter, ntime = ntime, timeforce = timeforce) %>%
      filter(max != "0")

  }


  print( as.data.frame(targets))


  # update targets
  self$targets <- targets %>% ungroup()

  # initialise filters

  filters_list= list()

  for(a in unique(targets$cmt))  filters_list[[a]] <- tibble()

  self$filters_pos_above  <- filters_list
  self$filters_pos_below  <- filters_list



  return(  message("Once any VP's or filters have been created, targets will be locked."))

})
