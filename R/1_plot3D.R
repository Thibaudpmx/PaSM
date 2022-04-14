## ---------------------------
## Script name: plot_3D
##
## Purpose of script: plot in 3 dimension
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
## Notes:  made for demonstrating the algorithm
## Not usefull in real use.
##
## ---------------------------

# x = expr(k2)
# y = expr(lambda0)
# z = expr(Vd)


VP_proj_creator$set("public", "plot_3D", function(x, y, z,  toaddneg = NULL,add_point =F ){

  x <- enexpr(x)
  y <- enexpr(y)
  z <- enexpr(z)

  neg_above <- rlang::invoke(self$filters_neg_above, .fn = bind_rows) #%>% filter(!!filtre)
  neg_below <- rlang::invoke(self$filters_neg_below, .fn = bind_rows) #%>% filter(!!filtre)
  pos_above <- rlang::invoke(self$filters_pos_above, .fn = bind_rows) #%>% filter(!!filtre)
  pos_below <- rlang::invoke(self$filters_pos_below, .fn = bind_rows) #%>% filter(!!filtre)


  allpoints <- bind_rows(self$poolVP %>% mutate(test = "Accepted"),
                         neg_above %>% mutate(test = "Rejected_Above"),
                         neg_below %>% mutate(test = "Rejected_Below")

  )


  pltly <-  plot_ly()%>%
    add_markers(type = "scatter3d",
                mode = "markers",
                data = allpoints,
                x = ~lambda0,
                y = ~k2,
                z = ~Vd,
                color = ~test,
                opacity = 1,
                colors = c('darkgreen', "orange", 'red'))





  # Add above

  # neg_above %>%
  #   slice(1)

  neg_above %>%
    mutate(a = pmap(list(k2, lambda0, Vd), function(k2, lambda0, Vd){

      filtrecube <- crossing(k2 = c(0, k2),
                             lambda0 = c(lambda0, 0.18),
                             Vd  = c(Vd, 50)) %>%
        slice(1,2,6,5,3,4,8,7)

      pltly <<- pltly %>%
        add_trace(type = 'mesh3d',
                  data = filtrecube,
                  x = ~lambda0,
                  y = ~k2,
                  z = ~Vd,
                  i = c(7, 0, 0, 0, 4, 4, 6, 1, 4, 0, 3, 6),
                  j = c(3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3),
                  k = c(0, 7, 2, 3, 6, 7, 1, 6, 5, 5, 7, 2),
                  facecolor = rep("orange", 12),
                  opacity = 0.4
        )

      "r"

    }))

  # Add below


  neg_below %>%
    mutate(a = pmap(list(k2, lambda0, Vd), function(k2, lambda0, Vd){

      filtrecube <- crossing(k2 = c(k2, 3),
                             lambda0 = c(0,  lambda0),
                             Vd  = c(0, Vd)) %>%
        slice(1,2,6,5,3,4,8,7)


      pltly <<- pltly %>%
        add_trace(type = 'mesh3d',
                  data = filtrecube,
                  x = ~lambda0,
                  y = ~k2,
                  z = ~Vd,
                  i = c(7, 0, 0, 0, 4, 4, 6, 1, 4, 0, 3, 6),
                  j = c(3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3),
                  k = c(0, 7, 2, 3, 6, 7, 1, 6, 5, 5, 7, 2),
                  facecolor = rep("red", 12),
                  opacity = 0.4
        )
      "r"

    }))



  print(pltly)
  "Done"

})
