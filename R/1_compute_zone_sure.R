## ---------------------------
## Script name: compute_zone_sure
##
## Purpose of script: Compute zones we are sure the VPs are accepted
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

VP_proj_creator$set("public", "compute_zone_sure", function(domain){

  # self$poolVP
  filtres <-   self$make_filters()

  above <- filter_reduc(df = self$poolVP %>% select(-simul), filtre = filtres[1],
                        param_increase = self$param_increase$tumVol, param_reduce = self$param_reduce$tumVol, direction = "above")
  below <- filter_reduc(df = self$poolVP %>% select(-simul), filtre = filtres[2],
                        param_increase = self$param_increase$tumVol, param_reduce = self$param_reduce$tumVol, direction = "below")


  param_unique <- map(self$param, function(x){



    temp <-  c(above[[x]], below[[x]] )  %>% unique



    maxcol <- paste0(x, "max")
    mincol <- paste0(x, "min")


    if(length(temp) == 1){

      df_temp <- tibble(max = temp, min = temp)

    }else{

      df_temp <-tibble(max =  temp) %>%
        # add_row(lambda0max = 0) %>%
        # add_row(lambda0max = 1) %>%
        arrange(max) %>%
        mutate(min = lag(max)) %>%
        slice(-1)
    }

    df_temp <- df_temp %>%
      select(min, max)

    names(df_temp) <- c(mincol, maxcol)

    df_temp
  })


  allsquares <- rlang::invoke(.fn = crossing, .args = param_unique )

  allsquares <-  allsquares %>%
    mutate(above = na_lgl, below = na_lgl)

  for(x in names(self$param_increase)[[1]]){

    filters <- self$make_filters(x) %>%
      map_chr(~ gsub("line\\$", "", .x))

    filters <- paste0("( ", filters, ")")

    other <- self$param[! self$param %in% c(self$param_increase[[x]], self$param_reduce[[x]])]

    for(a in c(self$param_increase[[x]], other)){

      filters[1] <- gsub(paste0(" ",a, " "), paste0(" ", a, "min "), filters[1])
      filters[2] <- gsub(paste0(" ",a, " "), paste0(" ", a, "max "), filters[2])

    }

    for(a in self$param_reduce[[x]]){

      filters[1] <- gsub(paste0(" ",a, " "), paste0(" ", a, "max "), filters[1])
      filters[2] <- gsub(paste0(" ",a, " "), paste0(" ", a, "min "), filters[2])

    }

    # ref <- above %>% slice(a)
    #
    #  allsquares %>%
    #   filter(!!parse_expr(filters[2]))

    #
    #
    # ref <- above %>% slice(a)
    #
    # # allsquares <- allsquares %>%
    #   # filter(!!parse_expr(filter_above))
    # allsquares %>%
    #   filter(k2max <= 0)
    #   filter(kemin >=0.6 & Vdmin >=0 & lambda0min >= 0.03 &  lambda1min >= 10 & k2max <= 0 & k1min ==0.5)

    filter_above <- filters[1]

    # above <- self$filters_neg_above %>% filter(cmt == x)



    if(nrow(above) > 0){
      for(a in 1:nrow(above)){

        ref <- above %>% slice(a)

        allsquares <-
          allsquares %>%
          mutate(above = if_else(!!parse_expr(filter_above), T, above))

      }
    }

    filter_below<- filters[2]
    # filters_neg_below <- self$filters_neg_below %>% filter(cmt == x)
    if(nrow(below) > 0){
      for(a in 1:nrow(below)){

        ref <- below %>% slice(a)

        allsquares <- allsquares %>%
          mutate(below = if_else(!!parse_expr(filter_below), T, below))

      }
    }
  }

  final <- allsquares %>% filter(above & below)
  self$zone_sure <- final
})
