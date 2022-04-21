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


  param_unique <- map(self$param[!self$param %in% c("w0", "k1")], function(x){



    temp <-  c(self$poolVP[[x]],0, Inf )  %>% unique




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
        slice(-1) #%>%
        # bind_rows(tibble(max =  temp, min = temp)) # if we accept lower dimensions (to heavy in practical)
    }

    df_temp <- df_temp %>%
      select(min, max)

    names(df_temp) <- c(mincol, maxcol)

    df_temp
  })


  allsquares <- rlang::invoke(.fn = crossing, .args = param_unique ) %>%
    mutate(w0min = 50 , w0max = 50, k1min = 0.5, k1max = 0.5)#%>%


  for(x in names(self$param_increase)[[1]]){

    allsquares[[paste0(x,"_above")]] <- na_lgl
    allsquares[[paste0(x,"_below")]] <- na_lgl

    filtres <- filters <- self$make_filters(x) %>%
      map_chr(~ gsub("line\\$", "", .x))

    # filtres

    above <- filter_reduc(df = self$poolVP %>% select(-simul), filtre = filtres[1],
                          param_increase = self$param_increase[[x]], param_reduce = self$param_reduce[[x]], direction = "above")
    below <- filter_reduc(df = self$poolVP %>% select(-simul), filtre = filtres[2],
                          param_increase = self$param_increase[[x]], param_reduce = self$param_reduce[[x]], direction = "below")


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

    allsquares$above <- na_lgl

    # allsquares %>% filter(above)
    if(nrow(above) > 0){
      for(a in 1:nrow(above)){

        print(a)
        ref <- above %>% slice(a)


        allsquares <-
          allsquares %>%
          mutate(above = if_else(!!parse_expr(filter_above), T, above))

      }
    }


    filter_below<- filters[2]

    allsquares$below <- na_lgl
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
#
self$poolVP %>%
  filter(rowid %in% c(12191, 6856))
#
#
#

allsquares %>%
  filter(k2min <= 0.25, k2max >=0.25, lambda0min <= 1.05, lambda0max >=1.05, kemin <= 0.3, kemax >=0.3, Vdmin <= 27, Vdmax >=27,
         lambda1min <= 6, lambda1max >=6 & k2max == 0.25 & lambda0min == 1.05 &  kemin == 0.3)


allsquares %>%
  filter(k2min == 0.25, k2max == 0.25, lambda0min <= 1.05, lambda0max >=1.05, kemin <= 0.3, kemax >=0.3, Vdmin <= 27, Vdmax >=27,
         lambda1min <= 6, lambda1max >=6)


allsquares %>%
  filter(k2min == 0.25, k2max == 0.25)
allsquares

allsquares %>%
  filter(k2min == 0.25)
  filter(k2min == 0.2, k2max == 0.4, lambda0min == 0.28, kemin == 0.9, Vdmin == 35, lambda1min == 4.8)
#
# self$poolVP %>%
#   filter(rowid %in% c(31741))
#
# # below %>%
#   filter(k2 <= 0.2, lambda0 >= 1.12, ke >= 1.1, Vd >=37,lambda1>=4.8 )
#
#
# self$poolVP %>% arrange(lambda1)
