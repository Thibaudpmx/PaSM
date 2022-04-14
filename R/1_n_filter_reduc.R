## ---------------------------
## Script name: n_filter_reduc
##
## Purpose of script: Reduce the number of filter after add_VP
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



VP_proj_creator$set("public", "n_filter_reduc", function(){




  filters <- self$make_filters() %>%
    gsub(pattern = "line\\$", replacement = "" )



  # Handle negative above

  filter_reduc(self$filters_neg_above, obj = self, direction = "above")

  self$filters_neg_above %>%
    select(-starts_with("rowid")) %>%
    # arrange(desc(k2), k2, lambda0, lambda1, Vd) %>%
    rowid_to_column() -> temp

  message("Reducing negative above filters")

  self$filters_neg_above <-  filter_reduc(self$filters_neg_above, obj = self, direction = "above")%>%
    mutate(PrimFilter = T)

  # if(nrow(temp) > 0 ){
  # for(a in 1:nrow(temp)){
  #   # print(a)
  #
  #   if(a %in% temp$rowid){
  #     ref <- temp %>%
  #       filter(rowid == a)
  #
  #     temp %>%
  #       mutate(test = !!parse_expr(filters[[1]]))  -> temp
  #
  #     temp$test[temp$rowid == a] <- F
  #     temp <- temp %>%
  #       filter(test == F)
  #   }
  #   self$filters_neg_above <- temp %>% mutate(PrimFilter = T) %>% select(-test)
  # }
  # }




  # Handle negative below

  # self$filters_neg_below %>%
  #   select(-starts_with("rowid")) %>%
  #   # arrange(desc(k2), k2, lambda0, lambda1, Vd) %>%
  #   rowid_to_column() -> temp

  message("Reducing negative below filters")

  self$filters_neg_below <-  filter_reduc(self$filters_neg_below, obj = self, direction = "below") %>%
    mutate(PrimFilter = T)

  message("Done")
  #   if(nrow(temp) > 0 ){
  #   for(a in 1:nrow(temp)){
  #     # print(a)
  #
  #     if(a %in% temp$rowid){
  #       ref <- temp %>%
  #         filter(rowid == a)
  #
  #       temp %>%
  #         mutate(test = !!parse_expr(filters[[2]]))  -> temp
  #
  #       temp$test[temp$rowid == a] <- F
  #       temp <- temp %>%
  #         filter(test == F)
  #     }
  #
  #     self$filters_neg_below <- temp %>% mutate(PrimFilter = T) %>% select(-test)
  #   }
  # }
  #



  #
  #   # for each cmtt
  # for(cmtt in unique((self$targets$cmt))){
  #
  # # filter pos above <low
  #
  #
  #
  # self$filters_pos_above[[cmtt]] %>%
  #   select(-starts_with("rowid")) %>%
  #   # arrange(desc(k2), k2, lambda0, lambda1, Vd) %>%
  #   rowid_to_column() -> temp
  #
  #
  # if(nrow(temp) > 0){
  #
  #   message("Reducing positive above filters")
  # for(a in 1:nrow(temp)){
  #   # print(a)
  #
  #   if(a %in% temp$rowid){
  #     ref <- temp %>%
  #       filter(rowid == a)
  #
  #     temp %>%
  #       mutate(test = !!parse_expr(filters[[1]]))  -> temp
  #
  #     temp$test[temp$rowid == a] <- F
  #     temp <- temp %>%
  #       filter(test == F)
  #   }
  # }
  #
  #   self$filters_pos_above[[cmtt]] <- temp %>% mutate(PrimFilter = T) %>% select(-test)
  #
  # }
  #
  #
  # # handle pos belw
  #
  # self$filters_pos_below[[cmtt]] %>%
  #   select(-starts_with("rowid")) %>%
  #   # arrange(desc(k2), k2, lambda0, lambda1, Vd) %>%
  #   rowid_to_column() -> temp
  #
  # message("Reducing positive above filters")
  #
  # if(nrow(temp) > 0){
  #   message("Reducing positive above filters")
  #
  # for(a in 1:nrow(temp)){
  #   # print(a)
  #
  #   if(a %in% temp$rowid){
  #     ref <- temp %>%
  #       filter(rowid == a)
  #
  #     temp %>%
  #       mutate(test = !!parse_expr(filters[[2]]))  -> temp
  #
  #     temp$test[temp$rowid == a] <- F
  #     temp <- temp %>%
  #       filter(test == F)
  #   }
  # }
  #   self$filters_pos_below[[cmtt]] <- temp %>% mutate(PrimFilter = T) %>% select(-test)
  # }
  #
  #
  #
  # }
  #

})
