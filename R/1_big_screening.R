## ---------------------------
## Script name: big_screening
##
## Purpose of script: Second main algorithm
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

VP_proj_creator$set("public", "big_screening", function(domain){


  if(!"param" %in% names(domain) & nrow(domain) == 2){

    domain <-   imap(domain, function(x,y){

      values <- domain[[y]]
      values <- values[order(values)]
      tibble(param = y, from = values[[1]], to = values[[2]]) %>%
        mutate(step = if_else(from == to, 0,1))

    }) %>%
      rlang::invoke(.fn = bind_rows)
  }
  nsize = 10^18

  param_fluct <- domain %>%
    filter(step != 0)

  # How many param for having 20000
  nperparam <- 200000^(1/nrow(param_fluct)) %>% ceiling()



  param_fluct %>%
    mutate(values = pmap(list(from, to , step), function(from, to , step){

      seqq <- seq(from, to, step)
      indic <- length(seqq) * 0:(nperparam-1)/(nperparam-1)
      indic[indic == 0 ] <- 1
      round(seqq[indic])

    })) %>%
    unnest() -> allparamcut

  allparamcut %>%
    group_split(param) %>%
    map(function(x){

      list()

      x %>%
        select(values) -> temp
      names(temp) <- unique(x$param)
      temp
    }) %>%
    rlang::invoke(.fn = crossing) -> toadd

  # fixed one
  domain %>%
    filter(step == 0) -> fixed

  if(nrow(fixed)>0){
    for(a in 1:nrow(fixed)){

      toadd[[fixed$param[[a]]]] <- fixed$from[[a]]
    }
  }

  nsize/nrow(toadd)

  VP_temp <- self#VP_proj_creator$new(sourcefile = "D:/these/Second_project/QSP/modeling_work/VT_simeoni/1_user_inputs/1_config_Lindner.r")
  VP_temp$filters_neg_above <- tibble()
  VP_temp$filters_neg_below <- tibble()

  VP_temp$add_VP(toadd)

  VP_temp$n_filter_reduc()
  if(nrow(VP_temp$poolVP) >0) return(VP_temp$poolVP)
  # VP_temp$plot_VP()
  group_split(VP_temp$filters_neg_above,rowid)[[1]] -> x
  nrowss <- nrow(VP_temp$filters_neg_above)
  whereVP <- group_split(VP_temp$filters_neg_above,rowid)[1:min(20,nrowss)] %>%
    map(function(x){

      y <- x
      for(a in unique(allparamcut$param)){

        temp <- allparamcut %>% filter(param == a)
        currentloc <- which(temp$values == y[[a]])

        if(currentloc == 1){
          value_if_increase <- 0
          value_if_decrease <- temp$values[[currentloc+1]]
        }else if (currentloc == nrow(temp)){

          value_if_increase <- temp$values[[currentloc-1]]
          value_if_decrease <-  max(temp$values)

        }else{
          value_if_increase <- temp$values[[currentloc-1]]
          value_if_decrease <- temp$values[[currentloc+1]]
        }


        y[[a]] <- case_when(a %in% self$param_increase$Pore ~ value_if_increase ,
                            a %in% self$param_reduce$Pore ~ value_if_decrease)

      }

      bind_rows(x, y)[, domain$param]

    })

  return(whereVP)

})
