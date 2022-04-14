## ---------------------------
## Script name: time_analyserr
##
## Purpose of script: analyse the time for various steps of main algorithm
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

VP_proj_creator$set("public", "time_analyserr", function(){

  timetoana <- self$timesaver

  timetoanaloop <- timetoana$poolVP_compteur %>%
    filter(!is.na(n))%>%
    select(-time, -nelim, -ninfo, -computmodel)

  # green filter

  timetoanaloop %>%
    distinct(n, wholegreenfilter , `TRUE-NA`,  `TRUE-TRUE`, `NA-TRUE`, time_add_above_fil, time_add_below_fil, time_addgreennofil)




  timetoanaloop
  names(timetoanaloop)


  # Impact filter neg

  timetoanaloop %>%
    select(timemodel, nremoved_above_fil, filter_neg_above, remneg_above_fil, nfilter_negab_bef,nfilters_negab_af ) %>%
    mutate(timesaved = nremoved_above_fil * timemodel / 2000)

  timetoanaloop %>%
    select(timemodel, nremoved_below_fil, filter_neg_below, remneg_below_fil, nfilter_negbel_bef,nfilters_negbel_af ) %>%
    mutate(timesaved = nremoved_below_fil * timemodel / 2000)

  timetoanaloop %>%
    gather("key", "value", filter_neg_above, remneg_above_fil)

  sum(timetoanaloop$remneg_above_fil)







})

