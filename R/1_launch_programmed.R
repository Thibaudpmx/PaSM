## ---------------------------
## Script name: launch_programmed
##
## Purpose of script: Do the main algorithm for each non usable value
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

VP_proj_creator$set("public", "launch_programmed", function(){

  while(nrow( self$action_programmed$fix_df) != 0){

    cat(blue(paste0("\n\n",nrow( self$action_programmed$fix_df), " remaining VP generations\n")))

    line <-  self$action_programmed$fix_df %>%
      slice(1)


    VP_df_temp <- self$action_programmed$VP_df %>%
      mutate(dummyforjoin = 1) %>%
      left_join(line %>% mutate(dummyforjoin = 1), by = "dummyforjoin" )

    self$add_VP(VP_df_temp, fillatend = F, reducefilteratend = F,use_green_filter = F, npersalve = 2000, time_compteur = F, pctActivGreen = 0.75)

    print(self)

    self$action_programmed$fix_df <- self$action_programmed$fix_df %>% slice(-1)



  }

})
