## ---------------------------
## Script name: plot_2D
##
## Purpose of script: plot in 2 dimension
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

VP_proj_creator$set("public", "plot_2D", function(x, y , cmt_green = "tumVol", toaddneg = NULL, plotMain = F, add_point =F , IDref = NULL, plotoreturn = 3){

  x <- enexpr(x)
  y <-  enexpr(y)



  # Apply filter to keep only on plan

  # if( is.null(IDref)) IDref <- sample(  self$poolVP$id, size = 1)
  #
  # line <- self$poolVP %>%
  #   filter(id == IDref)
  #
  # namesparam <- names(line)
  # namesparam <- namesparam[! namesparam %in% c("rowid", "id", "protocol", "simul", deparse(x), deparse(y))]
  # namesparam <- namesparam[!grepl("(_BU$)|(_AL$)", namesparam)]
  #
  # filtre <- line[, namesparam] %>%
  #   unlist() %>%
  #   imap_chr(~ paste0(.y, " == " ,.x)) %>%
  #   paste0(collapse = " & ") %>%
  #   parse_expr
  # end apply filter

  ytype <- names(self$filters_pos_above)[[1]]#Conc

  neg_above <- rlang::invoke(self$filters_neg_above, .fn = bind_rows) #%>% filter(!!filtre)
  neg_below <- rlang::invoke(self$filters_neg_below, .fn = bind_rows) #%>% filter(!!filtre)
  pos_above <- rlang::invoke(self$filters_pos_above, .fn = bind_rows) #%>% filter(!!filtre)
  pos_below <- rlang::invoke(self$filters_pos_below, .fn = bind_rows) #%>% filter(!!filtre)




  plot_dot <-

    ggplot()+
    geom_point(data = pos_above, aes(!!x, !!y), col = "darkgreen") +
    geom_point(data = pos_below, aes(!!x, !!y), col = "darkgreen") +
    geom_point(data = neg_below, aes(!!x, !!y, shape = cmt), col = "chocolate") +
    geom_point(data = neg_above , aes(!!x, !!y, shape = cmt), col = "red", alpha = 1)+
    # geom_point(data = self$filters_neg_below, aes(k2, lambda0), col = "red", alpha = 1)+
    theme_bw()+
    ggtitle( "VP tested")+
    theme(plot.title = element_text(hjust = 0.5)); plot_dot

  if(length(unique(self$targets$cmt) == 1))  plot_dot <- plot_dot + guides(shape = NULL)

  if(!is.null(toaddneg)){

    toaddneg %>%
      left_join(self$poolVP) %>%
      filter(is.na(rowid)) -> addneg

    plot_dot <- plot_dot +
      geom_point(data = addneg, aes(!!x, !!y), col = "red")+
      geom_point(data = self$poolVP, aes(!!x, !!y), col = "darkgreen")

  }


  if(add_point == T) plot_dot <-   plot_dot+
    geom_point(data = self$poolVP, aes(!!x, !!y), col = "darkgreen")

  rectangles_above <- tibble()
  rectangles_below <- tibble()

  for(a in unique(self$targets$cmt) ){


    xana <- case_when(deparse(x) %in% self$param_increase[[a]] ~ "inc",
                      deparse(x) %in% self$param_reduce[[a]] ~ "dec",
                      deparse(x) %in% self$param_no_impact[[a]] ~ "no_impact",
                      T ~ "None")


    yana <- case_when(deparse(y) %in% self$param_increase[[a]] ~ "inc",
                      deparse(y) %in% self$param_reduce[[a]] ~ "dec",
                      deparse(y) %in% self$param_no_impact[[a]] ~ "no_impact",
                      T ~ "None")


    rectangles_above  <- neg_above %>%
      filter(cmt == a) %>%
      mutate(xmin = case_when(xana %in% c("dec", "no_impact") ~ -Inf, T ~ !!x),
             xmax = case_when(xana %in% c("dec","None") ~ !!x, T ~Inf),
             ymin = case_when(yana %in% c("dec", "no_impact") ~ -Inf, T ~!!y),
             ymax = case_when(yana  %in% c("dec","None") ~ !!y, T ~Inf)) %>%
      bind_rows(rectangles_above)

    rectangles_below  <- neg_below %>%
      filter(cmt == a) %>%
      mutate(xmin = case_when(xana %in% c("inc","no_impact") ~ -Inf, T ~ !!x),
             xmax = case_when(xana %in% c("inc","None") ~ !!x, T ~Inf),
             ymin = case_when(yana %in% c("inc","no_impact") ~ -Inf, T ~!!y),
             ymax = case_when(yana %in% c("inc","None") ~ !!y, T ~Inf))%>%
      bind_rows(rectangles_below)




  }

  rectangles <- bind_rows(rectangles_above %>% mutate(where = "above"), rectangles_below %>% mutate(where = "below"))


  plot1 <-  plot_dot +
    geom_rect(data = rectangles, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,  fill = where, col = where), alpha = 0.3)+
    ggtitle( "zone rejection")+
    scale_color_manual(values = c("red", "chocolate"))+
    scale_fill_manual(values = c("red", "chocolate"));plot1

  if(plotoreturn == 1) return(plot1)
  # plot1 +
  #   geom_point(data = self$poolVP, aes(!!x, !!y), col = "darkgreen", alpha = 0.2)
  #



  rectangles_above_lower <-  pos_above %>%
    mutate(xmin = case_when(xana == "dec" ~ 0, T ~ !!x),
           xmax = case_when(xana == "dec" ~ !!x, T ~Inf),
           ymin = case_when(yana == "dec" ~ 0, T ~!!y),
           ymax = case_when(yana == "dec" ~ !!y, T ~Inf))

  rectangles_below_upper  <- pos_below%>%
    mutate(xmin = case_when(xana == "inc" ~ 0, T ~ !!x),
           xmax = case_when(xana == "inc" ~ !!x, T ~Inf),
           ymin = case_when(yana == "inc" ~ 0, T ~!!y),
           ymax = case_when(yana == "inc" ~ !!y, T ~Inf))





  plot2 <-   plot_dot +
    geom_rect(data = rectangles_below_upper, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "darkgreen",  col = "darkgreen")+
    ggtitle( "zone below upper limit")

  plot3 <-  plot_dot +
    geom_rect(data = rectangles_above_lower, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "darkgreen",  col = "darkgreen")+
    ggtitle( "zone above lower limit")

  # plot_dot +
  #   geom_rect(data = rectangles_above_lower, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "red",  col = "red")+
  #   geom_rect(data = rectangles_below_upper, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.3,  fill = "blue",  col = "blue")
  #   ggtitle( "zone above lower limit")

  if(plotMain == F) return(plot_grid(plot_dot, plot1, plot2, plot3))


  self$poolVP %>%
    select(!!x, !!y) -> temp


  crossing(a = unique(temp[[deparse(x)]]), b = unique(temp[[deparse(x)]])) %>%
    filter(b >= a) %>%
    mutate(height = map2(a, b, function(a,b){

      temp %>%
        filter(!!x>=a & !!x<=b) %>%
        group_by(!!x) %>%
        summarise(max = max(!!y), min = min(!!y)) -> temp2

      tibble( floor = max(temp2$min), ceiling = min(temp2$max))
      # temp2

    })) %>%
    unnest() %>%
    filter(ceiling >= floor) %>%
    arrange(desc(b)) %>%
    group_by(a, floor, ceiling) %>%
    slice(1) %>%
    rename(xmin = a, xmax = b, ymin = floor, ymax = ceiling)-> col_green
  #
  #
  plot4<-
    ggplot()+
    geom_rect(data = rectangles, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,  fill = where, col = where), alpha = 1)+
    scale_color_manual(values = c("red", "chocolate"))+
    scale_fill_manual(values = c("red", "chocolate"))+
    geom_rect(data = col_green, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 1,  fill = "darkgreen",  col = "darkgreen")+
    theme_bw()+
    ggtitle( "Total")+
    theme(plot.title = element_text(hjust = 0.5))


  if(!is.null(toaddneg)){


    plot4 <- plot4 +
      geom_point(data = addneg, aes(!!x, !!y), col = "red", alpha = 0)+
      geom_point(data = self$poolVP, aes(!!x, !!y), col = "darkgreen", alpha = 0)

  }

  plot_grid(plot4,plot_grid(plot_dot, plot1, plot2, plot3))

})
