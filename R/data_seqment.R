# data <- read.table("D:/Peccary_Annexe/Exemple_demo/DATA/Simeoni.txt", header = T, sep = ";", na.strings = ".") %>%
#   as_tibble
# YTYPEcol <-  expr(YTYPE )
# filter <- expr(EVID == 0)
# groups <- exprs(Dose, YTYPE)
# groups <- exprs(protocols, cmt)
#' @export
#'
#'
data_segment <- function(data , ..., ntime = 3, timeforce = NULL, filter = NULL, plot = T){

  filter <- enexpr(filter)
  groups <- enexprs(...)

  if(!is.null(filter)) data <- data %>% filter(!!filter)

data %>%
  group_by(!!!groups) %>%
  nest() %>%
  # slice(1) %>%
  # pull(data) -> x; x <- x[[1]]
  mutate(sampling = map(data, function(x){

    # x %>%
    #   group_by(time) %>%
    #   tally %>%
    #   arrange(desc(n))

    if(!is.null(timeforce)){

      times <- timeforce

    }else{
   x %>% distinct(time) %>%
      arrange(time) %>%
      pull(time) -> times

    times <- map_dbl(1:ntime, ~ times[(.x/ntime * length(times)) %>% round])

    }

    x %>%
      group_by(ID) %>%
      nest() %>%
      crossing(time = times) %>%
      # pull(data) -> y; y <- y[[3]]
      mutate(value = map2_dbl(data, time, function(y,z){

       test <-  y %>%
          filter(time == z)

       if(nrow(test) == 1 ) return(test$OBS)
       #
       timebef <- y$time[y$time < z]; timebef <- timebef[length(timebef)]

       timeaf <-  y$time[y$time > z]
      if(length(timeaf) == 0) return(NA)

        y$OBS[y$time %in% c(timebef, timeaf[[1]])] %>% median()
# return(2)
      })) %>%
      filter(!is.na(value)) %>%
      # filter(time == 50)
      group_by(time) %>%
      summarise(min = min(value), max = max(value))

  })) %>%
  select(-data) %>%
  unnest(sampling) -> targets

if(plot == T) print(data_segment_plot(data, targets))

targets


}

#' @export
#'
#'
data_segment_plot <- function(data, targets,  filter = NULL){


  data %>%
    filter(cmt %in% targets$cmt) %>%
    ggplot()+
    geom_line(aes(time, OBS, group = ID)) +
    facet_wrap(cmt~protocol, scales = "free")+
    geom_segment(data = targets,
                 aes(x = time, xend = time, y = min, yend = max), col ="red", size = 2)+
    scale_y_log10()

}
