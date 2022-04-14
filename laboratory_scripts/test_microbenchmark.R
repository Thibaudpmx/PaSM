library(microbenchmark)
# Test 1 (one row per protocol) ------------------------------------------------------------------

test <- crossing(a = rnorm(n = 100, mean = 0, 5),
                 b = rnorm(n = 100, mean = 0, 5),
                 c = rnorm(n = 1000, mean = 0, 5))

testproto <- letters[1:3]
# first lets test rowid_to_columns vs mutate(id = 1:nrow) -> rowid_to_columns champion
microbenchmark(test %>% rowid_to_column("id"),
               test %>% mutate(id = 1:nrow(test) ),  times = 30, unit = "ms")

# then, main test
test1 <- function(){
  test %>%
    rowid_to_column("id") %>% # attribute a id for each VP
    crossing(protocol = testproto ) %>% # each VP must have one row per protocol
    rowid_to_column("rowid")
}

test2 <- function(){

  temp <-  test %>%
    rowid_to_column("id")

  map(testproto, function(x){

    temp %>% mutate(protocol = x)

  }) %>%
    bind_rows() %>%
    rowid_to_column("rowid")

}


microbenchmark(test1(), test2(), times = 5) # test2 big winner
# Unit: milliseconds
# expr       min        lq      mean   median       uq      max neval cld
# test1() 2890.4191 2977.4620 3095.5487 3055.608 3104.879 3449.375     5   b
# test2()  658.7454  832.3338  959.5845 1083.546 1088.341 1134.956     5  a




# Test 2th  reducing filter -------------------------------------------------
# the idea is: can we arrange the values for bein in the right order

test <- tibble(A = rnorm(n = 2000, mean = 20,sd = 5), B = rnorm(n = 2000, mean = 20,sd = 5), C = rnorm(n = 2000, mean = 20,sd = 5))


test1 <- function(arrangee){
df2 <- test %>%
  select(-starts_with("iddummy")) %>%
  rowid_to_column("iddummy")

filtre <- "A <= ref$A & B <= ref$B & C >= ref$C"

if(arrangee == T ){
df2 <- df2 %>%
  arrange(desc(A), desc(B), C)
}

for(a in df2$iddummy){
  # print(a)

  if(a %in% df2$iddummy){
    ref <- df2 %>%
      filter(iddummy == a)

    df2 %>%
      mutate(test = !!parse_expr(filtre))  -> df2

    df2$test[df2$iddummy == a] <- F
    df2 <- df2 %>%
      filter(test == F)
  }
}
}


microbenchmark(test1(arrangee = F), test1(arrangee = T), times = 30)
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval cld
# test1() 573.7099 595.2016 617.3127 608.7725 623.6155 834.0578    30   b
# test2() 189.8674 198.3476 206.4515 202.1570 215.5486 235.0098    30  a


# Test 2pratique  reducing filter -------------------------------------------------

test1 <- function(arrangee){
df2 <- df %>%
  select(-starts_with("iddummy")) %>%
  rowid_to_column("iddummy")

filtre <- gsub("line\\$", "", filtre)

if(arrangee == T ){


list_arrange <- list()
list_arrange <- map(obj$param_increase[[a]],~ expr(desc(!!parse_expr(.x))))
list_arrange <- c(list_arrange, map(obj$param_reduce[[a]],~ expr(!!parse_expr(.x))))

df2 <- df2 %>%
  arrange(!!!list_arrange)

}



for(a in df2$iddummy){
  # print(a)

  if(a %in% df2$iddummy){
    ref <- df2 %>%
      filter(iddummy == a)

    df2 %>%
      mutate(test = !!parse_expr(filtre))  -> df2

    df2$test[df2$iddummy == a] <- F
    df2 <- df2 %>%
      filter(test == F)
  }
}
}

# a <- "tumVol"
microbenchmark(test1(arrangee = F), test1(arrangee = T), times = 30)
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval cld
# test1() 551.3989 560.4982 578.0551 566.1296 575.7230 790.6074    30   b
# test2() 145.2174 146.9940 155.2679 150.9399 156.9059 230.7469    30  a

# Largely worth it !!

# Now it has been implemented, let's see the impact:
microbenchmark(
filter_reduc(df = filters_neg_below ,obj = self, direction = "below", arrangeDf = T),
filter_reduc(df = filters_neg_below ,obj = self, direction = "below", arrangeDf = F),times = 30
)
# Results: 147 ms vs 350, so less impressive but still worth it

microbenchmark(
  filter_reduc(df = filters_neg_below ,obj = self, direction = "above", arrangeDf = T),
  filter_reduc(df = filters_neg_below ,obj = self, direction = "above", arrangeDf = F),times = 30
)
# Results: 32 ms vs 145 !!



# test3 the Filter application ------------------------------------------------------

# strategy number one: a loop for each filter

poolVP_id <- poolVP <- crossing(A = 1:100, B = 1:100, C = 1:100) %>% rowid_to_column("id")

sizeFilter <- 100

filters_neg_above_reduc <- tibble(A = sample(1:100, size = sizeFilter, replace = T),
                                  B = sample(1:100, size = sizeFilter, replace = T),
                                  C = sample(1:100, size = sizeFilter, replace = T))

all_param <- names(filters_neg_above_reduc)
filtre <- "A  <= ref$A & B <= ref$B & C <= ref$C"
filtre <- rep(filtre, 10) %>% paste0(collapse = " & ")


test1 <- function(nfilter = 100, foldrepeafi = 1){

  filters_neg_above_reduc <- tibble(A = sample(1:100, size = nfilter, replace = T),
                                    B = sample(1:100, size = nfilter, replace = T),
                                    C = sample(1:100, size = nfilter, replace = T))

  filtre <- "A  <= ref$A & B <= ref$B & C <= ref$C"
  filtre <- rep(filtre, foldrepeafi) %>% paste0(collapse = " & ")

  t0 <- Sys.time()
  poolVP2_id <- poolVP_id
poolVP2 <- poolVP


  for(a in 1:nrow(filters_neg_above_reduc)){

    ref <- filters_neg_above_reduc %>% slice(a)

    poolVP_id %>%
      mutate(test = !!parse_expr(filtre)) %>%
      filter(test == T) %>% pull(id) ->idtorem

  poolVP2 <- poolVP2 %>%
       filter(! id %in% idtorem)
  # if(length(idtorem) > 0)aeza

  }
difftime( Sys.time(), t0)

}

# strategy number two: a unique giant filter

test2 <- function(nfilter = 100, foldrepeafi = 1){

  filters_neg_above_reduc <- tibble(A = sample(1:100, size = nfilter, replace = T),
                                    B = sample(1:100, size = nfilter, replace = T),
                                    C = sample(1:100, size = nfilter, replace = T))

  filtre <- "A  <= ref$A & B <= ref$B & C <= ref$C"
  filtre <- rep(filtre, foldrepeafi) %>% paste0(collapse = " & ")

  t0 <- Sys.time()

  poolVP2_id <- poolVP_id
  poolVP2 <- poolVP


  temp <- filters_neg_above_reduc %>%
    mutate(filtre = filtre)

for(a in all_param){

  temp %>%
    mutate(filtre = map2_chr(filtre, !!parse_expr(a), function(filtre, x){

      gsub(paste0("ref\\$", a),x,  filtre)

    })) -> temp

}

  temp %>%
    pull(filtre) -> filtres


  filtre_line <- paste0("(", filtres, ")") %>% paste0(collapse = "|")


      poolVP_id %>%
        mutate(test = !!parse_expr(filtre_line)) %>%
        filter(test == T) %>% pull(id) ->idtorem

      poolVP2 <- poolVP2 %>%
        filter(! id %in% idtorem)


      difftime( Sys.time(), t0)
}


microbenchmark(test2(nfilter = 10),
               test2(nfilter = 100),
               test2(nfilter = 1000),times = 1)

microbenchmark(test2(nfilter = 100, foldrepeafi = 1),
               test2(nfilter = 100, foldrepeafi = 10),
               test2(nfilter = 100, foldrepeafi = 50),times = 1)

microbenchmark(test1(), test2(), times = 10, unit = "s")

# Unit: seconds
# expr       min        lq      mean    median        uq       max neval cld
# test1() 0.8546946 0.8763057 0.9578492 0.9113031 0.9971446 1.1882779    10   b
# test2() 0.1702156 0.1973067 0.2252590 0.2073278 0.2179361 0.4288857    10  a

filters_neg_above_reduc <- bind_rows(filters_neg_above_reduc, filters_neg_above_reduc,
                                     filters_neg_above_reduc, filters_neg_above_reduc,filters_neg_above_reduc,
                                     filters_neg_above_reduc,filters_neg_above_reduc,filters_neg_above_reduc)
microbenchmark(test1(), test2(), times = 1, unit = "s")


# Test 3 pratique ---------------------------------------------------------

poolVP2_id <- poolVP_id
poolVP2 <- poolVP


test1 <- function(){

  t0 <- Sys.time()

  poolVP2_id <- poolVP_id
  poolVP2 <- poolVP

  if(nrow(filters_neg_above_reduc) > 0){
    for(a in 1:nrow(filters_neg_above_reduc)){

      ref <- filters_neg_above_reduc %>% slice(a)

      poolVP2_id %>%
        mutate(test = !!parse_expr(filters[[ref$cmt]][["above"]])) %>%
        pull(id) -> idtorem

      poolVP2 <- poolVP2 %>%
        filter(! id %in% idtorem)

    }
  }

  difftime( Sys.time(), t0)

}

test2 <- function(){


  t0 <- Sys.time()

  poolVP2_id <- poolVP_id
  poolVP2 <- poolVP


  temp <- filters_neg_above_reduc %>%
    mutate(filtre = map_chr(cmt,  ~ filters[[.x]][["above"]]))

  for(a in all_param){

    temp %>%
      mutate(filtre = map2_chr(filtre, !!parse_expr(a), function(filtre, x){

        gsub(paste0("ref\\$", a),x,  filtre)

      })) -> temp

  }

  temp %>%
    pull(filtre) -> filtres


  filtre_line <- paste0("(", filtres, ")") %>% paste0(collapse = "|")


  poolVP_id %>%
    mutate(test = !!parse_expr(filtre_line)) %>%
    filter(test == T) %>% pull(id) ->idtorem

  poolVP2 <- poolVP2 %>%
    filter(! id %in% idtorem)


  difftime( Sys.time(), t0)
}


microbenchmark(test1(), test2(), times = 10, unit = "s")




