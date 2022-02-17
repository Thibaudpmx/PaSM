# Each time you make a modificaiton of the config file, do a VT_source before trying anything.

# Verification number 1 ---------------------------------------------------

# Just make sur RxoDE compiled your model correctly, with no error message

VT_source()
model_RxODE

# Verification number 2 ---------------------------------------------------
# Here we want to verify that simulations works

VT_source()

# without any input and param modif
simulations()

# with input but no param modif

add_events <- tibble(cmt = c("Venetoclax")) %>% # Any compartments
  mutate(time = 0, amt = 0)

simulations(add_events = add_events)

# with input and param modif

res <- simulations(add_events = add_events, ind_param = c( any_param = 0))

# I encourage you to verify with plots the dynamic of your model
# (Replace Venetoclax by any number of compartments)

res %>%
  gather("key", "value", Venetoclax) %>%
  ggplot()+
  geom_line(aes(time, value))+
  facet_wrap(~key, scales = "free")


# Verification number 3 ---------------------------------------------------
# just verify output is indeed TRUE or FALSE

simulations(add_events = add_events, returnSim = F)

