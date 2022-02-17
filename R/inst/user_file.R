
# Step 1 create or load project ---------------------------------------------------
create_VT_Project("#path_project")

# Step 2: fulfil the file

#path_file_config
#path_file_verif

# Step3: create celltheque

toadd <- crossing(param_of_choise = seq(1,10,0.1)) # what you want to add in your celltheque as individuals

add_events <- tibble(cmt = c("com_of_choice")) %>%  # Administration sampling, with "concX" being replaced by concentration
  mutate(time = 0, amt = "conc1")                   # mentionned in the configuration files

celltheque_produc(file.name = "test_1.RDS", toadd = toadd, add_events = add_events) # the celltheque production function!

# Step 3 : create or load a virtual tumor

cellt <- cellthequeLoad(drug = 1,return = 2, update = T)


# Step 4 : create or load a virtual tumor

VT <- createVirtualTumor(filter_data = Drug == 1 & conc4 == 0 & Cell_line_bin == 1)

# Step 5 : plot and optimize virtual tumor

VT_plot(VT)

VT_optim(VT)

VT_optim2(VT)


VT_plot(VT)

