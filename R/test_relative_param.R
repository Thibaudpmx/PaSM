# library(peccary)
# library(QSPVP)
# library(RxODE)
# library(progress)
# Step 1 create or load project ---------------------------------------------------
# create_VT_Project("D:/these/Second_project/QSP/modeling_work/VT_simeoni")
#
# values <- tibble(k1 = 10,
#                   k2 = 2,
#                   ke =1,
#                   lambda0 = 2,
#                   lambda1 = 100,
#                   Vd =  40
# )
#
#
#                   # protocols = c("dose0", "dose50", "dose100")
#
# model_RxODE <-  RxODE({
#   d/dt(Central) <- -ke * Central
#   Conc <- Central/Vd
#
#   tumVol <- X1 + X2 + X3 + X4
#   growth <- lambda0 * X1/((1 + (lambda0 * tumVol/lambda1)^psi)^(1/psi))
#   d/dt(X1) <- growth - X1 * Conc * k2
#   d/dt(X2) <- X1 * Conc * k2 - k1 * X2
#   d/dt(X3) <- k1 * (X2 - X3)
#   d/dt(X4) <- k1 * (X3 - X4)
# })
#
#
# model_RxODE <-  RxODE({
#   d/dt(Central) <- -ke * Central
#   Conc <- Central/Vd
#
#   tumVol <- X1 + X2 + X3 + X4
#
#   if(tumVol < lambda1/lambda0){
#
#   growth <- lambda0 * X1
#
#   }else{
#
#     growth <- lambda1
#   }
#
#
#   d/dt(X1) <- growth - X1 * Conc * k2
#   d/dt(X2) <- X1 * Conc * k2 - k1 * X2
#   d/dt(X3) <- k1 * (X2 - X3)
#   d/dt(X4) <- k1 * (X3 - X4)
# })


# goal : take a model, a number of parameter to try, to extract the value

find_relative <- function(output,  values = toadd, param = NULL, protocol = "dose50", sensitivity = 1E-6){

  # output  = expr(tumVol)
  output <- enexpr(output)

  if(is.null(param)) param <- names(values)

  ref <-  simulations( values, add_events = protocols[[protocol]]) %>%
      mutate(ref = !!output) %>%
      select(time, ref)

  tibble(param) %>%
    mutate(test = map_chr(param, function(x){

     values2 <- values
     values2[[x]] <- values2[[x]] * 2

     simx2 <- simulations( values2, add_events = protocols[[protocol]])

     for_analysis <- ref %>%
       left_join(simx2, by = "time") %>%
       mutate(test = !!output - ref) %>%
       mutate(test = if_else(abs(test) < sensitivity, 0, test)) %>%
       pull(test)

     lengthref <- length(for_analysis)

     # is it always higher?
     if(sum(for_analysis == 0) == lengthref) return("all0_indep_or_need_increase")
    if(sum(for_analysis >= 0) == lengthref) return("increase")
    if(sum(for_analysis <= 0) == lengthref) return("decrease")

     return("time_dep")

    }))
}


#
# values <- tibble(k1 = 0.5,
#                  k2 = 2,
#                  ke =1,
#                  lambda0 = 2,
#                  lambda1 = 4,
#                  Vd =  40
# )
# x <- "k1"
# values2 <- values
# values2[[x]] <- values2[[x]] * 10
#
# simulations(ind_param =  values, add_events = protocols[[protocol]]) %>%
#   mutate(test = "ref") %>%
#   bind_rows(
# simulations( values2, add_events = protocols[[protocol]]) %>%
#   mutate(test = "x2")
#   ) %>%
#   gather("key", "value", tumVol, X1, X2, X3) %>%
#   filter(time  < 20) %>%
#   ggplot()+
#
#   geom_line(aes(time, value, col = test))+
#   facet_wrap(~key, scales = "free")
#
#
# simulations(ind_param =  values, add_events = protocols[[protocol]]) %>%
#   mutate(test = "ref") %>%
#   bind_rows(
#     simulations( values2, add_events = protocols[[protocol]]) %>%
#       mutate(test = "x2")
#   ) %>%
#   # gather("key", "value", tumVol, X1, X2, X3) %>%
#   # filter(time  < 20) %>%
#   ggplot()+
#
#   geom_line(aes(time, tumVol, col = test))
# 2**(1:10)
#
#
#
# # Simeoni demo plot -------------------------------------------------------
#
# # ind_param <- tibble(k1 = c(0.1, 10), k2 = c(0.5, 0.5), ke = c(1, 1), lambda0 = c(1, 1), lambda1 = c(70, 70), psi = c(20, 20), Vd = c(5, 5), nameparset = c("k1 = 0.1", "k1 = 10"))
# #
# # parameters_df <- ind_param %>%
# #   rownames_to_column("id") %>%
# #   group_by(id, nameparset) %>%
# #   nest(.key = "parameter")
# #
# # times <- seq(0L, 20L, by = 0.1)
# #
# # events_df <- tibble(Proto = "1", cmt = "Central", time = "0", amt = "50", method = "add", ADM = "1", Perf = structure(1L, .Label = c("None", "rate", "time"), class = c("ordered", "factor")), evid = 1) %>%
# #   mutate(time = map(time, ~eval(parse_expr(.x)))) %>%
# #   unnest(time) %>%
# #   mutate(Proto = "") %>%
# #   bind_rows(crossing(time = seq(0L, 20L, by = 0.1), evid = 0, cmt = c("Central", "X1", "X2", "X3", "X4", "Conc", "tumVol", "growth"), Proto = "")) %>%
# #   arrange(time) %>%
# #   mutate(amt = as.double(amt), time = as.double(time)) %>%
# #
# #   group_by(Proto) %>%
# #   nest(.key = "events")
# #
# # states <- c(Central = 0, X1 = 50, X2 = 0, X3 = 0, X4 = 0)
# #
# # model <- RxODE({
# #   d/dt(Central) <- -ke * Central
# #   Conc <- Central/Vd
# #   tumVol <- X1 + X2 + X3 + X4
# #   growth <- lambda0 * X1/((1 + (lambda0 * tumVol/lambda1)^psi)^(1/psi))
# #   pctProlif <- X1 / tumVol
# #   d/dt(X1) <- growth - X1 * Conc * k2
# #   d/dt(X2) <- X1 * Conc * k2 - k1 * X2
# #   d/dt(X3) <- k1 * (X2 - X3)
# #   d/dt(X4) <- k1 * (X3 - X4)
# # })
# #
# #
# # simulations <- crossing(parameters_df, events_df) %>%
# #   mutate(simul = map2(parameter, events, function(parameter, events) {
# #     as_tibble(model$solve(parameter, events, states))
# #   })) %>%
# #   select(-parameter, -events) %>%
# #   unnest
# #
# #
# # breaks2 <- map(-40:40, ~1:9 * 10^.x) %>%
# #   reduce(c)
# #
# # labels2 <- as.character(breaks2)
# #
# # labels2[-seq(1, length(labels2), 9)] <- ""
# #
# # simulations %>%
# #   gather(-id, -nameparset, -Proto, -time, key = "key", value = "value") %>%
# #   mutate(color = paste0(nameparset, "\n", Proto, "\n") %>%
# #            gsub(pattern = "\n\n\n?", replacement = "\n")) %>%
# #   filter(key %in% c("tumVol")) %>%
# #   ggplot + geom_line(aes(time, value, col = fct_reorder(color, as.double(id))), size = 1.5) +
# #   labs(x = "Time", y = "", col = "") +
# #   theme_bw() +
# #   facet_wrap(~key, scales = "free") +
# #   scale_y_log10()+
# #   geom_hline(data = tibble(key = "tumVol", value = 70), aes(yintercept = value, lty = "threshold"))+
# #   scale_linetype_manual(values = 2)
# #   # # geom_vline(xintercept = 1.7, col = "red")+
# #   # geom_vline(xintercept = 5.4, col = "blue")
# #
# # simulations %>%
# #  filter(tumVol >= 68) %>%
# #   group_by(nameparset) %>%
# #   slice(1)
