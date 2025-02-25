load("../../data/surv_data.RData")
load("../../data/fun_data.RData")

fpc_obj <- AR1_PACE(fun_data, surv_data, nbasis=10, pve=0.9) # 90% leads to 3 FPC components
model <- AR1_FRAILTY(formula=Event(t_start, t_stop, id, status) ~ z1, data=surv_data, fpca_obj=fpc_obj, para0=c(0.5, 0.5))
plot_beta(model)
