data(simDat)

fit <- AR1_FRAILTY(formula=Recur(t_start %to% t_stop, id, event=status) ~ z1, sdat=sdat, fdat=fdat)

basesurv(fit)
plot(fit, which="basesurv")
