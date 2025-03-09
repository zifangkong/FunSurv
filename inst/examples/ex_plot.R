fit <- AR1_FRAILTY(formula=Recur(t_start %to% t_stop, id, event=status) ~ z1, sdat=sdat, fdat=fdat)
plot(fit, which="beta")
plot(fit, which="fpc")
plot(fit, which="basesurv")
