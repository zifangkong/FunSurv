data(simDat)
fit <- AR1_FRAILTY(Recur(t_start %to% t_stop, id, status) ~ z1,
                   sdat = sdat, fdat = fdat, iter.max = 30)
plot(fit, what = "beta")
plot(fit, what = "fpc")
plot(fit, what = "basesurv")
