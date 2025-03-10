data(simDat)

fit <- AR1_FRAILTY(Recur(t_start %to% t_stop, id, status) ~ z1,
                   sdat = sdat, fdat = fdat)

summary(fit)
