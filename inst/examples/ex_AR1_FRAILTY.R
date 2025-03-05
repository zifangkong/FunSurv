data(simDat)

fit <- AR1_FRAILTY(formula=Recur(t_start %to% t_stop, id, event=status) ~ z1, sdat=sdat, fdat=fdat)

str(basesurv(fit))
print(fit)
summary(fit)
plot(fit)
