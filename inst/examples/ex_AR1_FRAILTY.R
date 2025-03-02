data(simDat)

## Calculate event-specific FPC scores
## 90% leads to 3 FPC components
fpc_obj <- AR1_PACE(fdat, sdat, nbasis = 10, pve = 0.9) 

## Fit the functional regression with autoregressive frailty 
fit <- AR1_FRAILTY(formula = Events(t_start, t_stop, id, status) ~ z1, data = sdat, fpca_obj = fpc_obj, para0 = c(0.5, 0.5))

basesurv(fit)
plot(fit)


