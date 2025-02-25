## ex1
Events(1:5, 2:6, id=1:5, status=c(1,1,0,1,0)) 

## ex2
load("../../data/surv_data.RData")
Events(surv_data$t_start, surv_data$t_stop, id=surv_data$id, status=surv_data$id)
