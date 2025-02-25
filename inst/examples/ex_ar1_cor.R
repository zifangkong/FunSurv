## ex1
ar1_cor(n=5, rho=0.3)

## ex2
ar1_cor(n=5, rho=2) # return an error message that rho must be between 0 and 1.

### first derivative of the AR(1) structure with respect to rho
dar1_cor.drho(n=5, rho=0.3)
dar1_cor.drho(n=5, rho=3) # return an error message that rho must be between 0 and 1.
