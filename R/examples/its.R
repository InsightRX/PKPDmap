library(PKPDsim)
library(PKPDmap)

############################################################
## Example A
############################################################
modA <- new_ode_model(code = "dAdt[1] = -KA*A[1]; dAdt[2] = KA*A[1] -(CL/V)*A[2];", 
                      obs = list(cmt = 2, scale = "V"))
datA <- read.csv(file = "~/git/nonmem_examples/Hands_onA/acop.csv")
dosesA <- datA[datA$evid == 1,]
tdmA <- datA[datA$evid == 0 & datA$id < 101,] 
tdmA <- tdmA[,c("id", "time", "dv", "evid")]
colnames(tdmA) <- c("id", "t", "y", "evid")
parA <- list(KA = 2, CL = 20, V = 200)
regA <- new_regimen(amt = 10000, time = c(0), type = "bolus")
fitA <- run_its(parameters = parA, regimen = regA, model = modA, 
               omega = c(0.5, 0.1, 0.5), 
               fixed = "KA",
               err = list(add = 0.1, prop = 0.05),
               data = tdmA, max_iter = 10)

#library(ggplot2)
ggplot(data.frame(fitA$res), aes(x = ipred, y = dv)) + geom_point()
ggplot(data.frame(fitA$res), aes(x = pred, y = dv)) + geom_point()
fitA$parameters

############################################################
## Example B
############################################################
modB <- new_ode_model(code = "dAdt[1] = -(CL/V)*A[1];", obs = list(scale = "V/1000"))
cols <- c("ID", "TIME","DV", "AMT", "EVID", "CRCL", "AGE", "SEX")
datB <- read.table(file = "~/git/nonmem_examples/Hands_onB/pktab1", skip=1, header=T)[,cols]
colnames(datB) <- tolower(colnames(datB))
colnames(datB)[c(2:3)] <- c("t", "y")
dosesB <- datB[datB$evid == 1 & datB$id == 1,]
tdmB <- datB[datB$evid == 0 & datB$id < 101,] 
parB <- list(CL = 5, V = 50)
reg <- new_regimen(amt = 100, n = 2, time = c(0, 24, 48), type = "bolus")

fitB <- run_its(
  parameters = parB, 
  regimen = reg, 
  model = modB, 
  omega = c(0.5, 0.1, 0.5), 
  err = list(add = 0.1, prop = 0.05),
  data = tdmB,
  max_iter = 10)
