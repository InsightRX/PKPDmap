library(PKPDsim)
library(PKPDmap)
library(dplyr)

############################################################
## Example A
############################################################
modA <- new_ode_model(code = "dAdt[1] = -KA*A[1]; dAdt[2] = KA*A[1] -(CL/V)*A[2];", 
                      obs = list(cmt = 2, scale = "V"))
datA <- read.csv(file = "~/git/nonmem_examples/Hands_onA/acop.csv")
dosesA <- datA %>% filter(evid == 1)
tdmA <- datA %>% filter(evid == 0, id < 101) %>%
  select(id = id, t = time, y = dv, evid = evid) #%>%
#  filter(t %in% c(4, 19))
parA <- list(KA = 2, CL = 20, V = 200)
regA <- new_regimen(amt = 10000, time = c(0), type = "bolus")
fitA <- run_its(parameters = parA, regimen = regA, model = modA, 
               omega = c(0.5, 0.1, 0.5), 
               fixed = "KA",
               err = list(add = 0.1, prop = 0.05),
               data = tdmA, max_iter = 10)

#library(ggplot2)
fitA$res %>% data.frame() %>%
  ggplot(aes(x = ipred, y = dv)) + geom_point()
fitA$res %>% data.frame() %>%
  ggplot(aes(x = pred, y = dv)) + geom_point()
fitA$parameters

############################################################
## Example B
############################################################

modB <- new_ode_model(code = "dAdt[1] = -(CL/V)*A[1];", obs = list(scale = "V/1000"))
datB <- read.table(file = "~/git/nonmem_examples/Hands_onB/pktab1", skip=1, header=T) %>%
   dplyr::select(id = ID, t = TIME, y = DV, amt = AMT, evid = EVID, crcl = CRCL, age = AGE, sex = SEX)
dosesB <- datB %>% filter(evid == 1, id == 1)
tdmB <- datB %>% filter(evid == 0, id < 101) 
parB <- list(CL = 5, V = 50)
reg <- new_regimen(amt = 100, n = 2, time = c(0, 24, 48), type = "bolus")

fit <- run_its(parameters = par,regimen = reg,model = mod, 
        omega = c(0.5, 0.1, 0.5), 
        err = list(add = 0.1, prop = 0.05),
        data = tdm, max_iter = 10)
