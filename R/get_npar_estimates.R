#' Get non-parametric estimates given a set of support points, model, and data
#' 
#' @param parameters grid (data.frame) of support points
#' @param error error model
#' @param model `PKPDsim` model
#' @param regimen `PKPDsim` regimen
#' @param t_obs vector of observations
#' @param data vector of obsersved data
#' @export
get_npar_estimates <- function(parameters = NULL, 
                               error = list(prop = 0.05, add = 0.1), 
                               model = NULL,
                               regimen = NULL, 
                               t_obs = c(24),
                               data = NULL) {
  grid <- expand.grid(x = parameters[,1], y = parameters[,2])
  all <- c()
  ll <- c()
  for(i in 1:length(grid[,1])) {
    par_tmp <- list(CL = grid[i, 1], V = grid[i, 2])
    tmp <- sim(ode = model, parameters = par_tmp, regimen = regimen, t_obs = t_obs,
               only_obs = TRUE)$y
    all <- rbind(all, tmp)
    ll <- c(ll, get_likelihood_of_data(data = data, ipred = tmp, error = error))
  }
  tmp <- data.frame(cbind(grid, ll))
  tmp$CL_tmp <- tmp[,1] * tmp$ll 
  tmp$V_tmp  <- tmp[,2] * tmp$ll 
  CL_av <- sum(tmp$CL_tmp) / sum(tmp$ll)
  V_av <- sum(tmp$V_tmp) / sum(tmp$ll)
  return(list(parameters = list(CL_av, V_av), prob = tmp))
}
