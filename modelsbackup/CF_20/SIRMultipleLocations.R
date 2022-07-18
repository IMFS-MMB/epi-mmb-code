#Multiple Location Eichenbaum et al Model
#Coronavirus and Stay-at-Home Project
#This file specifically runs the code to compare the
#single and multiple location model
#Figure 9
#Note: Individual locations are referred to as counties throughout code
#------------------------------------------------------------------

#Clear Contents
#rm(list = ls())

#Data Packages
library(data.table)

#Equation Solver Packages
library(nleqslv)
library(nloptr)

#Graph Packages
library(ggplot2)
library(ggthemes)


#------------------------------------------------------------------
#------------------------------------------------------------------
#Set Parameters, Compute Steady State, Set Initial Guesses
#------------------------------------------------------------------
#------------------------------------------------------------------
currentdir = getwd()

#Parameters
A_prod = 10         #Productivity (10 makes PV_i 1% lower than PV_s)
theta = 9           #Disutility of labor (9 makes steady state hours 1/3)
phi_inf = 0.8       #Infected labor productivity
beta = 0.98^(1/52)  #Discount rate
pi_s1 = 0.0185      #Prob of infection through consumption
pi_s2 = 1.8496           #Prob of infection through labor2
pi_s3 = 0.2055        #Prob of infection exogenously
pi_d = 0.01*7/18    #Prob of dying
pi_r = 7/18 - pi_d  #Prob of recovering

#Compute Steady State Values
hours_rss = (1/theta)^(1/2)            #labor recovered (post-infection steady state)
c_rss = A_prod*hours_rss                #consumption recovered
u_rss = log(c_rss)-theta/2*hours_rss^2  #utility recovered
val_rss = 1/(1-beta)*u_rss             #PV utility recovered (Steady State)
val_rssConsUnits = val_rss*c_rss       #PV utility in cons. units (Urss*Marg.Util.Cons); value of life
hours_iss = (1/theta)^(1/2)            #labor infected
c_iss = phi_inf*A_prod*hours_iss       #consumption infected
u_iss = log(c_iss)-theta/2*hours_iss^2 #utility infected
val_iss = 1/(1 - beta*(1-pi_r-pi_d))*(u_iss + beta*pi_r*val_rss)  #PV utility infected (Steady State)
val_sss = val_rss #PV utility susceptible (Steady State)

#Initialize Model
num_periods = 100   #Event horizon
num_locations = 21  #Number of states
num_perloc = num_periods*num_locations
Pop_0 = rep(1/num_locations, num_locations)     #Initial Populations
tableinival <- read.table(paste(currentdir , "\\initval.txt",sep = ""))
iniI0 = tableinival[1,1]
#epsilon = rep(0.001, num_locations)             #Initial Infected
epsilon = rep(iniI0, num_locations)             #Initial Infected
hours_S = rep(hours_rss, num_perloc) #Hours Susceptible
hours_I = rep(hours_iss, num_perloc) #Hours Infected
hours_R = rep(hours_rss, num_perloc) #Hours Recovered
policy_guess = rep(0, num_perloc)    #Policy Variable


#------------------------------------------------------------------
#------------------------------------------------------------------
#Create Infection Transition Matrix and Function to Simulate SIR Model
#------------------------------------------------------------------
#------------------------------------------------------------------
#Transition Matrices (Cross-location, Over time)
#Always true, just here to collapse lines
if (TRUE) { 
  #Initializes a data.table with location being infected
  #County refers to virus destination
  dt_distance = data.table(County = rep(1:num_locations, num_locations))
  
  #Denotes which county infection originates and Population of origin county
  #cty_travel refers to virus origin
  dt_distance[, `:=`(cty_travel = 1:num_locations,
                     Pop_0 = Pop_0), by = County]
  
  #Set Travel cost across locations
  dt_distance[, tau_exogenous := 1 - (abs(County - cty_travel)/num_locations), 
              by = County]
  #Adjust for Population
  dt_distance[, tau_exoPop := tau_exogenous*Pop_0]
  
  #Sorts data.table by County being infected
  setkey(dt_distance, County)
  
  #Matrix for transmission across locations
  #Rows denote the county being infected (Home H), columns infection origin (i)
  mat_tranInfect = matrix(dt_distance$tau_exoPop, 
                          nrow = num_locations, byrow = TRUE)
  
  #Standard SIR transmission matrix over time
  #Equations 3-7 in paper, this also includes population
  #T_t column is a placeholder, the value is updated in 'fn_sim'
  mat_tranSIR = matrix(c(c(1, 0, 0, 0, 0, -1), #S_t column
                         c(0, 1-pi_r-pi_d, 0, 0, 0, 1), #I_t column
                         c(0, pi_r, 1, 0, 0, 0), #R_t column
                         c(0, pi_d, 0, 1, 0, 0), #D_t column
                         c(0, -pi_d, 0, 0, 1, 0), #Pop_t column
                         c(0, 0, 0, 0, 0, 1)), #T_t column
                       nrow = 6, ncol = 6)
}

#Function for simulating the virus transmission
#Inputs: Hours and Policy Guesses
#Output: List of simulation results
#Note: This 'list of vectors' is later converted to a 'data.table' for efficiency
fn_sim = function(hours_S, hours_I, hours_R, mu_c){
  #All inputs should be length num_locations*num_periods
  
  #This is set up so arguments are sorted first by period, then county
  #Thus all counties in a given period are located adjacent to each other
  County = rep(1:num_locations, num_periods)
  Pop_init = rep(Pop_0, num_periods)
  period = rep(1:num_periods, times = rep(num_locations, num_periods))
  hours_st = hours_S
  hours_it = hours_I
  hours_rt = hours_R
  policy_var = mu_c
  
  #Solve for Recovered Unknowns
  lambda_rbt = theta*hours_rt/A_prod #Constraint
  c_rt = 1/((1+policy_var)*lambda_rbt) #Consumption
  u_rt = log(c_rt) - (theta/2)*hours_rt^2 #Utility
  Gamma_t = (1+policy_var)*c_rt - A_prod*hours_rt #Transfers
  
  #Solve for Infected Unknowns
  lambda_ibt = theta*hours_it/(A_prod*phi_inf) #Constraint
  c_it = 1/((1+policy_var)*lambda_ibt) #Consumption
  u_it = log(c_it) - (theta/2)*hours_it^2 #Utility
  
  #Solve for Susceptible Unknowns
  c_st = (A_prod*hours_st + Gamma_t)/(1 + policy_var) #Consumption
  u_st = log(c_st) - (theta/2)*hours_st^2 #Utility
  
  
  #------------------------------------------------------------------
  #------------------------------------------------------------------
  #Simulate SIR Model
  #Results in the probability of infection in period t (tau_t)
  #This prob. is needed to compute the PV function for susceptibles
  #------------------------------------------------------------------
  #------------------------------------------------------------------
  
  #Pre-allocate Simulation
  mat_sir = matrix(c(rep(1-epsilon, num_periods), #S_t
                     rep(epsilon, num_periods), #I_t
                     rep(0, num_perloc), #R_t
                     rep(0, num_perloc), #D_t
                     Pop_init, #Pop_t
                     rep(epsilon, num_periods)), #T_t
                   nrow = num_perloc, ncol = 6)
  
  #Set Initial Values
  mat_lag = matrix(c(1-epsilon, #S_t
                     epsilon, #I_t
                     rep(0, num_locations), #R_t
                     rep(0, num_locations), #D_t
                     Pop_0, #Pop_t
                     epsilon), #T_t
                   nrow = num_locations, ncol = 6)
  
  #Newly Transmitted Infections from Virus Outbreak (T_t, Equation 3 in paper)
  #Note: This is the only step of the simulation that differs from ERT
  #Notice that mat_tranInfect is an NxN matrix, mat_lag[, j] is a Nx1 vector
  #This accounts for the virus transmission inherited from all locations
  #and thus results in a Nx1 vector of new transmissions
  period_temp = period == 1
  mat_lag[, 6] = mat_lag[,1]*(
    pi_s1*(c_st[period_temp])*(mat_lag[,2]*c_it[period_temp]) + 
      pi_s2*(hours_st[period_temp])*(mat_lag[, 2]*hours_it[period_temp]) + 
      pi_s3*(mat_tranInfect%*%mat_lag[, 2]))
  
  #Simulate SIR model
  for (i in 2:num_periods){
    #Calculate all variables other than T_t
    period_temp = period == i
    mat_temp = mat_lag%*%mat_tranSIR
    
    #Calculate T_t
    mat_temp[, 6] = mat_temp[,1]*(
      pi_s1*(c_st[period_temp])*(mat_temp[,2]*c_it[period_temp]) + 
        pi_s2*(hours_st[period_temp])*(mat_temp[, 2]*hours_it[period_temp]) + 
        pi_s3*(mat_tranInfect%*%mat_temp[, 2]))
    
    #Update Values
    mat_sir[period_temp, ] = mat_temp
    mat_lag = mat_temp
  }
  
  #Compute period probability of infection
  tau_t = mat_sir[, 6]/mat_sir[, 1]
  
  #------------------------------------------------------------------
  #------------------------------------------------------------------
  #Solve for Present Value Functions and Lagrange multiplier for infections
  #------------------------------------------------------------------
  #------------------------------------------------------------------
  
  #Initialize terminal conditions
  val_rt = u_rt + beta*val_rss
  val_it = u_it + beta*((1-pi_r-pi_d)*val_iss + pi_r*val_rss)
  val_st = u_st + beta*((1-tau_t)*val_sss + tau_t*val_iss)
  
  #Placeholders
  period_temp = period == num_periods
  val_rlead = val_rt[period_temp]
  val_ilead = val_it[period_temp]
  val_slead = val_st[period_temp]
  
  #Pre-allocate
  lambda_tau_t = rep(0, num_perloc)
  lambda_tau_t[period_temp] = beta*(val_iss - val_sss)
  
  #Recursively compute PV functions
  for (i in (num_periods-1):1){
    #Period placeholder
    period_temp = period == i
    tau_temp = tau_t[period_temp]
    
    #Calculate values
    val_rtemp = u_rt[period_temp] + beta*val_rlead 
    val_itemp = u_it[period_temp] + beta*((1-pi_r-pi_d)*val_ilead + pi_r*val_rlead) 
    val_stemp = u_st[period_temp] + beta*((1-tau_temp)*val_slead + 
                                            tau_temp*val_ilead)
    
    #Allocate values
    val_rt[period_temp] = val_rtemp
    val_it[period_temp] = val_itemp
    val_st[period_temp] = val_stemp
    lambda_tau_t[period_temp] = beta*(val_ilead - val_slead)
    
    #Update placeholders
    val_rlead = val_rtemp
    val_ilead = val_itemp
    val_slead = val_stemp
  }
  
  
  #------------------------------------------------------------------
  #------------------------------------------------------------------
  #Solve for susceptible lagrange multipliers first
  #------------------------------------------------------------------
  #------------------------------------------------------------------
  
  #Lagrange multiplier for budget constraint
  lambda_sbt = (1/c_st + lambda_tau_t*pi_s1*mat_sir[, 2]*c_it)/(1 + policy_var)
  
  #Compile Results
  dt_sim = list('period' = period, 'county' = County, 
                'S_t' = mat_sir[, 1], 'I_t' = mat_sir[, 2], 
                'R_t' = mat_sir[, 3], 'D_t' = mat_sir[, 4],
                'hours_st' = hours_st, 'hours_it' = hours_it, 
                'hours_rt' = hours_rt, 'policy_var' = policy_var, 
                'c_st' = c_st, 'c_it' = c_it, 'c_rt' = c_rt, 
                'val_st' = val_st, 'val_it' = val_it, 'val_rt' = val_rt, 
                'Gamma_t' = Gamma_t, 'lambda_sbt' = lambda_sbt, 
                'lambda_tau_t' = lambda_tau_t)
  return(dt_sim)
}

#------------------------------------------------------------------
#------------------------------------------------------------------
#Functions to solve for hours and policy
#------------------------------------------------------------------
#------------------------------------------------------------------

#Function for solving optimal labor decision
#Inputs: Hours and Policy Guesses
#Output: Minimization error for optimal labor decision
fn_min_labor = function(x, policy_guess){
  
  #Simulate model
  dt_temp = fn_sim(hours_S = x[1:(num_perloc)],
                   hours_I = x[(num_perloc+1):(2*num_perloc)],
                   hours_R = x[(2*num_perloc+1):(3*num_perloc)],
                   mu_c = policy_guess)
  
  #Pre-allocate error
  err_labor = rep(0, 3*num_perloc)
  
  #Note: 'Three' unknowns (hours_st, hours_it, hours_rt), so 'three' eqns.
  #Equation 1
  err_labor[1:(num_perloc)] = (1+dt_temp$policy_var)*dt_temp$c_it - 
    phi_inf*A_prod*dt_temp$hours_it - dt_temp$Gamma_t
  #Equation 2
  err_labor[(num_perloc+1):(2*num_perloc)] = 
    -dt_temp$Gamma_t*(dt_temp$S_t + dt_temp$I_t + dt_temp$R_t) + 
    dt_temp$policy_var*(dt_temp$S_t*dt_temp$c_st + dt_temp$I_t*dt_temp$c_it + 
                          dt_temp$R_t*dt_temp$c_rt)
  #Equation 3
  err_labor[(2*num_perloc+1):(3*num_perloc)] = 
    -theta*dt_temp$hours_st + 
    A_prod*dt_temp$lambda_sbt + 
    dt_temp$lambda_tau_t*pi_s2*dt_temp$I_t*dt_temp$hours_it
  
  return(err_labor)
}

#Function for solving optimal policy decision
#Inputs: Policy Guess (length: num_periods*num_locations)
#Output: Total Country Utility
fn_min_policy = function(x){
  
  #Solve for optimal hours given policy guess (x)
  hours_temp = nleqslv(c(hours_S_guess, hours_I_guess, hours_R_guess), 
                       fn_min_labor, jac = NULL, 
                       policy_guess = x,
                       method = 'Broyden', global = 'dbldog',
                       control = list(allowSingular = FALSE, 
                                      xtol = 1e-04, ftol = 1e-04))
  #Store labor decisions
  hours_S_guess <<- hours_temp$x[1:(num_perloc)]
  hours_I_guess <<- hours_temp$x[(num_perloc+1):(2*num_perloc)]
  hours_R_guess <<- hours_temp$x[(2*num_perloc+1):(3*num_perloc)]
  
  #Simulate model with optimal hours
  dt_temp = fn_sim(hours_S = hours_S_guess, 
                   hours_I = hours_I_guess, 
                   hours_R = hours_R_guess, 
                   mu_c = x)
  
  #Policy is to maximize total country utility
  err_U = -Pop_0%*%(dt_temp$S_t[1:num_locations]*dt_temp$val_st[1:num_locations] + 
              dt_temp$I_t[1:num_locations]*dt_temp$val_it[1:num_locations] + 
              dt_temp$R_t[1:num_locations]*dt_temp$val_rt[1:num_locations])

  return(err_U)
}

#Function for solving optimal national policy decision
#Inputs: Policy Guess (length: num_periods)
#Output: Total Country Utility
fn_min_policy_natl = function(x){
  
  #Expand same policy to all locations
  pol_natl = rep(x, each = num_locations)
  
  #Get optimal hours given policy
  hours_temp = nleqslv(c(hours_S_guess, hours_I_guess, hours_R_guess), 
                       fn_min_labor, jac = NULL, 
                       policy_guess = pol_natl,
                       method = 'Broyden', global = 'dbldog',
                       control = list(allowSingular = FALSE, 
                                      xtol = 1e-04, ftol = 1e-04))
  #Store labor decisions
  hours_S_guess <<- hours_temp$x[1:(num_perloc)]
  hours_I_guess <<- hours_temp$x[(num_perloc+1):(2*num_perloc)]
  hours_R_guess <<- hours_temp$x[(2*num_perloc+1):(3*num_perloc)]
  
  #Simulate model with optimal hours
  dt_temp = fn_sim(hours_S = hours_S_guess, 
                   hours_I = hours_I_guess, 
                   hours_R = hours_R_guess, 
                   mu_c = pol_natl)
  
  #Policy is to maximize total country utility
  err_U = -Pop_0%*%(dt_temp$S_t[1:num_locations]*dt_temp$val_st[1:num_locations] + 
              dt_temp$I_t[1:num_locations]*dt_temp$val_it[1:num_locations] + 
              dt_temp$R_t[1:num_locations]*dt_temp$val_rt[1:num_locations])

  return(err_U)
}


#------------------------------------------------------------------
#------------------------------------------------------------------
#Model Simulations
#Simulate 21-location vs 1-location model
#------------------------------------------------------------------
#------------------------------------------------------------------
#21-location Model (Always true, just here to collapse lines)
#Note: 21 locations will take ~10 minutes
if (TRUE){
#Solve for and store optimal hours given policy guess
opt_hours = nleqslv(c(rep(hours_rss, num_perloc), 
                      rep(hours_iss, num_perloc), 
                      rep(hours_rss, num_perloc)), 
                    fn_min_labor, jac = NULL, 
                    policy_guess = policy_guess,
                    method = 'Broyden',
                    global = 'dbldog',
                    control = list(allowSingular = FALSE, trace = 1,
                                   xtol = 1e-04, ftol = 1e-04))

#Simulate model with optimal hours and save results as data.table
dt_results = as.data.table(fn_sim(hours_S = opt_hours$x[1:(num_perloc)],
                                  hours_I = opt_hours$x[(num_perloc+1):(2*num_perloc)],
                                  hours_R = opt_hours$x[(2*num_perloc+1):(3*num_perloc)],
                                  mu_c = policy_guess))

#Sum consumption and hours across agents within location and period
dt_results[, `:=`(C_cty = S_t*c_st + I_t*c_it + R_t*c_rt, 
                  hours_cty = S_t*hours_st + I_t*hours_it + R_t*hours_rt)]

#Compute national totals
dt_results[, `:=`(S_Agg = sum(Pop_0*S_t), 
                  I_Agg = sum(Pop_0*I_t),
                  C_Agg = sum(Pop_0*C_cty),
                  D_Agg = sum(Pop_0*D_t),
                  R_Agg = sum(Pop_0*R_t),
                  hours_Agg = sum(Pop_0*hours_cty)), 
           by = period]
}

#Single-location
#Note: 1 locations will take <10 secs
if(TRUE){
  #Re-initialize Model
  num_periods = 100   #Event horizon
  num_locations = 1  #Number of states
  num_perloc = num_periods*num_locations
  Pop_0 = rep(1/num_locations, num_locations)     #Initial Populations
  epsilon = rep(0.001, num_locations)             #Initial Infected
  hours_S = rep(hours_rss, num_perloc) #Hours Susceptible
  hours_I = rep(hours_iss, num_perloc) #Hours Infected
  hours_R = rep(hours_rss, num_perloc) #Hours Recovered
  policy_guess = rep(0, num_perloc)    #Policy Variable
  
  #Transition matrix
  if (TRUE){ 
    #Initializes a data.table with location being infected
    #County refers to virus destination
    dt_distance = data.table(County = rep(1:num_locations, num_locations))
    
    #Denotes which county infection originates and Population of origin county
    #cty_travel refers to virus origin
    dt_distance[, `:=`(cty_travel = 1:num_locations,
                       Pop_0 = Pop_0), by = County]
    
    #Set Travel cost across locations
    dt_distance[, tau_exogenous := 1 - (abs(County - cty_travel)/num_locations), 
                by = County]
    #Adjust for Population
    dt_distance[, tau_exoPop := tau_exogenous*Pop_0]
    
    #Sorts data.table by County being infected
    setkey(dt_distance, County)
    
    #Rows denote the county being infected (Home H), columns infection origin (i)
    mat_tranInfect = matrix(dt_distance$tau_exoPop, 
                            nrow = num_locations, byrow = TRUE)
  }
  
  #Solve for and store optimal hours given policy guess
  opt_hours = nleqslv(c(rep(hours_rss, num_perloc), 
                        rep(hours_iss, num_perloc), 
                        rep(hours_rss, num_perloc)), 
                      fn_min_labor, jac = NULL, 
                      policy_guess = policy_guess,
                      method = 'Broyden',
                      global = 'dbldog',
                      control = list(allowSingular = FALSE, trace = 1,
                                     xtol = 1e-04, ftol = 1e-04))

  #Simulate model with optimal hours and save results as data.table
  dt_results_base = as.data.table(fn_sim(hours_S = opt_hours$x[1:(num_perloc)],
                                    hours_I = opt_hours$x[(num_perloc+1):(2*num_perloc)],
                                    hours_R = opt_hours$x[(2*num_perloc+1):(3*num_perloc)],
                                    mu_c = policy_guess))

  #Merge Single-location results
  dt_results_base[, `:=`(S_baseline = S_t,
                         I_baseline = I_t,
                         C_baseline = S_t*c_st + I_t*c_it + R_t*c_rt, 
                         hours_baseline = S_t*hours_st + I_t*hours_it + R_t*hours_rt)]
  dt_results = merge(dt_results, 
                     dt_results_base[, c('S_baseline', 'I_baseline', 
                                       'C_baseline', 'hours_baseline', 'period')], 
                   by = 'period')

}

paramnames = c("c_rss", "hours_rss")
paramvals = c(c_rss,hours_rss )
dt_params = data.frame(name = paramnames, values = paramvals )

#------------------------------------------------------------------
#------------------------------------------------------------------
#Make Graphs
#Susceptible, Infected, Consumptions, Hours
#------------------------------------------------------------------
#------------------------------------------------------------------
if(0){
  #Susceptible
  ggplot(dt_results) + 
    geom_line(aes(x = period, y = 100*S_t, group = county), 
              color = 'grey', alpha = 0.7) +
    geom_line(aes(x = period, y = 100*S_Agg, group = 1, color = 'Multiple Locations'), 
              size = 1.2) +
    geom_line(aes(x = period, y = 100*S_baseline, group = 1, color = 'Baseline'), 
              size = 1.2) +
    scale_color_stata() +
    scale_y_continuous(breaks = seq(40, 100, 10),
                       labels = seq(40, 100, 10),
                       limits = c(40, 100)) +
    labs(color = '', x = 'Weeks', y = '% of Initial Population') +
    theme(axis.text = element_text(size = 24),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 24)) +
    theme_hc()
  
  #Infected
  ggplot(dt_results) + 
    geom_line(aes(x = period, y = 100*I_t, group = county), 
              color = 'grey', alpha = 0.7) +
    geom_line(aes(x = period, y = 100*I_Agg, group = 1, color = 'Multiple Locations'), 
              size = 1.2) +
    geom_line(aes(x = period, y = 100*I_baseline, group = 1, color = 'Baseline'), 
              size = 1.2) +
    scale_color_stata() +
    scale_y_continuous(breaks = seq(0, 5, 1),
                       labels = seq(0, 5, 1),
                       limits = c(0, 5.2)) +
    labs(color = '', x = 'Weeks', y = '% of Initial Population') +
    theme(axis.text = element_text(size = 24),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 24)) +
    theme_hc()
  
  #Consumption
  ggplot(dt_results) + 
    geom_line(aes(x = period, y = 100*C_cty/c_rss - 100, group = county), 
              color = 'grey', alpha = 0.7) +
    geom_line(aes(x = period, y = 100*C_Agg/c_rss - 100, group = 1, 
                  color = 'Multiple Locations'), 
              size = 1.2) +
    geom_line(aes(x = period, y = 100*C_baseline/c_rss - 100, group = 1, 
                  color = 'Baseline'), 
              size = 1.2) +
    scale_color_stata() +
    scale_y_continuous(breaks = seq(-14, 0, 2), 
                       labels = seq(-14, 0, 2),
                       limits = c(-10, 0)) +
    labs(color = '', x = 'Weeks', y = '% Deviation from SS') +
    theme(axis.text = element_text(size = 24),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 24)) +
    theme_hc()
  
  #Hours
  ggplot(dt_results) + 
    geom_line(aes(x = period, y = 100*hours_cty/hours_rss - 100, group = county), 
              color = 'grey', alpha = 0.7) +
    geom_line(aes(x = period, y = 100*hours_Agg/hours_rss - 100, group = 1, 
                  color = 'Multiple Locations'), 
              size = 1.2) +
    geom_line(aes(x = period, y = 100*hours_baseline/hours_rss - 100, group = 1, 
                  color = 'Baseline'), 
              size = 1.2) +
    scale_color_stata() +
    scale_y_continuous(breaks = seq(-12, 0, 2), 
                       labels = seq(-12, 0, 2),
                       limits = c(-10, 0)) +
    labs(color = '', x = 'Weeks', y = '% Deviation from SS') +
    theme(axis.text = element_text(size = 24),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 24)) +
    theme_hc()
}
  

write.csv(dt_results,paste(currentdir , "\\dtresults.csv",sep = ""), row.names = FALSE)
write.csv(dt_params,paste(currentdir , "\\dtparams.csv",sep = ""), row.names = FALSE)

  