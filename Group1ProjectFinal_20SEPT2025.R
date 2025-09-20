########################################################################
# PROJECT BRIEF / TASK DESCRIPTION
# ----------------------------------------------------------------------
# At the onset of the rainy season in a village of fifteen thousand people
# in Sub-Saharan Africa, about eight mosquitoes buzz around each person.
# Twenty-five villagers are already infectious with a vector-borne disease,
# and about twelve percent of the mosquitoes carry the parasite.
#
# In humans, the average infectious period is about twenty-five days, after
# which individuals immediately return to the susceptible class. Mosquitoes
# live only around ten days, bite humans at a rate of 0.3 bites per day,
# seek hosts at a rate of 0.1, and each bite has a 0.45 probability of
# transmitting infection from mosquito to human. Records suggest that the
# transmission rate from infectious humans to mosquitoes is about 0.5.
# Once infected, mosquitoes remain infectious until they die. Both the human
# and mosquito populations are assumed to be closed during the period of
# analysis (no births/deaths/migration in humans; mosquito population kept
# roughly constant by replacement births).
#
# Tasks for Participants (as provided)
# 1. Write the flow diagram of the baseline model (i.e. without interventions),
#    its assumptions, model parameters and the simulation code.
# 2. Simulate the dynamics of this vector-borne disease over 365 days without
#    interventions.
#    a. Calculate the prevalence of the disease in the population at the end
#       of the simulation period.
#    b. Provide an interpretation of the computed prevalence.
# 3. Update the model structure: consider the following intervention scenarios
#    i. Scenario 1: Vector control — decreases the host seeking rate of
#       mosquitoes by 60%. Assume the vector control program has an annual
#       budget of 50,000 USD to implement the intervention.
#    ii. Scenario 2: Vaccination — moving a proportion of susceptible humans
#        into a vaccinated class. Immunity wanes on average after 30 days,
#        returning individuals to susceptibility. About 30% of the susceptible
#        population are vaccinated and the vaccine has an efficacy of 0.65.
#        Assume the cost of vaccination is 20 USD per person vaccinated.
#    iii. Scenario 3: Combined intervention — applying both vector control and
#         vaccination simultaneously, with the costs of both interventions
#         considered.
#    a. Write the flow diagram for each scenario to include the updated
#       information on the model structure and the proposed intervention.
#    b. State the model assumptions.
#    c. Give an example of an intervention that has a similar mode of action.
# 4. Obtain the model parameters of your updated model structure described in
#    the narrative.
# 5. Write the equations of the updated model structure: translate the flow
#    diagram into a system of differential equations.
# 6. Update the code: modify the R code to include the updated model structure.
# 7. Quantify intervention impact: run simulations for each intervention and
#    in combination and then calculate infections averted compared to the
#    baseline at the end of the simulation period. Interpret the results.
# 8. Cost-effectiveness: calculate the cost per case averted for each
#    intervention and in combination compared to the baseline scenario after
#    365 days. Which scenario gives better value for money? Interpret results.
#
# ----------------------------------------------------------------------
# SHORT WORKFLOW / SCRIPT GOALS
# ----------------------------------------------------------------------
# This R script implements and analyses a mathematical SIS–SI vector–host
# model for the scenario above and compares three intervention scenarios:
#   - Baseline (no intervention)
#   - Vector control (60% reduction in host-seeking)
#   - Continuous vaccination (target ~30% protected, VE=0.65, waning ~30 days)
#   - Combined (vector control + vaccination)
#
# Outputs produced by this script:
#  - Time-series for S, I (humans) and S, I (mosquitoes) over 365 days
#  - Cumulative incidence (humans and mosquitoes)
#  - Prevalence at day 365 and interpretation
#  - Infections averted vs baseline for each intervention
#  - Cost-per-case-averted (USD) and a cost-effectiveness table/plot
#  - Interactive plots for visualization (plotly)
#
# ----------------------------------------------------------------------
# REQUIRED LIBRARIES (and purpose)
# ----------------------------------------------------------------------
# deSolve   : solves systems of ordinary differential equations (ODEs)
# tidyverse : data wrangling and reshaping (dplyr, tidyr, tibble)
# plotly    : interactive plotting of time-series and cost-effectiveness
#
# ----------------------------------------------------------------------
# KEY PARAMETERS (used in the script)
# ----------------------------------------------------------------------
# Nh (human population)            = 15,000
# k (mosquitoes per human)         = 8
# Nm (total mosquitoes)            = Nh * k = 120,000
# Ih0 (initial infectious humans)  = 25
# Im0 (initial infectious mosq)    = 12% of Nm
# γ (human recovery rate)          = 1/25 per day
# μ_m (mosquito death rate)        = 1/10 per day
# biting rate                      = 0.3 bites per mosquito per day
# seeking                          = 0.1 host-seeking success
# P_mh (mosq->human prob)          = 0.45 per infectious bite
# β (human->mosq transmission)     = 0.5
#
# Intervention parameters:
#  - Vector control: 60% reduction in seeking; annual cost $50,000
#  - Vaccination: VE = 0.65; waning ω = 1/30 per day; ν calibrated to ~30%
#                 target protected fraction; cost per dose = $20

########################################################################
# SIS–SI Vector–Host Model — Continuous Vaccination
########################################################################

## ---------------------------------------------------------
## 0. Load libraries
## ---------------------------------------------------------
library(deSolve)    # ODE solver
library(tidyverse)  # Data wrangling
library(plotly)     # Interactive plotting

## =========================================================
## TASK 1: Baseline model - flow diagram, assumptions, params
## =========================================================

## 1A. Flow diagram (Baseline, *no* vaccination compartments)
## ---------------------------------------------------------
## Humans (SIS):
##   Sh  --(infection: λ_h * Sh)-->  Ih
##   Ih  --(recovery: γ * Ih)----->  Sh
##
## Mosquitoes (SI with replacement births):
##   Sm  --(infection: λ_m * Sm)-->  Im
##   Sm, Im --(death: μ_m * Sm/Im)--> removed
##   births = μ_m * Nm replace deaths so Nm ~ constant
##
## Trackers:
##   CH: cumulative human infections (∫ λ_h * Sh dt)
##   CM: cumulative mosquito infections (∫ λ_m * Sm dt)

## 1B. Assumptions (baseline)
## ---------------------------------------------------------
## - Human population Nh fixed over 1 year (no births/deaths/migration).
## - Mosquito population Nm maintained by replacement births = μ_m * Nm.
## - Humans: SIS (no lasting immunity).
## - Mosquitoes: SI (no recovery; infectious until death).
## - Frequency-dependent transmission via prevalence:
##       λ_h = α * (Im / Nh)
##       λ_m = β * (Ih / Nh)
## - Homogeneous mixing; deterministic ODE model.

## 1C. Parameters (baseline)
## ---------------------------------------------------------
Nh     <- 15000                    # total humans
k      <- 8                        # mosquitoes per human
Nm     <- Nh * k                   # total mosquitoes
Ih0    <- 25                       # initial infectious humans
Im0    <- 0.12 * Nm                # initial infectious mosquitoes (12%)
Sh0    <- Nh - Ih0                 # initial susceptible humans
Sm0    <- Nm - Im0                 # initial susceptible mosquitoes

biting  <- 0.3                     # bites per mosquito per day
seeking <- 0.1                     # host-seeking success rate
P_mh    <- 0.45                    # prob. human infected per infectious bite
alpha   <- biting * seeking * P_mh # mosquito → human transmission coefficient
beta    <- 0.5                     # human → mosquito transmission coefficient
gamma   <- 1/25                    # human recovery rate (mean inf. period 25 d)
mu_m    <- 1/10                    # mosquito death rate (mean lifespan 10 d)

params_base <- list(alpha = alpha, beta = beta, gamma = gamma, mu_m = mu_m, Nh = Nh, Nm = Nm)

state0_base <- c(
  Sh = Sh0,
  Ih = Ih0,
  Sm = Sm0,
  Im = Im0,
  CH = 0,   # cumulative human infections
  CM = 0    # cumulative mosquito infections
)

## 1D. Baseline ODE system (SIS–SI)
## ---------------------------------------------------------
sis_si <- function(t, state, p) {
  with(as.list(c(state, p)), {
    Nh_now <- Sh + Ih
    Nm_now <- Sm + Im
    
    # Forces of infection (FOI)
    lambda_h <- alpha * (Im / Nh_now)   # mosquito → human
    lambda_m <- beta  * (Ih / Nh_now)   # human → mosquito
    
    # Humans (SIS)
    dSh <- -lambda_h * Sh + gamma * Ih
    dIh <-  lambda_h * Sh - gamma * Ih
    
    # Mosquitoes (SI with replacement births)
    births <- mu_m * Nm_now
    dSm <- births - lambda_m * Sm - mu_m * Sm
    dIm <- lambda_m * Sm - mu_m * Im
    
    # Cumulative incidence
    dCH <- lambda_h * Sh
    dCM <- lambda_m * Sm
    
    list(c(dSh, dIh, dSm, dIm, dCH, dCM))
  })
}

## =========================================================
## TASK 2: Simulate baseline dynamics over 365 days
## =========================================================
times <- seq(0, 365, by = 1)  # daily time grid for 1 year

out_base <- as.data.frame(ode(y = state0_base, times = times, func = sis_si, parms = params_base))

## 2a. Prevalence at day 365 (humans)
prev_end_base <- tail(out_base$Ih, 1) / Nh
cat(sprintf("\n[Task 2.a] Baseline human prevalence at day 365: %.6f (%.4f%%)\n",
            prev_end_base, 100 * prev_end_base))

## 2b. Interpretation
if (prev_end_base < 0.01) {
  message("[Task 2.b] Interpretation: Very low endemicity at 1 year under baseline (<1%).")
} else if (prev_end_base < 0.10) {
  message("[Task 2.b] Interpretation: Low but persistent transmission at 1 year (1%–10%).")
} else {
  message("[Task 2.b] Interpretation: Substantial ongoing transmission at 1 year (≥10%).")
}

## Baseline plot 
p_base <- subplot(
  # Humans panel
  plot_ly(out_base, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h", line = list(color = "blue")) %>%
    add_lines(y = ~Ih, name = "I_h", line = list(color = "red")),
  # Mosquitoes panel
  plot_ly(out_base, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m", line = list(color = "green")) %>%
    add_lines(y = ~Im, name = "I_m", line = list(color = "orange")),
  nrows = 1, shareX = TRUE
) %>%
  layout(
    title = list(text = "<b>Baseline Dynamics</b>", x = 0.5),
    annotations = list(
      list(x = 0.20, y = 0.95, text = "Humans", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 14, color = "black")),
      list(x = 0.80, y = 0.95, text = "Mosquitoes", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 14, color = "black"))
    )
  )
p_base

## =========================================================
## TASK 3: Interventions (continuous vaccination in Scenarios 2 & 3)
## =========================================================

## -------------------------
## Scenario 1: Vector control (reduce host-seeking by 60%)
## -------------------------
## Flow diagram change:
## - Same compartments as baseline; only parameter α decreases
##   seeking_vec = seeking * (1 - 0.60), α_vec = biting * seeking_vec * P_mh
## Assumptions:
## - Acts by reducing mosquito–human contact (no direct change in mosquito lifespan here).
## Budget:
## - 50,000 USD per year (used in cost-effectiveness analysis).

seeking_vec <- seeking * (1 - 0.60)
alpha_vec   <- biting * seeking_vec * P_mh
params_vec  <- modifyList(params_base, list(alpha = alpha_vec))

state0_vec <- state0_base
out_vec <- as.data.frame(ode(y = state0_vec, times = times, func = sis_si, parms = params_vec))

p_vec <- subplot(
  plot_ly(out_vec, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h", line = list(color = "blue")) %>%
    add_lines(y = ~Ih, name = "I_h", line = list(color = "red")),
  plot_ly(out_vec, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m", line = list(color = "green")) %>%
    add_lines(y = ~Im, name = "I_m", line = list(color = "orange")),
  nrows = 1, shareX = TRUE
) %>%
  layout(
    title = list(text = "<b>Scenario 1: Vector Control (60% ↓ host-seeking)</b>", x = 0.5),
    annotations = list(
      list(x = 0.20, y = 0.95, text = "Humans", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 14, color = "black")),
      list(x = 0.80, y = 0.95, text = "Mosquitoes", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 14, color = "black"))
    )
  )
p_vec

## -------------------------
## Scenario 2: Vaccination (continuous, daily)
## -------------------------
## Vaccination is a *continuous flow* every day.
## We model protection with efficacy VE:
##  - A per-capita vaccination *attempt* rate ν is applied to susceptibles (S_h).
##  - Of those, a fraction VE become fully protected and move to V_h.
##  - The remaining (1-VE) do not gain protection and effectively stay in S_h.
##  - Immunity wanes at rate ω: V_h --ω--> S_h (mean 30 days).
##
## Choosing ν to target ~30% protected at any time (ignoring infection):
##  In the simplified S↔V system without infection, the protected fraction at equilibrium is
##      f = (VE * ν) / (VE * ν + ω)
##  We solve for ν given target f_target = 0.30:
##      VE*ν / (VE*ν + ω) = 0.30  ⇒  VE*ν = 0.30*(VE*ν + ω)
##      VE*ν = 0.30*VE*ν + 0.30*ω  ⇒  0.70*VE*ν = 0.30*ω
##      ν = (0.30/0.70) * (ω / VE)  ≈ 0.428571 * (ω / VE)
VE    <- 0.65
omega <- 1/30
nu    <- (0.30/0.70) * (omega / VE)   # continuous vaccination rate calibrated to ~30% protected

## ODEs (humans add V_h; mosquitoes unchanged):
##  Let λ_h = α * (Im / Nh), λ_m = β * (Ih / Nh), Nh = Sh + Ih + Vh
##  dSh = -λ_h*Sh + γ*Ih + ω*Vh - (VE*nu)*Sh
##  dIh =  λ_h*Sh - γ*Ih              # only Sh can be infected (Vh fully protected)
##  dVh = (VE*nu)*Sh - ω*Vh
##  dSm = births - λ_m*Sm - μ_m*Sm;  dIm = λ_m*Sm - μ_m*Im; births = μ_m*(Sm+Im)
##  Trackers:
##    dCH = λ_h*Sh
##    dCM = λ_m*Sm
##    dCV = (nu)*Sh        # cumulative vaccination *attempts* (doses administered), for costing

sis_si_vaccination_cont <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh_now <- Sh + Ih + Vh
    Nm_now <- Sm + Im
    
    lambda_h <- alpha * (Im / Nh_now)  # mosquito → human
    lambda_m <- beta  * (Ih / Nh_now)  # human → mosquito
    
    # Humans with continuous vaccination
    inflow_V <- VE * nu * Sh            # protected doses move S→V
    dSh <- -lambda_h * Sh + gamma * Ih + omega * Vh - inflow_V
    dIh <-  lambda_h * Sh - gamma * Ih
    dVh <-  inflow_V - omega * Vh
    
    # Mosquitoes (SI with replacement births)
    births <- mu_m * Nm_now
    dSm <- births - lambda_m * Sm - mu_m * Sm
    dIm <- lambda_m * Sm - mu_m * Im
    
    # Cumulative metrics
    dCH <- lambda_h * Sh
    dCM <- lambda_m * Sm
    dCV <- nu * Sh         # all vaccination attempts/doses (includes effective + ineffective)
    
    list(c(dSh, dIh, dVh, dSm, dIm, dCH, dCM, dCV))
  })
}

# Initial state (no protected at t=0; continuous vaccination starts immediately)
state0_vac <- c(
  Sh = Sh0,
  Ih = Ih0,
  Vh = 0,
  Sm = Sm0,
  Im = Im0,
  CH = 0,
  CM = 0,
  CV = 0   # cumulative vaccination doses administered
)

params_vac <- params_base
params_vac$nu <- nu
params_vac$VE <- VE
params_vac$omega <- omega

out_vac <- as.data.frame(ode(y = state0_vac, times = times,
                             func = sis_si_vaccination_cont, parms = params_vac))

# Plot (include V_h in humans panel, purple)
p_vac <- subplot(
  plot_ly(out_vac, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h", line = list(color = "blue")) %>%
    add_lines(y = ~Ih, name = "I_h", line = list(color = "red")) %>%
    add_lines(y = ~Vh, name = "V_h", line = list(color = "purple")),
  plot_ly(out_vac, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m", line = list(color = "green")) %>%
    add_lines(y = ~Im, name = "I_m", line = list(color = "orange")),
  nrows = 1, shareX = TRUE
) %>%
  layout(
    title = list(text = "<b>Scenario 2: Cont. Vac. (VE=0.65, waning=30d, ~30% protected)</b>", x = 0.5),
    annotations = list(
      list(x = 0.20, y = 0.95, text = "Humans", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 14, color = "black")),
      list(x = 0.80, y = 0.95, text = "Mosquitoes", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 14, color = "black"))
    )
  )
p_vac

## -------------------------
## Scenario 3: Combined (Vector control + Continuous vaccination)
## -------------------------
## Changes = Scenario 2 human system + α reduced as in Scenario 1
params_comb <- params_vac
params_comb$alpha <- alpha_vec

state0_comb <- state0_vac
out_comb <- as.data.frame(ode(y = state0_comb, times = times,
                              func = sis_si_vaccination_cont, parms = params_comb))

p_comb <- subplot(
  plot_ly(out_comb, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h", line = list(color = "blue")) %>%
    add_lines(y = ~Ih, name = "I_h", line = list(color = "red")) %>%
    add_lines(y = ~Vh, name = "V_h", line = list(color = "purple")),
  plot_ly(out_comb, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m", line = list(color = "green")) %>%
    add_lines(y = ~Im, name = "I_m", line = list(color = "orange")),
  nrows = 1, shareX = TRUE
) %>%
  layout(
    title = list(text = "<b>Scenario 3: Combined (Vector Control + Cont. Vac.)</b>", x = 0.5),
    annotations = list(
      list(x = 0.20, y = 0.95, text = "Humans", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 14, color = "black")),
      list(x = 0.80, y = 0.95, text = "Mosquitoes", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 14, color = "black"))
    )
  )
p_comb

## =========================================================
## TASK 4: Parameter summary (including vaccination/combined)
## =========================================================
params_table <- tibble(
  Parameter = c(
    "Nh","k","Nm","Ih0","Im0","Sh0","Sm0",
    "biting","seeking","P_mh","alpha (baseline)","beta","gamma","mu_m",
    "Vector control reduction","alpha (vector control)",
    "VE (vaccine efficacy)","omega (waning)","nu (vaccination rate)"
  ),
  Value = c(
    Nh,k,Nm,Ih0,Im0,Sh0,Sm0,
    biting,seeking,P_mh,alpha,beta,gamma,mu_m,
    0.60,alpha_vec,
    VE,omega,nu
  )
)
cat("\n[Task 4] Model parameters summary:\n")
print(params_table)

## =========================================================
## TASK 5: Equations (summary in comments)
## =========================================================
## Baseline (SIS–SI):
##  dSh = -α*(Im/Nh)*Sh + γ*Ih
##  dIh =  α*(Im/Nh)*Sh - γ*Ih
##  births = μ_m*(Sm+Im)
##  dSm = births - β*(Ih/Nh)*Sm - μ_m*Sm
##  dIm = β*(Ih/Nh)*Sm - μ_m*Im
##  dCH = α*(Im/Nh)*Sh ; dCM = β*(Ih/Nh)*Sm
##
## Continuous vaccination (VE, ν, ω; V_h protected):
##  Nh = Sh + Ih + Vh ; λ_h = α*(Im/Nh) ; λ_m = β*(Ih/Nh)
##  inflow_V = VE*ν*Sh
##  dSh = -λ_h*Sh + γ*Ih + ω*Vh - inflow_V
##  dIh =  λ_h*Sh - γ*Ih
##  dVh =  inflow_V - ω*Vh
##  (mosquito eqns as baseline)
##  dCH = λ_h*Sh ; dCM = λ_m*Sm ; dCV = ν*Sh  (doses/day)
##
cat("\n[Task 5] Equations summarised in comments above.\n")

## =========================================================
## TASK 6: Updated code — implemented above (ODE functions)
## =========================================================
cat("\n[Task 6] Updated code defined and executed for all scenarios.\n")

## =========================================================
## TASK 7: Quantify intervention impact (infections averted)
## =========================================================
## Use cumulative human infections (CH) at t=365.
CH_base_end <- tail(out_base$CH, 1)
CH_vec_end  <- tail(out_vec$CH, 1)
CH_vac_end  <- tail(out_vac$CH, 1)
CH_comb_end <- tail(out_comb$CH, 1)

cases_averted_vec  <- CH_base_end - CH_vec_end
cases_averted_vac  <- CH_base_end - CH_vac_end
cases_averted_comb <- CH_base_end - CH_comb_end

cat("\n[Task 7] Cumulative human infections (CH) at 365 days and averted cases vs baseline:\n")
cat(sprintf(" Baseline CH:        %.1f\n", CH_base_end))
cat(sprintf(" Vector control CH:  %.1f  -> Averted: %.1f\n", CH_vec_end,  cases_averted_vec))
cat(sprintf(" Vaccination CH:     %.1f  -> Averted: %.1f\n", CH_vac_end,  cases_averted_vac))
cat(sprintf(" Combined CH:        %.1f  -> Averted: %.1f\n", CH_comb_end, cases_averted_comb))

if (cases_averted_comb > max(cases_averted_vec, cases_averted_vac)) {
  message("[Task 7] Interpretation: Combined intervention averts the most infections over 1 year.")
} else if (max(cases_averted_vec, cases_averted_vac) > cases_averted_comb) {
  message("[Task 7] Interpretation: A single intervention averts more cases than combined — revisit assumptions or parameterisation.")
} else {
  message("[Task 7] Interpretation: Similar impacts across interventions under current parameters.")
}

## =========================================================
## TASK 8: Cost-effectiveness (cost per case averted)
## =========================================================
## Costs over 365 days:
##  - Vector control: fixed annual budget 50,000 USD
##  - Vaccination: 20 USD per *dose administered* (we integrate CV(t))
##  - Combined: sum of the two
cost_vector_control <- 50000
cost_per_vaccinee   <- 20

# Total doses administered = CV(t=365)
CV_vac_end  <- tail(out_vac$CV, 1)
CV_comb_end <- tail(out_comb$CV, 1)

cost_vaccination    <- CV_vac_end  * cost_per_vaccinee
cost_vaccination_cb <- CV_comb_end * cost_per_vaccinee
cost_combined       <- cost_vector_control + cost_vaccination_cb

# Helper to compute CPCA and avoid division by zero/negative
cpca <- function(cost, averted) ifelse(averted > 0, cost/averted, Inf)

ce_vec  <- cpca(cost_vector_control, cases_averted_vec)
ce_vac  <- cpca(cost_vaccination,   cases_averted_vac)
ce_comb <- cpca(cost_combined,      cases_averted_comb)

cat("\n[Task 8] Cost per case averted (USD):\n")
cat(sprintf(" Vector control:   %s (Cost = $%.0f)\n",
            ifelse(is.finite(ce_vec), sprintf('%.2f', ce_vec), "NA / not cost-effective"),
            cost_vector_control))
cat(sprintf(" Vaccination:      %s (Cost = $%.0f; Doses = %.1f)\n",
            ifelse(is.finite(ce_vac), sprintf('%.2f', ce_vac), "NA / not cost-effective"),
            cost_vaccination, CV_vac_end))
cat(sprintf(" Combined:         %s (Cost = $%.0f; Doses = %.1f + Vector control $%.0f)\n",
            ifelse(is.finite(ce_comb), sprintf('%.2f', ce_comb), "NA / not cost-effective"),
            cost_combined, CV_comb_end, cost_vector_control))

# CE table (sorted by CPCA)
ce_tbl <- tibble(
  Scenario = c("Vector control", "Vaccination", "Combined"),
  Cost_USD = c(cost_vector_control, cost_vaccination, cost_combined),
  Cases_Averted = c(cases_averted_vec, cases_averted_vac, cases_averted_comb),
  Cost_per_Case_Averted_USD = c(ce_vec, ce_vac, ce_comb)
) %>% arrange(Cost_per_Case_Averted_USD)

cat("\n[Task 8] Cost-effectiveness table (sorted by CPCA):\n")
print(ce_tbl, n = Inf)

best <- ce_tbl$Scenario[1]
message(sprintf("[Task 8] Interpretation: '%s' provides the lowest cost per case averted over 365 days (best value for money under current assumptions).", best))

## CE plane plot (single panel with bold title)
p_ce <- plot_ly(ce_tbl, x = ~Cases_Averted, y = ~Cost_USD, text = ~Scenario,
                type = "scatter", mode = "markers+text", textposition = "top center",
                marker = list(size = 12)) %>%
  layout(
    title = list(text = "<b>Cost-Effectiveness Plane: Cost vs Infections Averted (365 days)</b>", x = 0.5),
    xaxis = list(title = "Infections Averted (vs baseline)"),
    yaxis = list(title = "Cost (USD)")
  )
p_ce

## =========================================================
## Optional: Comparison plot of cumulative incidence CH
## =========================================================
comp_df <- tibble(
  time = times,
  Baseline = out_base$CH,
  `Vector control` = out_vec$CH,
  `Vaccination (continuous)` = out_vac$CH,
  `Combined (continuous vac + VC)` = out_comb$CH
) %>%
  pivot_longer(-time, names_to = "Scenario", values_to = "CH")

p_comp <- plot_ly(comp_df, x = ~time, y = ~CH, color = ~Scenario) %>%
  add_lines() %>%
  layout(title = list(text = "<b>Cumulative Human Incidence (CH) over 365 days</b>", x = 0.5),
         xaxis = list(title = "Time (days)"),
         yaxis = list(title = "Cumulative human infections"))
p_comp



