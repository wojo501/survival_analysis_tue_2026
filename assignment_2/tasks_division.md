# Option A — 4 Parts, One Person Each

## Part 1: Sokhail — Practical Context + Research Question

* Define the patient population.
* Define the event: time from discharge until first unplanned readmission within 90 days.
* Explain why survival analysis is appropriate.
* State the research question.
  * *Example question:* How does misspecifying the survival distribution affect estimation of the effect of remote monitoring on time to hospital readmission?

## Part 2: Wojciech — Survival Model + Model Assumptions

* Choose Weibull as the true survival model.
* Explain why Weibull makes sense for readmission risk.
* Decide whether hazard should decrease or increase over time.
* Specify shape/scale parameters.
* Describe the basic model:
  
  $$
  h(t \mid X) = h_0(t)\exp(\beta_1X_1 + \beta_2X_2 + \beta_3X_3)
  $$

* Briefly justify the model choice.

## Part 3: Michal — Covariates

* Define at least three covariates:
  * Age
  * Comorbidity score
  * Remote-monitoring treatment
* Give their distributions.
* Give expected directions of effects.
* Give true coefficient values for the simulation.
* Explain relationships between covariates if relevant.
* *Note:* The assignment explicitly asks for at least three covariates, their definitions, expected effects, and relationships where relevant.

## Part 4: Mikolaj — Censoring + Data-Generating Mechanism

* Explain administrative censoring at 90 days.
* Explain possible independent dropout.
* Explain whether censoring is independent.
* Set the sample size, e.g., $N=1000$.
* Specify how covariates are generated.
* Specify how survival times are generated.
* Specify how censoring times are generated.
* Define the observed time:
  
  $$
  Y_i = \min(T_i, C_i, 90)
  $$

* Define the event indicator:
  
  $$
  \delta_i = I(T_i \le C_i, \ T_i \le 90)
  $$

* *Note:* This person basically writes the mathematical simulation blueprint, which Part I explicitly requires.