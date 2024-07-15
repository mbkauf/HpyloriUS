Meeting with Fernando 8/23
1. Updates on coding
    a. SIS model working now 
        - start with no treatment
        - need antibiotic use
    b. Cleaned up code using linter in VSCode
    c. SIS ABR model coded up
    d. Updated race/ethnicity groups allowable
2. Waiting to hear back from Mercedes about antibiotic use
3. Questions:
    a. Calibrating WAIFW matrix lambda prime and lamba xi
        - get rid of lambda xi. 
        - more or less infectious
        - would have been a different multiplier
        - Manuel (expert at Stanford) and Julie Parsonette 
    b. Can you explain how we use birth cohorts?
        - pick any life table
        - calibrate to period prevalence
        - maybe use period before abr
        - use 1991 period? 
            - Ask Jorge for period, not birth cohort prevalence
            - period FOI
    c. How do we project cases?
        - bring to present time

Scenarios
- one time treatment for everyone. Assuming 100% eff. 
- then less effectiveness
- include test and treat
- probability of mutation/resitance (VOI)

Start generating validation figures
- similar to 2nd chapter
    - prevalence of resitance over time
        - by age groups
- send email to josh, jeremy, and jason about race mixing

- function aproxfun values for different time periods (step-wise or linear)
- ask Costanza for test statistics

- 

Meeting with Fernando 8/30
1. Tranmission matrix
    - Keeling and Rohani
2. Scenario analysis results
- compare FOI predicted vs empirical

- data nhanes
- predicted prevalence from catalyitic model
- model predicted prevalence

predited foi vs model predicted foi

reverse S and I dash line
- just infectious over time
- no pooling of infected
- vertical lines for treatment dates


Meeting with Fernando 9/25
- model predicted instead of fitted (foi graph)
- look for steady state (compare year 400 to NHANES)

Good morning. Helicobactor pylori or H pylori is an extremely common bacterial infection that is estimated to be prevalent in two thirds of the world’s population. It is also the strongest known risk factor for gastric cancer. While more common globally, H pylori is still prevalent in the United States, especially among racial and ethnic minorities.

Our objective is to build and validate a race-specific, age-structured, Susceptible-Infected-Susceptible dynamic transmission model of H pylori infection and treatment in the US. guide decision making related to h. pylori erradication policies using a ....

We developed an ordinary differential equation SIS model. We calibrated race and age specific transmission parameters to match model predicted force of infection to empirical estimates from other CISNET related work. 

small finding. 

If you’d like to learn more, come see me at poster number #. Thank you.


For the poster:
- Figure for relative impact of scenarios
- update text
- re-calibrate with different starting values
    - get to steady state before starting projections

Bigger picture:
- antiobiotic resistance


Model Name:

HP
HS
GC
G-CHS
HaS-GasCan-Sim

HelP-GC-HaS-Sim
Helicbactor pylori Gastric Cancer Ha

Meeting 1/22/2024
1) Seattle
2) Updates
    a) Probabilistic code
        - have code to run in parallel
        - triangular for mutation rate
            - switch to beta
        - Cholesky decomp for WAIFW
        - missing any other variables?
        - alphas uniform?
        - predicted mean and CrI
        - sensitivity and specificty of tests (correlations)
        - fit a functional form of antibiotic usage over time

    b) GPUs?
        - Run on Sherlock?

    c) Talk timeline 
        - Concerns about dissertation timeline

        look into current treatment guidelines
        cost of sensitivity test
        background 
        generate FOI over time accounting treatment
        changes in prevalence and bounds
        hexamaps - FOI over time by age

        Present to NCI

Meeting 2/5/2024
1) Talk FOI plot, why white appears to have FOI vary over time
    a) Scale of difference
    b) There was a small issue with the code when I sent, I haven't re-run the plots
2) Antibiotic use over time
    a) this is a really time consuming task, can I table for the time being?
3) IMIS
    a) How should we think about priors? Right now it's random uniform
        - Draw from NM results?
        - Random uniform?
        - Uniform from NM bounds
        - Log normals from NM results
        - exponentiate to avoid going outside the bounds
        - logit transformation
    b) Number of runs? 1000 from posterior (once it's working) Make iter 40
    - figure out how many to do overnight
    c) Run in parallel (ask Fernando for code)
4) PLot resistance over time
5) Plot step 2
6) Plot prevalence of H. pylori overtime with bounds from other studies
    - Chelsea/Mercedes has data in excel file
7) Funds for packaging for these generators
    - Plots of policies
    - Diagram of model

8) 1991, 2000, 2010, 2019

Harvard-Stanford Meeting 2/8/2024
1) My updates
    1) IMIS calibration for tranmission matrix parameters
    2) Validation targets
        - Chelsea/Mercedes - prevalence/resistance over time data?

- look for proportion of macrolides that are clarithromycin in US 
- plot highest likelihood parameters sets from IMIS
- DEoptim?


2/22/2024 Meeting
My updates:
1) Used webplot digitizer on clarithromycin rates from 1999-2012
    1) There appears to be a cyclical trend, but thought we could potentially fit curve to data
    2) Assume increase from 1993 to 1999. check 2 or 3 shapes. validate to targets
2) Hessian for RP
    1) The nearPD function results in a SEs that are far to large
    2) Calibrate log parameters and exponeniate later - testing know
    3) SIR with bounds around coefficient point estimates
    4) Save trace of GOF
    5) Plots that show convergence in appendix
3) Lower bound 1 - other alphas UB 
    1) Prior as dirchelet
    2) use density of prior to constrain optimization
    3) log likelihood and add log prior can optimize over log posterior
    4) Nick Menzies calibration paper (S1)
        - parameters 1-4 are dirchelet


My TODO list for today:
1) DONE - Calibrate log parameters
2) DONE - Save trace of GOF and parameter values
    - Also for assortative
3) Test probabilistic WAIFW with new parameters
    - If this didn't work, may need a quick SIR

applyfun()

either continue downward or constant
.2 over 5 years

3/18 Meeting
1) Calibration updates
    1) Use same v_inf for assort and RP
    2) Alphas are calibrated, but same issue that we have with RP betas, very large SEs
        1) I can try SIR or IMIS next?
    3) Show plots

treat group with high prevalence

IMIS 1000 resample, sample size of posterior in 1000, start with 10 iterations
