func_patient_char_gen <- function(age, # mean sd
                                  sex, #gender
                                  cd4_counts, #mean sd
                                  viral_load, #mean sd
                                  mortality_rate, # shape scale
                                  visits
                                  # lost_to_follow_up_rate
                                  ){
  
  # step 1: Create alist to fill
  temp_list = list(
    bl_age = 0,
    sex = 0,
    bl_cd4_count = 0,
    bl_vl = 0,
    trt_outcome = 0,
    years_of_follow_up = 0#,
       
    )
  # browser()
    # Step 2: generate results
    temp_list[['bl_age']] <- rlnorm(1, mean = age[1], sd = age[2])
    temp_list[['sex']] <- sample(c('F', 'M'), 1, prob = c(sex, abs(1- sex)))
    temp_list[['bl_cd4_count']] <- rlnorm(1, mean = cd4_counts[1], sd = cd4_counts[2])
    temp_list[['bl_vl']] <- rlnorm(1, mean = viral_load[1], sd = viral_load[2])
    # temp_list[['prob_mortality']] <- rweibull(1, shape = mortality_rate[1], scale = mortality_rate[1])
    prob_mortality <- rweibull(1, shape = mortality_rate[1], scale = mortality_rate[1])
    temp_list[['trt_outcome']] <- sample(c('D', 'A'), 1, prob = c(prob_mortality, abs(1-(prob_mortality))))
    temp_list[['years_of_follow_up']] <- rnorm(1, mean = visits[1], sd = visits[2])
  print(temp_list)
  return(temp_list)
}
x <- list()
for (i in 1:10) {
  x[[i]] <- func_patient_char_gen(age = c(3.1, .05), # meanlog sdlog
                                  sex = 0.52, #gender
                                  cd4_counts = c(4.2, 0.6), #meanlog sdlog
                                  viral_load = c(4.75, 2.16), #meanlog sdlog
                                  mortality_rate = c(1.2, 1) ,# shape scale
                        visits = c(5.8, 1.42)
                                  # lost_to_follow_up_rate
)
  # x=1
  # print(temp_list)
}

library(simstudy)
pt_dataset <- readRDS('data/africos_cohort.rds') # data edited from the AFRICOS_sample selection.R

# set.seed(11)
m <- .31; f <- .69 # Yapa et al. 2022
min_base_vl <- 5e3; max_base_vl <- 1e6
mean_cd4 <- 389.5 # Yapa et al. 2022
#Set preliminary values
l     <- 381.8;
u     <- 389.5;
n     <- 397.1;
alpha <- 0.05;
#Compute sample mean and SD
crit <- qt(alpha/2, df = n-1, lower.tail = FALSE);
mean_cd4 <- (l+u)/2;
SD_cd4   <- (u-l)*sqrt(n)/(2*crit)
vl_model1 <- nlme::lme(
  fixed = logvl ~ bs(time_t, 3),
  random = ~ 1 | id,
  data = pt_dataset %>%
    filter(!is.na(`VL Copies/mL`)) %>%
    mutate(logvl = log(`VL Copies/mL`, 10),
           time_t = `Duration of started ART (years)`,
           id = `SUBJECT ID (CHAR)`),
  na.action = na.exclude, control = lmeControl(opt = "optim")
)
b0 <- coef(summary(vl_model1))[[1]]
b1 <- coef(summary(vl_model1))[[2]]
b2 <- coef(summary(vl_model1))[[3]]
b3 <- coef(summary(vl_model1))[[4]]
m_var <- as.numeric(VarCorr(vl_model1)[[2]])

def <- defData(varname = "base_age", dist = "normal", formula = 26, variance = 7)
def <- defData(def,varname = "sex", dist = "categorical", formula = 'm;f')
def <- defData(def, varname = "base_CD4", dist = "normal", formula = mean_cd4, variance = SD_cd4)
def <- defData(def, varname = "base_vl", dist = "uniform", formula = 'min_base_vl;max_base_vl')
def <- defData(def, varname = "nCount", dist = "noZeroPoisson", formula = 6) # nCount defines the number of measurements for an individual
def <- defData(def, varname = "mInterval", dist = "gamma", formula = 3, variance = 1e-4) # mInterval specifies the average time between intervals for a subject (in quarters)
def <- defData(def, varname = "vInterval", dist = "nonrandom", formula = 1e-3) # vInterval specifies the variance of those interval times. If vInterval is set to 0 or is not defined, the interval for a subject is determined entirely by the mean interval.
set.seed(11)
dt <- genData(200, def)
# dt[id %in% c(8, 121)]
# The resulting longitudinal data for these two subjects can be inspected after 
# a call to addPeriods. Notice that no parameters need to be set since all information 
# resides in the data set itself:
dtPeriod <- addPeriods(dt) %>%
  mutate(time = time/4)
# dtPeriod[id %in% c(8, 121)]
# If a time-sensitive measurement is added to the data set …
def2 <- defDataAdd(varname = "vl", dist = "normal", formula = '..b0 + ..b1 * time + ..b2 * time + ..b3 * time', variance = m_var)
dtPeriod <- addColumns(def2, dtPeriod)


# https://cran.r-project.org/web/packages/simstudy/vignettes/longitudinal.html
# https://www.google.com/search?q=simulating+longitudinal+data+in+r&rlz=1C1GCEB_enZA917ZA917&oq=simulating+longitudianl+data&aqs=chrome.1.69i57j0i13i512j0i22i30i625j0i22i30j0i390l4.14006j0j7&sourceid=chrome&ie=UTF-8
