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



def <- defData(varname = "xbase", dist = "normal", formula = 20, variance = 3)
def <- defData(def, varname = "nCount", dist = "noZeroPoisson", formula = 6)
def <- defData(def, varname = "mInterval", dist = "gamma", formula = 30, variance = 0.01)
def <- defData(def, varname = "vInterval", dist = "nonrandom", formula = 0.07)
dt <- genData(200, def)
dt[id %in% c(8, 121)]

# https://cran.r-project.org/web/packages/simstudy/vignettes/longitudinal.html
# https://www.google.com/search?q=simulating+longitudinal+data+in+r&rlz=1C1GCEB_enZA917ZA917&oq=simulating+longitudianl+data&aqs=chrome.1.69i57j0i13i512j0i22i30i625j0i22i30j0i390l4.14006j0j7&sourceid=chrome&ie=UTF-8
