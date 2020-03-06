setwd("~/nBox/r/wuhanflu/sporadic_surveillance/yuemu")
set.seed(888)
clean_up_nicely = FALSE
source('code/functions.r')

for(scen in 1:15)
{
  cat('Scenario',scen,'\n')
  surveillance = read_data(paste0('data/surveillance_scenario_',scen,'.csv'))
  pars = read_pars('data/parameters_scenarios_1_4.csv')
  if(scen>=5)pars = read_pars('data/parameters_scenarios_5_6.csv')
  if(scen>=7)pars = read_pars('data/parameters_scenarios_7_15.csv')
  mcmcoutput = mcmc(pars,surveillance,10000,100)
  summariser(mcmcoutput,paste0('output/summary_scenario_',scen,'.txt'))
}
if(clean_up_nicely)rm(mcmcoutput,pars,surveillance,casestodate,casestoday,loglikelihood,mcmc,mh,modelmean,parsepar,read_data,read_pars,summariser,clean_up_nicely,scen)


