setwd("~/nBox/Work Folder/Ncovs/Estimating the size of an outbreak of COVID/")
set.seed(888)
clean_up_nicely = TRUE
source('functions_pneumonia_ILI.r')

surveillance = read_data('surveillance.csv')
pars = read_pars('parameters.csv')
mcmcoutput = mcmc(pars,surveillance,10000,100)
summariser(mcmcoutput,'summary.txt')

if(clean_up_nicely)rm(mcmcoutput,pars,surveillance,casestodate,casestoday,
                      loglikelihood,mcmc,mh,modelmean,parsepar,read_data,read_pars,summariser,clean_up_nicely)


