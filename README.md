# COVID-19-Outbreak-Size
These codes outline a simple Bayesian model designed to estimate the outbreak size during the exponential growth phase of the COVID-19 epidemic from one or two surveillance streams providing counts of cases meeting various criteria. We illustrate it through scenarios based upon virological surveillance from a network of influenza-like illness consultations in primary care and from pneumonia cases in hospitals, but the approach generalizes to other surveillance streams such as mortalities. 

functions_pneumonia.r:                     functions for pneumonias surveillance system;
functions_pneumonia_ILI.r:                 functions for a smattering of pneumonias and influenza-like illnesses surveillance system;
functions_previous_not_tested_pneumonia.r: functions for pneumonias surveillance system without testing pneumonias for previous days;
mcmc.r:                                    MCMC code to estimate the outbreak size;
parameters.csv:                            parameters involved in simulation;
surveillance.csv:                          fake data for the simulation.

