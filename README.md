# Thermal guild fidelity in a warming world: resolving specialization and flexibility in temperate lake fish

Welcome to our repository of all code and data associated with this otolith microchemistry work!  


## Overview

--------------------------------------------------
How consistent are the thermal niches of freshwater fishes? We used otolith δ18O and δ13C to quantify thermal habitat use and reconstruct average growth temperatures, finding two thermal guilds (cool- and cold-water) and a range of thermal specialization within the community. Across 15 lakes, brook trout maintained consistent growing temperatures that aligned with thermal tolerance ranges, but experienced warmer growing temperatures where oxythermal stress was higher. Overall, this work highlights the use of δ18O for visualizing in-situ thermal niches. 

## Objectives 
---------------------------------------------

* Quantify thermal guild variation between species, within a single lake. 

* Quantify thermal guild variation within a thermally sensitive species, across multiple lakes


## Repository Structure
--------------------------------------------------


This repository uses the following dependencies: 


Files are set up to follow the analysis workflow. You will need to generate to load `Step1_DataSetupFunctions.R` in each time you run a script. Additionally, bayesian estimates for niche area and overlap are generated with the same posterior distributions. You will need to generate and save the posterior in `Step3_ItnerspecificLML_SIBERSetup.R`, so that you can load into subsequent scripts. 


## Requirements 
------------------------------------------------------------


` 
library(tidyverse)
library(viridis)
library(SIBER)
library(factoextra)
library(lme4)
library(future.apply) ## Package for parallelizing lapply
library(progressr)
`

