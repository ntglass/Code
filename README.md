# Code
Examples of my coding!

Thanks so much for viewing! This document will be updated to reflect additional file uploads.

ThreshEx 1_Nick.R is an analysis of USDA collected plant coverage from 1996 to 2016, which is publicly available on the Jornada LTER website: https://lter.jornada.nmsu.edu/. The goal of the experiment was to test if grazing pressure combined with woody encroachment could explain the ecological state transitions occuring in drylands across the world: one from semiarid grasslands to shrub lands. The experiment was conducted by the USDA Agricultural Research Service and began in 1996 with intense grazing episodes occurring in 18 paddocks in Jornada del Muerto north of Las Cruces, New Mexico, USA. Approximately 80% of grass biomass was removed by cows that grazed either in the summer or the winter. 6 paddocks were not grazed. In addition, there was a shrub removal treatment: 12 paddocks had Honey Mesquite shrubs clipped and treated with herbicide, while the other 12 paddocks had shrubs left intact. The treatments ended in 2000, and USDA-ARS researchers returned to record plant cover and diversity in 2002, 2007, and 2016. The code I've provided downloads the publicly available data and transforms the plant counts into percent cover values. We then create a working dataset with the grazing and shrub treatments as categorical explanatory factors and grazing severity, post-treatment shrub cover, and post-treatment diversity as continous/integer explanatory variables. Response variables are the recovery of Black Grama, the once-dominant grass of the region, the productivity of Honey Mesquite shrubs described in Rain-Use Efficiency, and plant species richness and diversity across time. We use mixed-effects models, simpler regressions, and PERMANOVAs to examine the data multiple ways. We compare models using AIC where applicable, residuals, and rudimentary model visualizations. Once we find a good model, we extract the effect sizes and 95% confidence intervals. Finally, the code ends by visualizing the data and chosen models in some pretty figures that project our estimates of plant productivity and diversity into the future.
