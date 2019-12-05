# Simulating host parasite dynamics

## Installation

You can install this package by using the following commands:

```{r}
# install.packages(devtools)
devtools::install_github("Euphrasiologist/host_parasite_simulations")
library(HoPaSim)
```

I am completely new to modelling host-parasite dynamics. I just thought it would be interesting to do, as I have been investigating this with emprical data. The object here was to simulate an interaction between a host and a parasite, to look at dynamics primarily from the point of view of the parasite. The object was also to sort of 'tailor' this for plants. Plant parasites are interesting because they cannot move; perhaps that makes the dynamics easier to code! I am frightfully bad at this, and I am really interested in taking this somewhere. I will annotate the code as fully as I can and hopefully someone, somewhere can help make these models (1) work properly, (2) under a variety of scenarios and (3) make them more efficient! They are pretty slow.

This is also designed to be a package, so it should be eventually easy peasy to get these functions onto your computer.



Max Brown 03.12.19.
