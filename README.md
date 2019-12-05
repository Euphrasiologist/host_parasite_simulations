Simulating host parasite dynamics
================
Max R Brown
05 December, 2019

# Simulating host parasite dynamics

I am completely new to modelling host-parasite dynamics. I just thought
it would be interesting to do, as I have been investigating this with
emprical data. The object here was to simulate an interaction between a
host and a parasite, to look at dynamics primarily from the point of
view of the parasite. The object was also to sort of ‘tailor’ this for
plants. Plant parasites are interesting because they cannot move;
perhaps that makes the dynamics easier to code\! I am frightfully bad at
this, and I am really interested in taking this somewhere. I will
annotate the code as fully as I can and hopefully someone, somewhere can
help make these models (1) work properly, (2) under a variety of
scenarios and (3) make them more efficient\! They are pretty slow.

## Installation

You can install this package by using the following commands:

``` r
# install.packages(devtools)
# devtools::install_github("Euphrasiologist/host_parasite_simulations")
library(HoPaSim)
```

## Rationale

These functions are designed to model a simple host parasite interaction
with a number of assumptions. A host grid is generated, which is fixed
for all future generations (i.e. they are perennial hosts). A parasite
grid is then overlaid, which contains ephemeral annual plant parasites
which persist only for a single generation. Depending on the benefit a
host gifts the parasite, this will determine the number of progeny
seeded in the next generation. Explicitly, the assumptions are:

  - that hosts are fixed in position and perennial, and are unaffected
    by parasitism. This is probably not totally realistic, because in
    reality, plant parasites have the capability of being really quite
    virulent.
  - no effect on parasite of nearby hosts. There is no additional
    weighted matrix which accounts for the spatial effects of
    parasitism. In reality, we know that plant parasites can attach to
    multiple hosts at the same time. Maybe for a future implementation.
  - parasites seed randomly over the field. Seed disperal in parasitic
    plants vary. Some are dispersed by ants, others mechanically or by
    wind.
  - all parasites survive to end of generation if they land on a host,
    and selection is \> 0. In the wild, parasites may well die due to
    environmental conditions or pathologies of host defence etc.
  - parasites do not produce any offspring without a host, and all
    parasites with a host reproduce. This is related to the previous
    point, however, it is known for some species (especially facultative
    parasites) that seed can be produced even without a host plant
    species.

Given these (not unreasonable) assumptions we can create a model which
can simulate parasite population sizes over time measured in
generations. These can be for single host parasite systems, single
parasite and multiple hosts, or multiple parasite species and multiple
host species.

## Single host single parasite system

This is the most simple system to model. In the functions, the user can
determine the size of the field to play on, the selective advantage
measured in terms of offspring that survive to maturity in the next
generation if you land on the host, number of host individuals and
number of parasites. Number of generations can also be specified.

The output is extremely verbose, as each iteration of the spatial matrix
structure is saved \[maybe change this?\]. The interesting part is
probably the first slot in the output, accessed as shown below, which
saves population sizes over time

``` r
library(data.table); library(ggplot2)
sim1 <- HoPaSim::onehost_oneparasite(field.size = 25^2, s = 2, host.number = 400, parasite.number = 20, gens = 100)

sim1.2 <- setDT(sim1[[1]])

head(sim1.2)
#>    generation P.size
#> 1:          1     20
#> 2:          2     18
#> 3:          3     23
#> 4:          4     26
#> 5:          5     31
#> 6:          6     41
```

A plot can then be generated using the
data.

``` r
ggplot(sim1.2, aes(x = generation, y = P.size))+geom_line() + theme_bw() + xlab(label = "Generations") + ylab(label = "Population Size")
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

For these parameters under this model, the parasite sweeps rapidly to
parasitise every possible host in about 20 generations. It is capped at
the upper limit of 400 as this is the number of host plants I seeded the
field with.

This example can be easily replicated many times using the function
`onehost_oneparasite_dynamics()`, with an extra parameter argument,
reps, which can be used to specify the number of replicate simulations
to use.

``` r
# use same parameters as above
sim2 <- suppressWarnings(HoPaSim::onehost_oneparasite_dynamics(field.size = 25^2, s = 2, host.number = 400, parasite.number = 20, gens = 100, reps = 30))
head(sim2)
#>    generation P.size rep
#> 1:          1     20   1
#> 2:          2     12   1
#> 3:          3     15   1
#> 4:          4     19   1
#> 5:          5     26   1
#> 6:          6     31   1

ggplot(sim2, aes(x = generation, y = P.size))+geom_line(aes(group = rep)) + theme_bw() + xlab(label = "Generations") + ylab(label = "Population Size")
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

From this, it is easy to see that there is certainly variability
associated with this simulation, but they are all roughly doing a
similar thing. The driving engine for the simulations is the R function
`sample`, and this is causing the variability.
