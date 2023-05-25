# Mars Atmosphere History
## Running the code
To run the model from the command line just type:

`python3 main.py`

Or if you want to run in interactive mode, type
`ipython3`
to enter ipython3 and from ipython3 you can then run the model by typing:

`run -i main.py `

Note that you only need to type `ipython3` once to enter interactive python, but you may wish to run the model multiple times using `run -i main.py`

## Overview
This Python model is based on the paper by Kurokawa et al. (2018, [here](http://dx.doi.org/10.1016/j.icarus.2017.08.020)). 

The model starts 4.5 billion years ago with an initial atmospheric pressure (usually 1.5 bar), and solves for the evolution of the atmosphere up until the current day. The initial isotope abundance is assumed to be volcanic. Evolution over time is calculated by dividing the time period into smaller time-steps and incrementing variables over these time-steps. Typically, the time-steps are around 1 million years. 

Over a time-step the number of impacts are calculated using a crater chronology model. This number is scaled so that, over the whole 4.5 billion years, the total impactor mass equals some predefined value (normally 2 $\times$ 10<sup>21</sup> kg).

Each impact has a size distribution and velocity distribution individually sampled from separate size and velocity distributions using Monte Carlo sampling. For each impact, during a time interval, we calculate the loss of the atmosphere due to atmospheric eroson. We also assume a volatile fraction of the impactor and add this to the atmosphere. Depending on the impactor type (asteroid or comet) we add elemental trace abundances to the atmosphere. Not only do we add the elemental abundances, but also the isotope composition for the different sources (e.g. the ratios of isotopes of the same element).

Sputtering and photochemical escape for each element in the atmosphere is calculated using a look-up table for CO<sub>2</sub>, with other elements scaled by this value. The isotope composition is also affected by sputtering (lighter isotopes are more likely to escape). 

Interplanetary dust particles add a continuous flux of elements over time. They also alter the isotope composition (with sources based on observations). 

Volcanic degassing is another continuous source, which changes over time. It adds H<sub>2</sub>O, CO<sub>2</sub>, and N<sub>2</sub> to the atmosphere. However, it also affects the isotope composition. 

The scale height of the atmosphere affects atmospheric erosion and sputtering. It is calculated throughout the simulation based on the gravitational field of Mars and the atmospheric mass. 

By default the atmosphere is assumed to collapse when the pressure drops below 0.5 bar. Simulations by Forget et al. suggest that this threshold may change depending on the obliquity of Mars. 

## What does the model calculate?
It calculates the evolution of atmospheric pressure, trace element composition, and isotope ratio. In simulations where the atmosphere collapses it reproduces the current day composition and isotope ratio. 
