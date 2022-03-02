# photon-emission-correlation-spectroscopy

Supporting code for "Photon emission correlation spectroscopy as an analytical tool for quantum defects"
https://arxiv.org/abs/2111.01252

Repository contains an example of processing experimental time-tagged-time-resolved photon emission statistics data and an example of simulating photon emission statistics

`Photon Statistic Simulator/example.m` defines a 5-level toy model with spin dynamics and simulates the model's optical dynamics through `PhotonStatisticSim.m`.

`time-tagged-time-resolved/crosscorrelation_analysis_script_example.m` analyzes example hBN time-tagged-time-resolved photon arrival data through the processing script, ` TTTR_cross_correlation.m`, including a demonstration of partitioning data between different blinking states.
