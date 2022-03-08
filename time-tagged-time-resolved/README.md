# time-tagged-time-resolved

Directory for processing experimental time-tagged-time-resolved photon emission statistics data.

`crosscorrelation_partitions_example.m` analyzes and plots example hBN time-tagged-time-resolved photon arrival data including a demonstration of partitioning data between different blinking states.

`Data/PicoHarp/1iPicoData_072121/` is a directory containing the example data in the form of the raw time-tagged data in the `.ptu` format that is output by the time-correlated single-photon counter.
`Data/1iPicoData_072121.mat` is a file containing the metadata that points to the location and names of the `.ptu` files.

`tttr-functions` contains functions for processing tttr data. ` `tttr-functions/TTTR_cross_correlation.m` is the main processing function and the rest are supporting functions.

`Processed data` contains example data output from `crosscorrelation_partitions_example.m`



