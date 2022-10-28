# AgeingBrain
Files to build brain model during learning

These are the final files used to run code for "Ageing Brains" project

bipartite_graph2.m -> this file creates the bipartite graph and distributes edge weights with normal and uniform distribution. It generates figures and saves .csv files to be used in subsequent script running

RouletteSelection.m -> function that will be used in Ip_generation_threshold3.m

Ip_generation_threshold3.m -> file generates input pattern and finds threshold for them. It then simulates learning process in old and young brains, creating figures at every step

input_pattern_analysis.m -> file performs input pattern similarity as well as output pattern similarity for similar input patterns

output_pattern_similarity_KURF.m -> file performs output pattern similarity for young and old learning and similarity for identical output patterns
