# parsimonious_model
parsimonious model by M. Gram, D.J. Miller, 2006:
"Unsupervised Learning of Parsimonious Mixtures on Large Spaces With Integrated Feature and Component Selection"

As an example, first run generate_random_num to synthetically generate a test-case. Then set M_max to a reasonable starting point (>=5).


Run [MODEL,METRICS] = parsimonious(Data_whole, M_max);


It works in Octave, and it needs statistics and ndpar packages.
