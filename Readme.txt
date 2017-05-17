This repository contains all the code used for the simulations in the paper ***

To understand the code within the program it is strongly suggested to
read the paper as well as the Supplementary Info.

The files are:
community_construction_coex.py:
Contains all functions for the invariant community structure case
community_construction_repl.py:
Contains all functions for the invariant community structure case
The files community_construction* contain the same function names,
the functions with the same name do the same, one with invariant com. structure
and the other with changing community structure.

parameters of communities.py:
Generates communities, this programm was used to generate the files
"coex, com_para.p" and "repl, com_para.p"

coex, com_para.p:
Pickle file that contains 2*100'000 communities. Once with positive 
environmental change and once with negative. Communities with invariant com. str.
repl, com_para.p:
Pickle file that contains 4*2*100'000 communities. Four with positive 
environmental change and four with negative. Community structure changes,
2 for each of the cases: p=0.95, p=0.50, p is random, p = 0.05

help_functions.py:
Contains helpfunctions, mainly the plot command used by most plot_fig*  files

plot_fig*:
Generates the figure(s) indicated. Correspond to the "Figure *,pdf" files

Examples.txt:
Contains a list of simple examples. Complicated examples can be found in the
plot_fig* files