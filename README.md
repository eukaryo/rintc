# rintc
  
$ make

to compile the source and obtain the executable "rintc". I succeed to compile them in the following environment:

gcc/7.2.0 on CentOS Linux release 7.5.1804 .

then, 

$ ./rintc \<rna\> \<str\> \<W\> \<algo\>

where

rna: RNA sequence in upper case, e.g. CCCCAAAAGGGG

str: RNA structure enclosed in quotes, e.g. "((((....))))"

W: maximum-span constraint, e.g. 12

algo: name of the algorithm, e.g. RintCwithDFT



For other algos and replicating the experiments, see main.cpp. 
      
example:

./rintc CCCCAAAAGGGG "((((....))))" 12 RintCwithDFT

# How to reproduce
For reproduce the experiments perfomed in the paper, use the commands below:

$ ./rintc \<ComputationTimeExperimentW100\>

$ ./rintc \<ComputationTimeExperimentWN\>

$ ./rintc \<HeatResistanceExperiment\>

$ ./rintc \<AccuracyExperiment151\>

$ ./rintc \<RrnaAccuracyExperiment\>
