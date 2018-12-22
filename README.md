Evolution of MHC repertoire and number of copies:
==============================

This work is part of the project [**Evolution of the number of copies of MHC genes: testing the optimality hypothesis and exploring alternatives.**](https://sites.google.com/site/evobiolab/projects) supported by the Polish Science Centre (NCN).

This branch **MhcEvo2000** has the ability to be run parellel on multi-core CPUs using openMP with *#pragma omp parallel* preprocessor calls. **It aims at simulating the sexual selection component in MHC gene copy number variation**.

This is a new development - the program used for research and papers is in the branch **InDeapthInfection**.

Documentation:
-----------

Go from here to *./Doc/html/* and launch the *./Doc/html/index.html* file in your web browser.
The documentation was created with Doxygen.

Authors:
--------
Coded by **Piotr Bentkowski** - bentkowski.piotr@gmail.com  
in collaboration with **Jacek Radwan** - jradwan@amu.edu.pl

Evolutionary Biology Group  
Faculty of Biology  
Adam Mickiewicz University  
ul. Umultowska 89  
61-614 Poznań  
Poland  

www: https://sites.google.com/site/evobiolab/home  
PI : **Jacek Radwan** - jradwan@amu.edu.pl

How to compile:
-----------
The program was written in [C++14 standard](https://en.wikipedia.org/wiki/C%2B%2B14) so if you are using GCC, then version gcc 4.8 seems to be the minimum requirement (I had 5.4 and 6.2). This program has some serious dependencies on [C++ Boost Libraries](http://www.boost.org/). Should compile smoothly on most modern GNU/Linux distros with Boost Libs installed. Having [Scons build tool](http://www.scons.org/) might be useful too. Basic compilation works fine on Ubuntu 16.04 LTS with mentioned packages installed by running the command:
```bash
g++ -O3 -o MHC_model main.cpp Gene.cpp Antigen.cpp Host.cpp Pathogen.cpp H2Pinteraction.cpp Random.cpp Tagging_system.cpp Environment.cpp DataHandler.cpp -fopenmp -std=c++14
```

The code here can be also used as a toolbox for your research. You can stitch your own *main_yourown.cpp* file with your scenario and tailored procedures (we did so for our research) and compile it using Scons script which is part of this code bundle. To do so run:
```Bash
scons -Q scenario="Scenarios/main_default.cpp"
```
replacing *Scenarios/main_default.cpp* with the path to your own custom *main_xyz.cpp*. Use our *main_core.cpp* file as the base to develop your scenarios of MHC evolution dynamics. WARNING! Using our Scons script overwrites the *main.cpp* file in the main source directory - always store your scenarios somewhere else than the *main.cpp* file!

Sometimes your HPC Cluster is lame and old and it has fairly outdated compiler (e.g. gcc < 4.8). Then you can statically link the libraries on your fancy brand new PC running the latest Linux distro and send the no-dependencies executable to cluster. Compile like this:
```bash
g++ -static -O3 -o MHC_model main.cpp Gene.cpp Antigen.cpp Host.cpp Pathogen.cpp H2Pinteraction.cpp Random.cpp Tagging_system.cpp Environment.cpp DataHandler.cpp -fopenmp -std=c++14
```

Or run the Scons script:
```shell
scons -Q scenario="Scenarios/main_default.cpp" linking="static"
```
But it may throw some warnings… e.g. */usr/lib/gcc/x86_64-linux-gnu/5/libgomp.a(target.o): In function `gomp_target_init`: (.text+0xba): warning: Using 'dlopen' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking*.

How to run:
-----------

The program takes exactly 17 parameters. These are:

*  **00** - Program's name
*  **01** - Number of threads program will try to use on a multi-core CPU. Giving 0 will make the program use all CPU cores available.
*  **02** - Number of bits in a gene.
*  **03** - Number of bits in an antigen.
*  **04** - Host population size.
*  **05** - Pathogen population size.
*  **06** - Number of pathogen species.
*  **07** - Number of genes in one host chromosome (they have two chromosomes) when the model is being initialised.
*  **08** - Number of pathogen generations per one host generation.
*  **09** - Number of host generations (effective length of model run).
*  **10** - Probability of mutation in hosts ([0,1] range).
*  **11** - Probability of mutation in pathogens ([0,1] range).
*  **12** - The heterozygote advantage / lack of advantage mode. It has to be 10 for heterozygote advantage or 11 for lack of thereof.
*  **13** - Probability of deleting a gene in the host ([0,1] range).
*  **14** - Probability of duplicating a gene in the host ([0,1] range).
*  **15** - Maximal number of genes permitted in one host chromosome.
*  **16** - Number of sexual partners an individual checks out before selecting one for mating.
*  **17** - Alpha factor for the host fitness function ([0,1] range).
 
 Historically we also had one more argument, that can be added if you want:
 - 18 - Fraction of the antigen's bits which are forbidden from changing (a.k.a. *No Mutation Bits*).

This program can recognise simple errors in the argument list (a probability value out of [0,1] range, negative
values when only positive are allowed etc.), but will not recognise when they don't make a 'biological' sense.

The output and data visualisation:
-----------

Program produces a number of text files containing desired data. The file *InputParameters.json* contains the values of parameters fed to the model. The file *HostsGeneDivers.csv* contains the main statistics of the model, including the time stamp used in the visualisation process of the other files with data. This two files are the basics to analyse the simulation output. Other files are optional and they can be switched off by commenting out some lines in the main function (*main.cpp* file) what should speed up the program runtime.

All the files are:

*  ***InputParameters.json*** - contains run's parametrisation in JSON text format 
*  ***HostsGeneDivers.csv*** - contain basic statistics of the output for host genomes 
*  ***HostGenomesFile.[t].csv*** - contains all the genomes of all the cells for host population in time *t*. There can be more then one file like this for different time snapshots. 
*  ***PathoGenomesFile.[t].csv*** - contains all the genomes of all the cells for pathogen population in time *t*. There can be more then one file like this for different time snapshots. 
*  ***HostGeneNumbTotal.csv*** - total number of genes in each individual host in each time step (linked to *HostMHCsNumbUniq.csv*, that each column in *HostGeneNumbTotal.csv* is represents the same cells as in *HostMHCsNumbUniq.csv*). 
*  ***HostMHCsNumbUniq.csv*** - number of unique MHC alleles in each individual host in each time step (linked to*HostGeneNumbTotal.csv*, that each column in *HostGeneNumbTotal.csv* is represents the same cells as in *HostMHCsNumbUniq.csv*). 
*  ***NoMutationInPathoList.csv*** - list of conserved antigen sites. Each line contains one pathogen species, each number indicates the index of a "no-mutation" site in this species antigen. When empty it means there is no mutation restrictions. 
*  ***NumberOfMhcAfterMating.csv*** -number of the unique MHC types in each individual host in each time step after the mating procedure took place.
*  ***NumberOfMhcBeforeMating.csv*** - number of the unique MHC types in each individual host in each time step before the mating procedure took place.
*  ***NumberOfMhcInMother.csv*** - number of the unique MHC types in each individual host that is selecting a partner (a.k.a. "mother") in each time step during mating procedure. Each individual has a corresponding partner at the same index in the file *NumberOfMhcInFather.csv*. 
*  ***NumberOfMhcInFather.csv*** - number of the unique MHC types in each individual host that has been selected as a mating (a.k.a. "father") in each time step during mating procedure. Each individual has a corresponding partner at the same index in the file *NumberOfMhcInMother.csv*. 
*  ***PresentedPathogenNumbers.csv*** - number of presented pathogens by each individual in each time step. 


Visualisation is done using Python 3.6 scripts containing a a lot of calls to Numpy, Matplotlib and other scietific Python libraries. You may wish to consider using the [Python Anaconda](https://www.anaconda.com/download/) for your Pythonic endeavours. Visualisation and stats scripts can be found in *PyScripts* directory. 

To load the environment for output analysis run:
```shell
conda env create -f mhcEvoCondaEnv.yaml
```
file `mhcEvoCondaEnv.yaml` can be found in *PyScripts* directory.
