/**
 * @mainpage "Modelling MHC evolutionary dynamics"
 *
 * @authors Coded by <b>Piotr Bentkowski</b>: <a href="mailto:bentkowski.piotr@gmail.com
 * ?subject=MHC evolutionary dynamics - model code v1.0">bentkowski.piotr@gmail.com</a> \n
 * in collaboration with <b>Jacek Radwan</b>: <a href="mailto:jradwan@amu.edu.pl
 * ?subject=MHC evolutionary dynamics - model code v1.0">jradwan@amu.edu.pl</a> \n
 * \n
 * <a href="https://sites.google.com/site/evobiolab/home">Evolutionary Biology Group</a> \n
 * Faculty of Biology \n
 * Adam Mickiewicz University \n
 * ul. Umultowska 89 \n
 * 61-614 Poznań \n
 * Poland \n
 *
 *
 * @version 0.9.1
 *
 * @section intro Introduction
 * This program is an implementation of the evolution dynamics of MHC described
 * in a number of papers e.g.: \n
 * 1. Borghans, J. a M., Beltman, J. B., & De Boer, R. J. (2004). <i>MHC polymorphism under host-pathogen
 * coevolution.</i> Immunogenetics, 55(11), 732–9. http://doi.org/10.1007/s00251-003-0630-5 \n
 * 2. Ejsmond, M., Babik, W., & Radwan, J. (2010).<i> MHC allele frequency distributions under parasite-driven
 * selection: A simulation model.</i> BMC Evolutionary Biology, 10, 332. http://doi.org/10.1186/1471-2148-10-332 \n
 * 3. Ejsmond, M. J., & Radwan, J. (2009). <i>MHC diversity in bottlenecked populations: a simulation model.</i>
 * Conservation Genetics, 12(1), 129–137. http://doi.org/10.1007/s10592-009-9998-6 \n
 * 4. Ejsmond, M. J., Radwan, J., & Wilson, A. B. (2014). <i>Sexual selection and the evolutionary dynamics of
 * the major histocompatibility complex.</i> Proceedings of the Royal Society B: Biological Sciences, 281, 20141662.
 * http://dx.doi.org/10.1098/rspb.2014.1662 \n
 * 5. Ejsmond, M. J., & Radwan, J. (2015). <i>Red Queen Processes Drive Positive Selection on Major Histocompatibility Complex (MHC) Genes.</i> PLOS Comput. Biol. 11:e1004627. http://doi.org/10.1371/journal.pcbi.1004627 \n \n
 * <b>Read them before you start tweaking anything on your own.</b>
 *
 * @section Compilation
 * The program was written in <a href="https://en.wikipedia.org/wiki/C%2B%2B14">C++14 standard</a> so if
 * you are using GCC then version gcc 4.8 seems to be the minimum requirement (I have 5.4). This program has
 * some serious dependencies on <a href="http://www.boost.org/">C++ Boost Libraries</a>. Should compile
 * smoothly on most modern GNU/Linux distros with Boost Libs installed. Having <a href="http://www.scons.org/">Scons
 * build tool</a> might be useful too. Basic compilation works fine on Ubuntu 16.04 LTS with mentioned packages
 * installed by running the command: \n
 * \n
 * $<b>
 * g++ -O3 -o SCBuild/MHC_model main.cpp nlohmann/json.hpp Gene.cpp Antigen.cpp Host.cpp Pathogen.cpp H2Pinteraction.cpp Random.cpp Tagging_system.cpp Environment.cpp DataHandler.cpp -fopenmp -std=c++14
 * </b>\n \n
 * The code here can be also used like a toolbox for your research. You can stitch your own <i>main_yourown.cpp</i> file
 * with your scenario and custom procedures (we did so for our research) and compile it using Scons script
 * which is part of this code bundle. To do so run:\n
 * \n
 * $<b> scons -Q scenario="Scenarios/main_default.cpp" linking="dynamic"</b>\n \n
 * replacing "Scenarios/main_default.cpp" with the path to your own custom <i>main_xyz.cpp</i>. Use our
 * <i>main_core.cpp</i> file as the base to develop your own scenarios of MHC evolution dynamics. WARNING!
 * Using our Scons script overwrites the <i>main.cpp</i> file in the main source directory - always
 * store your scenarios somewhere else then the <i>main.cpp</i> file!\n
 * \n
 * Sometimes the HPC Cluster you have is lame and old and it has fairly outdated compiler (e.g. gcc < 4.8). Then you can statically link the libraries on your fancy brand new PC running latest Linux distro and send the no-dependencies executable to cluster. Compile like this: \n \n
 * $<b>
 * g++ -static -O3 -o MHC_model main.cpp Gene.cpp Antigen.cpp Host.cpp Pathogen.cpp H2Pinteraction.cpp RandomNumbs.cpp Tagging_system.cpp Environment.cpp DataHandler.cpp -std=c++14 </b>\n
 * \n
 * Or run the Scons script:\n
 * \n
 * $<b> scons -Q scenario="Scenarios/main_default.cpp" linking="static" </b>\n \n
 * But it may throw some warnings e.g.: \n
 * /usr/lib/gcc/x86_64-linux-gnu/5/libgomp.a(target.o): In function gomp_target_init :
 * (.text+0xba): warning: Using 'dlopen' in statically linked applications requires at runtime the shared libraries 
 * from the glibc version used for linking \n\n  
 *
 * @section Parameters
 * <b>In the most advanced scenario the program takes exactly 17 parameters. These are:</b> \n
 * <b> 00 </b>- Program's name \n
 * <b> 01 </b>- Number of threads program will try to use on a multi-core CPU. Giving 0 will make the program use all
 * CPU cores available. \n
 * <b> 02 </b>- Number of bits in a gene. \n
 * <b> 03 </b>- Number of bits in an antigen. \n
 * <b> 04 </b>- Host population size. \n
 * <b> 05 </b>- Pathogen population size.\n
 * <b> 06 </b>- Number of pathogen species.\n
 * <b> 07 </b>- Number of genes in one host chromosome (they have two chromosomes) when the model
 * is being initialised. \n
 * <b> 08 </b>- Number of pathogen generations per one host generation. \n
 * <b> 09 </b>- Number of host generations (effective length of model run). \n
 * <b> 10 </b>- Probability of mutation in hosts ([0,1] range). \n
 * <b> 11 </b>- Probability of mutation in pathogens ([0,1] range). \n
 * <b> 12 </b>- The heterozygote advantage / lack of advantage mode. It has to be
 * 10 for heterozygote advantage or 11 for lack of thereof. \n
 * <b> 13 </b>- Probability of deleting a gene in the host ([0,1] range). \n
 * <b> 14 </b>- Probability of duplicating a gene in the host ([0,1] range). \n
 * <b> 15 </b>- Maximal number of genes permitted in one host chromosome. \n
 * <b> 16 </b>- Number of sexual partners an individual checks out before selecting one for mating. \n
 * <b> 17 </b>- Alpha factor for the host fitness function ([0,1] range).\n
 * Historically we also had one more argument, that can be added if you want : \n
 *  18- Fraction of the antigen's bits which are forbidden from changing (a.k.a. "No Mutation Bits"). \n
 * \n
 * This program can recognise simple errors in the argument list (a probability value out of
 * [0,1] range, negative values when only positive are allowed etc.), but will not recognise
 * when they don't make a 'biological' sense.
 *
 * @section The output and data visualisation
 * Program produces a number of text files containing desired data. The file
 * <i>InputParameters.json</i> contains the values of parameters fed to the model. The
 * file <i>HostsGeneDivers.csv</i> contains the main statistics of the model,
 * including the time stamp used in the visualisation process of the other files
 * with data. This two files are the basics to analyse the simulation output. Other
 * files are optional and they can be switched off by commenting out some
 * lines in the main function (<i>main.cpp</i> file) what should speed up the program runtime.
 * \n\n
  *  All the files are:\n
 * <b><i>InputParameters.json </i></b> - contains run's parametrisation in JSON text format \n
 * <b><i>HostsGeneDivers.csv</i></b> - contain basic statistics of the output for host genomes \n
 * <b><i>HostGenomesFile.[t].csv</i></b> - contains all the genomes of all the cells for host population
 * in time <i>t</i>. There can be more then one file like this for different time snapshots. \n
 * <b><i>PathoGenomesFile.[t].csv</i></b> - contains all the genomes of all the cells for pathogen population
 * in time <i>t</i>. There can be more then one file like this for different time snapshots. \n
 * <b><i>HostGeneNumbTotal.csv </i></b> - total number of genes in each individual host in each
 * time step (linked to <i>HostMHCsNumbUniq.csv </i>, that each column in <i>HostGeneNumbTotal.csv </i>
 * is represents the same cells as in <i>HostMHCsNumbUniq.csv </i>). \n
 * <b><i>HostMHCsNumbUniq.csv </i></b> -  number of unique MHC alleles in each individual host in each
 * time step (linked to <i>HostGeneNumbTotal.csv </i>, that each column in <i>HostGeneNumbTotal.csv </i>
 * is represents the same cells as in <i>HostMHCsNumbUniq.csv </i>). \n
 * <b><i> NoMutationInPathoList.csv </i></b> - list of conserved antigen sites. Each line contains one pathogen
 * species, each number indicates the index of a "no-mutation" site in this species antigen. When empty it means
 * there is no mutation restrictions. \n
 * <b><i> NumberOfMhcAfterMating.csv </i></b> -number of the unique MHC types in each individual host in each
 * time step after the mating procedure took place. \n
 * <b><i> NumberOfMhcBeforeMating.csv </i></b> - number of the unique MHC types in each individual host in each
 * time step before the mating procedure took place. \n
 * <b><i> NumberOfMhcInMother.csv </i></b> - number of the unique MHC types in each individual host that is selecting
 * a partner (a.k.a. "mother") in each time step during mating procedure. Each individual has a corresponding partner 
 * at the same index in the file <i>NumberOfMhcInFather.csv</i>. \n
 * <b><i> NumberOfMhcInFather.csv </i></b> - number of the unique MHC types in each individual host that has been
 * selected as a mating (a.k.a. "father") in each time step during mating procedure. Each individual has a corresponding
 * partner at the same index in the file <i>NumberOfMhcInMother.csv</i>. \n
 * <b><i> PresentedPathogenNumbers.csv </i></b> - number of presented pathogens by each individual in each time step. \n

 * \n
 * Visualisation is done using Python 3.6 scripts containing a a lot of calls to Numpy, Matplotlib and other scietific
 * Python libraries. You may wish to consider using the
 * <a href="https://www.anaconda.com/download/">Python Anaconda project</a>. These scripts
 * can be found in <i>PyScripts</i> directory. \n
 * To load the environment for output analysis run: \n\n
 * $ <b>conda env create -f mhcEvoCondaEnv.yaml</b> \n
 *
 */
