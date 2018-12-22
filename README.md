Evolution of MHC repertoire and number of copies:
==============================

This work is part of the project [**Evolution of the number of copies of MHC genes: testing the optimality hypothesis and exploring alternatives.**](https://sites.google.com/site/evobiolab/projects) supported by the Polish Science Centre (NCN).

Submitted - waiting for reviews.

Documentation:
-----------

Go from here to *./Doc/html/* and launch the *./Doc/html/index.html* file in your web browser. The documentation was created with Doxygen.

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
The program was written in https://en.wikipedia.org/wiki/C%2B%2B14 C++14 standard so if you are using GCC then version gcc 4.8 seems to be the minimum requirement. This program has
some serious dependencies on http://www.boost.org/ C++ Boost Libraries. Should compile smoothly on most modern GNU/Linux distros with Boost Libs installed. Having http://www.scons.org/ Scons build tool might be useful too. Basic compilation works fine on Ubuntu 14.04 LTS as well as 16.04 LTS with mentioned packages installed by running the command:
```bash
g++ -O3 -o MHC_model main.cpp Gene.cpp Antigen.cpp Host.cpp Pathogen.cpp H2Pinteraction.cpp RandomNumbs.cpp Tagging_system.cpp Environment.cpp DataHandler.cpp -std=c++14
```

It should also compile with flag -std=c++11 only with minor warnings. The code here can be also used as a toolbox for your research. You can stitch your own *main_yourown.cpp* file with your scenario and tailored procedures (we did so for our research) and compile it using Scons script which is part of this code bundle. To do so run:
```Bash
scons -Q scenario="Scenarios/main_default.cpp"
```
replacing *Scenarios/main_default.cpp* with the path to your own custom *main_xyz.cpp*. Use our *main_core.cpp* file as the base to develop your scenarios of MHC evolution dynamics. WARNING!
Using our Scons script overwrites the *main.cpp* file in the main source directory - always
store your scenarios somewhere else than the *main.cpp* file!

Sometimes the HPC Cluster you have is lame and old and it has fairly outdated compiler (e.g. gcc < 4.8). Then you can statically link the libraries on your fancy brand new PC running "The Awesome Linux 3000" distro and send the no-dependencies executable to cluster. Compile like this:
```bash
g++ -static -O3 -o MHC_model main.cpp Gene.cpp Antigen.cpp Host.cpp Pathogen.cpp H2Pinteraction.cpp RandomNumbs.cpp Tagging_system.cpp Environment.cpp DataHandler.cpp -std=c++14
```

Or run the Scons script:
```shell
scons -Q scenario="Scenarios/main_default.cpp" linking="static"
```
