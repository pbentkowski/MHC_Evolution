#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
'''
This is SCons script to compile MHC evolution simulation model.
The main idea behind it is that in directory 'Scenarios' one
creates a main_sth_sth.cpp file which is build from blocks coded
in the source files. This script compiles the 'scenario' and
writes the executable to SCBuild directory.
One retrieves the executable from there and runs the show. Be sure
you complied the model on the architecture you gonna run in at.

Created on Thu May 28 13:46:49 2015
for Evolutionary Biology Group, Faculty of Biology
    Adam Mickiewicz University, Poznań, Poland
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
'''
from shutil import copyfile
from os import getcwd

cxxflaggs = "-std=c++1y"
cppath = '/usr/include/boost/'

# normal compilation
env_dynamic = Environment(CCFLAGS='-O3',
              CPPPATH=cppath,
              CXXFLAGS=cxxflaggs)

# debugging compilation
env_dbg = Environment(CCFLAGS='-g',
          CPPPATH='/usr/include/boost/',
          CXXFLAGS=cxxflaggs)

# static compilation
env_static = Environment(CCFLAGS='-O3 -static',
             CPPPATH=cppath,
             CXXFLAGS=cxxflaggs)

scenario = ARGUMENTS.get('scenario', 0)
try:
  mainn = str(scenario)
  #local_main = getcwd() + "/" + mainn.split("/")[-1]
  local_main = getcwd() + "/main.cpp"
  copyfile(mainn, local_main)
except:
  print "Using the default main.cpp to build the program."
  mainn = "Scenarios/main_default.cpp"
  local_main = getcwd() + "/main.cpp"
  copyfile(mainn, local_main)

OUTprog = "SCBuild/" + mainn.split("/")[-1].split(".")[0]
SRS = ['DataHarvester.cpp', 'Environment.cpp', 'Gene.cpp',
       'Antigen.cpp','H2Pinteraction.cpp', 'Host.cpp',
       'Pathogen.cpp','RandomNumbs.cpp', 'Tagging_system.cpp',
       local_main]
if True:
    t = env_static.Program(target=OUTprog, source=SRS)
else:
    t = env_dynamic.Program(target=OUTprog, source=SRS)
Default(t)
