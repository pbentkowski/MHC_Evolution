Python scripts for MHC Evo data visualization
----------------------------------------------

The scripts are written in Python 3.6. We recommend Anaconda scientific Python distribution for using them:
https://www.anaconda.com/download/

When you have Anaconda installed run:
```shell
conda env create -f mhcEvoCondaEnv.yaml
```
to get all the stuff you need to run the scripts (and unfortunately some more...)

#### What the scripts do

* ***animateMHCsHist.py*** - reads files with genome size histograms (`HostGeneNumbTotal_ChrOne.csv` and `HostMHCsNumbUniq_ChrOne.csv`) and a general host statistics file `HostsGeneDivers.csv` to render an animation how the size of the hosts' chromosomes and the number of MHC unique alleles in them evolve. You probably need a video codec like i.g. *ffmpeg* for *matplotlib* to be able
to create a MP4 animation clip.
* ***antigen_similiraty.py*** - checks the similarities in all the antigens in the hole pathogen population. Function `loadThePopulation(FILE)` loads the file created by the modelling framework into a handy data structure and the rest calculates various statistics.
* ***evolution_big_stats.py*** - iterates through directories and looks for the file with the Host population snapshot called `HostGenomesFile.XXXX.csv` and the `InputParameters.json` file with the parameters used it the run. It extracts information about the genes origin like ancestry tree and MRCA.
* ***evolution_mut_count.py*** - loads the file `HostGenomesFile.XXXX.csv` (final snapshot of the host population) and calculates how many mutation got fixated during the MHCs' evolution. Plots the histogram.
* ***file_len.py*** - just an utility function. Counts the number of lines in a text file.
* ***get_params_from_inputFiles.py*** - searches for `InputParameters.json` files, pulls out parameters from them and renders them in one line which can be feed as input to the model's program.
* ***host_heterozygoty_check.py*** - checks what percentage of host population has no MHC gene repetitions in their genomes.
* ***infection_vs_MHC_stats.py*** - **??????**
* ***K_succes_in_N.py*** - calculates the theoretical probability of [getting a run of K or more successes (heads) in a row in N Bernoulli trials (coin flips)?](http://www.askamathematician.com/2010/07/q-whats-the-chance-of-getting-a-run-of-k-successes-in-n-bernoulli-trials-why-use-approximations-when-the-exact-answer-is-known/)
* ***mhc_plot_selected_runs.py*** - plots the generic statistics referring to time evolution of MHCs in hosts. These are the number of MHC types, diversity of MHC types, host fitness etc.
* ***MHC_segregation_sex_select.py*** - calculates if MHC genes are correlated in co-occurrence in individuals as a result e.g. sex selection.
* ***MHC_similiraty.py*** - calculates and plots similarities between antigens and genes in pathogen and host populations. Does it separately for hosts and for pathogens.
* ***packed_plots_of_MHC_alleles.py*** - walks the directory tree looking for model runs which are characterized by same parametrisation as the template file provided by the user (e.g `Template.json`). Then process these results by fancy stats and plots that processed output on a nice graph.
* ***patho_bitgene_filter.py*** - filters the file with Pathogen population data to leave off only the bit string representation of the antigen removing the mutation history data.
* ***pathogen_spp_cooccur.py*** - **??????**
* ***plot_Chrom_size_last_shot.py*** - plots histogram of number of MHC alleles and all MHC genes in one chromosome at the end of simulation.
* ***sex_scenarios_comp.py*** - uses a post-processed file (you may need to edit it manually) and creates a nice boxplot based on the data in that file. Check the `datype` data type to see what kind of file you need. Script `packed_plots_of_MHC_alleles.py` may be useful in creating this file. E.g. is `Integr_16_1e5`.
* ***sex_selection_on_MHC_numb.py*** - uses the files `NumberOfMhcInMother.csv`, `NumberOfMhcInFather.csv`, `NumberOfMhcBeforeMating.csv` and `InputParameters.json` to analyse the strength of selection preference on partners' MHC type number depending on the sexual selection scenario used.
* ***testing_bitstrings.py*** -  **??????**
* ***transform_paramInput_2_json.py*** - iterates through the old simulation results and replaces the old `InputParameters.csv` file with a newer Json version of the parameter input file, that is easier to compare with templates in multi-simulation comparisons.

Ipython Notebooks mostly for plotting:

* ***plotting_for_publication.ipynb***
* ***Mutation_mode_compr.ipynb***  
* ***plotting_for_publication.ipynb***
* ***sex_selection.ipynb***
