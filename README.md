# CHIMERA_TERRESTRIAL_PLANETS
Welcome to the CHIMERA "terresrial" planet transmission retrieval code. This was the code used in Tremblay et al. 2020. Further description of the code/assumptions is in that publication.  Before you begin, we highly recommend that you first do the jupyter notebook tutorials in https://github.com/mrline/CHIMERA before using this.  There you will gain familiarity with retrieval concepts, program flow, and installation of pymultinest.  

Instructions for "use":
1. Pull this repository (all the folders maintaining the same directory structure)
2. The ABSCOEFF_TERRESTRIAL folder is empty.  Go to https://www.dropbox.com/sh/hql4yc4zjdomzqs/AAC8SCAMPZGBvfWOv4B6TR0va?dl=0 or the corrosponding zenodo (https://zenodo.org/record/3630884#.XjJWnhdKjOR) and download the correlated-K pickles. Such files are too large for github.
3. Go into the TEMPLATE_TRANSMISSION_TERRESTRIAL and open/run the make_spec.py routine. This will generate a "simulated" spectrum given various typical atmosphere parameters (e.g., gas mixing ratios, temperature, clouds, star/planet radii etc.).  See the comments in make_spec.py for more details. Feel free to play around with this routine. Upon running make_spec.py in command line, a plot of the spectrum with simulated noise will appear. Also a "data.pic" will be generated as well, which contains the wavelength grid, transit depth, and transit depth errors.
4.  To run the retrieval (again, we assume you are familiar with pymultinest and its api) on multiple cores (assuming you installed it correctly), type mpirun -np 4 python call_pymultinest.py.  With 1000 livepoints and 4 cores, it took just over 7 hours to do the 2-11 um spectrum at R=100 (this case, about 365,793 likelihood/model evaluations).  
5. Run the routine plot_PMN.py to produce the corner plot.
6. If you want to compute "detection significances" you'll need to modify call_pymultinest.py by removing a gas parameter (and fixing to -12 or something) and re-running.  You can compute the bayes factors by comparing the evidence from the full model with all gases and the one with one less gas. See Tremblay+2020 for what this means.

The corner plots for the R=50 and R=100 cases in Tremblay et al. 2020 can be found in the zip file (Corner_Plots_Tremblay_2020.zip)
