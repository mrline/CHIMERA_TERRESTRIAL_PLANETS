# CHIMERA_TERRESTRIAL_PLANETS
Welcome to the CHIMERA "terresrial" planet transmission retrieval code. This was the code used in Tremblay et al. 2020. Further description of the code/assumptions is in that publication.  Before you begin, we highly recommend that you first do the jupyter notebook tutorials in https://github.com/mrline/CHIMERA before using this.  There you will gain familiarity with retrieval concepts, program flow, and installation of pymultinest.  

Instructions for "use":
1. Pull this repository (all the folders maintaining the same directory structure)
2. The ABSCOEFF_TERRESTRIAL folder is empty.  Go to https://www.dropbox.com/sh/hql4yc4zjdomzqs/AAC8SCAMPZGBvfWOv4B6TR0va?dl=0 and download the correlated-K pickles.
3. Go into the TEMPLATE_TRANSMISSION_TERRESTRIAL and open/run the make_spec.py routine. This will generate a "simulated" spectrum given various typical atmosphere parameters (e.g., gas mixing ratios, temperature, clouds, star/planet radii etc.).  See the comments in make_spec.py for more details. Feel free to play around with this routine. Upon running make_spec.py in command line, a plot of the spectrum with simulated noise will appear. Also a "data.pic" will be generated as well, which contains the wavelength grid, transit depth, and transit depth errors.
4. 
