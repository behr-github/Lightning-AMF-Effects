# Lightning-AMF-Effects
This repository contains the code base that supports "Quantification of the effect of modeled lightning NO2 on UV-visible air mass factors"

The two .m files in the top level are the main programs that you should be aware of.

* `lnox_amf_paper_figs.m` calls the necessary functions to make the figures included in the paper. From here, you can follow the function call path to see what analysis functions were used.
* `match_wrf2aircraft.m` handles the matching of WRF data to DC3 data.

The namelists necessary to run WRF are stored in the WRF-Namelits subdirectory. The WPS namelist is common for all runs, there are three WRF
namelist.input files:

1. One with lightning off and no nudging
2. One with lightning on and no nudging
3. One with lightning on and nudging on

Not all permutations of mol NO/flash and flashrate are included; the only settings that need changed for that are the 
n\_ic and n\_cg values under &chem and the flashrate\_factor value under &phys.
