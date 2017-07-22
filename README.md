# Lightning-AMF-Effects
This repository contains the code base that supports "Quantification of the effect of modeled lightning NO2 on UV-visible air mass factors"
The two .m files in the top level are the main programs that you should be aware of.
`lnox_amf_paper_figs.m` calls the necessary functions to make the figures included in the paper. From here, you can follow the function call path to see what analysis functions were used.
`match_wrf2aircraft.m` handles the matching of WRF data to DC3 data.
