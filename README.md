# noble_gases_in_dust

The present tool provides two general models for gas implantation in dust, the Single Grain Model (SGM) and Multiple Grain Model (MGM) based on Verchovsky et al. 2003 PASA (https://ui.adsabs.harvard.edu/abs/2003PASA...20..329V/abstract). 
The output is the concentration of the gas implanted, based on the size of the grain population and the range of implantation depth (S). These data can be used to study isotopic anomalies in meteoritic materials for elements where we expect that the implanted component is dominating over the condensed component, like for noble gases.

The S grid obtained can be eventually mapped on an ion energy grid depending on the ion size and on the grain type. This conversion would required to user to map the obtained S distribution here with some particle interaction model for dust (a reference is given in the notebook). 

In comparison with the MGM model results, as example the data for noble gases Kr and Xe measured in presolar the KJ samples of SiC grains by Lewis et al. 1994 (https://ui.adsabs.harvard.edu/abs/1994GeCoA..58..471L/abstract). 

Some basic plotting is done in comparison with observations, and in general the results by Verchovsky et al. 2004 ApJ (https://ui.adsabs.harvard.edu/abs/2004ApJ...607..611V/abstract) are confirmed: Kr was implanted with a single low implantation range (or low energy) component in dust, while Kr shows two components, a high implantation range (or high energy) component and a low energy component.

Notice that for Xe we adopt a sampling different than the one adopted by Lewis et al. 1994, since we use all the KJ groups, and the LQ samples are not used. Therefore, the G-component and N-component slopes obtained are different. Of course an analogous approach could be used as in Lewis et al., but in the case of Xe it would imply to use slopes outside the grain size range (KJA, KJG and KJH) considered to fit the grains data with a linear regression. 
