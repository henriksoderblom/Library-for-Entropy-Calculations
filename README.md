# Library-for-Entropy-Calculations
A library of functions for calculating different properties with data obtained with the simulation tool UppASD. The functions can be used for data obtained in other ways, but the files should be correctly written (see CAUTION).

Tool for my master thesis "Comparison of Methods for Characterization of Magnetocaloric Materials" in physics. Code is commented but not documented yet, since it is still changing while I am working on the thesis. Documentation will be available after thesis is finished.


CAUTION: At the moment the functions are defined with specific file names. This will be changed to a better system. In the code it is commented which columns are being used for calculations, and what quantities those columns are representing in the calculations.

Info:
  Entropy_with_simulatedCv: Calculates entropy using specific heat capacity obtained from UppASD as C=(<E^2>-<E>^2)*k_B*(1/k_B*T)^2.
  
  Entropy_with_calculatedCv: Calculates entropy using specific heat calculated from dU/dT.
  
  Entropy_with_CanonicalEnsemble: Calculating entropy using the magnon density of states (canonical ensemble).

  plotEntropy: Plotting total entropies and entropy differences. 
  
  plotMagnonSpectrum: Plotting adiabatic and/or correlated magnon spectrum.
  
  plotAMSmDOS: Plotting the magnon density of states obtained from the adiabatic magnon spectrum. Integration over the whole function is also done to check value.
  
  
  
