# PMCML-Polarized-Monte-Carlo-for-Multiple-Layers
Polarized Monte Carlo for multiple layers(PMCML) is a light transport model used to model the polarized light propagation in a biological medium

This model studies the polarized light transport in a biological tissues with different layers and heterogeneties.
PMCML has been built upon the light scattering model for polarized light developed by _Jessica.et.al_. All biological tissues contains different types of cell, anisotropy and heterogeneties. To account for changes in the polarized in such heterogenetic medium, we have studied the light transport of polarized light to understand the variations in multilayered medium and how polarized light state changes in a multiple layered tissue having different optical properties. 

The boundary mismatch conditions and weight attenuation due to Fresnel's reflection has been implemented along with the light transport modelled using Monte Carlo method.

![boundary](https://user-images.githubusercontent.com/86607064/125151172-dadea100-e162-11eb-84fd-2f968b4493da.png)
The Fig 1. shows the 3 layer model we have considered with a slab geometry with scttering of polarized light modeled in Meridian plane geometry. The reflectiona and transmission beams from an interface is marked with the corresponding angles of incidence and refraction(transmission) taken for computation.

The Mueller matrix generated using this code has been analysed with different scattering, absorption coefficients, nanoparticles and wavelength.

**RESULTS**



**INSTRUCTIONS FOR COMPILING AND EXECUTING**
1. Load all necessary C header files and compiler
2. Load the folder of pmcml.zip
3. Use an IDE to generate the make file
4. 'iquv' executable will be created
5. Execute using ./iquv
6. 16 data files shall be generated 
7. Postprocessing of data files and computation of Mueller matrix can be done using these.


***NOTE: Users can define their own values for absorption, scatterer size, wavelength and incident Stokes vector, refractive indices and layer thickness with number of layers.



The postprocessing was performed using MATLAB R2020b. 

