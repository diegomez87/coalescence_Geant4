# coalescence_Geant4
Autors: Diego Gomez Coral: diegomez@fisica.unam.mx
        Eulogio Serradilla
        
Coalescence implementation in Geant4
In order to generate antideuterons in proton-proton 
collisions or proton-nucleus collisions a subroutine 
has been added into G4TheoFSGenerator.

Such addition simply simulate antideuteron production,
comparing momentum differences between antiprotons and 
antineutrons with a free parameter called coalescence momentum.
When differences are less than the value of the parameter
antideuterons are created. The value of the coalescence parameter
has been chosen to fit experimental data and has been set to a 
wide energy range.
