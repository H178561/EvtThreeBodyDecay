# EvtThreeBodyDecay

## Register the model in EvtGen

In `src/EvtGenModels/EvtModelReg.cpp' add the lines

```
#include "EvtGenModels/EvtThreeBodyDecays.hh"
```
and 
```
modelist.registerModel( new EvtThreeBodyDecays);
```

## usage in Particle Decay Generators

For the usage in pythia or other generators with EvtGen one has to call the model in the .dec File with the respecting .json File

```
Decay Lambda_c+
    1.00  p+ pi+ K- THREEBODYDECAYS lc2ppik-lhcb-2683025.json;
Enddecay
CDecay anti-Lambda_c-

End

```
