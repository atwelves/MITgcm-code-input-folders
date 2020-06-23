# MITgcm-code-input-folders
Self-shading and meltwater spreading control transition from light to iron limitation in an Antarctic coastal polynya

Code folders and input files for the experiments described in the above paper

no_melt: experiment with ice shelf melt suppressed
melt_pump: experiment with ice shelf melt active
gmw_iron: melt_pump plus an explicit source of iron in meltwater
fixed_zeu: gmw_iron minus phytoplankton self-shading
flat_iron: gmw_iron minus iron uptake
max_yield: gmw_iron minus iron uptake or phytoplankton self-shading

Sensitivity experiments are performed by changing files in data.obcs and data.exf

MITgcm checkpoint: 66k

https://github.com/MITgcm/MITgcm/releases/tag/checkpoint66k

All that is needed to perform these experiments is then the init and rbcs sub-folders, which will be made available upon request from andrew.twelves@ed.ac.uk
