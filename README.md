This repository contains code for simulating the flow-driven adiabatic inversion process underlying pCASL MRI. The scripts in this repository originate from two submissions: Magnetic Resonance in Medicine (MRM) and MethodsX.

MRM submission includes three main functions:
1. "MRM_Main1a_AdiabaticInversion_SingleVelocity_LabeledScan": A Bloch-equation–based simulation for the label scan of pCASL MRI, focusing on a single flow velocity.
2. "MRM_Main1b_AdiabaticInversion_SingleVelocity_ControlScan": A parallel simulation for the control scan, also focusing on a single flow velocity.
3. "MRM_Main2_AdiabaticInversion_MultiVelocity": Extending the previous two simulations to include variable flow velocities, allowing investigation of the dependence of labeling efficiency on flow velocity.

MethodsX submission includes one main function:
"MethodsX_Main1_BIC_Protocol_FlowPatternSim": This simulation consists of three subcomponents and uses the previously established relationship between flow velocity and labeling efficiency (Magn Reson Med, 2025; DOI: 10.1002/mrm.70086).
Part 1: Examining the effects of arterial radius and flow velocity on labeling efficiency. Arterial radius and flow velocity were varied independently to generate multiple combinations, and labeling efficiency was determined for each.
Part 2: Comparing cross-sectional labeling-efficiency distributions between laminar and plug flow profiles.
Part 3: Comparing global labeling efficiency between laminar and plug flows across a range of peak flow velocities.

The remaining files and functions are called by the main functions listed above.
