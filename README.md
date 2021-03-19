# FastElectrodeSegmentor

Module contructed for 3D Slicer using C++ and Python.

The algorithm is able to find and reconstruct deep electrodes from brain CT from SEEG patients using a 4D search structure. The method can reconstruct all founded electrodes from patient, to do it, method use MRI too and execute another steps like image registration, brain tissue segmentation, and tissue reconstruction, at the end, the results shows all brain reconstruction with electrodes inside an accurate position.

The method has dependencies as follows:

[SlicerElastix](https://github.com/lassoan/SlicerElastix)

[ModifiedqEntropySegmentation](https://github.com/MarcosSoares10/ModifiedqEntropySegmentation)

[ROBEXBrainExtraction](https://github.com/MarcosSoares10/ROBEXBrainExtraction)
