cmake_minimum_required(VERSION 3.13.4)

project(FastElectrodeSegmentor)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/FastElectrodeSegmentor")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Marcos Oliveira (University of São Paulo)")
set(EXTENSION_DESCRIPTION "CSIM group Module")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/FastElectrodeSegmentor.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/FastElectrodeSegmentor/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(FastElectrodeSegmentor)
add_subdirectory(Interface)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
