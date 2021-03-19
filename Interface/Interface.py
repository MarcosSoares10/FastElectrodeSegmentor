import os
import platform
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import SimpleITK as sitk
import sitkUtils
import subprocess
from datetime import datetime
import shutil
import math
import numpy as np


#
# Interface
#

class Interface(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Interface" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
        
    self.parent.contributors = ["John Doe (AnyWare Corp.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
It performs a simple thresholding on the input volume and optionally captures a screenshot.
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.
    self.parentPath = os.path.dirname(parent.path)

#
# InterfaceWidget
#

class InterfaceWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the CT to the algorithm." )
    parametersFormLayout.addRow("Input CT Volume: ", self.inputSelector)

        #
    # input volume selector
    #
    self.inputSelectorMRI = slicer.qMRMLNodeComboBox()
    self.inputSelectorMRI.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelectorMRI.selectNodeUponCreation = True
    self.inputSelectorMRI.addEnabled = False
    self.inputSelectorMRI.removeEnabled = False
    self.inputSelectorMRI.noneEnabled = True
    self.inputSelectorMRI.showHidden = False
    self.inputSelectorMRI.showChildNodeTypes = False
    self.inputSelectorMRI.setMRMLScene( slicer.mrmlScene )
    self.inputSelectorMRI.setToolTip( "Pick the MRI to the algorithm." )
    parametersFormLayout.addRow("Input MRI Volume: ", self.inputSelectorMRI)

    #
    # output volume selector
    #
    '''self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = True
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output to the algorithm." )
    parametersFormLayout.addRow("Output Volume: ", self.outputSelector)'''

    #
    # threshold value
    #
    #self.imageThresholdSliderWidget = ctk.ctkSliderWidget()
    #self.imageThresholdSliderWidget.singleStep = 1
    #self.imageThresholdSliderWidget.minimum = 0
    #self.imageThresholdSliderWidget.maximum = 255
    #self.imageThresholdSliderWidget.value = 1
    #self.imageThresholdSliderWidget.setToolTip("Set threshold value for computing the output image. Voxels that have intensities lower than this value will set to zero.")
    #parametersFormLayout.addRow("Image threshold", self.imageThresholdSliderWidget)

    #
    # check box to trigger taking screen shots for later use in tutorials
    #
    #self.enableScreenshotsFlagCheckBox = qt.QCheckBox()
    #self.enableScreenshotsFlagCheckBox.checked = 0
    #self.enableScreenshotsFlagCheckBox.setToolTip("If checked, take screen shots for tutorials. Use Save Data to write them to disk.")
    #parametersFormLayout.addRow("Enable Screenshots", self.enableScreenshotsFlagCheckBox)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = True
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    #self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    #self.inputSelectorMRI.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
   

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    #self.onSelect()

  def cleanup(self):
    pass

  #def onSelect(self):
  #  self.applyButton.enabled = self.inputSelector.currentNode() and self.inputSelectorMRI.currentNode() and self.outputSelector.currentNode()

  def onApplyButton(self):
    logic = InterfaceLogic()
    #imageThreshold = self.imageThresholdSliderWidget.value
    imageThreshold = 240
    logic.run(self.inputSelector.currentNode(), self.inputSelectorMRI.currentNode(),imageThreshold)

#
# InterfaceLogic
#

class InterfaceLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True

  '''def isValidInputOutputData(self, inputVolumeNode, inputMRIVolumeNode):
    """Validates if the output is not the same as input
    """
    if not inputVolumeNode:
      slicer.util.errorDisplay('isValidInputOutputData failed: no CT input volume node defined') #logging.debug
      return False
    if not inputMRIVolumeNode:
      slicer.util.errorDisplay('isValidInputOutputData failed: no MRI input volume node defined')
      return False
    if inputVolumeNode.GetID()==outputVolumeNode.GetID() and inputMRIVolumeNode.GetID()==outputVolumeNode.GetID():
      slicer.util.errorDisplay('isValidInputOutputData failed: one of inputs and output volume is the same. Create a new volume for output to avoid this error.')
      return False
    return True'''
  


  def makeModel(self,volumeNode, labelNumber, smoothValue,namenode):
    #the volume need to present must valid
    if not volumeNode:
      return


    # set up the model maker node
    parameters = {}
    parameters['Name'] = "EditorModel"
    parameters["InputVolume"] = volumeNode.GetID()
    parameters['FilterType'] = "Sinc"

    # build only the currently selected model.
    parameters['Labels'] = labelNumber #using default label, that can further improved by using label that user selected
    parameters["StartLabel"] = -1
    parameters["EndLabel"] = -1

    parameters['GenerateAll'] = False
    parameters["JointSmoothing"] = False
    parameters["SplitNormals"] = True
    parameters["PointNormals"] = True
    parameters["SkipUnNamed"] = True
    parameters["Decimate"] = 0.25
    parameters["Smooth"] = smoothValue #defaul smooth parameter

    #
    # output
    # - make a new hierarchy node if needed
    #
    numNodes = slicer.mrmlScene.GetNumberOfNodesByClass( "vtkMRMLModelHierarchyNode" )
    outHierarchy = None
    for n in range(numNodes):
      node = slicer.mrmlScene.GetNthNodeByClass( n, "vtkMRMLModelHierarchyNode" )
      if node.GetName() == "Editor Models":
        outHierarchy = node
        break

    if not outHierarchy:
      outHierarchy = slicer.vtkMRMLModelHierarchyNode()
      outHierarchy.SetScene( slicer.mrmlScene )
      outHierarchy.SetName( namenode )
      slicer.mrmlScene.AddNode( outHierarchy )

    parameters["ModelSceneFile"] = outHierarchy

    modelMaker = slicer.modules.modelmaker

    #
    # run the task (in the background)
    # - use the GUI to provide progress feedback
    # - use the GUI's Logic to invoke the task
    # - model will show up when the processing is finished
    #
    slicer.cli.run(modelMaker, None, parameters)

  def remove(self, path):
    """ param <path> could either be relative or absolute. """
    if os.path.isfile(path) or os.path.islink(path):
        os.remove(path)  # remove the file
    elif os.path.isdir(path):
        shutil.rmtree(path)  # remove dir and all contains
    else:
        raise ValueError("file {} is not a file or dir.".format(path))

  def saveNode(self, node, filename, properties={}):
        """Save 'node' data into 'filename'.
    It is the user responsibility to provide the appropriate file extension.
    User has also the possibility to overwrite the fileType internally retrieved using
    method 'qSlicerCoreIOManager::fileWriterFileType(vtkObject*)'. This can be done
    by specifying a 'fileType'attribute to the optional 'properties' dictionary.
    """
        from slicer import app
        properties["nodeID"] = node.GetID()
        properties["fileName"] = filename
        if hasattr(properties, "fileType"):
            filetype = properties["fileType"]
        else:
            filetype = app.coreIOManager().fileWriterFileType(node)
        return app.coreIOManager().saveNodes(filetype, properties)
  
  def checkifelectrodeexists(self, arrayelectrodes,target):
    check = False
    for electrode in arrayelectrodes:
        if np.all(electrode == target):
            check = True
            break
        else:
            check = False
    return check

  def distance(self, source, target):  
    distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(source, target)]))
    print("Euclidean distance from x to y: ",distance)
   
  def run(self, inputVolume, inputMRIVolume, Thresholdvalue, enableScreenshots=0):
    """
    Run the actual algorithm
    """

    #if not self.isValidInputOutputData(inputVolume, inputMRIVolume, outputVolume):
    #slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
    #  return False

    logging.info('Processing started')

     
    #volmn = sitkUtils.PullVolumeFromSlicer(inputVolume.GetID())
    #logging.info('A')
    #auximg=sitk.Cast(sitk.RescaleIntensity(volmn), sitk.sitkUInt8)
    #logging.info('B')
    #volumesLogic = slicer.modules.volumes.logic()
    #rescaledinputimage = volumesLogic.CloneVolume(slicer.mrmlScene, outputVolume, 'Rescaled Input Image')
    #thresholdedVolumeNode = volumesLogic.CloneVolume(slicer.mrmlScene, outputVolume, 'Thresholded Image')
    #print auximg.GetWidth()
    #print auximg.GetHeight()
    #print auximg.GetDepth()
    #print auximg.GetOrigin()
    #print auximg.GetSpacing()
    #print auximg.GetDirection()
    

    
    #sitkUtils.PushVolumeToSlicer(auximg, rescaledinputimage.GetID())
    
    #Corregistro CT MRI
    #ROBEX
    #SEGMENTAR ELETRODOS
    #SEGMENTAR TECIDOS CEREBRO
    #RECONSTRUCAO 3D
    


  ###################################################################################################################### 
    #os.path.abspath(__file__)
    #os.path.dirname(os.path.abspath(__file__))

    now = datetime.now()
    distinctsubpathname = str(now.strftime("%m%d%Y_%H%M%S"))

    #Remove Directories
    if os.path.exists(os.path.dirname(__file__) + "/Tmp"):
      self.remove(os.path.dirname(__file__) + "/Tmp")
    
    os.mkdir(os.path.dirname(__file__) + "/Tmp")
    os.mkdir(os.path.dirname(__file__) + "/Tmp/"+distinctsubpathname)
    os.mkdir(os.path.dirname(__file__) + "/Tmp/"+distinctsubpathname+"/transform")
    os.mkdir(os.path.dirname(__file__) + "/Tmp/"+distinctsubpathname+"/result")

    


    self.ElastixExecutable = os.path.dirname(__file__)+"/Elastix/bin/elastix"
    self.TransformixExecutable = os.path.dirname(__file__)+"/Elastix/bin/transformix"
    self.fixedInputFile = os.path.dirname(__file__) + "/Tmp/"+distinctsubpathname+"/fixedtmp.nii"
    self.movingInputFile = os.path.dirname(__file__) + "/Tmp/"+distinctsubpathname+"/movingtmp.nii"
    self.transformpath = os.path.dirname(__file__) + "/Tmp/"+distinctsubpathname+"/transform"
    self.resultelastixpath = os.path.dirname(__file__) + "/Tmp/"+distinctsubpathname+"/result"
    self.parameterfile = os.path.dirname(__file__) + "/Resources/Parameters_Rigid.txt"
    self.transformparameterfile = self.transformpath + "/TransformParameters.0.txt"
    self.resourcespath = os.path.dirname(__file__) + "/Resources"
    self.electrodefilepath = os.path.dirname(__file__) + "/Tmp/"+distinctsubpathname

    volumesLogic = slicer.modules.volumes.logic()


    

    

    ############################################DESATIVEDFORTEST######################################################

    if inputMRIVolume is None:
        namemriaux= "MRI_ATLAS_"+str(now.strftime("%m%d%Y_%H%M%S"))
        slicer.util.loadVolume(self.resourcespath+"/atlasImage.mha", {'name': namemriaux})
        self.parameterfile2 = os.path.dirname(__file__) + "/Resources/Parameters_BSpline.txt"
        self.transformparameterfile = self.transformpath + "/TransformParameters.1.txt"
        properties = {'useCompression': 0}
        slicer.util.saveNode(inputVolume, self.fixedInputFile, properties)
        slicer.util.saveNode(slicer.util.getNode(namemriaux), self.movingInputFile, properties)
        subprocess.call(r'"'+self.ElastixExecutable+'" -f "'+self.fixedInputFile+'" -m "'+self.movingInputFile+'" -out "'+self.transformpath+'" -p "'+self.parameterfile+'" -p "'+self.parameterfile2+'"', shell=True)
        subprocess.call(r'"'+self.TransformixExecutable+'" -tp "'+self.transformparameterfile+'" -out "'+self.resultelastixpath+'" -in "'+self.movingInputFile+'" -def all', shell=True)
    else:
        properties = {'useCompression': 0}
        slicer.util.saveNode(inputVolume, self.fixedInputFile, properties)
        slicer.util.saveNode(inputMRIVolume, self.movingInputFile, properties)
        subprocess.call(r'"'+self.ElastixExecutable+'" -f "'+self.fixedInputFile+'" -m "'+self.movingInputFile+'" -out "'+self.transformpath+'" -p "'+self.parameterfile+'"', shell=True)
        subprocess.call(r'"'+self.TransformixExecutable+'" -tp "'+self.transformparameterfile+'" -out "'+self.resultelastixpath+'" -in "'+self.movingInputFile+'"', shell=True)
    
    #slicer.mrmlScene.RemoveNode(inputVolume)
    #slicer.mrmlScene.RemoveNode(inputMRIVolume)

    outputTransformNode = slicer.vtkMRMLTransformNode()
    outputTransformNode.SetName("transform_"+str(now.strftime("%m%d%Y_%H%M%S")))
    slicer.mrmlScene.AddNode(outputTransformNode)

    outputTransformPath = os.path.join(self.resultelastixpath, "deformationField.mhd")
    [success, loadedOutputTransformNode] = slicer.util.loadTransform(outputTransformPath, returnNode = True)
    if success:
      if loadedOutputTransformNode.GetReadAsTransformToParent():
        outputTransformNode.SetAndObserveTransformToParent(loadedOutputTransformNode.GetTransformToParent())
      else:
        outputTransformNode.SetAndObserveTransformFromParent(loadedOutputTransformNode.GetTransformFromParent())
      #slicer.mrmlScene.RemoveNode(loadedOutputTransformNode)

   
       
        
    
    nameaux= "Brain_Reg_"+str(now.strftime("%m%d%Y_%H%M%S"))
    slicer.util.loadVolume(self.resultelastixpath+"/result.mhd", {'name': nameaux})
    


    outputVolume = volumesLogic.CloneVolume(slicer.mrmlScene, inputVolume, "Aux_Crop_"+str(now.strftime("%m%d%Y_%H%M%S")))
    brainmask = volumesLogic.CreateAndAddLabelVolume( slicer.mrmlScene, inputVolume, "brainmask" )
    path2files = os.path.dirname(__file__)
    
    cliParams={}
    cliParams["inputVolume"] = slicer.util.getNode(nameaux).GetID()
    cliParams["outputVolume"] = outputVolume.GetID()
    cliParams["brainMask"] = brainmask.GetID()
    if platform.system() is "Windows":
      cliParams["datPath"] = path2files + '\\Resources\\dat\\'
      cliParams["refVolsPath"] = path2files + '\\Resources\\ref_vols'
    else:
      cliParams["datPath"] = path2files + '/Resources/dat/'
      cliParams["refVolsPath"] = path2files + '/Resources/ref_vols'

    slicer.cli.run(slicer.modules.robexbrainextractioncore, None, cliParams, wait_for_completion=True)
    BrainCropVolume = volumesLogic.CloneVolume(slicer.mrmlScene, inputVolume, "Brain_Crop_"+str(now.strftime("%m%d%Y_%H%M%S")))
    
    cliParams={}
    cliParams["InputVolume"] = outputVolume.GetID()
    cliParams["OutputVolume"] = BrainCropVolume.GetID()
    slicer.cli.run(slicer.modules.castscalarvolume, None, cliParams, wait_for_completion=True)
    #slicer.mrmlScene.RemoveNode(slicer.util.getNode("Aux_Crop_"+str(now.strftime("%m%d%Y_%H%M%S"))))

    
    PreprocessesCT = volumesLogic.CloneVolume(slicer.mrmlScene, inputVolume, "Preprocessed_CT_"+str(now.strftime("%m%d%Y_%H%M%S")))
    cliParams = {'InputVolume': inputVolume.GetID(), 'MaskVolume': brainmask.GetID(), 'OutputVolume': PreprocessesCT.GetID()}
    cliNode = slicer.cli.run(slicer.modules.maskscalarvolume, None, cliParams, wait_for_completion=True)
    
    ############################################DESATIVEDFORTEST######################################################

    electrodeslabelimage = volumesLogic.CreateAndAddLabelVolume( slicer.mrmlScene, inputVolume, "electrodes" )
    #cliParams = {'inputVolume': inputVolume.GetID(), 'outputVolume': electrodeslabelimage.GetID(),'threshold' : Thresholdvalue,'datPath' : self.electrodefilepath}
    
    cliParams = {'inputVolume': PreprocessesCT.GetID(), 'outputVolume': electrodeslabelimage.GetID(),'threshold' : Thresholdvalue,'datPath' : self.electrodefilepath}
    
    cliNode = slicer.cli.run(slicer.modules.fastelectrodesegmentor, None, cliParams, wait_for_completion=True)
    #cliNode = slicer.cli.run(slicer.modules.fastelectrodesegmentorsquare, None, cliParams, wait_for_completion=True)
    #cliNode = slicer.cli.run(slicer.modules.fastelectrodesegmentorcircle, None, cliParams, wait_for_completion=True)

    #now = datetime.now()
    #nameaux= "BRAIN_"+str(now.strftime("%m%d%Y_%H%M%S"))
    #slicer.util.loadVolume(self.tmpOutputFile, {'name': nameaux})
    
    Brain3D = volumesLogic.CreateAndAddLabelVolume( slicer.mrmlScene, BrainCropVolume, "Brain3D" )
    cliParams = {'inputVolume': BrainCropVolume.GetID(), 'outputVolume': Brain3D.GetID()}
    cliNode = slicer.cli.run(slicer.modules.modifiedqentropysegmentation, None, cliParams, wait_for_completion=True)

        
    
    #Rendering Electrodes
    self.makeModel(Brain3D, "1,2,3", 5,"BRAIN")
    
    
    #self.makeModel(Brain3D, 1, 5,"CSF")
    #self.makeModel(Brain3D, 2, 5,"GM")
    #self.makeModel(Brain3D, 3, 5,"WM")
    
    #self.makeModel(electrodeslabelimage, 307, 5, "Electrodes")

    with open(self.electrodefilepath+"/electrodescordinatesFile.txt", 'r') as filehandle:
    
      columnarray = []
      rowarray = []
      slicenumarray = []
      electrodes = []
      for line in filehandle:
          currentPlace = line[:-1]
          
          if  len(currentPlace.strip()) != 0:
              aux =currentPlace.split(",")
              columnarray.append(int(aux[0]))
              rowarray.append(int(aux[1]))
              slicenumarray.append(int(aux[2]))
              
          if  len(currentPlace.strip()) == 0:
              targetarray = ([round(np.median(columnarray)),round(np.median(rowarray)),round(np.median(slicenumarray))])
              if (self.checkifelectrodeexists(electrodes,targetarray)==False):
                  electrodes.append(targetarray)            
              columnarray=[]
              rowarray=[]
              slicenumarray=[]
          
          #i = i+1
              #print("Blank line")
                  #AuxValues.append(currentPlace)
      #print(i)

    filehandle.close()

    nameelectrodesfiducial= "Electrodes_"+str(now.strftime("%m%d%Y_%H%M%S"))
    mlogic = slicer.modules.markups.logic()   
    

    fidNode = slicer.util.getNode(mlogic.AddNewFiducialNode(nameelectrodesfiducial))
    volumeNode = inputVolume

    electrodenum = 0

    for electrode in electrodes:
      electrodenum = electrodenum+1
      # Get position of highest voxel value
      point_Ijk = [electrode[0], electrode[1], electrode[2]]

      # Get physical coordinates from voxel coordinates
      volumeIjkToRas = vtk.vtkMatrix4x4()
      volumeNode.GetIJKToRASMatrix(volumeIjkToRas)
      point_VolumeRas = [0, 0, 0, 1]
      volumeIjkToRas.MultiplyPoint(np.append(point_Ijk,1.0), point_VolumeRas)

      # If volume node is transformed, apply that transform to get volume's RAS coordinates
      transformVolumeRasToRas = vtk.vtkGeneralTransform()
      slicer.vtkMRMLTransformNode.GetTransformBetweenNodes(volumeNode.GetParentTransformNode(), None, transformVolumeRasToRas)
      point_Ras = transformVolumeRasToRas.TransformPoint(point_VolumeRas[0:3])

      # Add a markup at the computed position and print its coordinates
      fidNode.AddFiducial(point_Ras[0], point_Ras[1], point_Ras[2], "Contact_"+str(electrodenum))











    ############################################DESATIVEDFORTEST######################################################
    
    #elastix.exe -f ../Tmp/fixedtmp.nii -m ../Tmp/movingtmp.nii -out transform -p Parameters_Rigid.txt
    #transformix.exe -tp transform/TransformParameters.0.txt -out result -in ../Tmp/movingtmp.nii
    #print(cliParams)
    #registration = volumesLogic.CloneVolume(slicer.mrmlScene, inputMRIVolume, 'Registration')
    #self.tmpRobexInputFile = "..\\Tmp\\result\\result.mhd"
    #self.tmpRobexOutputFile = "..\\Tmp\\outputtmp.nii"
    #self.tmpOutputFile = os.path.dirname(__file__)+"/Tmp/outputtmp.nii"
    
    #self.saveNode(slicer.util.getNode(nameaux), self.tmpOutputFile)
    #slicer.util.saveNode(slicer.util.getNode(nameaux), self.tmpOutputFile, properties)
    #slicer.mrmlScene.RemoveNode(slicer.util.getNode(nameaux))
    
    #slicer.mrmlScene.RemoveNode(inputVolume)
    #slicer.mrmlScene.RemoveNode(inputMRIVolume)

    #self.RobexExecutable = os.path.dirname(__file__)+"\ROBEX"

    #os.chdir(self.RobexExecutable)

    #os.system("ROBEX.exe "+self.tmpRobexInputFile+" "+self.tmpRobexOutputFile)
    
    
    #now = datetime.now()
    #nameaux= "BRAIN_"+str(now.strftime("%m%d%Y_%H%M%S"))
    #slicer.util.loadVolume(self.tmpRobexOutputFile, {'name': nameaux})
    
    
    #logging.info(os.path.dirname(os.path.abspath(__file__))+"\ROBEX\ROBEX.exe")

    #os.system("ROBEX.exe ..\Tmp\inputtmp.nii ..\Tmp\outputtmp.nii")
    
    # Capture screenshot
    #if enableScreenshots:
    #  self.takeScreenshot('InterfaceTest-Start','MyScreenshot',-1)

    logging.info('Processing completed')

    return True


class InterfaceTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_Interface1()

  def test_Interface1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import SampleData
    SampleData.downloadFromURL(
      nodeNames='FA',
      fileNames='FA.nrrd',
      uris='http://slicer.kitware.com/midas3/download?items=5767')
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = InterfaceLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
