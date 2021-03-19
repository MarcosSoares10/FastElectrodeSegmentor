#include "itkImageFileWriter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkPluginUtilities.h"
#include "FastElectrodeSegmentorCLP.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkImage.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <ctime>
#include "itkFlatStructuringElement.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkImageIterator.h"
#include "itkImageRegionIterator.h"


// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

int*** hypercube(int size, int radius) {

        int*** hcube = new int** [size];

        for (int i = 0; i < size; i++)
        {
            hcube[i] = new int* [size];
            for (int in = 0; in < size; in++)
            {
                hcube[i][in] = new int[size];
            }
        }


        int x0 = int(size / 2);
        int y0 = int(size / 2);
        int z0 = int(size / 2);

        for (int a = 0; a < size; a++)
        {
            for (int b = 0; b < size; b++)
            {
                for (int c = 0; c < size; c++)
                {
                    hcube[a][b][c] = 0;
                }
            }
        }

        for (int x = x0 - radius; x < x0 + radius + 1; x++)
        {
            for (int y = y0 - radius; y < y0 + radius + 1; y++)
            {
                for (int z = z0 - radius; z < z0 + radius + 1; z++)
                {
                    hcube[x][y][z] = 255;
                }
            }
        }
        return hcube;
    }

int*** electrode(int size, int radius) {

        int*** elec = new int** [size];

        for (int i = 0; i < size; i++)
        {
            elec[i] = new int* [size];
            for (int in = 0; in < size; in++)
            {
                elec[i][in] = new int[size];
            }
        }


        int x0 = int(size / 2);
        int y0 = int(size / 2);
        int z0 = int(size / 2);

        for (int a = 0; a < size; a++)
        {
            for (int b = 0; b < size; b++)
            {
                for (int c = 0; c < size; c++)
                {
                    elec[a][b][c] = 0;
                }
            }
        }

        for (int x = x0 - radius; x < x0 + radius + 1; x++)
        {
            for (int y = y0 - radius; y < y0 + radius + 1; y++)
            {
                for (int z = z0 - radius; z < z0 + radius + 1; z++)
                {
                    int deb = radius - int(x0 - x) - int(y0 - y) - int(z0 - z);
                    if (deb >= 0) {
                        elec[x][y][z] = 255;
                    }
                }
            }
        }
        return elec;
    }


int*** CreateMovingMinibatch(int size) {
        int*** movingcube = new int** [size];
        for (int i = 0; i < size; i++)
        {
            movingcube[i] = new int* [size];
            for (int in = 0; in < size; in++)
            {
                movingcube[i][in] = new int[size];
            }

        }

        for (int x = 0; x < size; x++)
        {
            for (int y = 0; y < size; y++)
            {
                for (int z = 0; z < size; z++)
                {
                    movingcube[x][y][z] = 0;
                }

            }

        }

        return movingcube;
    }


int*** CreateMainMatrix(int X, int Y, int Z) {
        int*** mainmatrix = new int** [X];

        for (int i = 0; i < X; i++)
        {
            mainmatrix[i] = new int* [Y];

            for (int j = 0; j < Y; j++)
                mainmatrix[i][j] = new int[Z];
        }


        for (int a = 0; a < X; a++)
        {
            for (int b = 0; b < Y; b++)
            {
                for (int c = 0; c < Z; c++)
                {
                    mainmatrix[a][b][c] = 0;
                }

            }

        }

        return mainmatrix;
    }

int CountPixels(int*** elementtocount, int size) {
        int countpixels = 0;
        for (int Sx = 0; Sx < size; Sx++)
        {
            for (int Sy = 0; Sy < size; Sy++)
            {
                for (int Sz = 0; Sz < size; Sz++)
                {
                    if (elementtocount[Sx][Sy][Sz] == 255) {
                        countpixels++;
                    }
                }

            }

        }
        return countpixels;
    }


bool Compare(int*** movingbatchelement, int*** hypercube, int size, int minimumpixels) {

            int totalzeros = (size * size * size) - CountPixels(hypercube, size);
            int targetpixelsintomask = 0;
            if (CountPixels(movingbatchelement, size) >= minimumpixels) {
                for (int Sx = 0; Sx < size; Sx++)
                {
                    for (int Sy = 0; Sy < size; Sy++)
                    {
                        for (int Sz = 0; Sz < size; Sz++)
                        {
                            if (hypercube[Sx][Sy][Sz] == 0 && movingbatchelement[Sx][Sy][Sz] == 0) {
                                targetpixelsintomask++;
                            }
                        }

                    }

                }
            }
            if (targetpixelsintomask == totalzeros)
            {
                return true;
            }
            else {
                return false;
            }

        }

int*** PopulateMinibatch(int*** movingminibatch, int*** mainmatrix, int sizehypercube, int z, int x, int y){
    for (int Sx = 0; Sx < sizehypercube; Sx++) {
        for (int Sy = 0; Sy < sizehypercube; Sy++) {
            for (int Sz = 0; Sz < sizehypercube; Sz++){
                 movingminibatch[Sx][Sy][Sz] = mainmatrix[Sx + z][Sy + x][Sz + y];
                                                       }
                                                   }
                                                }
    return movingminibatch;
}

void deleteauxMatrix(int*** matrix, int X, int Y) {
            // deallocate memory
            for (int i = 0; i < X; i++)
            {
                for (int j = 0; j < Y; j++)
                    delete[] matrix[i][j];

                delete[] matrix[i];
            }

            delete[] matrix;
        }


template <typename TPixel>
int DoIt( int argc, char * argv[], TPixel )
{
  PARSE_ARGS;

//  std::clock_t starttime;
//  double duration;

  typedef TPixel InputPixelType;
  typedef TPixel OutputPixelType;
  typedef TPixel RescaledPixelType;
  typedef TPixel ThresholdedPixelType;
  typedef TPixel DilatedPixelType;
  typedef TPixel ClosingPixelType;



  const unsigned int Dimension = 3;

  typedef itk::Image<InputPixelType, Dimension> InputImageType;
  typedef itk::Image<OutputPixelType, Dimension> OutputImageType;
  typedef itk::Image<ThresholdedPixelType, Dimension> ThresholdedImageType;
  typedef itk::Image<RescaledPixelType, Dimension> RescaledImageType;
  typedef itk::Image<DilatedPixelType, Dimension> DilatedImageType;
  typedef itk::Image<ClosingPixelType, Dimension> ClosingImageType;
  typedef itk::ImageFileReader<InputImageType>  ReaderType;


  std::ofstream electrodescordinatesFile;
  electrodescordinatesFile.open(datPath+"/electrodescordinatesFile.txt");

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputVolume.c_str());
  reader->Update();
  std::cout << "Rescale Image" << std::endl;
  typedef itk::RescaleIntensityImageFilter<InputImageType, RescaledImageType> RescaleType;
  typename RescaleType::Pointer rescaler = RescaleType::New();
  rescaler->SetInput(reader->GetOutput());
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);
  rescaler->Update();

  typename RescaledImageType::Pointer res = rescaler->GetOutput();
  typename RescaledImageType::RegionType region = res->GetLargestPossibleRegion();
  typename RescaledImageType::SizeType size = region.GetSize();
  int width = size[1];
  int height = size[0];
  int depth = size[2];
  typename RescaledImageType::IndexType start;
  start.Fill(0);
  res->SetRegions(region);
  region.SetIndex(start);
  region.SetSize(size);
  res->Allocate();

  typename OutputImageType::Pointer out = OutputImageType::New();
  out->CopyInformation(rescaler->GetOutput());
  out->SetRegions(rescaler->GetOutput()->GetRequestedRegion());
  typename OutputImageType::RegionType regionoutimg = out->GetLargestPossibleRegion();
  typename OutputImageType::SizeType sizeout = regionoutimg.GetSize();
  typename OutputImageType::IndexType startout;
  startout.Fill(0);
  regionoutimg.SetIndex(startout);
  regionoutimg.SetSize(sizeout);
  out->Allocate(); //Allocate memory for image pixel data
  out->FillBuffer(0.0);

  std::cout << "Loop Into Image" << std::endl;
  int sizehypercube = 9;
  int radiushypercube = 2;
  int*** cube = hypercube(sizehypercube, radiushypercube);
  //int*** sphere = electrode(sizehypercube, radiushypercube);
  int*** movingminibatch = CreateMovingMinibatch(sizehypercube);
  //int*** auxminibatch = CreateMovingMinibatch(sizehypercube);
  int*** mainmatrix = CreateMainMatrix(depth, height, width);
  int minimunpixelsacceptable = 5;



  TPixel level=0;
          for (int x = 0; x < depth; x++) {
              for (int y = 0; y < height ; y++) {
                  for (int z = 0; z < width; z++) {
                      start[0] = z;
                      start[1] = y;
                      start[2] = x;
                      level = res->GetPixel(start);
                      if (int(level) >= threshold) {
                      mainmatrix[x][y][z] = 255;
                      }else{ mainmatrix[x][y][z] = 0; }
                  }
              }
          }

 int reducesearch = 10;

 for (int z = 0; z < depth - sizehypercube; z++) {
     for (int x = reducesearch; x < height  - sizehypercube - reducesearch; x++) {
          for (int y = reducesearch; y < width - sizehypercube - reducesearch; y++) {
              /*ToDo: Testing improviments*/
              movingminibatch = PopulateMinibatch(movingminibatch, mainmatrix, sizehypercube, z, x, y);
                        if (Compare(movingminibatch, cube, sizehypercube, minimunpixelsacceptable) == true) {
                              for (int a = 0; a < sizehypercube; a++) {
                                  for (int b = 0; b < sizehypercube; b++) {
                                      for (int c = 0; c < sizehypercube; c++) {
                                          startout[0] = c + y;
                                          startout[1] = b + x;
                                          startout[2] = a + z;
                                          if (movingminibatch[a][b][c] == 255) {
                                              movingminibatch[a][b][c] = 307;
                                              electrodescordinatesFile << (c + y) << "," << (b + x) << "," << (a + z) << "\n";
                                                                               }
                                          out->SetPixel(startout, movingminibatch[a][b][c]);
                                                                                }
                                                                          }
                                                                      }
                                  electrodescordinatesFile << "\n";
                                                                                                            }
                                                                                           }
                                                                                   }
                                                  }
          electrodescordinatesFile.close();

          deleteauxMatrix(mainmatrix, depth, height);
          deleteauxMatrix(cube, sizehypercube, sizehypercube);
          deleteauxMatrix(movingminibatch, sizehypercube, sizehypercube);

          typedef itk::ImageFileWriter<OutputImageType> WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetFileName(outputVolume.c_str());
          writer->SetInput(out); //out rescaler->GetOutput()
          writer->SetUseCompression(1);
          writer->Update();
//          duration = (std::clock() - starttime) / (double)CLOCKS_PER_SEC;
//          std::cout << "Duration: " << duration << '\n';




  return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(inputVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0) );
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<signed char>(0) );
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0) );
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0) );
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cerr << "Unknown input image pixel component type: ";
        std::cerr << itk::ImageIOBase::GetComponentTypeAsString( componentType );
        std::cerr << std::endl;
        return EXIT_FAILURE;
        break;
      }
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
