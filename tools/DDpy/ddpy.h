/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ddpy_h_
#define model_ddpy_h_

#define _MODEL_NON_SINGULAR_DD_ 1 // 0 classical theory, 1 Cai's regularization method, 2 Lazar's regularization method

#include <string>
#include <tuple>
#include <list>
#include <map>
#include <stdlib.h> // EXIT_SUCCESS, EXIT_FAILURE

#include <filesystem>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
//#include <pybind11/smart_holder.h>

#include <Polycrystal.h>
//#include <SimplicialMesh.h>
//#include <DDconfigIO.h>
//#include <Grain.h>
#include <DefectiveCrystal.h>
#include <DislocationDynamicsBase.h>
//#include <DislocationSegment.h>
//#include <DislocationDynamicsModule.h>
#include <UniformExternalLoadController.h>
#include <IDreader.h>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <MicrostructureGeneratorBaseInMem.h>
#include <MicrostructureGeneratorInMem.h>
#include <MicrostructureSpecification.h>

namespace py = pybind11;
namespace ddpy
{
class DDInterface
{
   private:
      typedef model::DefectiveCrystal<3,0> DefectiveCrystalType;
      typedef model::DislocationDynamicsBase<3> DislocationDynamicsBaseType;
      std::string dddFolderPath;
      std::unique_ptr<DislocationDynamicsBaseType> ddBase;
      std::unique_ptr<DefectiveCrystalType> DC;

      // simulation box parameters required before SimplicialMesh
      std::vector<double> boxBounds; // xlo xhi ylo yhi zlo zhi
      double boxSkew;
      double boxVolume;  // [\AA^{2}]
      double burgersMagnitude;
      int solidSolutionNoiseMode;
      int stackingFaultNoiseMode;
      std::vector<std::string> acceptableLattices;//({"bcc","fcc"});
      std::vector<std::string> acceptableMaterials;//({"Cu","Fe_320"});

      bool debugFlag;

      void resetStaticIDs();
      void readddBase();

   public:
      typedef Eigen::Matrix<double,3,1> VectorDim;
      typedef Eigen::Matrix<double,3,3> MatrixDim;
      DDInterface( const std::string& dddFolderPathIn):
          dddFolderPath( dddFolderPathIn)
          , boxBounds( std::vector<double>({0,100,0,100,0,100}))
          , boxSkew( 0.0)
          , boxVolume( 0.0)
          , burgersMagnitude( -1.0)
          , solidSolutionNoiseMode( 0) //# 0=no noise, 1= read noise, 2=compute noise
          , stackingFaultNoiseMode( 0)
          , acceptableLattices( std::vector<std::string>({"bcc","fcc"}))
          , acceptableMaterials( std::vector<std::string>({"Fe_320","Cu"}))
          , debugFlag( true)
      {
         resetStaticIDs();
      };

      std::list<std::shared_ptr<model::MicrostructureSpecification>> microstructureSpecifications;

      //std::list<std::tuple< size_t, size_t, double>>
      std::map<std::tuple< size_t, size_t>, double>
         getResolvedShearStresses();
      //std::list<std::tuple< size_t, size_t, double>>
      std::map<std::tuple< size_t, size_t>, double>
         getResolvedShearStrains();
      //std::list<std::tuple< size_t, size_t, double>>
      std::map<std::tuple< size_t, size_t>, double>
         getPlasticStrains();
      //std::map<size_t, double>
      std::map<std::tuple<size_t,size_t>, double>
         getDensityPerSlipSystem();
      ////std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
      //std::map< std::pair<size_t,size_t>, VectorDim>
      //   getSlipSystemNormals() const;
      //std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
      //std::map< std::pair<size_t,size_t>, Eigen::Matrix<double,3,1> >
      //py::dict
      //   getSlipSystemBurgersVectors() const;

      std::string getFolderPath(){ return dddFolderPath;}
      void setEndingStep( const long int& endingStep);
      size_t getCurrentStep();
      void setCurrentStep( const long int& step);
      void setOutputFrequency( const long int& delta);
      void setBoxBounds(
            const double& xlo, // [m]
            const double& xhi,
            const double& ylo,
            const double& yhi,
            const double& zlo,
            const double& zhi
            );
      void runGlideSteps( size_t Nsteps);
      double getBurgersMagnitude();
      void readBurgersMagnitude( const std::string& materialPath);
      void readDefectiveCrystal();
      void regenerateMicrostructure();
      void writeConfigToTxt();
      void clearMicrostructureSpecifications();
      void specifyDipoles(
         const std::string& tag,
         const std::vector<int>& slipSystemIDs,
         const std::vector<int>& exitFaceIDs,
         const pybind11::array_t<double,
                  pybind11::array::c_style | pybind11::array::forcecast>&
                  pointsIn,
         const std::vector<double>& heights,
         const std::vector<int>& nodes,
         const std::vector<double>& glideSteps
         );

      void specifyLoops(
         const std::string& tag,
         const std::vector<int>& slipSystemIDs,
         const std::vector<double>& loopRadii,
         const std::vector<long int>& loopSegmentCounts, //periodicLoopSides
         const pybind11::array_t<double,
                  pybind11::array::c_style | pybind11::array::forcecast>&
                  loopCenters
         );
      void specifyLoopDensity(
               const std::string& tag,
               const double& periodicLoopTargetDensityIn,
               const long int& periodicLoopSegmentCountIn,
               const double& periodicLoopRadiusDistributionMeanIn,
               const double& periodicLoopRadiusDistributionStdIn
         );
      void specifyLoopDensitiesPerSlipSystem(
               const std::string& tag,
               const std::map<int, double>& loopDensitiesPerSlipSystemIn,
               const long int& periodicLoopSegmentCountIn,
               const double& periodicLoopRadiusDistributionMeanIn,
               const double& periodicLoopRadiusDistributionStdIn
         );
      void specifyPrismaticLoops(
         const std::string& tag,
         const std::vector<int>& slipSystemIDs,
         const std::vector<double>& prismaticLoopRadii,
         const pybind11::array_t<double,
                  pybind11::array::c_style | pybind11::array::forcecast>&
                  prismaticLoopCenters
         );
      void specifyPrismaticLoopDensity(
               const std::string& tag,
               const double& prismaticLoopTargetDensityIn,
               const double& prismaticLoopRadiusDistributionMeanIn,
               const double& prismaticLoopRadiusDistributionStdIn
         );
      void specifyPrismaticLoopDensitiesPerSlipSystem(
               const std::string& tag,
               const std::map<int, double>& prismaticLoopDensitiesPerSlipSystemIn,
               const double& prismaticLoopRadiusDistributionMeanIn,
               const double& prismaticLoopRadiusDistributionStdIn
         );
      void specifyDipoleDensity(
               const std::string& tag,
               const double& periodicDipoleTargetDensityIn
         );


      void setOutputPath( const std::string& outputPath);

      void generateMicrostructure();
      void regeneratePolycrystalFile(
            const py::array_t<double,
               py::array::c_style | py::array::forcecast>
               c2g,
               const std::string& lattice,
               const std::string& material,
               const std::string& meshFilePath
            );

      //template <>
      //void setExternalLoad<model::UniformExternalLoadController<typename DefectiveCrystalType>>(
      void setExternalLoad(
         std::optional< const pybind11::array_t<double,
            pybind11::array::c_style | pybind11::array::forcecast>>&
               ExternalStress0In, // 3x3 matrix
         std::optional< const pybind11::array_t<double,
            pybind11::array::c_style | pybind11::array::forcecast>>&
               ExternalStressRateIn, // 3x3 matrix
         std::optional< const pybind11::array_t<double,
            pybind11::array::c_style | pybind11::array::forcecast>>&
               ExternalStrain0In, // 3x3 matrix
         std::optional< const pybind11::array_t<double,
            pybind11::array::c_style | pybind11::array::forcecast>>&
               ExternalStrainRateIn, // 3x3 matrix
         std::optional< const pybind11::array_t<double,
            pybind11::array::c_style | pybind11::array::forcecast>>&
               MachineStiffnessRatio // Voigt format 11 22 33 12 23 13
            );
      //~DDInterface()
      //{
      //   return;
      //}

};

//   //std::map< size_t, VectorDim> printSlipSystemNormals()
//   //std::map< std::pair<size_t,size_t>, model::ReciprocalLatticeVector<3>>
//   // printSlipSystemNormals()
//   //std::map< std::pair<size_t,size_t>, VectorDim> printSlipSystemNormals()
//   //std::map< std::pair<size_t,size_t>, std::string> printSlipSystemNormals()



//}; // class DefectiveCrystalInterface

} // namespace ddpy

#endif
