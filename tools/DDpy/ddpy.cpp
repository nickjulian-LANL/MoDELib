/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2023 by Nicholas H. Julian <njulian@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ddpy_cpp_
#define model_ddpy_cpp_

#include <ddpy.h>

std::map<std::tuple< size_t, size_t>, double>
ddpy::DDInterface::getResolvedShearStresses()
{
   //Eigen::Matrix< double, 3, 3> stress;
   MatrixDim stress;
   std::map<std::tuple<size_t,size_t>,double> rss;
   if ( DC == nullptr)
   {
      std::cout << "error: getResolvedShearStresses(), "
        << " DefectiveCrystal not yet initialized" << std::endl;
      return rss;
   }
   if ( DC->externalLoadController != nullptr)
   {
      stress = DC->externalLoadController->stress( VectorDim::Zero());
   }
   else
   {
      std::cout << "error: getResolvedShearStresses(), "
        << " DC->externalLoadController isn't instantiated" << std::endl;
      return rss;
   }

   for ( const auto& grain : DC->poly.grains)
   {
      rss.clear();
      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
      { // loop over slip system
         std::tuple<size_t,size_t> key( grain.second.grainID, ss->sID);
         rss[ key] = (stress * (ss->unitNormal)).dot(ss->unitSlip);
      } // loop over slip system
   }
   return rss;
}

std::map<std::tuple< size_t, size_t>, double>
ddpy::DDInterface::getResolvedShearStrains()
{
   MatrixDim strain;
   std::map<std::tuple< size_t, size_t>, double> rss;
   if ( DC == nullptr)
   {
      std::cout << "error: getResolvedShearStrains(), "
        << " DefectiveCrystal not yet initialized" << std::endl;
      return rss;
   }
   if ( DC->externalLoadController != nullptr)
   {
      strain = DC->externalLoadController->strain( VectorDim::Zero());
   }
   else
   {
      std::cout << "error: getResolvedShearStrains(), "
        << " DC->externalLoadController == nullptr" << std::endl;
      return rss;
   }
   for ( const auto& grain : DC->poly.grains)
   {
      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
      { // loop over slip system
         std::tuple<size_t,size_t> key( grain.second.grainID, ss->sID);
         rss[ key] = (strain * (ss->unitNormal)).dot(ss->unitSlip);
      } // loop over slip system
   }
   return rss;
}

std::map<std::tuple<size_t,size_t>, double>
ddpy::DDInterface::getDensityPerSlipSystem()
{
   std::map<std::tuple<size_t,size_t>, double> densityPerSlipSystem;
   //std::map<size_t, double> densityPerSlipSystem;
   //densityPerSlipSystem[1] = 2.0; // debug
   //densityPerSlipSystem[3] = 0.1; // debug

   size_t ssID;
   size_t grainID;
   for ( const auto& loop : DC->DN->loops())
   {
      if ( loop.second.lock()->slipSystem() != nullptr)
      {
         grainID = loop.second.lock()->grain.grainID;
         ssID = loop.second.lock()->slipSystem()->sID;
         std::tuple<size_t,size_t> key( grainID, ssID);
         for ( const auto& loopLink : loop.second.lock()->loopLinks())
         {
            if ( loopLink->networkLink())
            {
               if ( ! loopLink->networkLink()->hasZeroBurgers())
               {
                  if (( !loopLink->networkLink()->isBoundarySegment())
                        &&(!loopLink->networkLink()->isGrainBoundarySegment())
                     )
                  {
                     // if ssID isn't in the map yet, instantiate w/value 0
                     densityPerSlipSystem.try_emplace( key, 0.0);
                     densityPerSlipSystem[ key ]
                        += loopLink->networkLink()->chord().norm()
                           /(
                            loopLink->networkLink()->loopLinks().size()
                            * ddBase->mesh.volume()
                            * std::pow( ddBase->poly.b_SI, 2)
                            // mesh.volume lacks |b|^3
                            // chord().norm() probably lacks |b| factor
                            // leaving two |b| factors in denominator
                            );
                  }
               }
            }
         }
      }
      else
      {
         //std::cout << "error: loop.second.lock()->slipSystem(): "
         //   << loop.second.lock()->slipSystem() << std::endl;
      }
   }
   return densityPerSlipSystem;
}

//std::list<std::tuple< size_t, size_t, double>>
std::map<std::tuple< size_t, size_t>, double>
ddpy::DDInterface::getPlasticStrains()
{
   //std::list<std::tuple< size_t, size_t, double>> plasticStrains;
   std::map<std::tuple< size_t, size_t>, double> plasticStrains;
   if ( DC == nullptr)
   {
      std::cout << "error: getPlasticStrains(), "
        << " DefectiveCrystal not yet initialized" << std::endl;
      return plasticStrains;
   }
   std::map<std::pair<int,int>,double> sspd(
         DC->DN->slipSystemPlasticDistortion());
   // [[grain ID, slip system ID], strain value in that slip system]

   //std::cout << "sspd.size: " << sspd.size() << std::endl;

   // copy into list of tuples to be returned
   for (const auto& itr : sspd)
   {
       plasticStrains[
             std::tuple( itr.first.first, itr.first.second)
         ] = itr.second;
             //.push_back(std::tuple( itr.first.first, itr.first.second, itr.second));
   }

   //for ( const auto& itr : plasticStrains)
   //{
   //   std::cout << "(grain,slipSystem,strain): " << "("
   //      << std::get<0>(itr) << ","
   //      << std::get<1>(itr) << ","
   //      << std::get<2>(itr) << ")" << std::endl;
   //}

   return plasticStrains;
}

//std::map<std::pair< size_t, size_t>, ddpy::DDInterface::VectorDim>
//   ddpy::DDInterface::getSlipSystemNormals() const
//{
//   size_t grainCount = 0;
//   size_t slipSystemCount = 0;
//   //std::map< std::pair<size_t,size_t>, model::ReciprocalLatticeVector<3>> normals;
//   std::map< std::pair<size_t,size_t>, VectorDim> normals;
//   //std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
//      normals;
//   if ( DC == NULL)
//   {
//      std::cout << "error: getSlipSystemNormals(), "
//        << " DefectiveCrystal not yet initialized" << std::endl;
//      return normals;
//   }
//   //std::map< std::pair<size_t,size_t>, std::string> normals;
//   // ((grain number, slip system number), plane normal)
//   for ( const auto& grain : DC->poly.grains)
//   {
//      std::cout << "grain " << grainCount << std::endl;
//      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
//      { // loop over slip system
//         //std::cout << "slip system " << slipSystemCount
//         //   << " plane normal:" << std::endl
//         //   << " (" << std::endl << ss->unitNormal
//         //   << ")" << std::endl;
//         normals.push_back(
//               std::tuple( grainCount, slipSystemCount, ss->unitNormal)
//               );
//         ++slipSystemCount;
//      }
//      ++grainCount;
//   }
//   for ( const auto& nn : normals)
//   {
//      std::cout << "grain " << std::get<0>( nn)
//         << " slip system " << std::get<1>( nn)
//         << std::endl << std::get<2>( nn) << std::endl;
//   }
//   return normals;
//}

//std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
//std::map< std::pair<size_t,size_t>, Eigen::Matrix<double,3,1> >
//   ddpy::DDInterface::getSlipSystemBurgersVectors() const
//{
//   size_t grainCount = 0;
//   size_t slipSystemCount = 0;
//   //std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
//   std::map< std::pair<size_t,size_t>, Eigen::Matrix<double,3,1>>
//      burgersVectors;
//   if ( DC == NULL)
//   {
//      std::cout << "error: getSlipSystemBurgersVectors(), "
//        << " DefectiveCrystal not yet initialized" << std::endl;
//      return burgersVectors;
//   }
//   for ( const auto& grain : DC->poly.grains)
//   {
//      ++grainCount;
//      std::cout << "grain " << grainCount << std::endl;
//      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
//      { // loop over slip system
//         ++slipSystemCount;
//         burgersVectors.emplace(
//               std::make_pair(
//                  std::make_pair( grainCount, slipSystemCount),
//                  ss->unitSlip)
//               );
//         //std::cout << "slip system " << slipSystemCount
//         //   << " burgers vector:" << std::endl
//         //   << " (" << std::endl << ss->unitSlip
//         //   << ")" << std::endl;
//      }
//   }
//   return burgersVectors;
//}

double ddpy::DDInterface::getBurgersMagnitude()
{
   if ( DC == nullptr)
   {
      std::cout << "error, getBurgersMagnitude: "
         << " DefectiveCrystal not yet instantiated." << std::endl;
      return -1;
   }
   if ( DC->DN == nullptr)
   {
      std::cout << "error, getBurgersMagnitude: "
         << " DefectiveCrystal::DefectiveNetwork not yet instantiated."
         << std::endl;
      return -1;
   }
   return DC->DN->poly.b_SI/1e-10;
}

void ddpy::DDInterface::readBurgersMagnitude(
      const std::string& materialPath)
{
   burgersMagnitude = model::TextFileParser( materialPath).readScalar<double>("b_SI",true);
   std::cout << "burgersMagnitude: " << burgersMagnitude << std::endl; // debug
   return;
}


size_t ddpy::DDInterface::getCurrentStep()
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot getCurrentStep until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return 0;
   }
   return DC->simulationParameters.runID;
}

void ddpy::DDInterface::setCurrentStep( const long int& step)
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot setCurrentStep until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }

   DC->simulationParameters.runID = step;
   DC->externalLoadController = DC->getExternalLoadController(
            DC->simulationParameters, *DC, step
            );
   DC->simulationParameters.manageRestart();
   return;
}

//template <>
//void ddpy::DDInterface::setExternalLoad<model::UniformExternalLoadController<typename DefectiveCrystalType>>(
void ddpy::DDInterface::setExternalLoad(
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
            MachineStiffnessRatioIn // Voigt format 11 22 33 12 23 13
      )
{
   if ( ddBase == nullptr) readddBase();
   if ( DC == nullptr)
   {
      std::cout << "error: cannot assign external load until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }

   if ( DC->externalLoadController == nullptr)
   {
      // instantiate a load controller by reading a file
      //  (e.g. uniformExternalLoadController.txt)
      DC->externalLoadController = DC->getExternalLoadController(
              DC->simulationParameters,
              *DC,
              DC->simulationParameters.runID
              );
   }

   if ( DC->simulationParameters.externalLoadControllerName
         != "UniformExternalLoadController")
   {
      std::cout << "error ddpy::DDInterface::setExternalLoad(): "
        << "DC->simulationParameters.externalLoadControllerName: "
        << DC->simulationParameters.externalLoadControllerName
        << " != UniformExternalLoadController ."
        << "Only UniformExternalLoadController is supported by "
        << "ddpy::DDInterface" << std::endl;
      return;
   }

   // Adapting calculations from UniformExternalLoadController constructor

   // read and assign ExternalStressRate if it was specified
   if ( ExternalStressRateIn.has_value())
   {
      auto ExternalStressRateInBuf = ExternalStressRateIn.value().unchecked<2>();
      for ( ssize_t ii=0; ii < ExternalStressRateIn.value().shape()[0]; ++ii)
         for ( ssize_t jj=0; jj < ExternalStressRateIn.value().shape()[1]; ++jj)
         {
            DC->externalLoadController->ExternalStressRate( ii, jj)
               = ExternalStressRateInBuf( ii, jj);
         }
      assert((
         DC->externalLoadController->ExternalStressRate
            - DC->externalLoadController->ExternalStressRate.transpose()
             ).norm()<DBL_EPSILON && "ExternalStressRate is not symmetric."
            );
   }

   // read and assign ExternalStrainRate if it was specified
   if ( ExternalStrainRateIn.has_value())
   {
      auto ExternalStrainRateInBuf = ExternalStrainRateIn.value().unchecked<2>();
      for ( ssize_t ii=0; ii < ExternalStrainRateIn.value().shape()[0]; ++ii)
         for ( ssize_t jj=0; jj < ExternalStrainRateIn.value().shape()[1]; ++jj)
         {
            DC->externalLoadController->ExternalStrainRate( ii, jj)
               = ExternalStrainRateInBuf( ii, jj);
         }
      assert((
               DC->externalLoadController->ExternalStrainRate
                - DC->externalLoadController->ExternalStrainRate.transpose()
              ).norm()
              < DBL_EPSILON && "ExternalStrainRate is not symmetric."
            );
   }

   // read and assign MachineStiffnessRatio if it was specified
   if ( MachineStiffnessRatioIn.has_value())
   {
      auto MachineStiffnessRatioInBuf = MachineStiffnessRatioIn.value().unchecked<1>();
      for ( ssize_t ii=0; ii < MachineStiffnessRatioIn.value().shape()[0]; ++ii)
      {
         DC->externalLoadController->MachineStiffnessRatio( ii)
               = MachineStiffnessRatioInBuf( ii);
      }
      DC->externalLoadController->MachineStiffnessRatio.block(0,0,1,DefectiveCrystalType::dim)
         =DC->externalLoadController->MachineStiffnessRatio.block(0,0,1,DefectiveCrystalType::dim)
            * (2+DC->externalLoadController->lambda);

      static constexpr int voigtSize = DefectiveCrystalType::dim*(DefectiveCrystalType::dim+1)/2;
      Eigen::Matrix<double,voigtSize,voigtSize>  Cinv=Eigen::Matrix<double,voigtSize,voigtSize>::Identity();
      double nu( DC->DN->poly.nu);
      double nu_use( DC->externalLoadController->nu_use);
      Cinv.block(0,0,DefectiveCrystalType::dim,DefectiveCrystalType::dim)
         <<
             (nu_use)/nu, -(nu_use),       -(nu_use),
            -(nu_use),   (nu_use)/nu,     -(nu_use),
            -(nu_use),   -(nu_use),       (nu_use)/nu;

      Eigen::Matrix<double,voigtSize,voigtSize>  machinestiffness
         =
         DC->externalLoadController->MachineStiffnessRatio.asDiagonal();
      DC->externalLoadController->stressmultimachinestiffness
         = (Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness*Cinv).inverse();
      DC->externalLoadController->strainmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness*Cinv).inverse()*machinestiffness;
   }

   // ExternalStrain is either assigned the value of ExternalStrain0 or
   //  left untouched.
   if ( ExternalStrain0In.has_value())
   {
      auto ExternalStrain0InBuf = ExternalStrain0In.value().unchecked<2>();
      for ( ssize_t ii=0; ii < ExternalStrain0In.value().shape()[0]; ++ii)
         for ( ssize_t jj=0; jj < ExternalStrain0In.value().shape()[1]; ++jj)
         {
            DC->externalLoadController->ExternalStrain( ii, jj)
               = ExternalStrain0InBuf( ii, jj);
         }
      assert(
               ( DC->externalLoadController->ExternalStrain
                  - DC->externalLoadController->ExternalStrain.transpose()
               ).norm()
               < DBL_EPSILON && "ExternalStrain0 is not symmetric."
            );
   }

   // If ExternalStress0In is specified, then assign it to ExternalStress
   //  and calculate plasticStrain as the difference between ExternalStrain
   //  and elasticstrain(ExternalStress,nu_use).
   // The following mimicks how ExternalStress and plasticStrain are
   //  assigned when reading the state of the external load from F/F_0.txt.
   if ( ExternalStress0In.has_value())
   {
      // TODO: this block has an error
      auto ExternalStress0InBuf = ExternalStress0In.value().unchecked<2>();
      for ( ssize_t ii=0; ii < ExternalStress0In.value().shape()[0]; ++ii)
         for ( ssize_t jj=0; jj < ExternalStress0In.value().shape()[1]; ++jj)
         {
            DC->externalLoadController->ExternalStress( ii, jj)
               = ExternalStress0InBuf( ii, jj);
         }

      MatrixDim pdr( DC->DN->plasticDistortion());
      DC->externalLoadController->plasticStrain
            = ( pdr +pdr.transpose())*0.5;
      MatrixDim dstrain(
         DC->externalLoadController->ExternalStrain // possibly new
         - DC->externalLoadController->plasticStrain
         );

      DC->externalLoadController->ExternalStress
         = DC->externalLoadController->stressconsidermachinestiffness(
               dstrain,
               DC->externalLoadController->ExternalStress
               );

      assert(
               (
                  DC->externalLoadController->ExternalStress
                   - DC->externalLoadController->ExternalStress.transpose()
               ).norm()
               < DBL_EPSILON && "ExternalStress0 is not symmetric."
            );
   }
   else
   {  // If ExternalStress0In is not specified
      //  then mimick assignments used by the load constructor when
      //  F/F_0.txt cannot be read.

      if (
            ( ExternalStrain0In.has_value())
            || ( MachineStiffnessRatioIn.has_value())
         )
      {  // If ExternalStrain or stiffness ratio were passed to this
         //  function, then the ExternalStress and plasticStrain would
         //  need to be updated.
         MatrixDim pdr( DC->DN->plasticDistortion());
         DC->externalLoadController->plasticStrain
            = ( pdr +pdr.transpose())*0.5;

         MatrixDim dstrain(
            DC->externalLoadController->ExternalStrain // possibly new
            - DC->externalLoadController->plasticStrain
            );

         MatrixDim S_stress( DC->externalLoadController->ExternalStress);
         DC->externalLoadController->ExternalStress
            = DC->externalLoadController->stressconsidermachinestiffness(
                  dstrain,
                  S_stress
                  );
      }
      // If none of ExternalStrain, stiffness ratio, or ExternalStress
      //  were passed to this function, the values previously assigned by
      //  the constructor should remain.
   }
   return;
}

void ddpy::DDInterface::setOutputFrequency( const long int& outputFrequency)
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot assign output frequency until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }
   if ( DC->DN == nullptr)
   {
      std::cout << "error: cannot assign output frequency until "
         << "DefectiveNetwork is instantiated" << std::endl;
      return;
   }

   int outputFrequencyInt = static_cast<int>( outputFrequency);

   if ( outputFrequencyInt <= 0)
   {
      std::cout << "error: static_cast<int>( outputFrequency) "
         << outputFrequencyInt << " <= 0" << std::endl;
      return;
   }

   DC->DN->outputFrequency = outputFrequencyInt;
   return;
}

void ddpy::DDInterface::setEndingStep( const long int& endingStep)
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot assign ending step until "
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }
   if ( endingStep <= DC->simulationParameters.Nsteps)
   {
      std::cout << "warning: assigning endingStep=" << endingStep
         << " <= currentStep:" << DC->simulationParameters.Nsteps
         << std::endl;
   }
   DC->simulationParameters.Nsteps = endingStep;
   return;
}

void ddpy::DDInterface::readddBase()
{
   //resetStaticIDs();
   if ( ddBase != nullptr) resetStaticIDs(); // TODO: is this necessary?
   ddBase =  std::unique_ptr<DislocationDynamicsBaseType>(
         new DislocationDynamicsBaseType( dddFolderPath)
         );
   return;
}

void ddpy::DDInterface::setBoxBounds(
      const double& xlo,
      const double& xhi,
      const double& ylo,
      const double& yhi,
      const double& zlo,
      const double& zhi
      )
{
   if ( burgersMagnitude <= 0)
   {
      std::cout << "error, setBoxBounds: burgersMagnitude <= 0\n"
         << "try calling readBurgersMagnitude( <material file path>) first"
         << std::endl;
      return;
   }
   boxBounds[0] = xlo * 1e-10 / burgersMagnitude;
   boxBounds[1] = xhi * 1e-10 / burgersMagnitude;
   boxBounds[2] = ylo * 1e-10 / burgersMagnitude;
   boxBounds[3] = yhi * 1e-10 / burgersMagnitude;
   boxBounds[4] = zlo * 1e-10 / burgersMagnitude;
   boxBounds[5] = zhi * 1e-10 / burgersMagnitude;
   boxVolume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo); // [\AA^{2}]

   return;
}

void ddpy::DDInterface::readDefectiveCrystal()
{
   //if ( DC != NULL) delete DC;
   if ( ddBase == nullptr) readddBase();
   DC = std::unique_ptr<DefectiveCrystalType>(
         new DefectiveCrystalType( *ddBase)
         );
   return;
}

void ddpy::DDInterface::generateMicrostructure()
{
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "generating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << std::endl; // debug

   if ( microstructureSpecifications.size() == 0)
   {
      std::cout << "error: generateMicrostructure() requires defects to be specified via respective functions" << std::endl;
      return;
   }
   for (const auto& spec : microstructureSpecifications) // debug
   { // debug
      std::cout << "spec->tag: " << spec->tag << std::endl; // debug
      std::cout << "spec->microstructureType: " // debug
         << spec->microstructureType<< std::endl; // debug
      std::cout << "spec->style: " << spec->style<< std::endl; // debug
   } // debug
   model::MicrostructureGeneratorInMem mg( *ddBase, microstructureSpecifications);
   return;
}


void ddpy::DDInterface::writeConfigToTxt()
{
   if ( DC == nullptr)
   {
      std::cout << "error: cannot write configuration to file until"
         << "DefectiveCrystal is instantiated" << std::endl;
      return;
   }
   if ( DC->DN == nullptr)
   {
      std::cout << "error: cannot write configuration to file until "
         << "DefectiveNetwork is instantiated" << std::endl;
      return;
   }
   std::cout << "attempting to write a configuration" << std::endl; // debug
   DC->DN->io().output( DC->simulationParameters.runID);
   return;
}

void ddpy::DDInterface::runGlideSteps( size_t Nsteps)
{
   if ( DC == nullptr) readDefectiveCrystal();

   setEndingStep( getCurrentStep() + Nsteps);
   DC->runGlideSteps();
   return;
}

//void ddpy::DDInterface::generateLoopDensity()
//void ddpy::DDInterface::generateDipoleDensity()
//void ddpy::DDInterface::generateDipole()
//{
//   // TODO:  emulate MicrostructureGenerator.cpp:57-108
//   // Can a periodicDipoleGenerator be instantiated without reading
//   //   microstructure parameters from a file? No.
//   // maybe create alternatives to the following which don't read from files:
//   //     MicrostructureGeneratorBase, PeriodicDipoleGenerator, ...
//   //
//   // PeriodicDipoleGenerator inherits from MicrostructureGeneratorBase
//   //
//   // MicrostructureGeneratorBase instantiates a TextFileParser and reads
//   //  type
//   //  style
//   //  tag
//   //
//   // PeriodicDipoleGenerator::generateDensity uses parser to read a file
//   //
//   // I don't want to instantiate a MicrostructureGenerator
//   //  because it will attempt to read the files listed in traitsIO.microstructureFile
//   //TODO: empty traitsIO.microstructureFile
//
//   // things to replace from TextFileParser parser found in periodicDipoleGenerator::generateIndividual:
//   const std::vector<int> periodicDipoleExitFaceIDs(this->parser.readArray<int>("periodicDipoleExitFaceIDs",true));
//   const Eigen::Matrix<double,Eigen::Dynamic,dim> periodicDipolePoints(this->parser.readMatrix<double>("periodicDipolePoints",periodicDipoleSlipSystemIDs.size(),dim,true));
//   const std::vector<double> periodicDipoleHeights(this->parser.readArray<double>("periodicDipoleHeights",true));
//   const std::vector<int> periodicDipoleNodes(this->parser.readArray<int>("periodicDipoleNodes",true));
//   const std::vector<double> periodicDipoleGlideSteps(this->parser.readArray<double>("periodicDipoleGlideSteps",true));
//
//   return;
//}

//void ddpy::DDInterface::regenerateMicrostructure()
//{
//   if ( ddBase == nullptr) readddBase();
//   //
//   // instantiate a MicrostructureGenerator
//   model::MicrostructureGeneratorInMem mg( *ddBase); // requires MicrostructureSpecifications
//   if ( debugFlag)
//      std::cout << "finished call to MicrostructureGeneratorInMem" << std::endl;
//
//   setCurrentStep(0); // DefectiveCrystal will use runID to read some things
//
//   // instantiate a DefectiveCrystalType DC( ddBase)
//   DC = std::unique_ptr<DefectiveCrystalType>(
//         new DefectiveCrystalType( *ddBase)
//         );
//   if ( debugFlag)
//      std::cout << "finished call to DefectiveCrystalType" << std::endl;
//
//   return;
//}

void ddpy::DDInterface::resetStaticIDs()
{
   model::StaticID<model::EshelbyInclusionBase<3>>::set_count(0);
   model::StaticID<model::Lattice<3>>::set_count(0);
   model::StaticID<model::DislocationMobilityBase>::set_count(0);
   model::StaticID<model::PeriodicPlaneNode<3>>::set_count(0);
   model::StaticID<model::PeriodicPlanePatch<3>>::set_count(0);
   model::StaticID<model::PeriodicPatchBoundary<3>>::set_count(0);
   //model::StaticID<model::Derived>::set_count(0);
   model::StaticID<model::SlipSystem>::set_count(0);
   //model::StaticID<model::LoopPathClipperNode>::set_count(0);
   model::StaticID<model::PlanarMeshFace<3>>::set_count(0);
   //model::StaticID<model::LoopPathClipperNode>::set_count(0);
   model::StaticID<model::NetworkNode<model::DislocationNode<3,0>>>::set_count(0);
   model::StaticID<model::NetworkLink<model::DislocationSegment<3,0>>>::set_count(0);
   model::StaticID<model::Loop<model::DislocationLoop<3,0>>>::set_count(0);
   model::StaticID<model::LoopNode<model::DislocationLoopNode<3,0>>>::set_count(0);
   model::StaticID<model::MeshPlane<3>>::set_count(0);
   model::StaticID<model::SimplexBase<3,0>>::set_count(0);
   model::StaticID<model::SimplexBase<3,1>>::set_count(0);
   model::StaticID<model::SimplexBase<3,2>>::set_count(0);
   model::StaticID<model::SimplexBase<3,3>>::set_count(0);
   //model::StaticID<model::SequentialOutputFile<prefix,auto>>::set_count(0);
   return;
}

void ddpy::DDInterface::specifyLoops(
         const std::string& tag,
         const std::vector<int>& periodicLoopSlipSystemIDsIn,
         const std::vector<double>& periodicLoopRadii,
         const std::vector<long int>& loopSegmentCountsIn,
         const pybind11::array_t<double,
                  pybind11::array::c_style | pybind11::array::forcecast>&
                  loopCentersIn
      )
{
   if ( ddBase == nullptr) readddBase();

   // copy Nx3 numpy array of points in 3-D to Eigen equivalents
   Eigen::Matrix<double, Eigen::Dynamic, 3> periodicLoopCenters;
   periodicLoopCenters.resize( loopCentersIn.shape()[0], loopCentersIn.shape()[1]);
   auto loopCentersInBuf = loopCentersIn.unchecked<2>();
   for ( ssize_t ii=0; ii < loopCentersIn.shape()[0]; ++ii)
      for ( ssize_t jj=0; jj < loopCentersIn.shape()[1]; ++jj)
         periodicLoopCenters( ii, jj) = loopCentersInBuf( ii, jj);

   // cast vectors of python integers to vector<int> (python int isn't int)
   std::vector<int> periodicLoopSlipSystemIDs( periodicLoopSlipSystemIDsIn.size(), 0);
   for ( size_t ii=0; ii < periodicLoopSlipSystemIDsIn.size(); ++ii)
      periodicLoopSlipSystemIDs[ ii] = static_cast<int>( periodicLoopSlipSystemIDsIn[ ii]);

   std::vector<long int> periodicLoopSides( loopSegmentCountsIn.size(), 0);
   for ( size_t ii=0; ii < loopSegmentCountsIn.size(); ++ii)
      periodicLoopSides[ ii] = static_cast<long int>( loopSegmentCountsIn[ ii]);

   // numpy loopRadii input is automatically seen as vector<double>

   // dummy instantiations
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   std::vector<double> prismaticLoopRadii;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PeriodicLoop", // type
               "individual", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, //periodicLoopTargetDensityIn, // double
               0, //periodicLoopSegmentCountIn, // long int
               0.0, //periodicLoopRadiusDistributionMeanIn // double
               0.0, //periodicLoopRadiusDistributionStdIn // double
               0.0, //periodicDipoleTargetDensityIn, // double
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, // prismaticLoopRadiiMean
               1.0 // prismaticLoopRadiiStd
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug

   return;
}

void ddpy::DDInterface::specifyLoopDensitiesPerSlipSystem(
         const std::string& tag,
         const std::map<int, double>& periodicLoopTargetDensitiesPerSlipSystemIn,
         const long int& periodicLoopSegmentCountIn,
         const double& periodicLoopRadiusDistributionMeanIn,
         const double& periodicLoopRadiusDistributionStdIn
      )
{
   if ( ddBase == nullptr) readddBase();
   // instantiate empty placeholders
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   std::vector<long int> loopSegmentCountsIn;
   //pybind11::array_t<double,
   //         pybind11::array::c_style | pybind11::array::forcecast>
   //         periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic, 3> periodicLoopCenters;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PeriodicLoop", // type
               "densitiesPerSlipSystem", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, // periodicLoopTargetDensityIn,
               periodicLoopSegmentCountIn,
               periodicLoopRadiusDistributionMeanIn,
               periodicLoopRadiusDistributionStdIn,
               0.0, //periodicDipoleTargetDensityIn,
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystemIn,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, // prismaticLoopRadiiMean
               1.0 // prismaticLoopRadiiStd
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug


   return;
}

void ddpy::DDInterface::specifyLoopDensity(
         const std::string& tag,
         const double& periodicLoopTargetDensityIn,
         const long int& periodicLoopSegmentCountIn,
         const double& periodicLoopRadiusDistributionMeanIn,
         const double& periodicLoopRadiusDistributionStdIn
      )
{
   if ( ddBase == nullptr) readddBase();

   // instantiate empty placeholders
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   std::vector<long int> loopSegmentCountsIn;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   //pybind11::array_t<double,
   //         pybind11::array::c_style | pybind11::array::forcecast>
   //         periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic, 3> periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PeriodicLoop", // type
               "density", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               periodicLoopTargetDensityIn,
               periodicLoopSegmentCountIn,
               periodicLoopRadiusDistributionMeanIn,
               periodicLoopRadiusDistributionStdIn,
               0.0, //periodicDipoleTargetDensityIn,
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, // prismaticLoopRadiiMean
               1.0 // prismaticLoopRadiiStd
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug

   return;
}
/**********************************************************************/
void ddpy::DDInterface::specifyPrismaticLoops(
         const std::string& tag,
         const std::vector<int>& prismaticLoopSlipSystemIDsIn,
         const std::vector<double>& prismaticLoopRadiiIn,
         const pybind11::array_t<double,
                  pybind11::array::c_style | pybind11::array::forcecast>&
                  prismaticLoopCentersIn
      )
{
   if ( ddBase == nullptr) readddBase();

   // copy Nx3 numpy array of points in 3-D to Eigen equivalents
   Eigen::Matrix<double, Eigen::Dynamic, 3> prismaticLoopCenters;
   prismaticLoopCenters.resize( prismaticLoopCentersIn.shape()[0], prismaticLoopCentersIn.shape()[1]);
   auto prismaticLoopCentersInBuf = prismaticLoopCentersIn.unchecked<2>();
   for ( ssize_t ii=0; ii < prismaticLoopCentersIn.shape()[0]; ++ii)
      for ( ssize_t jj=0; jj < prismaticLoopCentersIn.shape()[1]; ++jj)
         prismaticLoopCenters( ii, jj) = prismaticLoopCentersInBuf( ii, jj);

   // cast vectors of python integers to vector<int> (python int isn't int)
   std::vector<int> prismaticLoopSlipSystemIDs( prismaticLoopSlipSystemIDsIn.size(), 0);
   for ( size_t ii=0; ii < prismaticLoopSlipSystemIDsIn.size(); ++ii)
      prismaticLoopSlipSystemIDs[ ii] = static_cast<int>( prismaticLoopSlipSystemIDsIn[ ii]);

   // numpy loopRadii input is automatically seen as vector<double>

   // dummy instantiations
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<long int> periodicLoopSides;
   std::vector<double> periodicLoopRadii;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   Eigen::Matrix<double, Eigen::Dynamic,3> periodicLoopCenters;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PrismaticLoop", // type
               "individual", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, //periodicLoopTargetDensityIn, // double
               0, //periodicLoopSegmentCountIn, // long int
               0.0, //periodicLoopRadiusDistributionMeanIn // double
               0.0, //periodicLoopRadiusDistributionStdIn // double
               0.0, //periodicDipoleTargetDensityIn, // double
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadiiIn,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, // prismaticLoopRadiiMean
               1.0 // prismaticLoopRadiiStd
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug

   return;
}

void ddpy::DDInterface::specifyPrismaticLoopDensitiesPerSlipSystem(
         const std::string& tag,
         const std::map<int, double>& prismaticLoopTargetDensitiesPerSlipSystemIn,
         const double& prismaticLoopRadiiMeanIn,
         const double& prismaticLoopRadiiStdIn
      )
{
   if ( ddBase == nullptr) readddBase();
   // instantiate empty placeholders
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   std::vector<long int> loopSegmentCountsIn;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   //pybind11::array_t<double,
   //         pybind11::array::c_style | pybind11::array::forcecast>
   //         periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic, 3> periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic, 3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PrismaticLoop", // type
               "densitiesPerSlipSystem", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, // periodicLoopTargetDensityIn,
               0, //periodicLoopSegmentCountIn,
               1.0, //periodicLoopRadiusDistributionMeanIn,
               1.0, //periodicLoopRadiusDistributionStdIn,
               0.0, //periodicDipoleTargetDensityIn,
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopTargetDensitiesPerSlipSystemIn,
               prismaticLoopRadiiMeanIn,
               prismaticLoopRadiiStdIn
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug


   return;
}

void ddpy::DDInterface::specifyPrismaticLoopDensity(
         const std::string& tag,
         const double& prismaticLoopTargetDensityIn,
         const double& prismaticLoopRadiiMeanIn,
         const double& prismaticLoopRadiiStd
      )
{
   if ( ddBase == nullptr) readddBase();

   // instantiate empty placeholders
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> periodicDipoleExitFaceIDs;
   Eigen::Matrix< double, Eigen::Dynamic,3> periodicDipolePoints;
   std::vector<double> periodicDipoleHeights;
   std::vector<int> periodicDipoleNodes;
   std::vector<double> periodicDipoleGlideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   std::vector<long int> loopSegmentCountsIn;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   //pybind11::array_t<double,
   //         pybind11::array::c_style | pybind11::array::forcecast>
   //         periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic, 3> periodicLoopCenters;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
            new model::MicrostructureSpecification(
               "PrismaticLoop", // type
               "density", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs,
               periodicDipoleExitFaceIDs,
               periodicDipolePoints,
               periodicDipoleHeights,
               periodicDipoleNodes,
               periodicDipoleGlideSteps,
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, //periodicLoopTargetDensityIn,
               0, // periodicLoopSegmentCountIn,
               0.0, //periodicLoopRadiusDistributionMeanIn,
               0.0, //periodicLoopRadiusDistributionStdIn,
               0.0, //periodicDipoleTargetDensityIn,
               prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopTargetDensitiesPerSlipSystem,
               prismaticLoopRadiiMeanIn,
               prismaticLoopRadiiStd
         ));
   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug

   return;
}

/**********************************************************************/
void ddpy::DDInterface::specifyDipoleDensity(
               const std::string& tag,
               const double& periodicDipoleTargetDensityIn
      )
{
   if ( ddBase == nullptr) readddBase();

   // dummy instantiations. TODO: find another way to implement this class
   Eigen::Matrix<double, Eigen::Dynamic, 3> points;
   std::vector<double> heights;
   std::vector<int> periodicDipoleSlipSystemIDs;
   std::vector<int> exitFaceIDs;//( exitFaceIDsIn.size(), 0);
   std::vector<int> nodes;//( nodesIn.size(), 0);
   std::vector<double> glideSteps;
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   Eigen::Matrix<double,Eigen::Dynamic,3> periodicLoopCenters;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
         std::shared_ptr<model::MicrostructureSpecification>(
            //new model::PeriodicDipoleIndividualSpecification(
            new model::MicrostructureSpecification(
               "PeriodicDipole", // type
               "density", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDs, // std::vector<int>
               // dipoles
               exitFaceIDs, // std::vector<int>
               points, // Eigen::Matrix
               heights, // std::vector<double>
               nodes, // std::vector<int>
               glideSteps, // std::vector<double>
               // loops
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, //periodicLoopTargetDensityIn,
               0, //periodicLoopSegmentCountIn,
               0.0, //periodicLoopRadiusDistributionMeanIn
               0.0, //periodicLoopRadiusDistributionStdIn
               periodicDipoleTargetDensityIn,
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, //prismaticLoopRadiiMean // not used
               1.0 //prismaticLoopRadiiStd // not used
               )
            )
         );

   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug

   return;
}


void ddpy::DDInterface::specifyDipoles(
      const std::string& tag,
      const std::vector<int>& periodicDipoleSlipSystemIDsIn,
      const std::vector<int>& exitFaceIDsIn,
      const pybind11::array_t<double,
               pybind11::array::c_style | pybind11::array::forcecast>&
               pointsIn,
      const std::vector<double>& heights,
      const std::vector<int>& nodesIn,
      const std::vector<double>& glideSteps
      )
{
   if ( ddBase == nullptr) readddBase();

   std::cout << "pointsIn.shape(): " << pointsIn.shape() << std::endl; // debug
   std::cout << "periodicDipoleSlipSystemIDsIn.size(): " << periodicDipoleSlipSystemIDsIn.size() << std::endl; // debug
   std::cout << "exitFaceIDsIn.size(): " << exitFaceIDsIn.size() << std::endl; // debug

   Eigen::Matrix<double, Eigen::Dynamic, 3> points;
   points.resize( pointsIn.shape()[0], pointsIn.shape()[1]);
   auto pointsInBuf = pointsIn.unchecked<2>();
   for ( ssize_t ii=0; ii < pointsIn.shape()[0]; ++ii)
      for ( ssize_t jj=0; jj < pointsIn.shape()[1]; ++jj)
         points( ii, jj) = pointsInBuf( ii, jj);
   std::cout << "points.size(): " << points.size() << std::endl; // debug

   std::vector<int> periodicDipoleSlipSystemIDs( periodicDipoleSlipSystemIDsIn.size(), 0);
   for ( size_t ii=0; ii < periodicDipoleSlipSystemIDsIn.size(); ++ii)
      periodicDipoleSlipSystemIDs[ ii] = static_cast<int>( periodicDipoleSlipSystemIDsIn[ ii]);

   std::vector<int> exitFaceIDs( exitFaceIDsIn.size(), 0);
   for ( size_t ii=0; ii < exitFaceIDsIn.size(); ++ii)
      exitFaceIDs[ ii] = static_cast<int>( exitFaceIDsIn[ ii]);

   std::vector<int> nodes( nodesIn.size(), 0);
   for ( size_t ii=0; ii < nodesIn.size(); ++ii)
      nodes[ ii] = static_cast<int>( nodesIn[ ii]);

   // dummy instantiations. TODO: find another way to implement this class
   std::vector<int> periodicLoopSlipSystemIDs;
   std::vector<double> periodicLoopRadii;
   std::vector<long int> periodicLoopSides;
   Eigen::Matrix<double,Eigen::Dynamic,3> periodicLoopCenters;
   std::map<int, double> periodicLoopTargetDensitiesPerSlipSystem;
   std::map<int, double> prismaticLoopTargetDensitiesPerSlipSystem;
   Eigen::Matrix<double, Eigen::Dynamic,3> prismaticLoopCenters;
   std::vector<double> prismaticLoopRadii;
   std::vector<int> prismaticLoopSlipSystemIDs;

   microstructureSpecifications.emplace_back(
         std::shared_ptr<model::MicrostructureSpecification>(
            //new model::PeriodicDipoleIndividualSpecification(
            new model::MicrostructureSpecification(
               "PeriodicDipole", // type
               "individual", // style
               tag, // tag // std::string
               *ddBase,
               periodicDipoleSlipSystemIDsIn, // std::vector<int>
               // dipoles
               exitFaceIDs, // std::vector<int>
               points, // Eigen::Matrix
               heights, // std::vector<double>
               nodes, // std::vector<int>
               glideSteps, // std::vector<double>
               // loops
               periodicLoopSlipSystemIDs,
               periodicLoopRadii,
               periodicLoopSides,
               periodicLoopCenters,
               0.0, //periodicLoopTargetDensityIn,
               0, //periodicLoopSegmentCountIn,
               0.0, //periodicLoopRadiusDistributionMeanIn
               0.0, //periodicLoopRadiusDistributionStdIn
               0.0, //periodicDipoleTargetDensityIn,
               0.0, // prismaticLoopTargetDensityIn,
               periodicLoopTargetDensitiesPerSlipSystem,
               prismaticLoopSlipSystemIDs,
               prismaticLoopCenters,
               prismaticLoopRadii,
               prismaticLoopTargetDensitiesPerSlipSystem,
               1.0, //prismaticLoopRadiiMean // not used
               1.0 // prismaticLoopRadiiStd
               )
            )
         );

   for ( const auto& spec : microstructureSpecifications) // debug
      std::cout << "instantiating a MicrostructureGeneratorInMem with spec: " << spec->microstructureType << " " << spec->style << std::endl; // debug
   return;
}

void ddpy::DDInterface::setOutputPath( const std::string& outputPath)
{ // TODO: change return value to allow error checking
   // if the path doesn't exist, try to make it
   //if ( ! std::filesystem::is_directory( outputPath))
   //{
   //   // TODO: create the folders: outputPath
   //}
   //if (! std::filesystem::is_directory( outputPath + "/evl"))
   //{
   //   try
   //   {
   //      // TODO: create the folders: {outputPath}/evl
   //   }
   //   catch
   //   {
   //      return false;
   //   }
   //}
   //if (! std::filesystem::is_directory( outputPath + "/F"))
   //{
   //   try
   //   {
   //      // TODO: create the folders: {outputPath}/F
   //   }
   //   catch
   //   {
   //      return false;
   //   }
   //}
   //if (! std::filesystem::is_directory( outputPath + "/evl"))
   //{
   //   std::cout << "not found: " << outputPath + "/evl" << std::endl;
   //}
   //else
   //{
   //   std::cout << "found: " << outputPath + "/evl" << std::endl;
   //}
   //if ( DC == nullptr)
   //{
   //   std::cout << "error: defective crystal not yet read" << std::endl;
   //   return;
   //}
   if ( ddBase == nullptr) readddBase();
   std::cout << "prior evlFolder: "
      << ddBase->simulationParameters.traitsIO.evlFolder << std::endl; // debug
   std::cout << "assigning " << outputPath + "/evl" << " to evlFolder" << std::endl;
   ddBase->simulationParameters.traitsIO.evlFolder = outputPath + "/evl";
   ddBase->simulationParameters.traitsIO.auxFolder = outputPath + "/evl";
   ddBase->simulationParameters.traitsIO.fFolder = outputPath + "/F";
   ddBase->simulationParameters.traitsIO.fFile = outputPath + "/F/F_0.txt";
   ddBase->simulationParameters.traitsIO.flabFile = outputPath + "/F/F_labels.txt";
   std::cout << "results:" << std::endl
      << " ddBase->simulationParameters.traitsIO.evlFolder: " << ddBase->simulationParameters.traitsIO.evlFolder
      << std::endl
      << " ddBase->simulationParameters.traitsIO.auxFolder: "
      << ddBase->simulationParameters.traitsIO.auxFolder
      << std::endl
      << " ddBase->simulationParameters.traitsIO.fFolder: "
      << ddBase->simulationParameters.traitsIO.fFolder
      << std::endl;

   // output folders are at:
   //  ddBase->simulationParameters.traitsIO.evlFolder
   //  ddBase->simulationParameters.traitsIO.auxFolder
   //  ddBase->simulationParameters.traitsIO.fFolder
   return;
}

void ddpy::DDInterface::regeneratePolycrystalFile(
            const pybind11::array_t<double,
               pybind11::array::c_style | pybind11::array::forcecast>
               c2gIn,
            const std::string& lattice,
            const std::string& material,
            const std::string& meshFilePath
            )
{
   MatrixDim c2g;
   auto c2gNp = c2gIn.unchecked<2>(); // for reading c2g input np.array
   if (! ((c2gNp.shape(0) == 3) && (c2gNp.shape(1) == 3)))
   {
      std::cout << "error: regeneratePolycrystalFile( C2G) requires C2G "
         << "to be a 3x3 numpy array consisting of rows of normalized "
         << "basis vectors of the crystal lattice." << std::endl;
      return;
   }
   bool latticeIsAcceptable = false;
   bool materialIsAcceptable = false;
   for ( const auto& lat : acceptableLattices)
   {
      if ( lattice == std::string(lat)) latticeIsAcceptable = true;
   }
   for ( const auto& mat : acceptableMaterials)
   {
      if ( material == std::string(mat)) materialIsAcceptable = true;
   }
   if ( ! ( latticeIsAcceptable && materialIsAcceptable))
   {
      std::cout << "error: unacceptable lattice or material: "
         << lattice << ", " << material << std::endl;
      std::cout << "  acceptable lattice or materials are: ";
      for ( const auto& mat : acceptableMaterials) std::cout << mat << ", ";
      for ( const auto& lat : acceptableLattices) std::cout << lat << ", ";
   }
   std::cout << "lattice: " << lattice << ", material: " << material << std::endl; // debug

   c2g << c2gNp(0,0), c2gNp(0,1), c2gNp(0,2),
         c2gNp(1,0), c2gNp(1,1), c2gNp(1,2),
         c2gNp(2,0), c2gNp(2,1), c2gNp(2,2);
   //std::cout << "c2g:\n" << c2g << std::endl; // debug
   std::string outputFilePath = dddFolderPath + "/inputFiles/polycrystal.txt";

   // TODO: detect and move any existing polycrystal.txt file
   std::ofstream outputFile( outputFilePath,
              std::ofstream::out | std::ofstream::trunc);

   //std::cout << "lammps boundaries: " <<  std::endl;// debug
   //for ( const auto& bd : boxBounds) std::cout << bd << ", "; // debug
   //std::cout << std::endl; // debug

   double deltaX, deltaY, deltaZ;
   deltaX = abs( boxBounds[1] - boxBounds[0]);
   deltaY = abs( boxBounds[3] - boxBounds[2]);
   deltaZ = abs( boxBounds[5] - boxBounds[4]);

   MatrixDim AA; // AA scales the mesh: y=A(x-x0), where x is the input mesh
   if (! lattice.compare("bcc")) // .compare() returns 0 if they're equal
   {
      AA << -1.0, 1.0, 1.0,
         1.0, -1.0, 1.0,
         1.0, 1.0, -1.0;
      AA /= sqrt(3.0);
   }
   else if (! lattice.compare("fcc"))
   {
      AA << 0.0, 1.0, 1.0,
          1.0, 0.0, 1.0,
          1.0, 1.0, 0.0;
      AA /= sqrt(2.0);
   }
   else
   {
      std::cout << "error: lattice type not recognized while trying to "
         << "create matrix A" << std::endl;
      return;
   }

   AA = c2g * AA;

   MatrixDim f12, f31, f23;
   f12 << 1.0, boxSkew, 0.0,
       0.0, 1.0, 0.0,
       0.0, 0.0, 1.0;
   f31 << 1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       boxSkew, 0.0, 1.0;
   f23 << 1.0, 0.0, 0.0,
       0.0, 1.0, boxSkew,
       0.0, 0.0, 1.0;
   MatrixDim deformingMatrix;
   deformingMatrix = f12 * (f23 * f31);

   MatrixDim scalingMatrix;
   scalingMatrix << deltaX, 0, 0,
      0, deltaY, 0,
      0, 0, deltaZ;
   deformingMatrix = ( AA.inverse()) * (deformingMatrix * scalingMatrix);
   deformingMatrix = deformingMatrix.array().round();
   AA = AA * deformingMatrix;

   double x0x, x0y, x0z; // boxBounds required to identify x0
   x0x = boxBounds[0]/deltaX; // boxBounds[0] is xlo, deltaX = xhi -xlo
   x0y = boxBounds[2]/deltaY; // boxBounds[2] is ylo
   x0z = boxBounds[4]/deltaZ; // boxBounds[4] is zlo
   VectorDim x0;
   x0 << x0x, x0y, x0z; // scaling of the mesh is y=A(x-x0)

   outputFile << "materialFile=" << material + ".txt;" << std::endl;
   outputFile << "enablePartials=0;"
      << std::endl;
   outputFile << "absoluteTemperature = 300; # [K] simulation temperature"
      << std::endl;
   outputFile << "meshFile=" << meshFilePath << ";"
      << std::endl;
   outputFile << "C2G1="
        << std::setw(22) << std::setprecision(15) << c2g(0,0)
        << std::setw(22) << std::setprecision(15) << c2g(0,1)
        << std::setw(22) << std::setprecision(15) << c2g(0,2)
        << std::endl
        << std::setw(22) << std::setprecision(15) << c2g(1,0)
        << std::setw(22) << std::setprecision(15) << c2g(1,1)
        << std::setw(22) << std::setprecision(15) << c2g(1,2)
        << std::endl
        << std::setw(22) << std::setprecision(15) << c2g(2,0)
        << std::setw(22) << std::setprecision(15) << c2g(2,1)
        << std::setw(22) << std::setprecision(15) << c2g(2,2)
        << ";" << std::endl;
   //outputFile <<  c2g << ";" <<  std::endl; // precision is too low
   outputFile << std::endl;
   outputFile << "A="
        << std::setw(22) << std::setprecision(15) << AA(0,0)
        << std::setw(22) << std::setprecision(15) << AA(0,1)
        << std::setw(22) << std::setprecision(15) << AA(0,2)
        << std::endl
        << std::setw(22) << std::setprecision(15) << AA(1,0)
        << std::setw(22) << std::setprecision(15) << AA(1,1)
        << std::setw(22) << std::setprecision(15) << AA(1,2)
        << std::endl
        << std::setw(22) << std::setprecision(15) << AA(2,0)
        << std::setw(22) << std::setprecision(15) << AA(2,1)
        << std::setw(22) << std::setprecision(15) << AA(2,2)
        << ";" << std::endl;
   outputFile << std::endl << std::endl;
   outputFile << "x0="
      << std::setw(21) << std::setprecision(15) << x0(0)
      << std::setw(21) << std::setprecision(15) << x0(1)
      << std::setw(21) << std::setprecision(15) << x0(2)
      << ";" << std::endl;
   outputFile << "periodicFaceIDs= 0 1 2 3 4 5 " << ";" << std::endl;
   outputFile << std::endl;

   outputFile << "solidSolutionNoiseMode=" << solidSolutionNoiseMode
      << "; # 0=no noise, 1= read noise, 2=compute noise" << std::endl;
   outputFile << "stackingFaultNoiseMode=" << stackingFaultNoiseMode
      << ";" << std::endl;
   outputFile << std::endl;
   //outputFile << "solidSolutionNoiseFile_xz=" << ";" << std::endl;
   //outputFile << "solidSolutionNoiseFile_yz=" << ";" << std::endl;

   return;
}

void ddpy::DDInterface::clearMicrostructureSpecifications()
{
   microstructureSpecifications.clear();
   return;
}

PYBIND11_MODULE( ddpy, m) {
   namespace py = pybind11;
   m.doc() = "TODO: revise m.doc() in src/modelib2py.cpp";
   py::class_<ddpy::DDInterface>( m, "DDInterface")
      .def( py::init(
               // using a lambda function that returns an instantiation
               [](const std::string& modelibFolderPath)
               {
                  return std::unique_ptr< ddpy::DDInterface>(
                        new ddpy::DDInterface( modelibFolderPath)
                        );
               }), py::arg("modelibFolderPath").none(false)
            )
      .def("getCurrentStep",
            &ddpy::DDInterface::getCurrentStep
          )
      .def("setCurrentStep",
            &ddpy::DDInterface::setCurrentStep,
            py::arg("currentStep").none(false)
          )
      .def("setOutputFrequency",
            &ddpy::DDInterface::setOutputFrequency,
            py::arg("outputFrequency").none(false)
          )
      .def("runGlideSteps",
            &ddpy::DDInterface::runGlideSteps,
            py::arg("Nsteps").none(false)
          )
      .def("getResolvedShearStresses",
            &ddpy::DDInterface::getResolvedShearStresses
          )
      .def("getResolvedShearStrains",
            &ddpy::DDInterface::getResolvedShearStrains
          )
      .def("getPlasticStrains",
            &ddpy::DDInterface::getPlasticStrains
          )
      .def("getBurgersMagnitude",
            &ddpy::DDInterface::getBurgersMagnitude
          )
      .def("getDensityPerSlipSystem",
            &ddpy::DDInterface::getDensityPerSlipSystem
          )
      .def("readBurgersMagnitude",
            &ddpy::DDInterface::readBurgersMagnitude,
            py::arg("materialFilePath").none(false)
          )
      .def("regeneratePolycrystalFile",
            &ddpy::DDInterface::regeneratePolycrystalFile,
            py::kw_only(),
            py::arg( "c2g").none(false),
            py::arg( "lattice").none(false),
            py::arg( "material").none(false),
            py::arg( "meshFilePath").none(false)
          )
      .def("specifyLoops",
            &ddpy::DDInterface::specifyLoops,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "slipSystemIDs").none(false),
            py::arg( "loopRadii").none(false),
            py::arg( "loopSegmentCounts").none(false),
            py::arg( "loopCenters").none(false)
          )
      .def("specifyDipoles",
            &ddpy::DDInterface::specifyDipoles,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "slipSystemIDs").none(false),
            py::arg( "exitFaceIDs").none(false),
            py::arg( "points").none(false),
            py::arg( "heights").none(false),
            py::arg( "dipoleNodes").none(false),
            py::arg( "dipoleGlideSteps").none(false)
          )
      .def("specifyLoopDensity",
            &ddpy::DDInterface::specifyLoopDensity,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "loopDensity").none(false),
            py::arg( "loopSegmentCount").none(false),
            py::arg( "loopRadiusMean").none(false),
            py::arg( "loopRadiusStd").none(false)
          )
      .def("specifyLoopDensitiesPerSlipSystem",
            &ddpy::DDInterface::specifyLoopDensitiesPerSlipSystem,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "loopDensitiesPerSlipSystem").none(false),
            py::arg( "loopSegmentCount").none(false), // keep
            py::arg( "loopRadiusMean").none(false), // keep
            py::arg( "loopRadiusStd").none(false) // keep
          )
      .def("specifyPrismaticLoops", // TODO: implememnt
            &ddpy::DDInterface::specifyPrismaticLoops,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "slipSystemIDs").none(false),
            py::arg( "prismaticLoopRadii").none(false),
            py::arg( "prismaticLoopCenters").none(false)
          )
      .def("specifyPrismaticLoopDensity", // TODO: implememnt
            &ddpy::DDInterface::specifyPrismaticLoopDensity,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "prismaticLoopDensity").none(false),
            py::arg( "prismaticLoopRadiusMean").none(false),
            py::arg( "prismaticLoopRadiusStd").none(false)
          )
      .def("specifyPrismaticLoopDensitiesPerSlipSystem", // TODO: implememnt
            &ddpy::DDInterface::specifyPrismaticLoopDensitiesPerSlipSystem,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "prismaticLoopDensitiesPerSlipSystem").none(false),
            py::arg( "prismaticLoopRadiusMean").none(false), // keep
            py::arg( "prismaticLoopRadiusStd").none(false) // keep
          )
      .def("specifyDipoleDensity",
            &ddpy::DDInterface::specifyDipoleDensity,
            py::kw_only(),
            py::arg( "tag").none(false),
            py::arg( "dipoleDensity").none(false)
          )
      .def("generateMicrostructure",
         &ddpy::DDInterface::generateMicrostructure
         )
      .def("writeConfigToTxt",
            &ddpy::DDInterface::writeConfigToTxt
          )
      .def("readDefectiveCrystal",
            &ddpy::DDInterface::readDefectiveCrystal
          )
      .def("setOutputPath",
            &ddpy::DDInterface::setOutputPath,
            py::arg( "outputPath").none(false)
          )
      .def("setBoxBounds", // prerequisite: readBurgersMagnitude()
            &ddpy::DDInterface::setBoxBounds,
            py::arg( "xlo").none(false), // [m]
            py::arg( "xhi").none(false),
            py::arg( "ylo").none(false),
            py::arg( "yhi").none(false),
            py::arg( "zlo").none(false),
            py::arg( "zhi").none(false)
      )
      .def("clearMicrostructureSpecifications",
            &ddpy::DDInterface::clearMicrostructureSpecifications
      )
      .def("setExternalLoad", // prerequisite: readBurgersMagnitude()
            &ddpy::DDInterface::setExternalLoad,
            py::kw_only(),
            py::arg( "stress") = py::none(), //.none(true), // 3x3 matrix
            py::arg( "stressRate") = py::none(), // 3x3 matrix
            py::arg( "strain") = py::none(), // 3x3 matrix
            py::arg( "strainRate") = py::none(), // 3x3 matrix
            py::arg( "stiffnessRatio") = py::none() // Voigt format 11 22 33 12 23 13
          )
      //.def("getSlipSystemNormals"
      //      &ddpy::DDInterface::getSlipSystemNormals
      //    )
      //.def("getSlipSystemBurgersVectors"
      //      &ddpy::DDInterface::getSlipSystemBurgersVectors
      //    )
      ;
   //py::class_<ddpy::DDInterface>( m, "DDInterface")
   //   .def(
   //         py::init(
   //            // using a lambda function that returns an instantiation
   //            [](const std::string& folderName)
   //            {
   //               model::DislocationDynamicsBase ddBase( folderName);
   //               return std::unique_ptr< model::DefectiveCrystal>(
   //                     new model::DefectiveCrystal( ddBase)
   //                     );
   //            }
   //            ), py::arg("dddFolderPath")
   //         )
   //   .def("getCurrentStep",
   //         &ddpy::DDInterface::getCurrentStep
   //         )
   //   .def("setCurrentStep",
   //         &ddpy::DDInterface::setCurrentStep,
   //         py::arg("currentStep")
   //       )
   //   .def("runGlideSteps",
   //         &ddpy::DDInterface::runGlideSteps
   //       );
   //   // TODO: expose any other class members to pybind here
   //   //.def("printDisplacements", &ddpy::DDInterface::printDisplacements)
   //   ;
} //  PYBIND11_MODULE

#endif
