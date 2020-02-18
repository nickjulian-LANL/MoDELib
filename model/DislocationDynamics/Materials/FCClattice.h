/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_FCClattice_H_
#define model_FCClattice_H_

#include <memory>
#include <vector>
#include <Eigen/Dense>

#include <LatticeMath.h>
#include <SlipSystem.h>
#include <DislocatedMaterialBase.h>
#include <DislocationMobilityFCC.h>
#include <RationalLatticeDirection.h>

namespace model
{
    
    template<int dim>
    struct FCClattice
    {
        
    };
    
    template<>
    struct FCClattice<3> : public Lattice<3>
    {
        static constexpr int dim=3;
        static constexpr auto name="FCC";
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        static constexpr bool enable111planes=true;
        static constexpr bool enable110planes=false;
        static bool enablePartials;
        
        FCClattice(const MatrixDim& Q) :
        /* init */ Lattice<dim>(getLatticeBasis(),Q)
        {
            
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,dim> getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
          * in units of the crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << 0.0, 1.0, 1.0,
            /*   */ 1.0, 0.0, 1.0,
            /*   */ 1.0, 1.0, 0.0;
            
            return temp/sqrt(2.0);
        }
        
        /**********************************************************************/
        static std::vector<LatticePlaneBase> reciprocalPlaneNormals(const Lattice<dim>& lat)
        
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            
            
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat); // [011]
            LatticeVectorType a2(VectorDimI(0,1,0),lat); // [101]
            LatticeVectorType a3(VectorDimI(0,0,1),lat); // [110]
            
            std::vector<LatticePlaneBase> temp;
            
            if(enable111planes)
            {// {111} planes
                temp.emplace_back(a1,a3);           // is (-1, 1,-1) in cartesian
                temp.emplace_back(a3,a2);           // is ( 1,-1,-1) in cartesian
                temp.emplace_back(a2,a1);           // is (-1,-1, 1) in cartesian
                temp.emplace_back(a1-a3,a2-a3);     // is ( 1, 1, 1) in cartesian
            }
            
            if(enable110planes)
            {// {110} planes
                temp.emplace_back(a1+a2-a3,a3);
                temp.emplace_back(a1+a2-a3,a1-a2);
                temp.emplace_back(a1+a3-a2,a2);
                temp.emplace_back(a1+a3-a2,a1-a3);
                temp.emplace_back(a2+a3-a1,a1);
                temp.emplace_back(a2+a3-a1,a2-a3);
            }
            
            
            return temp;
        }

        /**********************************************************************/
//        static std::vector<std::shared_ptr<SlipSystem>> slipSystems(const DislocatedMaterialBase& materialBase,
        static std::vector<std::shared_ptr<SlipSystem>> slipSystems(const DislocatedMaterial<dim,Isotropic>& material,
                                                                    const Lattice<dim>& lat)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip systems of the Hexagonal lattice
          */
            typedef Eigen::Matrix<long int,3,1> VectorDimI;
            
            
            typedef LatticeVector<3> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat);
            LatticeVectorType a2(VectorDimI(0,1,0),lat);
            LatticeVectorType a3(VectorDimI(0,0,1),lat);
            
//            std::shared_ptr<DislocationMobilityBase> fccMobility(new DislocationMobilityFCC(materialBase));
            const std::shared_ptr<DislocationMobilityBase>& fccMobility(material.mobilities.at("fcc"));
            
            std::vector<std::shared_ptr<SlipSystem>> temp;
            
            if(enable111planes)
            {// <110>{111}
                if(enablePartials)
                {
                    
                    const Eigen::Matrix<double,2,2> A((Eigen::Matrix<double,2,2>()<<1.0,0.5,
                                                                                   0.0,sqrt(3.0)/2).finished());
                    
                    const Eigen::Matrix<double,2,1> N((Eigen::Matrix<double,2,1>()<<2,2).finished());
                    const Eigen::Matrix<double,2,1> D(Eigen::Matrix<double,2,1>::Ones());
                    const Eigen::Matrix<double,Eigen::Dynamic,3> f((Eigen::Matrix<double,Eigen::Dynamic,3>()<<0.0, 0.0, 0.0, // value at origin
                                                                    /*                                     */ 0.5 sqrt(3.0)/6 material.ISF;    % ISF
                                                                    0.25 sqrt(3)/4 182; % USF
                                                                    0.5 0 182
                                                                    ).finished());
                    const Eigen::Matrix<double,Eigen::Dynamic,5> df(TextFileParser("input.txt").readMatrixCols<double>("df",5,true));
                    
                    
                    std::shared_ptr<PeriodicLatticeInterpolant<2>> gammaSurface(new PeriodicLatticeInterpolant<2>(A,N.template cast<size_t>(),D.template cast<size_t>(),f,df));
                    
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1+a3)*(+1)),fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1+a3)*(-1)),fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a3)*(+1)),fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a1*2-a3)*(-1)),fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a1)*(+1)),fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, RationalLatticeDirection<3>(Rational(1,3),(a3*2-a1)*(-1)),fccMobility,nullptr));               // is (-1, 1,-1) in cartesian

                }
                else
                {
                    temp.emplace_back(new SlipSystem(a1,a3, a1,fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3,a1*(-1),fccMobility,nullptr));           // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3, a3,fccMobility,nullptr));               // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3,a3*(-1),fccMobility,nullptr));           // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3,a1-a3,fccMobility,nullptr));             // is (-1, 1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a1,a3,a3-a1,fccMobility,nullptr));             // is (-1, 1,-1) in cartesian
                    
                    temp.emplace_back(new SlipSystem(a3,a2, a3,fccMobility,nullptr));               // is ( 1,-1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2,a3*(-1),fccMobility,nullptr));           // is ( 1,-1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2, a2,fccMobility,nullptr));               // is ( 1,-1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2,a2*(-1),fccMobility,nullptr));           // is ( 1,-1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2,a3-a2,fccMobility,nullptr));             // is ( 1,-1,-1) in cartesian
                    temp.emplace_back(new SlipSystem(a3,a2,a2-a3,fccMobility,nullptr));             // is ( 1,-1,-1) in cartesian
                    
                    temp.emplace_back(new SlipSystem(a2,a1, a2,fccMobility,nullptr));               // is (-1,-1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1,a2*(-1),fccMobility,nullptr));           // is (-1,-1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1, a1,fccMobility,nullptr));               // is (-1,-1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1,a1*(-1),fccMobility,nullptr));           // is (-1,-1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1,a2-a1,fccMobility,nullptr));             // is (-1,-1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a2,a1,a1-a2,fccMobility,nullptr));             // is (-1,-1, 1) in cartesian
                    
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a1-a3,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a3-a1,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3,a2-a3,fccMobility,nullptr));       // is ( 1, 1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3,a3-a2,fccMobility,nullptr));       // is ( 1, 1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a1-a2,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
                    temp.emplace_back(new SlipSystem(a1-a3,a2-a3, a2-a1,fccMobility,nullptr));      // is ( 1, 1, 1) in cartesian
                }
            }
            

            
            if(enable110planes)
            {// <110>{110}
                temp.emplace_back(new SlipSystem(a1+a2-a3,a3, a3,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a1+a2-a3,a3, a3*(-1),fccMobility,nullptr));
                
                temp.emplace_back(new SlipSystem(a1+a2-a3,a1-a2,a1-a2,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a1+a2-a3,a1-a2,a2-a1,fccMobility,nullptr));
                
                temp.emplace_back(new SlipSystem(a1+a3-a2,a2, a2,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a1+a3-a2,a2, a2*(-1),fccMobility,nullptr));
                
                temp.emplace_back(new SlipSystem(a1+a3-a2,a1-a3,a1-a3,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a1+a3-a2,a1-a3,a3-a1,fccMobility,nullptr));
                
                temp.emplace_back(new SlipSystem(a2+a3-a1,a1, a1,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a2+a3-a1,a1, a1*(-1),fccMobility,nullptr));
                
                temp.emplace_back(new SlipSystem(a2+a3-a1,a2-a3,a2-a3,fccMobility,nullptr));
                temp.emplace_back(new SlipSystem(a2+a3-a1,a2-a3,a3-a2,fccMobility,nullptr));
            }
            
            return temp;
        }
        
        
        
    };
    
    bool FCClattice<3>::enablePartials=false;
    
} // namespace model
#endif

