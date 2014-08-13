/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _MODEL_ParticleSystemSerial_h_
#define _MODEL_ParticleSystemSerial_h_


#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif


#include <memory> // for auto_ptr
#include <time.h>
#include <iomanip> // std::scientific
#include <utility> // std::pair

#include <model/ParticleInteraction/PointSource.h>
#include <model/ParticleInteraction/FieldPoint.h>
#include <model/ParticleInteraction/ParticleSystemBase.h>
#include <model/SpaceDecomposition/SpatialCellObserver.h>


//#include <model/Threads/ParallelFor.h> // alterative to openmp using std::thread

namespace model {
    
    
    
    
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief Serial Version (openMP only)
     */
    template <typename _ParticleType >
    class ParticleSystemSerial : public ParticleSystemBase<_ParticleType>
    {
        typedef ParticleSystemBase<_ParticleType> ParticleSystemBaseType;
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
        typedef ParticleSystemSerial<_ParticleType> ParticleSystemSerialType;
        
    public:
        
        typedef typename ParticleSystemBaseType::SpatialCellType SpatialCellType;
        typedef typename ParticleSystemBaseType::ParticleContainerType ParticleContainerType;

        
        /**********************************************************************/
        template <typename FieldType>
        void computeNeighborField(const bool& useMultipole)
        {
            // Nearest-neighbor interaction
            //! -1 loop over all particles in the ParticleSystem (parallelized in OpenMP)
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (unsigned int k=0; k<this->size();++k)
            {
                //! -2 loop over neighbor cells of current particle
                for (typename SpatialCellType::CellMapType::const_iterator cIter =this->operator[](k).neighborCellsBegin();
                     /*                                                 */ cIter!=this->operator[](k).neighborCellsEnd();
                     /*                                               */ ++cIter)
                {
                    //! -3 loop over particles in the current neighbor cell
                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter =cIter->second->particleBegin();
                        /*                                                           */ qIter!=cIter->second->particleEnd();
                        /*                                                         */ ++qIter)
                    {
                        *static_cast<FieldPointBase<ParticleType,FieldType>* const>(&this->operator[](k)) += FieldType::compute(**qIter,this->operator[](k));
                    }
                }
                
                // Non-nearest-neighbor interaction
                if(useMultipole)
                {
                    typename SpatialCellType::CellMapType farCells(this->operator[](k).template farCells<_ParticleType>());
                    *static_cast<FieldPointBase<ParticleType,FieldType>* const>(&this->operator[](k)) += FieldType::multipole(this->operator[](k),farCells);
                }
            }
        }
        

        
        /**********************************************************************/
        template <typename OtherParticleType, typename FieldType>
        void computeField(std::deque<OtherParticleType>& fpDeq, const bool& useMultipole) const
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (unsigned int k=0; k<fpDeq.size();++k)
            {
                // Nearest-neighbor interaction
                typename SpatialCellType::CellMapType neighborCells(fpDeq[k].template neighborCells<_ParticleType>());
                
                //! -2 loop over neighbor cells of current particle
                for (typename SpatialCellType::CellMapType::const_iterator cIter =neighborCells.begin();
                     /*                                                 */ cIter!=neighborCells.end();
                     /*                                               */ ++cIter)
                {
                    //! -3 loop over particles in the current neighbor cell
                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter =cIter->second->particleBegin();
                        /*                                                           */ qIter!=cIter->second->particleEnd();
                        /*                                                         */ ++qIter)
                    {
                        //fpDeq[k].template field<FieldType>() += FieldType::compute(**qIter,fpDeq[k]);
                        *static_cast<FieldPointBase<OtherParticleType,FieldType>* const>(&fpDeq[k]) += FieldType::compute(**qIter,fpDeq[k]);
                    }
                }
                
                // Non-nearest-neighbor interaction
                if(useMultipole)
                {
                    typename SpatialCellType::CellMapType farCells(fpDeq[k].template farCells<_ParticleType>());
//
//                    fpDeq[k].template field<FieldType>() += FieldType::multipole(fpDeq[k],farCells);
                    *static_cast<FieldPointBase<OtherParticleType,FieldType>* const>(&fpDeq[k]) += FieldType::multipole(fpDeq[k],farCells);
                }
                
            }
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator<< (T& os, const ParticleContainerType& pS)
        {/*! Operator<< use ParticleType-specific operator<<
          */
            for (typename ParticleContainerType::const_iterator pIter=pS.begin();pIter!=pS.end();++pIter)
            {
                os<<(*pIter);
            }
            return os;
        }
        
    };

} // end namespace
#endif


//        /**********************************************************************/
//        template <typename OtherParticleType, typename FieldType>
//        void computeField(std::deque<OtherParticleType* const>& fpDeq) const
//        {
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//            for (unsigned int k=0; k<fpDeq.size();++k)
//            {
//                typename SpatialCellType::CellMapType neighborCells(fpDeq[k]->template neighborCells<_ParticleType>());
//
//                //! -2 loop over neighbor cells of current particle
//                for (typename SpatialCellType::CellMapType::const_iterator cIter =neighborCells.begin();
//                     /*                                                 */ cIter!=neighborCells.end();
//                     /*                                               */ ++cIter)
//                {
//                    //! -3 loop over particles in the current neighbor cell
//                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter =cIter->second->particleBegin();
//                        /*                                                           */ qIter!=cIter->second->particleEnd();
//                        /*                                                         */ ++qIter)
//                    {
//                        fpDeq[k]->template field<FieldType>() += FieldType::compute(**qIter,*fpDeq[k]);
//                    }
//                }
//            }
//
//        }




//            ParallelFor().runRange(0,this->size(),this,&ParticleSystemSerial<_ParticleType>::neighborFieldIteration<FieldType>); // alternative to openmp

//        /**********************************************************************/
//        template <typename FieldType>
//        void neighborFieldIteration(const size_t& _begin, const size_t& _end)
//        {
//            for (unsigned int k=_begin; k<_end;++k)
//            {
//                //! -2 loop over neighbor cells of current particle
//                for (typename SpatialCellType::CellMapType::const_iterator cIter =this->operator[](k).neighborCellsBegin();
//                     /*                                                 */ cIter!=this->operator[](k).neighborCellsEnd();
//                     /*                                               */ ++cIter)
//                {
//                    //! -3 loop over particles in the current neighbor cell
//                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter =cIter->second->particleBegin();
//                        /*                                                           */ qIter!=cIter->second->particleEnd();
//                        /*                                                         */ ++qIter)
//                    {
//                        *static_cast<FieldPointBase<ParticleType,FieldType>* const>(&this->operator[](k)) += FieldType::compute(**qIter,this->operator[](k));
//                    }
//                }
//            }
//        }

