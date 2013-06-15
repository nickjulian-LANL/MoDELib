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
#else
#define omp_get_thread_num() 0
//#define omp_get_max_threads() 1
#endif


#include <memory> // for auto_ptr
#include <time.h>
#include <iomanip> // std::scientific
#include <utility> // std::pair
#include <vector> // std::vector

#include <model/Utilities/SequentialOutputFile.h>

#include <model/SpaceDecomposition/SpatialCellProperty.h>

#include <model/ParticleInteraction/ParticleSystemBase.h>
#include <model/ParticleInteraction/SystemProperties.h>




namespace model {
    
    
    
    
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief Serial Version
     */
    template <typename _ParticleType, typename UserSystemProperties = SystemProperties<> >
    class ParticleSystemSerial :
    /* inheritance */  public ParticleSystemBase<_ParticleType,UserSystemProperties>
    {
        typedef ParticleSystemBase<_ParticleType,UserSystemProperties> ParticleSystemBaseType;
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
        //        typedef typename ParticleSystemBaseType::PositionType PositionType;
        typedef typename ParticleSystemBaseType::SpatialCellType SpatialCellType;
        typedef typename ParticleSystemBaseType::ParticleContainerType ParticleContainerType;
//        typedef typename ParticleSystemBaseType::SpatialCellObserverType SpatialCellObserverType;
        
        
        
        
        //        void loopOverNeighCells(typename ParticleContainerType::iterator& pIter)
        //        {
        //            for (typename SpatialCellType::CellMapType::const_iterator cIter=pIter->second->neighborCellsBegin();cIter!=pIter->second->neighborCellsEnd();++cIter)
        //            {
        //                // loop over particles in the neighbor cell
        //                for(typename SpatialCellType::ParticleContainerType::const_iterator qIter=cIter->second->particleBegin();qIter!=cIter->second->particleEnd();++qIter) // loop over neighbor particles
        //                {
        //                    InteractionType(*pIter->second,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
        //                }
        //            }
        //        }
        
        
    public:

        
        /**********************************************************************/
        template <typename InteractionType>
        void computeNeighborInteraction()
        {
            
//            double t0(clock());
#ifdef _OPENMP
            // compute the iterator range of each thread
            typedef typename ParticleContainerType::iterator ParticleContainerIteratorType;
            const int max_threads(omp_get_max_threads());
            std::vector<std::pair<ParticleContainerIteratorType,ParticleContainerIteratorType> > chunks;
            chunks.reserve(max_threads);
            size_t chunk_size= this->size() / max_threads;
            ParticleContainerIteratorType cur_iter = this->begin();
            for(int i = 0; i < max_threads - 1; ++i)
            {
                ParticleContainerIteratorType last_iter = cur_iter;
                std::advance(cur_iter, chunk_size);
                chunks.push_back(std::make_pair(last_iter, cur_iter));
            }
            chunks.push_back(std::make_pair(cur_iter, this->end()));
            
#pragma omp parallel for shared(chunks)
            // each thread picks one value of i
            for (int i=0;i<max_threads;++i)
            {
                // loop over particles assigned to i-th thread
                for (typename ParticleContainerType::iterator pIter=chunks[i].first;pIter!=chunks[i].second;++pIter)
                {
                    // loop over neighbor cells
                    for (typename SpatialCellType::CellMapType::const_iterator cIter=pIter->second->neighborCellsBegin();cIter!=pIter->second->neighborCellsEnd();++cIter)
                    {
                        // loop over particles in the neighbor cell
                        for(typename SpatialCellType::ParticleContainerType::const_iterator qIter=cIter->second->particleBegin();qIter!=cIter->second->particleEnd();++qIter) // loop over neighbor particles
                        {
                            InteractionType(*pIter->second,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                        }
                    }
                }
            }
#else
            // loop over all particles
            for (typename ParticleContainerType::iterator pIter=this->begin();pIter!=this->end();++pIter)
            {
                
//                pIter->second->computeStress();
                
                // loop over neighbor cells
                for (typename SpatialCellType::CellMapType::const_iterator cIter=pIter->second->neighborCellsBegin();cIter!=pIter->second->neighborCellsEnd();++cIter)
                {
                    // loop over particles in the neighbor cell
                    for(typename SpatialCellType::ParticleContainerType::const_iterator qIter=cIter->second->particleBegin();qIter!=cIter->second->particleEnd();++qIter) // loop over neighbor particles
                    {
                        InteractionType(*pIter->second,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
//                        InteractionType::compute(*pIter->second,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                    }
                }
            }
#endif
//            std::cout<<std::scientific<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]."<<std::endl;
        }
        
        
        //        /**********************************************************************/
        //        template <typename CellProperty>
        //        void computeCellProperty() const
        //        {
        ////#ifdef _OPENMP
        //
        //
        ////#else
        //            for (typename SpatialCellType::CellMapType::const_iterator cIter=SpatialCellObserverType::cellBegin();cIter!=SpatialCellObserverType::cellEnd();++cIter)
        //            {
        //                CellProperty(*cIter->second);
        //            }
        ////#endif
        //        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator<< (T& os, const ParticleContainerType& pS)
        {/*! Operator<< use ParticleType-specific operator<<
          */
            for (typename ParticleContainerType::const_iterator pIter=pS.begin();pIter!=pS.end();++pIter)
            {
                os<<(*pIter->second)<<std::endl;
            }
            return os;
        }
        
    };
    
    
    
    
    
} // end namespace
#endif

