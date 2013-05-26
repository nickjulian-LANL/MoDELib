/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * PIL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _ParticleSystem_h
#define _ParticleSystem_h

#include <mpi.h>

#include <metis.h> // partitioner

#include <memory> // for auto_ptr
#include <utility> // for std::pair
#include <map> 
#include <vector>
#include <deque> 
#include <boost/ptr_container/ptr_map.hpp> // TO BE CHANGED WITH ACTUAL MPI IMPLEMENTATION

#include <Eigen/Core>

#include <model/PIL/PilMPI.h>
#include <model/PIL/SystemProperties.h>
#include <model/PIL/SpatialCells/SpatialCellObserver.h>


#include <model/Utilities/SequentialOutputFile.h>


namespace pil {
    
    template <typename _ParticleType, typename UserSystemProperties = SystemProperties<> >
    class ParticleSystem :
    /* inheritance */  public PilMPI, 
    /* inheritance */  protected boost::ptr_map<size_t,_ParticleType>, // TO BE CHANGED WITH ACTUAL MPI IMPLEMENTATION
    /* inheritance */  private UserSystemProperties
    {
        //enum {dim=_ParticleType::dim};
        typedef SpatialCell<_ParticleType,_ParticleType::dim> SpatialCellType;
        typedef boost::ptr_map<size_t,_ParticleType> ParticleContainerType;
        typedef SpatialCellObserver<_ParticleType,_ParticleType::dim> SpatialCellObserverType;
        
        typedef std::deque<SpatialCell<_ParticleType,_ParticleType::dim>*> AssignedCellContainerType;
        AssignedCellContainerType assignedCells;

        
        int particleRankOffset;
        int assignedParticleSize;
        std::vector<int> particleRankOffsetVector;
        std::vector<int> particleSizeVector;
        /*****************************************/
        void assignReorderedIDs()
        {
            
            // Find Offset
            particleRankOffset=0;
            MPI_Exscan(&assignedParticleSize,&particleRankOffset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            //if (this->mpiRank==0) particleRankOffset=0; // make sure that on processor 0 particleRankOffset is 0

            MPI_Allgather(&assignedParticleSize,1,MPI_INT,&particleSizeVector[0],1,MPI_INT,MPI_COMM_WORLD);

            
            
  //          std::cout<<"Processor "<<this->mpiRank<<" :  particleRankOffset="<<particleRankOffset<<std::endl;
            
            int rID(0);
            for (typename AssignedCellContainerType::const_iterator cIter=assignedCells.begin();cIter!=assignedCells.end();++cIter){ // loop over cells assigned to this process
                for (typename SpatialCellType::ParticleContainerType::const_iterator pIter=(*cIter)->particleBegin();pIter!=(*cIter)->particleEnd();++pIter) // loop over particles in the current cell
                {
                    (*pIter)->rID=rID+particleRankOffset;
                    rID++;
                }
            }
            
            
        //MPI_Exscan(&assignedParticleSize,&particleRankOffset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        MPI_Allgather(&particleRankOffset,1,MPI_INT,&particleRankOffsetVector[0],1,MPI_INT,MPI_COMM_WORLD);
        
//        if (this->mpiRank==0){
//            for (int k=0;k<particleRankOffsetVector.size();++k)
//            {
//                std::cout<<particleRankOffsetVector[k]<<" "<<std::endl;
//            }
//        }
        
        }
        

        
    public:
        
    //    const int mpiRank;
    //    const int nProcs;

        
        typedef _ParticleType ParticleType; // make ParticleType available outside the class
        typedef typename ParticleType::PositionType PositionType;
        
        
        /*****************************************/
        ParticleSystem(const double& cellSize=1.0)
        /* init list */ 
        {
//            int mpiRank_temp, nProcs_temp;
//            MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
//            MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
            SpatialCellObserver<_ParticleType,_ParticleType::dim>::cellSize=cellSize;
            particleRankOffsetVector.resize(this->nProcs);
            particleSizeVector.resize(this->nProcs);
        }
        
        /*****************************************/
        template <typename ...AdditionalConstructorTypes>
        int addParticle(const PositionType& p, const AdditionalConstructorTypes&... args)
        {/*! Constructor with particle position and arbitrary number of
          *  additional parameters.
          */
            // 1- Generate a new ParticleType objects using std::auto_ptr
            std::auto_ptr<ParticleType> pN (new ParticleType(p,args...) );
            // 2- Obtain the sID (StaticID) of the new particle
            const size_t sID(pN->sID);
            // 3- Insert the particle in the ParticleSystem using sID as the key. Assert succesful insertion
            assert(this->insert(sID , pN ).second && "CANNOT INSERT PARTICLE IN PARTICLE SYSTEM.");
            
            
            //assignedCells
            
            
            // 4- Returns the sID of the particle.
            return sID;
        }
        
        
        
        /*****************************************/
        //template <typename InteractionType>
        void partionCells()
        {
            typedef typename SpatialCellObserverType::CellIdType CellIdType;
            typedef std::map<CellIdType,const int,model::CompareVectorsByComponent<int,_ParticleType::dim> > ijk2idMapType;
            //            typedef std::map<Eigen::Matrix<int,_ParticleType::dim,1>,const int,CompareVectorsByComponent<int,_ParticleType::dim> > ijk2idMapType;
            ijk2idMapType  ijk2idMap;;
            typedef std::map<int,const CellIdType> id2ijkMapType;
            //            typedef std::map<int,Eigen::Matrix<int,_ParticleType::dim,1> > id2ijkMapType;
            id2ijkMapType id2ijkMap;;
            
            
            std::vector<int> vWeights(SpatialCellObserverType::size(),1); // vector containing the weight of each cell

            
            int id(0);
            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter=SpatialCellObserverType::begin();cIter!=SpatialCellObserverType::end();++cIter)
            {
                ijk2idMap.insert(std::pair<Eigen::Matrix<int,_ParticleType::dim,1>,const int>(cIter->second->cellID,id));
                id2ijkMap.insert(std::pair<int,Eigen::Matrix<int,_ParticleType::dim,1> >(id,cIter->second->cellID));
                vWeights[id]=cIter->second->n2Weight();
                id++;
            }
            
            
            
            // - Populate vector of neighbors and offsets  
            std::vector<int> vNeighbors;
            std::vector<int> vOffsets;
            
            int offSet(0);
            for (typename SpatialCellObserverType::CellMapType::const_iterator cIter=SpatialCellObserverType::begin();cIter!=SpatialCellObserverType::end();++cIter)
            {
                vOffsets.push_back(offSet);
                for (typename SpatialCellObserverType::CellMapType::const_iterator nIter=cIter->second->neighborCellsBegin();nIter!=cIter->second->neighborCellsEnd();++nIter)
                {
                    if (cIter->second->cellID!=nIter->second->cellID)
                    {
                        typename ijk2idMapType::const_iterator mIter(ijk2idMap.find(nIter->second->cellID.transpose()));
                        vNeighbors.push_back(mIter->second);
                        offSet++;
                    }
                }
            }
            vOffsets.push_back(offSet); // push back an extra offset (required by Metis)
            


//            int mpiRank, nProcs;
//            MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
//            MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
            
            // partition into the same number of MPI processes
            int numPart = this->nProcs;
            //int numPart = PilMPI::nProcs;
            

            // Create a vector to store the output of the partitioner
            std::vector<int> partitionerRankVector(SpatialCellObserverType::size(),0); //
            
            if (numPart > 1)
            {
                int metisOptions[METIS_NOPTIONS];
                int numVert = SpatialCellObserverType::size(); // number of nodes
                int ncon = 1;  // number of weights per node
                int objval;
                float ubvec = 1.001; // load imbalance tolerance
                METIS_SetDefaultOptions(metisOptions);
                METIS_PartGraphKway(&numVert,&ncon,&vOffsets[0],&vNeighbors[0],
                                    &vWeights[0],NULL,NULL,&numPart,NULL,
                                    &ubvec,metisOptions,&objval,&partitionerRankVector[0]);
            }
            
            
            
            // reset
            assignedParticleSize=0;
            assignedCells.clear();
            // Assign ranks to cells
            for (unsigned int i=0; i<SpatialCellObserverType::size(); ++i)
            {
                typename id2ijkMapType::const_iterator idIter(id2ijkMap.find(i));
                const CellIdType cellID(idIter->second);
                
                typename SpatialCellObserverType::isCellType isC(SpatialCellObserverType::isCell(cellID));
                assert(isC.first && "CELL MUST EXIST!");
                isC.second->assignedRank=partitionerRankVector[i];
                //SpatialCellObserverType::getExistingCellByID(cellID).assignedRank=partitionerRankVector[i];
                
                if (this->mpiRank==partitionerRankVector[i])
                //if (PilMPI::mpiRank==partitionerRankVector[i])
                {
                    assignedCells.push_back(isC.second);
                    assignedParticleSize+=isC.second->size();
                }
            }
    
    
    
        }
        
        /*****************************************/
        template <typename InteractionType>
        void computeInteraction()
        {/*! Compute nearest-neighbor particle interaction according to InteractionType
          */
            
            assignReorderedIDs();
            
            
            InteractionType::resultVector.resize(3*this->size(),0.0);
            
            for (typename AssignedCellContainerType::const_iterator cIter=assignedCells.begin();cIter!=assignedCells.end();++cIter){ // loop over cells assigned to this process
                for (typename SpatialCellType::ParticleContainerType::const_iterator pIter=(*cIter)->particleBegin();pIter!=(*cIter)->particleEnd();++pIter) // loop over particles in the current cell
                {
                    for (typename SpatialCellType::CellMapType::const_iterator nIter=(*cIter)->neighborCellsBegin();nIter!=(*cIter)->neighborCellsEnd();++nIter) // loop over neighbor cells
                    {
                        for(typename SpatialCellType::ParticleContainerType::const_iterator qIter=nIter->second->particleBegin();qIter!=nIter->second->particleEnd();++qIter) // loop over neighbor particles
                        {
                            InteractionType(**pIter,**qIter);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
                        }
                    }
                }
            }
            

            //InteractionType::resultVerctor;
            
            //particleRankOffsetVector
            std::vector<int> interactionSizeVector(this->nProcs);            
            std::vector<int> interactionRankOffsetVector(this->nProcs);
            
            for (int k=0;k<interactionSizeVector.size();++k)
            {
                interactionSizeVector[k]=particleSizeVector[k]*InteractionType::DataPerParticle;
                interactionRankOffsetVector[k]=particleRankOffsetVector[k]*InteractionType::DataPerParticle;
                
            }
            
            
            MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
                           &InteractionType::resultVector[0],&interactionSizeVector[0],&interactionRankOffsetVector[0],MPI_DOUBLE,MPI_COMM_WORLD);
//            &InteractionType::resultVector[0],&interactionSizeVector[0],&interactionRankOffsetVector[0],InteractionType::ResultType,MPI_COMM_WORLD);

            
            if (this->mpiRank==0){
                for (int k=0;k<InteractionType::resultVector.size()/3;++k)
                {
                    std::cout<<InteractionType::resultVector[3*k]<<" "<<InteractionType::resultVector[3*k+1]<<" "<<InteractionType::resultVector[3*k+2]<<"\n";
                }
            
            }
            
            
        }
        
        
        /*****************************************/
        template <typename InteractionType>
        typename InteractionType::ResultType getInteractionResult(const int& pID)
        {
            return InteractionType::getResult(*(this->find(pID)->second));
        }
        
        /*****************************************/
        template <typename InteractionType>
        void resetInteraction()
        {/*! Compute full particle interaction according to InteractionType
          */
            for (typename ParticleContainerType::iterator iterJ=this->begin();iterJ!=this->end();++iterJ){ // FOR SYMMETRIC INTERACTION iterJ STARTS AT iterI
                InteractionType::reset(*iterJ->second);  // the constructor of InteractionType actually computes the binary interaction between *iterI and *iterJ
            }
        }
        
        /*****************************************/
        template <typename Property>
        const typename Property::PropertyType& getProperty() const
        {/*! Compute full particle interaction according to InteractionType
          */
            return *static_cast<const Property*>(this);
        }
        
        
        /* cells ****************************************/
        SpatialCellObserver<_ParticleType,_ParticleType::dim> cells()
        {
            return SpatialCellObserver<_ParticleType,_ParticleType::dim>();
        }
        
        /*****************************************/
        template<char P='P',bool autodelete=true>
        void MPIoutput()
        {
 //           int mpiRank;
 //           MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
            if (this->mpiRank == 0)
            //if (PilMPI::mpiRank == 0)
            {
                model::SequentialOutputFile<P,autodelete> pFile;
                pFile<<*this<<std::endl;
                
                model::SequentialOutputFile<'C',1> cFile;
                cFile<<cells()<<std::endl;
            }
        }
        
        /*****************************************/
        template <class T>
        friend T& operator<< (T& os, const ParticleSystem<ParticleType>& pS)
        {/*! Operator << use ParticleType specific << operator
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