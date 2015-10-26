/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNETWORKIO_H_
#define model_DISLOCATIONNETWORKIO_H_

#include <string>
#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>
#include <model/Utilities/UniqueOutputFile.h>
#include <model/Utilities/SequentialOutputFile.h>
#include <model/Utilities/SequentialBinFile.h>

#include <model/Utilities/TerminalColors.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/MPI/MPIcout.h>
#include <model/LatticeMath/LatticeMath.h>
#include <model/LatticeMath/LatticeMath.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename DislocationNetworkType>
    struct DislocationNetworkIO
    {
        
        enum {dim=DislocationNetworkType::dim};
        
        //    public:
        typedef typename DislocationNetworkType::VectorDimD VectorDimD;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        typedef typename DislocationNetworkType::NetworkNodeContainerType NetworkNodeContainerType;
        typedef typename DislocationNetworkType::GlidePlaneObserverType GlidePlaneObserverType;
        typedef typename DislocationNetworkType::SpatialCellObserverType SpatialCellObserverType;
        typedef typename SpatialCellObserverType::CellMapType CellMapType;
        typedef typename DislocationNetworkType::BvpSolverType::FiniteElementType FiniteElementType;
        typedef typename DislocationNetworkType::BvpSolverType::TrialFunctionType TrialFunctionType;
        typedef LatticeVector<dim> LatticeVectorType;
        
        enum {NdofXnode=NodeType::NdofXnode};
        
        static int  outputFrequency;
        static bool outputBinary;
        static bool outputGlidePlanes;
        static bool outputSpatialCells;
        static bool outputPKforce;
        static bool outputElasticEnergy;
        static bool outputMeshDisplacement;
        static bool outputDislocationLength;
        static bool outputPlasticDistortionRate;
        static bool outputQuadratureParticles;
        
        
        /* readVertices *******************************************************/
        static void readVertices(DislocationNetworkType& DN, const unsigned int& fileID)
        {/*! Reads file V/V_0.txt and creates DislocationNodes
          */
            typedef VertexReader<'V',9,double> VertexReaderType;
            VertexReaderType  vReader;	// sID,Px,Py,Pz,Tx,Ty,Tz,snID
            if (vReader.isGood(fileID,false)) // bin file exists
            {
                vReader.read(fileID,false);
                
            }
            else
            {
                if (vReader.isGood(fileID,true)) // txt file exists
                {
                    vReader.read(fileID,true);
                    
                }
                else
                {
                    assert(0 && "UNABLE TO READ VERTEX FILE V/V_x (x is the requested file ID).");
                    
                }
            }
            
            unsigned int kk(1);
            for (VertexReaderType::iterator vIter=vReader.begin();vIter!=vReader.end();++vIter)
            {
                const size_t nodeIDinFile(vIter->first);
                NodeType::set_count(nodeIDinFile);
                model::cout << "\r \r" << "Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<vReader.size()<<")"<<std::flush;
                
                LatticeVectorType L(VectorDimD(vIter->second.template segment<NdofXnode>(0)));
                const size_t nodeID(DN.insertVertex(L).first->first);
                assert(nodeID==nodeIDinFile);
                kk++;
            }
            model::cout<<std::endl;
        }
        
        /* readEdges **********************************************************/
        static void readEdges(DislocationNetworkType& DN, const unsigned int& fileID)
        {/*! Reads file E/E_0.txt and creates DislocationSegments
          */
            typedef EdgeReader  <'E',11,double>	EdgeReaderType;
            EdgeReaderType    eReader;	// sourceID,sinkID,Bx,By,Bz,Nx,Ny,Nz
            if (eReader.isGood(fileID,false)) // bin file exists
            {
                eReader.read(fileID,false);
                
            }
            else
            {
                if (eReader.isGood(fileID,true)) // txt file exists
                {
                    eReader.read(fileID,true);
                    
                }
                else
                {
                    assert(0 && "UNABLE TO READ EDGE FILE E/E_x (x is the requested file ID).");
                }
            }
            
            unsigned int kk(1);
            for (EdgeReaderType::iterator eIter=eReader.begin();eIter!=eReader.end();++eIter)
            {
                VectorDimD B(eIter->second.template segment<dim>(0  ).transpose()); // Burgers vector
                VectorDimD N(eIter->second.template segment<dim>(dim).transpose()); // Glide plane normal
                const size_t sourceID(eIter->first.first );
                const size_t   sinkID(eIter->first.second);
                model::cout << "\r \r" << "Creating DislocationSegment "<<sourceID<<"->"<<sinkID<<" ("<<kk<<" of "<<eReader.size()<<")              "<<std::flush;
                const bool success=DN.connect(sourceID,sinkID,LatticeVectorType(B));
                assert(success && "UNABLE TO CREATE CURRENT DISLOCATION SEGMENT.");
                kk++;
            }
            model::cout<<std::endl;
        }
        
        /* outputTXT **********************************************************/
        static void output(const DislocationNetworkType& DN, const unsigned int& runID)
        {/*! Outputs DislocationNetwork data to the following files (x is the runID):
          * ./E/E_x.txt (DislocationSegment(s) are always outputted)
          * ./V/V_x.txt (DislocationNode(s) are always outputted)
          * ./C/C_x.txt (DislocationCell(s) only if outputSpatialCells==true)
          * ./G/G_x.txt (GlidePlane(s) only if outputGlidePlanes==true)
          * ./P/P_x.txt (PK forces only if outputPKforce==true)
          * ./D/D_x.txt (mesh displacement only if outputMeshDisplacement==true)
          */
            model::cout<<"		Writing to "<<std::flush;
            
            //! 1- Outputs the Edge informations to file E_*.txt where * is the current simulation step
            if (outputBinary)
            {
                typedef std::pair<std::pair<int,int>,Eigen::Matrix<double,1,9> > BinEdgeType;
                SequentialBinFile<'E',BinEdgeType>::set_count(runID);
                SequentialBinFile<'E',BinEdgeType>::set_increment(outputFrequency);
                SequentialBinFile<'E',BinEdgeType> binEdgeFile;
                for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
                {
                    Eigen::Matrix<double,1,9> temp( (Eigen::Matrix<double,1,9>()<< linkIter->second.flow.cartesian().transpose(),
                                                     /*                                                          */ linkIter->second.glidePlaneNormal.transpose(),
                                                     /*                                                          */ linkIter->second.sourceTfactor,
                                                     /*                                                          */ linkIter->second.sinkTfactor,
                                                     /*                                                          */ linkIter->second.pSN()->sID).finished());
                    binEdgeFile.write(std::make_pair(linkIter->first,temp));
                }
                model::cout<<" E/E_"<<binEdgeFile.sID<<".bin"<<std::flush;
            }
            else
            {
                SequentialOutputFile<'E',1>::set_count(runID); // edgeFile;
                SequentialOutputFile<'E',1>::set_increment(outputFrequency); // edgeFile;
                SequentialOutputFile<'E',1> edgeFile;
                //edgeFile << *(const NetworkLinkContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
                for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
                {
                    edgeFile<< linkIter->second<<"\n";
                }
                model::cout<<" E/E_"<<edgeFile.sID<<".txt"<<std::flush;
            }
            
            //! 2- Outputs the Vertex informations to file V_*.txt where * is the current simulation step
            if (outputBinary)
            {
                typedef Eigen::Matrix<double,1,8> VertexDataType;
                typedef std::pair<int, VertexDataType> BinVertexType;
                SequentialBinFile<'V',BinVertexType>::set_count(runID);
                SequentialBinFile<'V',BinVertexType>::set_increment(outputFrequency);
                SequentialBinFile<'V',BinVertexType> binVertexFile;
                //                for (typename NetworkNodeContainerType::const_iterator nodeIter=DN.nodeBegin();nodeIter!=DN.nodeEnd();++nodeIter)
                for (const auto& node : DN.nodes())
                {
                    VertexDataType temp( (VertexDataType()<< node.second.get_P().transpose(),
                                          /*                                    */ node.second.get_T().transpose(),
                                          /*                                    */ node.second.pSN()->sID,
                                          /*                                    */ (node.second.meshLocation()==onMeshBoundary)).finished());
                    binVertexFile.write(std::make_pair(node.first,temp));
                }
                model::cout<<" V/V_"<<binVertexFile.sID<<".bin"<<std::flush;
            }
            else
            {
                SequentialOutputFile<'V',1>::set_count(runID); // vertexFile;
                SequentialOutputFile<'V',1>::set_increment(outputFrequency); // vertexFile;
                SequentialOutputFile<'V',1> vertexFile;
                //vertexFile << *(const NetworkNodeContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
                for (const auto& node : DN.nodes())
                {
                    vertexFile << (node.second)<<"\n";
                }
                model::cout<<", V/V_"<<vertexFile.sID<<".txt"<<std::flush;
            }
            
            if(outputSpatialCells)
            {
                //! 3- Outputs the nearest neighbor Cell structures to file C_*.txt where * is the current simulation step
                SequentialOutputFile<'C',1>::set_count(runID); // Cell_file;
                SequentialOutputFile<'C',1>::set_increment(outputFrequency); // Cell_file;
                SequentialOutputFile<'C',1> Cell_file;
                //              SpatialCellObserverType SPC;
                int cID(0);
                //for (typename CellMapType::const_iterator cellIter=SpatialCellObserverType::cellBegin();cellIter!=SpatialCellObserverType::cellEnd();++cellIter)
                for (const auto& cell : SpatialCellObserverType::cells())
                {
                    Cell_file<<cID<<"\t"<<cell.second->cellID.transpose()<<"\t"<<SpatialCellObserverType::cellSize()<<"\n";
                    ++cID;
                }
                model::cout<<", C/C_"<<Cell_file.sID<<std::flush;
            }
            
            if(outputGlidePlanes)
            {
                //! 4- Outputs the glide planes
                SequentialOutputFile<'G',1>::set_count(runID); // GlidePlanes_file;
                SequentialOutputFile<'G',1>::set_increment(outputFrequency); // GlidePlanes_file;
                SequentialOutputFile<'G',1> glide_file;
                glide_file << *dynamic_cast<const GlidePlaneObserverType*>(&DN);
                model::cout<<", G/G_"<<glide_file.sID<<std::flush;
            }
            
            if(outputPKforce)
            {
                SequentialOutputFile<'P',1>::set_count(runID); // Edges_file;
                SequentialOutputFile<'P',1>::set_increment(outputFrequency); // Edges_file;
                SequentialOutputFile<'P',1> p_file;
                int ll=0;
                for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
                {
                    const int qOrder(linkIter->second.rgauss.cols());
                    for (int q=0;q<qOrder;++q)
                    {
                        //                        p_file << ll*qOrder+q<<" "<< linkIter->second.rgauss.col(q).transpose()<<" "<<linkIter->second.pkGauss.col(q).transpose()<<"\n";
                        p_file << ll<<" "<< linkIter->second.rgauss.col(q).transpose()<<" "<<linkIter->second.pkGauss.col(q).transpose()<<"\n";
                        ll++;
                    }
                    //                    ll++;
                }
                model::cout<<", P/P_"<<p_file.sID<<std::flush;
            }
            
            if(outputElasticEnergy)
            {
                typedef typename DislocationNetworkType::DislocationParticleType::ElasticEnergy ElasticEnergy;
                SequentialOutputFile<'W',1>::set_count(runID);
                SequentialOutputFile<'W',1>::set_increment(outputFrequency);
                SequentialOutputFile<'W',1> w_file; //energy_file
                int ll=0;
                for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
                {
                    const int qOrder(linkIter->second.rgauss.cols());
                    for (size_t q=0;q<linkIter->second.quadratureParticleContainer.size();++q)
                    {
                        w_file << ll*qOrder+q<<" "<< linkIter->second.rgauss.col(q).transpose()<<" "<< linkIter->second.quadratureParticleContainer[q]->template field<ElasticEnergy>()<<"\n";
                    }
                    ll++;
                }
                model::cout<<", W/W_"<<w_file.sID<<std::flush;
            }
            
            if (DN.shared.use_bvp && outputMeshDisplacement && !(DN.runningID()%DN.shared.use_bvp))
            {
                model::SequentialOutputFile<'U',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'U',1>::set_increment(outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'U',true> d_file;
                d_file<<DN.shared.bvpSolver.displacement();
                
                
                //                for(const auto& node : DN.shared.bvpSolver.finiteElement().nodes())
                //                {
                //
                //                }
                
                /**************************************************************************/
                // Output displacement and stress on external mesh faces
                model::SequentialOutputFile<'S',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'S',1>::set_increment(outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'S',true> s_file;
                s_file<<DN.shared.bvpSolver.stress();
                
                model::cout<<", D/D_"<<d_file.sID<<"(FINISH HERE)"<<std::flush;
            }
            
            if (outputQuadratureParticles)
            {
                model::SequentialOutputFile<'Q',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'Q',1>::set_increment(outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'Q',true> q_file;
                for (const auto& particle : DN.particles())
                {
                    q_file<<particle<<"\n";
                }
                
                
            }
            
            // Output to F file
            UniqueOutputFile<'F'> f_file;
            model::cout<<" F/F_0.txt"<<std::flush;
            f_file<< DN.runningID()<<" "<<DN.get_totalTime()<<" "<<DN.get_dt()<<" ";
            
            if(outputDislocationLength)
            {
                const auto length=DN.networkLength();
                f_file<<length.first<<" "<<length.second<<" ";
            }
            
            if(outputPlasticDistortionRate)
            {
                Eigen::Matrix<double,dim,dim> pDR(DN.plasticDistortionRate());
                f_file<<pDR.row(0)<<" "<<pDR.row(1)<<" "<<pDR.row(2)<<" ";
            }
            
#ifdef userOutputFile
#include userOutputFile
#endif
            
            f_file<<std::endl;
            
        }
        
    };
    
    // Declare static data
    template <typename DislocationNetworkType>
    int DislocationNetworkIO<DislocationNetworkType>::outputFrequency=1;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputBinary=false;
    
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputGlidePlanes=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputSpatialCells=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputPKforce=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputElasticEnergy=true;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputMeshDisplacement=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputDislocationLength=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortionRate=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputQuadratureParticles=false;
    
    
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif

