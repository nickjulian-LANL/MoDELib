/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDauxIO_H_
#define model_DDauxIO_H_

#include <vector>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <cfloat>      // std::ifstream
#include <DDbaseIO.h>
#include <GlidePlaneBoundaryIO.h>
#include <GlidePlaneFactory.h>
#include <PeriodicGlidePlane.h>

namespace model
{
    
    
    template <int dim>
    struct DDauxIO : public DDbaseIO
    /*            */,private std::vector<GlidePlaneBoundaryIO<dim>>
    /*            */,private std::vector<PeriodicPlanePatchIO<dim>>
    {
        
        DDauxIO(const std::string& suffix="") :
        /* init */ DDbaseIO("evl","ddAux",suffix)
        {
            
        }
        
        /**********************************************************************/
        void setGlidePlaneBoundaries(const GlidePlaneFactory<dim>& dn)
        {
            glidePlanesBoundaries().clear();
            for(const auto& pair : dn.glidePlanes())
            {
                const auto glidePlane(pair.second.lock());
                if(glidePlane)
                {
                    for(const auto& seg : glidePlane->meshIntersections)
                    {
                        glidePlanesBoundaries().emplace_back(glidePlane->sID,*seg);
                    }
                }
            }
        }
        
        /**********************************************************************/
        void addPeriodicGlidePlane(const PeriodicGlidePlane<dim>& pgp)
        {
            for(const auto& patch : pgp.patches())
            {
                periodicGlidePlanePatches().emplace_back(*patch.second);
            }
//            std::cout<<"added "<<pgp.patches().size()<<" patches, now size="<<periodicGlidePlanePatches().size()<<std::endl;
        }
        
        /**********************************************************************/
        const std::vector<GlidePlaneBoundaryIO<dim>>& glidePlanesBoundaries() const
        {
            return *this;
        }
        
        std::vector<GlidePlaneBoundaryIO<dim>>& glidePlanesBoundaries()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::vector<PeriodicPlanePatchIO<dim>>& periodicGlidePlanePatches() const
        {
            return *this;
        }
        
        std::vector<PeriodicPlanePatchIO<dim>>& periodicGlidePlanePatches()
        {
            return *this;
        }
        
        /**********************************************************************/
        void writeTxt(const long int& runID)
        {
            
            const std::string filename(this->getTxtFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                // Write header
                file<<glidePlanesBoundaries().size()<<"\n";
                file<<periodicGlidePlanePatches().size()<<"\n";

                // Write Nodes
                for(const auto& gpBnd : glidePlanesBoundaries())
                {
                    file<<gpBnd<<"\n";
                }
                for(const auto& patch : periodicGlidePlanePatches())
                {
                    file<<patch<<"\n";
                }

                file.close();
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
        }
        
        /**********************************************************************/
        void writeBin(const size_t& runID)
        {
            
            const std::string filename(this->getBinFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                // Write header
                binWrite(file,glidePlanesBoundaries().size());
                binWrite(file,periodicGlidePlanePatches().size());

                // Write Nodes
                for(const auto& gpb : glidePlanesBoundaries())
                {
                    binWrite(file,gpb);
                }
                for(const auto& patch : periodicGlidePlanePatches())
                {
                    binWrite(file,patch);
                }
                
                file.close();
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
        }
        
        
        /**********************************************************************/
        void readBin(const size_t& runID)
        {
            const std::string filename(this->getBinFilename(runID));
            
            std::ifstream infile (filename.c_str(), std::ios::in|std::ios::binary);
            if(infile.is_open())
            {
                const auto t0=std::chrono::system_clock::now();
                model::cout<<"reading "<<filename<<std::flush;

                // Read header
                size_t sizeGP;
                infile.read (reinterpret_cast<char*>(&sizeGP), 1*sizeof(sizeGP));
                size_t sizePPP;
                infile.read (reinterpret_cast<char*>(&sizePPP), 1*sizeof(sizePPP));

                // Read body
                glidePlanesBoundaries().resize(sizeGP);
                infile.read (reinterpret_cast<char*>(glidePlanesBoundaries().data()),glidePlanesBoundaries().size()*sizeof(GlidePlaneBoundaryIO<dim>));
                periodicGlidePlanePatches().resize(sizePPP);
                infile.read (reinterpret_cast<char*>(periodicGlidePlanePatches().data()),periodicGlidePlanePatches().size()*sizeof(PeriodicPlanePatchIO<dim>));
                
                infile.close();
                printLog(t0);
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
            
        }

        
        /**********************************************************************/
        void readTxt(const size_t& runID)
        {
            const std::string filename(this->getTxtFilename(runID));
            
            std::ifstream infile (filename.c_str(), std::ios::in);
            if(infile.is_open())
            {

                const auto t0=std::chrono::system_clock::now();
                model::cout<<"reading "<<filename<<std::flush;
                std::string line;
                std::stringstream ss;

                size_t sizeGP;
                std::getline(infile, line);
                ss<<line;
                ss >> sizeGP;
                ss.clear();

                size_t sizePPP;
                std::getline(infile, line);
                ss<<line;
                ss >> sizePPP;
                ss.clear();

                glidePlanesBoundaries().clear();
                glidePlanesBoundaries().reserve(sizeGP);
                for(size_t k=0;k<sizeGP;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    glidePlanesBoundaries().emplace_back(ss);
                    ss.clear();
                }
                
                periodicGlidePlanePatches().clear();
                periodicGlidePlanePatches().reserve(sizePPP);
                for(size_t k=0;k<sizePPP;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    periodicGlidePlanePatches().emplace_back(ss);
                    ss.clear();
                }
                
                infile.close();
                printLog(t0);
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
            
        }
        
        void printLog(const std::chrono::time_point<std::chrono::system_clock>& t0) const
        {
            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            model::cout<<"  "<<glidePlanesBoundaries().size()<<" glidePlanesBoundaries "<<std::endl;
            model::cout<<"  "<<periodicGlidePlanePatches().size()<<" periodicPlanePatches"<<std::endl;
        }

        /**********************************************************************/
        void read(const size_t& runID)
        {
            if(isBinGood(runID))
            {
                readBin(runID);
            }
            else
            {
                if(isTxtGood(runID))
                {
                    readTxt(runID);
                }
                else
                {
                    std::cout<<"COULD NOT FIND INPUT FILEs evl/evl_"<<runID<<".bin or evl/evl_"<<runID<<".txt"<<std::endl;
                    assert(0 && "COULD NOT FIND INPUT FILEs.");
                }
            }
        }
        
    };
    
}
#endif


//        /**********************************************************************/
//        void writeTxt(const size_t& runID,
//                             const std::vector<GlidePlaneBoundaryIO<dim>> gpBoundaries)
//        {
//
//            const std::string filename(this->getTxtFilename(runID));
//            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
////            std::cout<<"Writing to "<<filename<<std::endl;
//            if(file.is_open())
//            {
//                // Write header
//                file<<gpBoundaries.size()<<"\n";
//
//                // Write Nodes
//                for(const auto& gpBnd : gpBoundaries)
//                {
//                    file<<gpBnd<<"\n";
//                }
//
//                file.close();
//            }
//            else
//            {
//                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
//                assert(false && "CANNOT OPEN FILE.");
//            }
//        }

//        /**********************************************************************/
//        template<typename DislocationNodeType>
//        void writeBin(const DislocationNodeType& dn,
//                             const long int& runID)
//        {
//
//            const std::string filename(this->getBinFilename(runID));
//            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
//            if(file.is_open())
//            {
//                // Write header
//                size_t nGP(0);
//                for(const auto& glidePlane : dn.glidePlanes())
//                {
//                    nGP+=glidePlane.second->meshIntersections.size();
//                }
//                binWrite(file,nGP);
//
//                // Write GlidePlaneBoundaries
//                for(const auto& glidePlane : dn.glidePlanes())
//                {
//                    for(const auto& seg : glidePlane->second->meshIntersections)
//                    {
//                        binWrite(file,GlidePlaneBoundaryIO<dim>(seg));
//                    }
//                }
//
//                file.close();
//            }
//            else
//            {
//                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
//                assert(false && "CANNOT OPEN FILE.");
//            }
//        }