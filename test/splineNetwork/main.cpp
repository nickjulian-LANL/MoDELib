/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>

//#include <Eigen/Dense>
#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/StaticID.h>


const int corder=1;

namespace model
{
    struct Dnetwork;
    struct Dnode;
//    struct DLnode;
    struct Dlink;
    struct Dloop;
    
    template<>
    struct TypeTraits<Dnetwork>
    {
        typedef Dnetwork LoopNetworkType;
        typedef Dnode  NodeType;
//        typedef DLnode LoopNodeType;
        typedef Dlink LinkType;
        typedef Dloop LoopType;
        typedef double FlowType;
        
        static constexpr FlowType zeroFlow=0.0;
    };
    
    template<>
    struct TypeTraits<Dnode> : public TypeTraits<Dnetwork>
    {
        
    };
    
    template<>
    struct TypeTraits<Dlink> : public TypeTraits<Dnetwork>
    {
        
    };
    
    template<>
    struct TypeTraits<Dloop> : public TypeTraits<Dnetwork>
    {
        
    };
    
//    template<>
//    struct TypeTraits<DLnode> : public TypeTraits<Dnetwork>
//    {
//        
//    };
}


#include <model/LoopNetwork/LoopNetwork.h>
#include <model/Geometry/Splines/SplineNode.h>
#include <model/Geometry/Splines/SplineSegment.h>
#include <model/IO/SequentialOutputFile.h>

//#include <Loop.h>



namespace model
{
    struct Dnetwork : public LoopNetwork<Dnetwork>{};
    
    struct Dnode : public SplineNode<Dnode,3,corder,CatmullRom>
    {
    
        Dnode(const Eigen::Vector3d& P ) : SplineNode<Dnode,3,corder,CatmullRom>(P){}
        
//    Node(const int& a)
//        {
//        }
        
    };
    
//    struct DLnode : public LoopNode<DLnode>
//    {
//        
//        //    Node(const int& a)
//        //        {
//        //        }
//        
//    };
    
    struct Dlink : public SplineSegment<Dlink,3,corder>//NetworkLink<Dlink>
    {
    
        Dlink(const std::shared_ptr<Dnode>& Ni,
              const std::shared_ptr<Dnode>& Nj) : SplineSegment<Dlink,3,corder>(Ni,Nj){}
        
    };
    
    
    struct Dloop : public Loop<Dloop>
    {
    
        Dloop(const Dnetwork& net, const double& flow) : Loop<Dloop>(net,flow){}
    };


}

using namespace model;


int main()
{

    Dnetwork DN;
    
    DN.insertDanglingNode((Eigen::Vector3d()<<-1.0,-1.0,0.0).finished());
    DN.insertDanglingNode((Eigen::Vector3d()<< 0.0,-1.0,0.0).finished());
    DN.insertDanglingNode((Eigen::Vector3d()<< 0.0, 0.0,0.0).finished());
    DN.insertDanglingNode((Eigen::Vector3d()<<-1.0, 0.0,0.0).finished());
    std::vector<size_t> loop0={0,1,2,3};

    
    DN.insertDanglingNode((Eigen::Vector3d()<<-1.0, 0.0,0.0).finished());
    DN.insertDanglingNode((Eigen::Vector3d()<< 0.0, 0.0,0.0).finished());
    DN.insertDanglingNode((Eigen::Vector3d()<< 0.0, 1.0,0.0).finished());
    DN.insertDanglingNode((Eigen::Vector3d()<<-1.0, 1.0,0.0).finished());
    std::vector<size_t> loop1={4,5,6,7};
    
    DN.insertDanglingNode((Eigen::Vector3d()<< 0.0, 0.0,0.0).finished());
    DN.insertDanglingNode((Eigen::Vector3d()<< 1.0, 0.0,0.0).finished());
    DN.insertDanglingNode((Eigen::Vector3d()<< 1.0, 1.0,0.0).finished());
    DN.insertDanglingNode((Eigen::Vector3d()<< 0.0, 1.0,0.0).finished());
    std::vector<size_t> loop2={8,9,10,11};
    
    
//    std::vector<size_t> loop0={0,1,2,3,4,5,6,7};
//        std::vector<size_t> loop1={1,2,12,13,7,8,14,15};
//    std::vector<size_t> loop2={8,7,16,17,18,19};

    //  std::vector<size_t> loop1={};

    DN.insertLoop(loop0,0.2);
    DN.insertLoop(loop1,0.2);
    DN.insertLoop(loop2,0.2);
    DN.clearDanglingNodes();
    DN.printLoops();

    
    int np=100;
    double du=1.0/(np-1);
    
    SequentialOutputFile<'P',true> file0;
    
    for(const auto& link : DN.networkLinks())
    {
        for(int i=0;i<np;++i)
        {
            file0<<link.second->get_r(i*du).transpose()<<std::endl;
        }
    }
    
    
    DN.contractSecond(3,4);
    DN.contractSecond(2,5);
    DN.contractSecond(2,8);
    DN.contractSecond(6,11);
    
    SequentialOutputFile<'P',true> file1;
    
    for(const auto& link : DN.networkLinks())
    {
        for(int i=0;i<np;++i)
        {
            file1<<link.second->get_r(i*du).transpose()<<std::endl;
        }
    }

    DN.node(0).second->set_P((Eigen::Vector3d()<<-1.0,-2.0,0.0).finished());
    
    
    SequentialOutputFile<'P',true> file2;
    
    for(const auto& link : DN.networkLinks())
    {
        for(int i=0;i<np;++i)
        {
            file2<<link.second->get_r(i*du).transpose()<<std::endl;
        }
    }
    
    SequentialOutputFile<'P',true> file3;
    
    for(const auto& link : DN.networkLinks())
    {
        const auto pair=link.second->hermite2posMatrix();

        if((pair.first*pair.second-link.second->hermiteDofs()).norm()>FLT_EPSILON)
        {
            std::cout<<pair.first<<std::endl<<std::endl;
            std::cout<<pair.second<<std::endl<<std::endl;
            std::cout<<pair.first*pair.second<<std::endl<<std::endl;
            std::cout<<link.second->hermiteDofs()<<std::endl<<std::endl;

            assert(0);
            //std::cout<<"WARNING"<<std::endl;
        }
        

        for(int i=0;i<np;++i)
        {
            file3<<link.second->sf(i*du)*pair.first*pair.second<<std::endl<<std::endl;
        }
    }
    
    //
//    DN.printLoops();
//
// //   DN.contract(2,12);
////    DN.printLoops();
//
////    std::vector<size_t> loop1={0,3,4,5};
////    std::vector<size_t> loop2={6,7};
////    
////    
////    DN.insertLoop(loop0,0.2);
////    DN.insertLoop(loop1,4.5);
////    DN.insertLoop(loop2,2.3);
////
////
//    
    std::cout<<"# loops="<<DN.loops().size()<<std::endl;
    std::cout<<"# nodes="<<DN.nodes().size()<<std::endl;
    std::cout<<"# loopLinks="<<DN.loopLinks().size()<<std::endl;
    std::cout<<"# networkLinks="<<DN.links().size()<<std::endl;
//    std::cout<<"# networkComponents="<<DN.networkComponents().size()<<std::endl;

////
////    std::cout<<"Expanding"<<std::endl;
////    DN.expand(3,4);
////    
////    DN.contract(0,1);
////    DN.contract(0,2);
////    
////    DN.contract(0,5);
////    DN.contract(0,4);
////    DN.contract(0,8);
////    
////    std::cout<<"# nodes="<<DN.nodes().size()<<std::endl;
////    std::cout<<"# loopLinks="<<DN.loopLinks().size()<<std::endl;
////    std::cout<<"# networkLinks="<<DN.links().size()<<std::endl;
////    
////    DN.checkLoops();
//    
//    
////    for(const auto& node : DN.nodes())
////    {
////        std::cout<<"node "<<node.second->sID<<" neighbor-size="<<node.second->loopLinks().size()<<std::endl;
////    }
////
////    DN.nodes();


    

    return 0;
    
}