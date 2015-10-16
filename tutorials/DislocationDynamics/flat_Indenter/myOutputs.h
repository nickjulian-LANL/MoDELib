/******************************************************************************/
// COMPUTATION OF DISPLACEMENT
double disp = 0.0;

// Obtain a list of nodes under the punch
const size_t nodeListID=DN.shared.bvpSolver.finiteElement().template createNodeList<FlatPunch>();

// Sum FEM displacement of those nodes
for(auto node : DN.shared.bvpSolver.finiteElement().nodeList(nodeListID))
{
    const Eigen::Matrix<double,dim,1> nodeDisp = DN.shared.bvpSolver.displacement().template segment<dim>(dim*node->gID);
    disp += nodeDisp(2);
}

// Compute dislocation displacement for those nodes
typedef BoundaryDisplacementPoint<DislocationNetworkType> FieldPointType;
typedef typename FieldPointType::DisplacementField DisplacementField;
std::deque<FieldPointType> fieldPoints; // the container of field points
for (auto node : DN.shared.bvpSolver.finiteElement().nodeList(nodeListID)) // range-based for loop (C++11)
{
    fieldPoints.emplace_back(*node);
}
DN.template computeField<FieldPointType,DisplacementField>(fieldPoints);

// Sum dislocation displacement to disp
for(auto node : fieldPoints)
{
    Eigen::Matrix<double,dim,1> nodeDisp = node.template field<DisplacementField>(); // line integral part of displacement
    
    if (DN.shared.use_virtualSegments) // solid-angle jump of virtual segments
    {
        for(const auto& segment : DN.links())
        {
            segment.second.addToSolidAngleJump(node.P,node.S,nodeDisp);
        }
    }
        
    disp += nodeDisp(2);
}

double avgdisp = disp/DN.shared.bvpSolver.finiteElement().nodeList(nodeListID).size();


/******************************************************************************/
// COMPUTATION OF LOAD
auto boundaryUnderPunch = FlatPunch::boundary(DN.shared.bvpSolver.finiteElement());
Eigen::Matrix<double,3,1> force(Eigen::Matrix<double,3,1>::Zero());
typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
boundaryUnderPunch.integrate(&DN.shared.bvpSolver,force,&BvpSolverType::bvpTraction); // integrate the bvp correction
boundaryUnderPunch.integrate(&DN.shared.bvpSolver,force,&BvpSolverType::ddTraction,DN); // integrate the dd traction


/******************************************************************************/
// PLASTIC STRAIN RATE
Eigen::Matrix<double,dim,dim> pSR(DN.plasticStrainRate());
std::pair<double,double> length=DN.networkLength(); // first is total immobilized length, second is total mobile length

/******************************************************************************/
// OUTPUT TO F_0.txt
UniqueOutputFile<'F'> f_file;
model::cout<<", F/F_0"<<std::flush;
f_file<< DN.runningID()<<" "<<DN.get_dt()<<"  "<<length.first<<" "<<length.second<<" "<<avgdisp<<"  "<<force(2)<<"  "<<pSR(0,0)<<"  "<<pSR(0,1)<<"  "<<pSR(0,2)<<"  "<<pSR(1,1)<<"  "<<pSR(1,2) <<"  "<<pSR(2,2)<<"  "<<std::endl;



