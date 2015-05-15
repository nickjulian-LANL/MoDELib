// BOUNDARY CONDITIONS
// Fix nodes at z=0
const size_t nodeList_0=finiteElement().template createNodeList<AtXmin<2>>();
Fix fix;
addDirichletCondition(nodeList_0,fix,{1,1,1},DN); // fix u1 u2 u3

// Displace top nodes using flat punch
if(DN.runningID()>=TensionTorsioner::relaxSteps)
{
    const size_t nodeList_1=finiteElement().template createNodeList<TensionTorsioner>();
    TensionTorsioner tt(finiteElement());
    tt.updateDisplacement(DN);
    
    if(TensionTorsioner::apply_tension &&  TensionTorsioner::apply_torsion)
    {
        addDirichletCondition(nodeList_1,tt,{1,1,1},DN); // prescribe u1, u2, u3
    }
    else if (TensionTorsioner::apply_tension && ! TensionTorsioner::apply_torsion)
    {
        addDirichletCondition(nodeList_1,tt,{0,0,1},DN); // prescribe only u3
    }
    else if (!TensionTorsioner::apply_tension &&  TensionTorsioner::apply_torsion)
    {
        addDirichletCondition(nodeList_1,tt,{1,1,0},DN); // prescribe only u1 & u2
    }
    else
    {
        std::cout<<"TensionTorsioner IS NOT APPLYING ANY LOADS"<<std::endl;
    }
}

// CALL SOLVER, passing the negative dislocation traction vector and the solution guess
displacement()=solve(-dislocationTraction.globalVector(),displacement());