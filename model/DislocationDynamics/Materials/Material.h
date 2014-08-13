/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_MATERIAL_H_
#define model_MATERIAL_H_

#include <cmath>
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/MPI/MPIcout.h> // defines mode::cout



namespace model {

    /**************************************************************************/
    /**************************************************************************/
    template<typename SymmetryType>
    class Material { };
    
    /**************************************************************************/
    /**************************************************************************/
    template<>
    class Material<Isotropic> {
        
        
        static int selectedMaterial;
        
        
        /**********************************************************************/
        template<int Z>
        static void select()
        {
            typedef PeriodicElement<Z,Isotropic> IM;
            
            mu=1.0;     // mu is used to normalized stress
            b =1.0;     // b is used to normalized length
            cs =1.0;    // shear wave speed is used to normalized velocity
            
            //            B =1.0;
//            Binv=1.0/B;
//            rho=IM::mu*std::pow(IM::b/IM::B,2)*IM::rho;  //! rho* = mu*(b/B)^2 * rho
            //        cs=IM::B/IM::b*std::pow(IM::rho*IM::mu,-0.5);
 //           cs=sqrt(mu/rho); // sc*=sqrt(mu*/rho*)
            nu=IM::nu;
            lambda=2.0*mu*nu/(1.0-2.0*nu);
            C1=1.0-nu;
            C2=1.0/(4.0*M_PI*C1);
            C3=1.0-2.0*nu;
            C4=0.5*C2;
            kT=1.38e-23/IM::mu/std::pow(IM::b,3)*T; // [-]
            
            dH0=PeriodicElement<Z,Isotropic>::dH0/IM::mu/std::pow(IM::b,3); // [-]
            p=PeriodicElement<Z,Isotropic>::p; // [-]
            q=PeriodicElement<Z,Isotropic>::q; // [-]
            A=PeriodicElement<Z,Isotropic>::A / IM::mu / (IM::b/sqrt(IM::mu/IM::rho));   // [1/K]
            tauP=PeriodicElement<Z,Isotropic>::tauP / IM::mu ;  // [-]
            
            std::string magentaColor    = "\033[0;35m";   // a magenta color
            std::string defaultColor    = "\033[0m";	   // the default color for the console
            model::cout<<magentaColor<<"Material is now: "<<IM::name<<defaultColor<<std::endl;
            
        }
        
        
    public:
        
        enum{Al=13, Ni=28, Cu=29, W=74, Fe=26};
        
        static double mu;
        static double b;
        static double B;
        static double Binv;
        static double rho;
        static double cs;
        static double nu;
        static double lambda;
        static double C1;
        static double C2;
        static double C3;
        static double C4;
        static double T;
        static double kT;
//        static double kb;
        static double tauIII;
        static double vAct;
        
        static double p;
        static double q;
        static double A;
//        static double Ta;
        static double tauP;

        static Eigen::Matrix<double,Eigen::Dynamic,2> dH0;

        /**********************************************************************/
        static void select(const unsigned int& Z)
        {/*!\param[in] Z the atomic number of the element to be selected
          *
          * Selects the element Z as the current material. All material properties
          * are updated accordingly.
          */
        
            switch (Z) {
                case Al:
                    selectedMaterial=Al;
                    select<Al>();
                    break;
                case Ni:
                    selectedMaterial=Ni;
                    select<Ni>();
                    break;
                case Cu:
                    selectedMaterial=Cu;
                    select<Cu>();
                    break;
                case W:
                    selectedMaterial=W;
                    select<W>();
                    break;
                case Fe:
                    selectedMaterial=Fe;
                    select<Fe>();
                    break;
                    
                default:
                    assert(0 && "Material not implemented.");
                    break;
            }
                        
        }
        
        /**********************************************************************/
        template <int dim>
        static void rotateCrystal(const Eigen::Matrix<double,dim,dim>& C2G)
        {/*!\param[in] C2G the crystal-to-global rotation matrix
          *
          * Rotates the default plane normals in 
          * PeriodicElement<Z,Isotropic>::CrystalStructure, where Z is the current 
          * selected material.
          */
            
            switch (selectedMaterial)
            {
                case Al:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Al,Isotropic>::CrystalStructure>(C2G);
                    break;
                case Ni:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Ni,Isotropic>::CrystalStructure>(C2G);
                    break;
                case Cu:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Cu,Isotropic>::CrystalStructure>(C2G);
                    break;
                case W:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<W,Isotropic>::CrystalStructure>(C2G);
                    break;
                case Fe:
                    CrystalOrientation<dim>::template rotate<PeriodicElement<Fe,Isotropic>::CrystalStructure>(C2G);
                    break;
                default:
                    assert(0 && "Material not implemented.");
                    break;
            }
        }
 
    };
    
    
    //double Material<Isotropic>::a2=1.0;  // square of core size a
    int Material<Isotropic>::selectedMaterial=29;
    double Material<Isotropic>::mu=1.0;
    double Material<Isotropic>::b=1.0;
    double Material<Isotropic>::B=1.0;
    double Material<Isotropic>::Binv=1.0;
    double Material<Isotropic>::rho=1.0;
    double Material<Isotropic>::cs=1.0;
    double Material<Isotropic>::nu=0.34;
    double Material<Isotropic>::lambda=2.0*1.0*0.34/(1.0-2.0*0.34);
    double Material<Isotropic>::C1=1.0-0.34;    // 1-nu
    double Material<Isotropic>::C2=1.0/(4.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
    double Material<Isotropic>::C3=1.0-2.0*0.34; // 1-2*nu
    double Material<Isotropic>::C4=1.0/(8.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
	double Material<Isotropic>::T=300.0;  // Temperature [K]
	//double Material<Isotropic>::kb=1.38e-23/48.0e9/pow(0.2556e-9,3);  // boltzmann constant
	double Material<Isotropic>::kT=1.38e-23/48.0e9/pow(0.2556e-9,3)*300.0;  // boltzmann constant * Temperature
	double Material<Isotropic>::tauIII=0.667e-3;  // critical resovled shear stress in stage III
	double Material<Isotropic>::vAct=300.0;  // Activation volume [b^3]

    double Material<Isotropic>::p=PeriodicElement<29,Isotropic>::p;
    double Material<Isotropic>::q=PeriodicElement<29,Isotropic>::q;
    double Material<Isotropic>::A=PeriodicElement<29,Isotropic>::A;
    double Material<Isotropic>::tauP=PeriodicElement<29,Isotropic>::tauP;
    
    Eigen::Matrix<double,Eigen::Dynamic,2> Material<Isotropic>::dH0=PeriodicElement<29,Isotropic>::dH0;

    
    
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif
