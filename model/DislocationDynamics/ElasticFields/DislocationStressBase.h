/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationStressBase_h_
#define _model_DislocationStress_h_

#include <Eigen/StdVector>
#include <Eigen/StdDeque>
#include <tuple>
#include <FieldBase.h>
#include <FieldPoint.h>
//#include <StressStraight.h>
#include <StressStraight.h>

//#include <BoundaryDislocationNetwork.h>


namespace model
{
    
    /*!\brief Class template that implements the calculation of the stress field
     * generated by a dislocation. Depending on the value of the macro
     * _MODEL_NON_SINGULAR_DD_, different dislocation theories are used for the
     * calculation.
     */
    template<short unsigned int _dim>
    struct DislocationStressBase
    /* inheritance */ : public FieldBase<double,_dim,_dim>
    {
        
        constexpr static int dim=_dim;
        typedef DislocationStress<_dim> DislocationStressType;
        typedef FieldBase<double,_dim,_dim> FieldBaseType;
        typedef typename FieldBaseType::MatrixType MatrixType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        typedef Material<dim,Isotropic> MaterialType;
        
        static bool use_multipole;
        
        //! Dislocation core size
        static  double a;
        
        //! Dislocation core size squared
        static  double a2;
        
        //! dim x dim identity matrix
        static const Eigen::Matrix<double,_dim,_dim> I;
        
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            
            a=TextFileParser(fileName).readScalar<double>("coreSize",true);
            assert(a>0.0 && "coreSize MUST BE > 0.");
            a2=a*a;
        }
        
        /* axialVector ******************************************/
        static VectorDim axialVector(const MatrixDim& M)
        {
            return (VectorDim()<< M(1,2)-M(2,1), M(2,0)-M(0,2), M(0,1)-M(1,0) ).finished();
        }
        
        /* skewMatrix ******************************************/
        static MatrixDim skewMatrix(const VectorDim& v)
        {
            return (MatrixDim()<<  0.0,  v(2), -v(1),
                    /*                */ -v(2),   0.0,  v(0),
                    /*                */  v(1), -v(0),   0.0).finished();
        }
        
        /**********************************************************************/
        static MatrixType get(const MatrixType& temp)
        {
            return MaterialType::C2 * (temp+temp.transpose());
        }
        
        /**********************************************************************/
        template <typename ParticleType>
        static MatrixType addSourceContribution(const ParticleType& field,
                                                const std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>>& ssdeq)
        {
			MatrixType temp(MatrixType::Zero());
			for (const auto& sStaight: ssdeq)
			{
                temp +=sStaight.nonSymmStress(field.P);
	        }
	        return temp;
        }
        
        /**********************************************************************/
        template <typename ParticleType>
        static MatrixType addSourceContribution(const ParticleType&)
        {
            return MatrixType::Zero();
        }
        
        /**********************************************************************/
        template <typename ParticleType, typename CellContainerType>
        static MatrixType multipole(const ParticleType& field,const CellContainerType& farCells)
        {/*!@param[in] field  the FieldPoint at which stress is computed
          * @param[in] farCells container of SpatialCell(s) that are not neighbors of field
          *\returns the Cauchy stress contribution of the farCells on field.
          *
          * \f[
          * \begin{align}
          * \sigma_{ij}=\oint_\mathcal{L}S_{ijkl}R(\mathbf x-\mathbf x')b'_k\xi'_l\ dL'
          * \end{align}
          * \f]
          * Now let \f$\mathbf{x'}=\mathbf{x}^c-\tilde{\mathbf{x}}\f$, then
          * \f[
          * \begin{align}
          * \sigma_{ij}=\oint_\mathcal{L}S_{ijkl}R(\mathbf x-\mathbf x^c+\tilde{\mathbf{x}})b'_k\xi'_l\ dL'
          * \end{align}
          * \f]
          * Now assume
          * \f[
          * \begin{align}
          * |\tilde{\mathbf{x}}|\ll |\mathbf x-\mathbf x^c|
          * \end{align}
          * \f]
          * then, expanding \f$R\f$ for \f$\tilde{\mathbf{x}}\rightarrow0\f$ yields
          * \f[
          * \begin{align}
          * R(\mathbf x-\mathbf x^c+\tilde{\mathbf{x}})\approx R(\mathbf x-\mathbf x^c)+\partial_pR(\mathbf x-\mathbf x^c)\tilde{x}_p+\ldots
          * \end{align}
          * \f]
          * and
          * \f[
          * \begin{align}
          * \sigma_{ij}=S_{ijkl}R(\mathbf x-\mathbf x^c)\oint_\mathcal{L}b'_k\xi'_l\ dL'+S_{ijkl}\partial_pR(\mathbf x-\mathbf x^c)\oint_\mathcal{L}\tilde{x}_pb'_k\xi'_l+\ldots
          * \end{align}
          * \f]
          * First order expansion:
          * \f[
          * \begin{align}
          * \mathbf\sigma^0(\mathbf x)=\frac{\mu}{4\pi(1-\nu)R^2}&\left\{
          * (1-\nu)\left[\hat{\mathbf\xi}'\otimes(\mathbf b\times\hat{\mathbf R})+(\mathbf b\times\hat{\mathbf R})\otimes\hat{\mathbf\xi}'\right]
          * +\left[(\hat{\mathbf \xi}'\times\mathbf b)\otimes\hat{\mathbf R}+\hat{\mathbf R}\otimes(\hat{\mathbf \xi}'\times\mathbf b)\right] \right. \nonumber\\
          * &\left.+\hat{\mathbf R}\cdot(\mathbf b\times\hat{\mathbf\xi}')\left[3\hat{\mathbf R}\otimes\hat{\mathbf R}+\mathbf I \right]
          * \right\}
          * \end{align}
          *\f]
          * Letting
          * \f[
          * \begin{align}
          * (\mathbf b\times\hat{\mathbf R})\otimes\hat{\mathbf\xi})_{ij}=\epsilon_{ikm}b_kR_m\xi_j=\epsilon_{ikm}R_m\alpha_{kj}=\mathbf{S}\cdot\mathbf{\alpha}
          * \end{align}
          * \f]
          * \f[
          * \begin{align}
          * \mathbf{S}=\left[\begin{array}{ccc}
          * 0&R_3&-R_2\\
          * -R_3&0&R_1\\
          * R_2&-R_1&0
          * \end{array}\right]
          * \end{align}
          * \f]
          * \f[
          * \begin{align}
          * \mathbf b\times\hat{\mathbf \xi}'=\epsilon_{ikm}b_k\xi_m=a_i
          * \end{align}
          * \f]
          * \f[
          * \begin{align}
          * \mathbf{a}=\left[\begin{array}{c}
          * \alpha_{23}-\alpha_{32}\\
          * \alpha_{31}-\alpha_{13}\\
          * \alpha_{12}-\alpha_{21}
          * \end{array}\right]
          * \end{align}
          * \f]
          * we obtain
          * \f[
          * \begin{align}
          * \mathbf\sigma^0(\mathbf x)=\frac{\mu}{4\pi(1-\nu)R^2}&\left\{
          * (1-\nu)\left[(\mathbf{S}\cdot\mathbf{\alpha})^T+\mathbf{S}\cdot\mathbf{\alpha}\right]
          * -\left[\mathbf{a}\otimes\hat{\mathbf R}+\hat{\mathbf R}\otimes\mathbf{a}\right] \right. \nonumber\\
          * &\left.+(\hat{\mathbf R}\cdot\mathbf a)\left[3\hat{\mathbf R}\otimes\hat{\mathbf R}+\mathbf I \right]
          * \right\}
          * \end{align}
          * \f]
          */
            MatrixType temp(MatrixType::Zero());
            for(auto cell : farCells)
            {
#ifdef _MODEL_ENABLE_CELL_VERTEX_ALPHA_TENSORS_
                for(size_t v=0;v<Pow<2,dim>::value;++v)
                {
                    VectorDim R(field.P-cell.second->vertices().col(v));
                    const double R2(R.squaredNorm());
                    R/=sqrt(R2); // normalize R;
                    const MatrixDim& alpha(std::get<0>(*cell.second)[v]);
                    const VectorDim a(axialVector(alpha));
                    const MatrixDim S(skewMatrix(R));
                    temp += (MaterialType::C1*S*alpha-a*R.transpose()+0.5*R.dot(a)*(3.0*R*R.transpose()+I))/R2;
                    
                }
#else
                VectorDim R(field.P-cell.second->center);
                const double R2(R.squaredNorm());
                R/=sqrt(R2); // normalize R;
                const MatrixDim& alpha(std::get<0>(*cell.second));
                const VectorDim a(axialVector(alpha));
                const MatrixDim S(skewMatrix(R));
                temp += (MaterialType::C1*S*alpha-a*R.transpose()+0.5*R.dot(a)*(3.0*R*R.transpose()+I))/R2;
#endif
            }
            return temp;
        }
        
    };
    
    /**************************************************************************/
    // Static data members
    
    template<short unsigned int _dim>
    bool DislocationStressBase<_dim>::use_multipole=true;
    
    //! Dislocation core size
    template<short unsigned int _dim>
    double DislocationStressBase<_dim>::a=1.0;
    
    //! Dislocation core size squared
    template<short unsigned int _dim>
    double DislocationStressBase<_dim>::a2=1.0;
    
    //! Identity matrix
    template<short unsigned int _dim>
    const Eigen::Matrix<double,_dim,_dim> DislocationStressBase<_dim>::I=Eigen::Matrix<double,_dim,_dim>::Identity();
    
}	// close namespace
#endif



