#ifndef ADMMOGLASSOFUSEDTALL_H
#define ADMMOGLASSOFUSEDTALL_H

#include "FADMMBaseTwoPenalties.h"
#include "Linalg/BlasWrapper.h"
#include "Spectra/SymEigsSolver.h"
#include "ADMMMatOp.h"
#include "utils.h"

using Rcpp::IntegerVector;
using Rcpp::CharacterVector;

using namespace Spectra;


// minimize  1/2 * ||y - X * beta||^2 + lambda * ||beta||_1
//
// In ADMM form,
//   minimize f(x) + g(z)
//   s.t. x - z = 0
//
// x => beta
// z => -X * beta
// A => X
// b => y
// f(x) => 1/2 * ||Ax - b||^2
// g(z) => lambda * ||z||_1
class ADMMogLassoFusedTall: public FADMMBaseTwo<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>
{
protected:
    typedef float Scalar;
    typedef double Double;
    typedef Eigen::Matrix<Double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Double, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<const Matrix> MapMat;
    typedef Eigen::Map<const Vector> MapVec;
    typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMatR;
    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;
    typedef const Eigen::Ref<const Vector> ConstGenericVector;
    typedef Eigen::SparseMatrix<double> SpMat;
    typedef Eigen::SparseMatrix<int, Eigen::RowMajor> SpMatIntR;
    typedef Eigen::SparseVector<double> SparseVector;
    typedef Eigen::LLT<Matrix> LLT;

    MapMat datX;                  // data matrix
    MapVec datY;                  // response vector
    const SpMatR C;               // pointer to C matrix
    const SpMat D;                // pointer to D matrix
    //const MapVec D;             // pointer D vector
    // VectorXd D;

    int nobs;                     // number of observations
    int nvars;                    // number of variables
    int M;                        // length of nu (total size of all groups)
    int ngroups;                  // number of groups
    int nfused;                   // number of fused penalties

    Vector XY;                    // X'Y
    MatrixXd XX;                  // X'X
    SparseMatrix<double,Eigen::ColMajor> CCol;
    VectorXd CC;                  // C'C diagonal
    SpMat DD;                     // pointer to D matrix
    VectorXd Cbeta;               // C * beta
    VectorXd Dbeta;               // D * beta

    VectorXd group_weights;       // group weight multipliers
    VectorXd fused_weights;       // fused lasso weight multipliers
    CharacterVector family;       // model family (gaussian, binomial, or Cox PH)
    IntegerVector group_idx;      // indices of groups

    double newton_tol;            // tolerance for newton iterations
    int newton_maxit;             // max # iterations for newton-raphson
    bool dynamic_rho;

    Scalar lambda0;               // minimum lambda to make coefficients all zero

    //Eigen::DiagonalMatrix<double, Eigen::Dynamic> one_over_D_diag; // diag(1/D)

    double avg_group_weights;
    Scalar lambda;                // L1 penalty
    Scalar lambdaF;               // Fused penalty
    bool rho_unspecified;         // was rho unspecified? if so, we must set it
    LLT solver;                   // matrix factorization
    VectorXd savedEigs;           // saved eigenvalues


    virtual void block_soft_threshold(VectorXd &gammavec, VectorXd &d,
                                      const double &lam, const double &step_size)
    {
        // This thresholding function is for the most
        // basic overlapping group penalty, the
        // l1/l2 norm penalty, ie
        //     lambda * sqrt(beta_1 ^ 2 + beta_2 ^ 2 + ...)

        // d is the vector to be thresholded
        // gammavec is the vector to be written to

        int itrs = 0;

        for (int g = 0; g < ngroups; ++g)
        {
            double ds_norm = (d.segment(group_idx(g), group_idx(g+1) - group_idx(g))).norm();
            double thresh_factor = std::max(0.0, 1 - step_size * lam * group_weights(g) / (ds_norm) );

            for (int gr = group_idx(g); gr < group_idx(g+1); ++gr)
            {
                gammavec(itrs) = thresh_factor * d(gr);
                ++itrs;
            }
        }
    }


    // x -> Ax
    void A_mult (Vector &res, Vector &beta)  { res.swap(beta); }
    // y -> A'y
    void At_mult(Vector &res, Vector &nu)  { res.swap(nu); }
    // z -> Bz
    void B_mult (Vector &res, Vector &gamma) { res = -gamma; }
    // ||c||_2
    double c_norm() { return 0.0; }


    virtual void soft_threshold(Vector &res, const Vector &vec, const double &penalty)
    {
        int v_size = vec.size();
        res.setZero();

        //const double *ptr = vec.data();
        for(int i = 0; i < v_size; i++)
        {
            double pen_fact = penalty * std::abs(fused_weights(i));

            // if(ptr[i] > pen_fact)
            //     res(i) = ptr[i] - pen_fact;
            // else if(ptr[i] < -pen_fact)
            //     res(i) = ptr[i] + pen_fact;
            
            if(vec(i) > pen_fact)
                res(i) = vec(i) - pen_fact;
            else if(vec(i) < -pen_fact)
                res(i) = vec(i) + pen_fact;
        }
    }

    static void soft_threshold(SparseVector &res, const Vector &vec, const double &penalty)
    {
        int v_size = vec.size();
        res.setZero();
        res.reserve(v_size);

        const double *ptr = vec.data();
        for(int i = 0; i < v_size; i++)
        {
            if(ptr[i] > penalty)
                res.insertBack(i) = ptr[i] - penalty;
            else if(ptr[i] < -penalty)
                res.insertBack(i) = ptr[i] + penalty;
        }
    }

    void next_beta(Vector &res)
    {
        Vector rhs =   XY - CCol.adjoint() * adj_nu.head(M) - D.adjoint() * adj_nu.tail(nfused);
        
        rhs       += rho * (CCol.adjoint() * adj_gamma      + D.adjoint() * adj_eta);

        res.noalias() = solver.solve(rhs);
    }

    virtual void next_gamma(Vector &res)
    {
        Cbeta      = CCol * main_beta;
        Vector vec = Cbeta + adj_nu.head(M) / rho;
   
        block_soft_threshold(res, vec, lambda, 1 / rho);
    }

    virtual void next_eta(Vector &res)
    {
        Dbeta      = D * main_beta;
        Vector vec = Dbeta + adj_nu.tail(nfused) / rho;
        soft_threshold(res, vec, lambdaF / rho);
    }

    void next_residual(Vector &res)
    {
        res.head(M)  = Cbeta;
        res.head(M) -= aux_gamma;
        
        res.tail(nfused)  = Dbeta;
        res.tail(nfused) -= aux_eta;
    }
    void rho_changed_action() {}
    void update_rho() {}


    // Calculate ||v1 - v2||^2 when v1 and v2 are sparse
    static double diff_squared_norm(const SparseVector &v1, const SparseVector &v2)
    {
        const int n1 = v1.nonZeros(), n2 = v2.nonZeros();
        const double *v1_val = v1.valuePtr(), *v2_val = v2.valuePtr();
        const int *v1_ind = v1.innerIndexPtr(), *v2_ind = v2.innerIndexPtr();

        Scalar r = 0.0;
        int i1 = 0, i2 = 0;
        while(i1 < n1 && i2 < n2)
        {
            if(v1_ind[i1] == v2_ind[i2])
            {
                Scalar val = v1_val[i1] - v2_val[i2];
                r += val * val;
                i1++;
                i2++;
            } else if(v1_ind[i1] < v2_ind[i2]) {
                r += v1_val[i1] * v1_val[i1];
                i1++;
            } else {
                r += v2_val[i2] * v2_val[i2];
                i2++;
            }
        }
        while(i1 < n1)
        {
            r += v1_val[i1] * v1_val[i1];
            i1++;
        }
        while(i2 < n2)
        {
            r += v2_val[i2] * v2_val[i2];
            i2++;
        }

        return r;
    }

    // Faster computation of epsilons and residuals
    double compute_eps_primal()
    {
        double r = std::max(Dbeta.norm(), std::max(Cbeta.norm(), aux_gamma.norm()));
        return r * eps_rel + std::sqrt(double(dim_dual)) * eps_abs;
    }


public:
    ADMMogLassoFusedTall(const Eigen::Ref<const MatrixXd>  &datX_, //ConstGenericMatrix &datX_,
                         ConstGenericVector &datY_,
                         const SpMatR &C_,// const VectorXd &D_,
                         const SpMat &D_,
                         const int nobs_, 
                         const int nvars_, 
                         const int M_,
                         const int ngroups_, 
                         const int nfused_,
                         Rcpp::CharacterVector family_,
                         VectorXd group_weights_,
                         VectorXd fused_weights_,
                         Rcpp::IntegerVector group_idx_,
                         const bool dynamic_rho_,
                         const double newton_tol_ = 1e-5,
                         const int newton_maxit_ = 100,
                         const double eps_abs_ = 1e-6,
                         const double eps_rel_ = 1e-6) :
    FADMMBaseTwo<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>
             (datX_.cols(), C_.rows(), D_.rows(), C_.rows() + D_.rows(),
              eps_abs_, eps_rel_),
              datX(datX_.data(), datX_.rows(), datX_.cols()),
              datY(datY_.data(), datY_.size()),
              C(C_),
              D(D_),
              nobs(nobs_),
              nvars(nvars_),
              M(M_),
              ngroups(ngroups_),
              nfused(nfused_),
              XY(datX.transpose() * datY),
              XX(XtX(datX)),
              CCol(Eigen::SparseMatrix<double>(M_, nvars_)),
              CC(nvars_),
              DD(XtX(D_)),
              Cbeta(C_.rows()),
              Dbeta(D_.rows()),
              group_weights(group_weights_),
              fused_weights(fused_weights_),
              family(family_),
              group_idx(group_idx_),
              newton_tol(newton_tol_),
              newton_maxit(newton_maxit_),
              dynamic_rho(dynamic_rho_),
              lambda0(XY.cwiseAbs().maxCoeff())
    { }

    double get_lambda_zero() const { return lambda0; }

    // init() is a cold start for the first lambda
    void init(double lambda_, double rho_, double lambdaF_)
    {
        main_beta.setZero();
        aux_gamma.setZero();
        aux_eta.setZero();
        dual_nu.setZero();

        adj_gamma.setZero();
        adj_eta.setZero();
        adj_nu.setZero();

        lambda  = lambda_;
        lambdaF = lambdaF_;
        rho     = rho_;

        // store ColMajor version of C
        CCol = C;
        

        // create vector CC, whose elements are the number of times
        // each variable is in a group
        for (int k=0; k < CCol.outerSize(); ++k)
        {
            double tmp_val = 0;
            for (SparseMatrix<double>::InnerIterator it(CCol,k); it; ++it)
            {
                tmp_val += it.value();
            }
            CC(k) = tmp_val;
        }
        

        //Linalg::cross_prod_lower(XX, datX);

        if(rho <= 0)
        {
            rho_unspecified = true;
            // Construct matrix operation object using the wrapper class DenseSymMatProd
            DenseSymMatProd<double> op(XX);
            
            // Construct eigen solver object, requesting the largest three eigenvalues
            SymEigsSolver<DenseSymMatProd<double>> eigs(op, 2, 5);
            
            // Initialize and compute
            eigs.init();
            int nconv = eigs.compute(SortRule::LargestAlge);
            
            // Retrieve results
            Vector evals;
            if(eigs.info() == CompInfo::Successful)
                evals = eigs.eigenvalues();
            
            savedEigs = evals;

            //float lam_fact = lambda;
            rho = std::pow(evals[0], 1.0 / 3.0) * std::pow(lambda, 2.0 / 3.0);
            /*
            if (lam_fact < savedEigs[1])
            {
                rho = std::sqrt(savedEigs[1] * std::pow(lam_fact * 4, 1.35));
            } else if (lam_fact * 0.25 > savedEigs[0])
            {
                rho = std::sqrt(savedEigs[1] * std::pow(lam_fact * 0.25, 1.35));
            } else
            {
                rho = std::pow(lam_fact, 1.05);
            }
             */
        } else {
            rho_unspecified = false;
        }

        //XX.diagonal().array() += rho;

        MatrixXd matToSolve(XX);
        matToSolve.diagonal() += rho * CC;
        matToSolve += rho * DD;
        

        // precompute LLT decomposition of (X'X + rho * (D'D + C'C))
        solver.compute(matToSolve.selfadjointView<Eigen::Lower>());

        eps_primal = 0.0;
        eps_dual = 0.0;
        resid_primal = 1e30;
        resid_dual = 1e30;

        adj_a = 1.0;
        adj_c = 1e30;

        rho_changed_action();
    }
    // when computing for the next lambda, we can use the
    // current main_beta, aux_gamma, dual_nu and rho as initial values
    void init_warm(double lambda_, double lambdaF_)
    {

        main_beta.setZero();
        aux_gamma.setZero();
        aux_eta.setZero();
        dual_nu.setZero();

        adj_gamma.setZero();
        adj_eta.setZero();
        adj_nu.setZero();

        lambda  = lambda_;
        lambdaF = lambdaF_;

        if (rho_unspecified)
        {
            //float lam_fact = lambda;
            rho = std::pow(savedEigs[0], 1.0 / 3.0) * std::pow(lambda, 2.0 / 3.0);
            /*
            if (lam_fact < savedEigs[1])
            {
                rho = std::sqrt(savedEigs[1] * std::pow(lam_fact * 4, 1.35));
            } else if (lam_fact * 0.25 > savedEigs[0])
            {
                rho = std::sqrt(savedEigs[1] * std::pow(lam_fact * 0.25, 1.35));
            } else
            {
                rho = std::pow(lam_fact, 1.05);
            }
             */

        }
        MatrixXd matToSolve(XX);
        matToSolve.diagonal() += rho * CC;
        matToSolve += rho * DD;

        // precompute LLT decomposition of (X'X + rho * (C'C + D'D))
        solver.compute(matToSolve.selfadjointView<Eigen::Lower>());

        eps_primal = 0.0;
        eps_dual = 0.0;
        resid_primal = 1e30;
        resid_dual = 1e30;

        adj_a = 1.0;
        adj_c = 1e30;
        rho_changed_action();
    }

    virtual double get_loss()
    {
        double loss;
        loss = (datY - datX * main_beta).array().square().sum();
        return loss;
    }

    virtual void update_adaptive_group_weights(VectorXd &weights_)
    {
        group_weights = weights_;
    }

    virtual void update_adaptive_fused_weights(VectorXd &weights_)
    {
        fused_weights = weights_;
    }

    // this function actually returns beta
    virtual VectorXd get_gamma()
    {
        VectorXd beta_return(nvars);
        for (int k=0; k < CCol.outerSize(); ++k)
        {
            int rowidx = 0;
            bool current_zero = false;
            bool already_idx = false;
            for (SparseMatrix<double>::InnerIterator it(CCol,k); it; ++it)
            {

                if (aux_gamma(it.row()) == 0.0 && !current_zero)
                {
                    rowidx = it.row();
                    current_zero = true;
                } else if (!current_zero && !already_idx)
                {
                    rowidx = it.row();
                    already_idx = true;
                }


            }
            beta_return(k) = aux_gamma(rowidx);
        }
        return beta_return;
    }

    // this function returns gamma
    virtual VectorXd get_aux_gamma()
    {
        return aux_gamma;
    }

    virtual MatrixXd get_hessian()
    {
        return XX;
    }

    virtual VectorXd get_xty()
    {
        return XY;
    }
    
    virtual int get_nselected(VectorXd &beta_vector)
    {
        int len = beta_vector.size();
        int nselected = 0;
        for (int j = 0; j < len; ++j)
        {
            if (beta_vector(j) != 0)
            {
                nselected++;    
            }
        }
        return nselected;
    }

};



#endif // ADMMOGLASSOFUSEDTALL_H
