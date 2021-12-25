#ifndef ADMMOGLASSOCOXPHFUSEDTALL_H
#define ADMMOGLASSOCOXPHFUSEDTALL_H

#include "FADMMBaseTwoPenalties.h"
#include "Linalg/BlasWrapper.h"
#include "Spectra/SymEigsSolver.h"
#include "ADMMMatOp.h"
#include "utils.h"
#include "iostream"

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
class ADMMogLassoCoxPHFusedTall: public FADMMBaseTwo<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>
{
protected:
  typedef float Scalar;
  typedef double Double;
  typedef Eigen::Matrix<Double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Matrix<Double, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<int, Eigen::Dynamic, 1> Vectori;
  typedef Eigen::Map<const Matrix> MapMat;
  typedef Eigen::Map<const Vector> MapVec;
  typedef Eigen::Map<const Vectori> MapVeci;
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMatR;
  typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;
  typedef const Eigen::Ref<const Vector> ConstGenericVector;
  typedef const Eigen::Ref<const Vectori> ConstGenericVectori;
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::SparseMatrix<int, Eigen::RowMajor> SpMatIntR;
  typedef Eigen::SparseVector<double> SparseVector;
  typedef Eigen::LLT<Matrix> LLT;

  MapMat datX;                  // data matrix
  MapVec datY;                  // response vector
  MapVeci delta;                // indicator of event
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

  VectorXi CIndex;              // index of data with event
  LLT solver;                   // matrix factorization
  bool rho_unspecified;         // was rho unspecified? if so, we must set it
  Scalar lambda;                // L1 penalty
  Scalar lambdaF;               // Fused penalty
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
    block_soft_threshold(res, vec, lambda, 1/rho);
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

  void rho_changed_action()
  {
    MatrixXd matToSolve(XX);
    matToSolve.diagonal() += rho * CC;
    matToSolve += rho * DD;

    // precompute LLT decomposition of (X'X + rho * D'D)
    solver.compute(matToSolve.selfadjointView<Eigen::Lower>());
  }
  void update_rho() {}

  void compute_rho()
  {
    if (rho_unspecified)
    {
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

      float lam_fact = lambda;
      //rho = std::pow(evals[0], 1.0 / 3) * std::pow(lambda, 2.0 / 3);
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

    }
  }



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
    ADMMogLassoCoxPHFusedTall(ConstGenericMatrix &datX_,
                              ConstGenericVector &datY_,
                              ConstGenericVectori &delta_,
                              const SpMatR &C_,// const VectorXd &D_,
                              const SpMat &D_,
                              int nobs_, int nvars_, int M_,
                              int ngroups_, int nfused_,
                              Rcpp::CharacterVector family_,
                              VectorXd group_weights_,
                              VectorXd fused_weights_,
                              Rcpp::IntegerVector group_idx_,
                              bool dynamic_rho_,
                              double newton_tol_ = 1e-5,
                              int newton_maxit_ = 100,
                              double eps_abs_ = 1e-6,
                              double eps_rel_ = 1e-6) :
    FADMMBaseTwo<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>
    (datX_.cols(), C_.rows(), D_.rows(), C_.rows() + D_.rows(),
     eps_abs_, eps_rel_),
     datX(datX_.data(), datX_.rows(), datX_.cols()),
     datY(datY_.data(), datY_.size()),
     delta(delta_.data(), delta_.size()),
     C(C_),
     D(D_),
     nobs(nobs_),
     nvars(nvars_),
     M(M_),
     ngroups(ngroups_),
     nfused(nfused_),
     XY(datX.transpose() * datY),
     XX(datX_.cols(), datX_.cols()),
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
     lambda0(set_lambda_zero())
     { }


    double set_lambda_zero() {
        int nevent = delta.sum();
        double lambda00 = 0;
        VectorXd XY_lambda0(nvars);
        getCindex();
        XY_lambda0.fill(0);
        for (int i = 0; i <= nevent-1; i++){
            XY_lambda0 = XY_lambda0 + datX.row(CIndex(i)).transpose() - (datX.bottomRows(nobs-CIndex(i)).colwise().sum()/(nobs-CIndex(i))).transpose();
        }
        lambda00 = XY_lambda0.cwiseAbs().maxCoeff();
        return lambda00;
    }

    double get_lambda_zero() const{
        return lambda0;
    }


    // get the compact event index vector
    void getCindex()
    {
        int nevent = delta.sum();
        VectorXi cindex(nevent);
        nevent = 0;
        for (int i = 0; i < delta.size(); ++i)
        {
            if (delta(i)==1)
            {
                cindex(nevent) = i;
                nevent += 1;
            }
        }
        CIndex = cindex;
    }



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
        rho = rho_;


        // store ColMajor version of C
        CCol = C;

        // create vector CC, whose elements are the number of times
        // each variable is in a group
        for (int k = 0; k < CCol.outerSize(); ++k)
        {
            double tmp_val = 0;
            for (SparseMatrix<double>::InnerIterator it(CCol,k); it; ++it)
            {
                tmp_val += it.value();
            }
            CC(k) = tmp_val;
        }

        if(rho <= 0)
        {
            rho_unspecified = true;
        } else {
            rho_unspecified = false;
        }


        eps_primal = 0.0;
        eps_dual = 0.0;
        resid_primal = 1e30;
        resid_dual = 1e30;

        adj_a = 1.0;
        adj_c = 1e30;

    }
    // when computing for the next lambda, we can use the
    // current main_beta, aux_gamma, dual_nu and rho as initial values
    void init_warm(double lambda_, double lambdaF_)
    {
        lambda = lambda_;
        lambdaF = lambdaF_;

        eps_primal = 0.0;
        eps_dual = 0.0;
        resid_primal = 1e30;
        resid_dual = 1e30;

        adj_a = 1.0;
        adj_c = 1e30;
    }

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

    virtual int solve(int maxit)
    {

        VectorXd beta_prev;


        int i;
        int j;
        getCindex();
        for ( i = 0; i < newton_maxit; ++i)
        {
            int eventnn = delta.sum();
            int nobs = delta.size();
            VectorXd W;
            VectorXd z;
            VectorXd expeta;
            VectorXd risketa(eventnn);
            VectorXd factor1(nobs);
            VectorXd factor2(nobs);

            risketa.fill(0);
            factor1.fill(0);
            factor2.fill(0);


            beta_prev = main_beta;

            //Cindex:=index of only event;delta = indicator of death

            expeta = ((datX * main_beta).array()).exp();


            //calculate risk set

            for (int index = 0; index < eventnn - 1; ++index)
            {
                risketa(index) = expeta.segment(CIndex(index),CIndex(index+1)-CIndex(index)).sum();
            }

            for (int iter = CIndex(eventnn-1); iter<= expeta.size()-1; ++iter)
            {
                risketa(eventnn-1) += expeta[iter];
            }

            for (int iter = eventnn - 1; iter >= 1; --iter)
            {
                risketa(iter-1) +=risketa(iter);
            }

            // calculate factor1

            (factor1.segment(0,CIndex(1))).fill(1/(risketa(0)));
            (factor2.segment(0,CIndex(1))).fill(1/(risketa(0)*risketa(0)));
            for (int index = 1; index < eventnn-1; ++index)
            {
                (factor1.segment((CIndex(index)),CIndex(index+1)-CIndex(index))).fill(factor1(CIndex(index)-1)+ (1/(risketa(index))));
                (factor2.segment((CIndex(index)),CIndex(index+1)-CIndex(index))).fill(factor2(CIndex(index)-1)+ (1/(risketa(index)*risketa(index))));
            }
            (factor1.tail(nobs-CIndex(eventnn-1))).fill(factor1(CIndex(eventnn-1)-1)+ 1/(risketa(eventnn-1)));
            (factor2.tail(nobs-CIndex(eventnn-1))).fill(factor2(CIndex(eventnn-1)-1)+ 1/(risketa(eventnn-1)*risketa(eventnn-1)));

            //std::cout << "expeta:" << expeta << std::endl;
            //std::cout << "delta:" << delta << std::endl;
            //std::cout << "risketa:" << risketa << std::endl;
            //std::cout << "factor1: " << factor1 << std::endl;
            //std::cout << "factor2: " << factor2 << std::endl;



            // calculate Jacobian
            W = factor1.array()*expeta.array()-factor2.array()*expeta.array().square();
            /*      std::cout << "W:" << W << std::endl;
            std::cout << "datX1:" << datX.block(0,0,31,1) << std::endl;
            std::cout << "datX2:" << datX.block(0,50,31,1) << std::endl;
            std::cout << "datX3:" << datX.block(0,100,31,1) << std::endl;*/



            // calculate Gradient
            z = (datX * main_beta).array() + (delta.cast<double>().array()-factor1.array()*expeta.array())*(W.array().inverse());
            // make sure no weights are too small

            //prob.fill(0.5);

            for (int kk = 0; kk < W.size(); ++kk)
            {
                if (W(kk) < 1e-5)
                {
                    W(kk) = 1e-5;
                }
            }


            // compute X'WX
            XX = XtWX(datX, W);


            // compute X'Wz
            XY = datX.adjoint() * (z.array()*W.array()).matrix();

            // compute rho after X'WX is computed
            compute_rho();
            //std::cout << "rho: " << rho << std::endl;

            // reset LDLT solver with new XX
            rho_changed_action();

            if (i > 0)
            {
                // reset values that need to be reset
                // for ADMM
                // and keep lambda the same
                init_warm(lambda, lambdaF);
            }

            for(j = 0; j < maxit; ++j)
            {
                old_gamma = aux_gamma;
                old_eta   = aux_eta;
                // old_nu = dual_nu;
                std::copy(dual_nu.data(), dual_nu.data() + dim_dual, old_nu.data());

                //std::cout << "beta_pre: " << beta_prev.sum() << std::endl;
                //std::cout << "main_beta: " << main_beta.sum() << std::endl;

                update_beta();
                update_eta();
                update_gamma();
                update_nu();
                //std::cout << "eps_primal: " << eps_primal << std::endl;
                //std::cout << "eps_dual: " << eps_dual << std::endl;
                //std::cout << "adj_nu: " << adj_nu << std::endl;
                //std::cout << "adj_gamma: " << adj_gamma << std::endl;
                // print_row(i);

                if(converged())
                    break;

                double old_c = adj_c;
                adj_c = compute_resid_combined();

                if(adj_c < 0.999 * old_c)
                {
                    double old_a = adj_a;
                    adj_a = 0.5 + 0.5 * std::sqrt(1 + 4.0 * old_a * old_a);
                    double ratio = (old_a - 1.0) / adj_a;
                    adj_gamma        = (1 + ratio) * aux_gamma - ratio * old_gamma;
                    adj_eta          = (1 + ratio) * aux_eta   - ratio * old_eta;
                    adj_nu.noalias() = (1 + ratio) * dual_nu   - ratio * old_nu;
                } else {
                    adj_a = 1.0;
                    adj_gamma = old_gamma;
                    adj_eta   = old_eta;
                    // adj_nu = old_nu;
                    std::copy(old_nu.data(), old_nu.data() + dim_dual, adj_nu.data());
                    adj_c = old_c / 0.999;
                }
                // only update rho after a few iterations and after every 40 iterations.
                // too many updates makes it slow.
                if(i > 5 && i % 2500 == 0)
                    update_rho();
            }

            //VectorXd dx = beta_prev - main_beta;
            //if (std::abs(XY.adjoint() * dx) < newton_tol)
            if (stopRule(main_beta, beta_prev, newton_tol))
            {
                //std::cout << "iters:\n" << i+1 << std::endl;
                break;
            }

        }
        // print_footer();

        return i + 1;
    }

    virtual void update_adaptive_group_weights(VectorXd &weights_)
    {
        group_weights = weights_;
    }

    virtual void update_adaptive_fused_weights(VectorXd &weights_)
    {
        fused_weights = weights_;
    }

    // this function returns gamma
    virtual VectorXd get_aux_gamma()
    {
        return aux_gamma;
    }

    // this function returns XTX or XTWX
    virtual MatrixXd get_hessian()
    {
        return XX;
    }

    virtual VectorXd get_xty()
    {
        return XY;
    }

};


#endif // ADMMOGLASSOCOXPHFUSEDTALL_H
