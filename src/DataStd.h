#ifndef DATASTD_H
#define DATASTD_H

#include <RcppEigen.h>
#include <Eigen/Core>

#ifdef __AVX__
#include "Linalg/AVX.h"
#endif
using Eigen::MatrixXd;

template <typename Double = double>
class DataStd
{
private:
    typedef Eigen::SparseVector<double> SparseVector;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Array <double, Eigen::Dynamic, 1> Array;
    typedef const Eigen::Ref<const Vector> ConstGenericVector;
    typedef Eigen::Ref<Array> ArrayRef;
    
    typedef Eigen::MappedSparseMatrix<double> MSpMat;
    typedef Eigen::SparseVector<float> SpVecf;
    typedef Eigen::SparseMatrix<float> SpMatf;
    typedef Eigen::SparseVector<double> SpVec;
    typedef Eigen::SparseMatrix<double> SpMat;
    

    // flag - 0: standardize = FALSE, intercept = FALSE
    //             directly fit model
    // flag - 1: standardize = TRUE, intercept = FALSE
    //             scale x and y by their standard deviation
    // flag - 2: standardize = FALSE, intercept = TRUE
    //             center x, standardize y
    // flag - 3: standardize = TRUE, intercept = TRUE
    //             standardize x and y
    const int flag;

    const int n;
    const int p;

    double meanY;
    double scaleY;
    Array  meanX;
    Array  scaleX;

    static double sd_n(ConstGenericVector &v)
    {
#ifdef __AVX__
        double s, ss;
        const int vsize = v.size();
        get_ss_avx<double>(v.data(), v.size(), s, ss);
        s /= vsize;
        return std::sqrt(ss / vsize - s * s);
#else
        double mean = v.mean();
        Vector v_centered = v.array() - mean;

        return v_centered.norm() / std::sqrt(double(v.size()));
#endif
    }

    // spvec -> spvec / arr, elementwise
    static void elementwise_quot(SparseVector &spvec, Array &arr)
    {
        for(typename SparseVector::InnerIterator iter(spvec); iter; ++iter)
        {
            iter.valueRef() /= arr[iter.index()];
        }
    }

    // inner product of spvec and arr
    static double sparse_inner_product(SparseVector &spvec, Array &arr)
    {
        double res = 0.0;
        for(typename SparseVector::InnerIterator iter(spvec); iter; ++iter)
        {
            res += iter.value() * arr[iter.index()];
        }
        return res;
    }

public:
    DataStd(int n_, int p_, bool standardize, bool intercept) :
        flag(int(standardize) + 2 * int(intercept)),
        n(n_),
        p(p_),
        meanY(0.0),
        scaleY(1.0)
    {
        if(flag == 3 || flag == 2)
            meanX.resize(p);
        if(flag == 3 || flag == 1)
            scaleX.resize(p);
    }

    void standardize(MatrixXd &X, Vector &Y)
    {
        double n_invsqrt = 1.0 / std::sqrt(Double(n));

        // standardize Y
        switch(flag)
        {
            case 1:
                scaleY = sd_n(Y);
                Y.array() /= scaleY;
                break;
            case 2:
            case 3:
                meanY = Y.mean();
                Y.array() -= meanY;
                scaleY = Y.norm() * n_invsqrt;
                Y.array() /= scaleY;
                break;
            default:
                break;
        }

        // standardize X
        switch(flag)
        {
            case 1:
                for(int i = 0; i < p; i++)
                {
                    scaleX[i] = sd_n(X.col(i));
                    X.col(i).array() *= (1.0 / scaleX[i]);
                }
                break;
            case 2:
                for(int i = 0; i < p; i++)
                {
                    meanX[i] = X.col(i).mean();
                    X.col(i).array() -= meanX[i];
                }
                break;
            case 3:
                for(int i = 0; i < p; i++)
                {
                    /*meanX[i] = X.col(i).mean();
                    X.col(i).array() -= meanX[i];
                    scaleX[i] = X.col(i).norm() * n_invsqrt;
                    X.col(i).array() /= scaleX[i];*/
    #ifdef __AVX__
                    double *begin = &X(0, i);
                    double s, ss;
                    get_ss_avx<double>(begin, n, s, ss);
                    meanX[i] = s / n;
                    scaleX[i] = std::sqrt(ss - s * s / n) * n_invsqrt;
                    standardize_vec_avx<double>(begin, n, meanX[i], 1.0 / scaleX[i]);
    #else
                    // double *begin = &X(0, i);
                    // double *end = begin + n;
                    // meanX[i] = X.col(i).mean();
                    // std::transform(begin, end, begin, std::bind2nd(std::minus<double>(), meanX[i]));
                    // scaleX[i] = X.col(i).norm() * n_invsqrt;
                    // std::transform(begin, end, begin, std::bind2nd(std::multiplies<double>(), 1.0 / scaleX[i]));
                    meanX[i] = X.col(i).mean();
                    X.col(i).array() -= meanX[i];
                    scaleX[i] = X.col(i).norm() * n_invsqrt;
                    X.col(i).array() /= scaleX[i];
    #endif
                }
                break;
            default:
                break;
        }
    }
    
    
    
    
    
    void standardize(const MSpMat &X, Vector &Y)
    {
        double n_invsqrt = 1.0 / std::sqrt(Double(n));
        
        // standardize Y
        switch(flag)
        {
        case 1:
            scaleY = sd_n(Y);
            Y.array() /= scaleY;
            break;
        // case 2:
        // case 3:
        //     meanY = Y.mean();
        //     Y.array() -= meanY;
        //     scaleY = Y.norm() * n_invsqrt;
        //     Y.array() /= scaleY;
        //     break;
        default:
            break;
        }
        
        // standardize X
//         switch(flag)
//         {
//         case 1:
//             for(int i = 0; i < p; i++)
//             {
//                 //scaleX[i] = sd_n(X.col(i));
//                 //X.col(i).array() *= (1.0 / scaleX[i]);
//                 
//                     
//             }
//             // loop over nonzero elements of X
//             for (int k=0; k<X.outerSize(); ++k)
//                 for (SparseMatrix<double>::InnerIterator it(X,k); it; ++it)
//                 {
//                     it.value() *= (1.0 / scaleX[it.col()]);
//                     //it.row();   // row index
//                     //it.col();   // col index (here it is equal to k)
//                     //it.index(); // inner index, here it is equal to it.row()
//                 }
//             break;
// //         case 2:
// //             for(int i = 0; i < p; i++)
// //             {
// //                 meanX[i] = X.col(i).mean();
// //                 X.col(i).array() -= meanX[i];
// //             }
// //             break;
// //         case 3:
// //             for(int i = 0; i < p; i++)
// //             {
// //                 /*meanX[i] = X.col(i).mean();
// //                  X.col(i).array() -= meanX[i];
// //                  scaleX[i] = X.col(i).norm() * n_invsqrt;
// //                  X.col(i).array() /= scaleX[i];*/
// // #ifdef __AVX__
// //                 double *begin = &X(0, i);
// //                 double s, ss;
// //                 get_ss_avx<double>(begin, n, s, ss);
// //                 meanX[i] = s / n;
// //                 scaleX[i] = std::sqrt(ss - s * s / n) * n_invsqrt;
// //                 standardize_vec_avx<double>(begin, n, meanX[i], 1.0 / scaleX[i]);
// // #else
// //                 double *begin = &X(0, i);
// //                 double *end = begin + n;
// //                 meanX[i] = X.col(i).mean();
// //                 std::transform(begin, end, begin, std::bind2nd(std::minus<double>(), meanX[i]));
// //                 scaleX[i] = X.col(i).norm() * n_invsqrt;
// //                 std::transform(begin, end, begin, std::bind2nd(std::multiplies<double>(), 1.0 / scaleX[i]));
// // #endif
// //             }
// //             break;
//         default:
//             break;
//         }
    }
    
    
    

    void recover(double &beta0, ArrayRef coef)
    {
        switch(flag)
        {
            case 0:
                beta0 = 0;
                break;
            case 1:
                beta0 = 0;
                coef /= scaleX;
                coef *= scaleY;
                break;
            case 2:
                coef *= scaleY;
                beta0 = meanY - (coef * meanX).sum();
                break;
            case 3:
                coef /= scaleX;
                coef *= scaleY;
                beta0 = meanY - (coef * meanX).sum();
                break;
            default:
                break;
        }
    }

    void recover(double &beta0, SparseVector &coef)
    {
        switch(flag)
        {
            case 0:
                beta0 = 0;
                break;
            case 1:
                beta0 = 0;
                elementwise_quot(coef, scaleX);
                coef *= scaleY;
                break;
            case 2:
                coef *= scaleY;
                beta0 = meanY - sparse_inner_product(coef, meanX);
                break;
            case 3:
                elementwise_quot(coef, scaleX);
                coef *= scaleY;
                beta0 = meanY - sparse_inner_product(coef, meanX);
                break;
            default:
                break;
        }
    }

    double get_scaleY() { return scaleY; }
};



#endif // DATASTD_H
