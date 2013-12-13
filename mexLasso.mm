
/* Software SPAMS v2.1 - Copyright 2009-2011 Julien Mairal
 *
 * This file is part of SPAMS.
 *
 * SPAMS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SPAMS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SPAMS.  If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file
 *                toolbox decomp
 *
 *                by Julien Mairal
 *                julien.mairal@inria.fr
 *
 *                File mexLasso.h
 * \brief mex-file, function mexLasso
 * Usage: alpha = mexLasso(X,D,param);
 * Usage: alpha = mexLasso(X,G,DtR,param);
 * */

#include "Dhog.h"
#include "decomp.h"

template <typename T>
void lassoSolver(int n, int M, double * input, double * output) {
    
    //Arg 1: Input Matrix prhs[0]
    //Arg 2: dhog prhs[1]
    //Arg 3: struct with three parameters (lambda, mode, pos). For us: [lambda, 2,True] prhs[2]
    //
    //
    
    //Parsing matrix dimension
    
    int K=dhog_y;
    int nD=dhog_x;
    
    
    T lambda = 0.02;
    T lambda2 = 0;
    int L = K;
    int numThreads = -1;
    bool pos = true;
    bool verbose=true;
    bool ols = false;
    bool cholesky = false;
    constraint_type mode;
    mode= PENALTY;
    
    if (L > n && !(mode == PENALTY && isZero(lambda) && !pos && lambda2 > 0)) {
        if (verbose)
            printf("L is changed to %d\n",n);
        L=n;
    }
    if (L > K) {
        if (verbose)
            printf("L is changed to %d\n",K);
        L=K;
    }
    
    Matrix<T> X(input,n,M);
    Matrix<T> D(dhog,n,K);
    SpMatrix<T> alpha;
    
    if (cholesky) lasso<T>(X,D,alpha,L,lambda,lambda2,mode,pos,ols,numThreads,NULL);
    
    else lasso2<T>(X,D,alpha,L,lambda,lambda2,mode,pos,numThreads,NULL);
    
    
}


