//
//  lasso.h
//  DetectMe
//
//  Created by a on 05/12/13.
//  Copyright (c) 2013 Josep Marc Mingot Hidalgo. All rights reserved.
//

#include <iostream>
#ifndef _lasso_h
#include "Dhog.h"
#endif
#include "decomp.h"
class SolverOfLassoProblem{
    public:
    template <typename T>
    void lassoSolver(int n, int M, double * input, double * output){
        
        //Arg 1: Input Matrix prhs[0]
        //Arg 2: dhog prhs[1]
        //Arg 3: struct with three parameters (lambda, mode, pos). For us: [lambda, 2,True] prhs[2]
        //
        //
        
        //Parsing matrix dimension
        
        int dhog_x=800;
        int dhog_y=1024;
        
        int K=dhog_y;
        int nD=dhog_x;
        
        
        T lambda = 0.02;
        T lambda2 = 0;
        int L = K;
        int numThreads = -1;
        bool pos = true;
        bool verbose=true;
        bool ols = false;
        bool cholesky = true;
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
        
        Matrix<T> Xtrans(input,M,n);
        Matrix <T> X;
        Xtrans.transpose(X);
        
        Matrix<T> Dtrans((double *) dhog,K,n);
        Matrix <T> D;
        Dtrans.transpose(D);
        
        SpMatrix<T> alpha;
        Matrix<T> outputMatrix;
        
        if (cholesky) lasso<T>(X,D,alpha,L,lambda,lambda2,mode,pos,ols,numThreads,NULL);
        
        else lasso2<T>(X,D,alpha,L,lambda,lambda2,mode,pos,numThreads,NULL);
        
        alpha.toFull(outputMatrix);
        for(int i=0;i<1024;i++){
            for(int j=0; j<M; j++){
                output[M*i+j]=outputMatrix(i,j);
            }
        }
        
    }
};