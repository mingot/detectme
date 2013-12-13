//
//  UIImage+iHOG.m
//  DetectMe
//
//  Created by a on 03/12/13.
//  Copyright (c) 2013 Josep Marc Mingot Hidalgo. All rights reserved.
//

#import "UIImage+iHOG.h"
#include <stdlib.h>
#import "UIImage+Resize.h"
#import "UIImage+Rotation.h"
#import "FilterConstant.h"
#import "dgray.h"
#import "lasso.hpp"




@implementation UIImage (iHOG)

- (int) getPos: (int) x ypos: (int) y zpos: (int) z xsize: (int) sx ysize: (int) sy zsize: (int) sz

{
    return x*sy*sz+y*sz+z;
}



- (UIImage *) convertToIHOG
{
    
    NSDate * start = [NSDate date];
    //Constants
    int n=1000000;
    int k=1024;
    int iters=1000;
    int lambda=0.02;
    int nx=5;
    int ny=5;
    int sbin=8;
 
    SolverOfLassoProblem * sol;
    sol = new SolverOfLassoProblem();

    //Compute HOG features
    HogFeature * hog = [self obtainHogFeatures];
    //[self printMatlabFeatures: hog];
    int par = 5;
    int fx = hog.numBlocksX;
    int fy = hog.numBlocksY;
    int fz = hog.numFeaturesPerBlock+1;
    double *featuresNoPad = hog.features;
    double *features = (double*) calloc((fy+2*par)*(fx+2*par)*fz,sizeof(double));
    int c = 0;
    int ihog=0,jhog=0;
    double *windows = (double*) calloc(ny*nx*fz*(fy+2*par-ny+1)*(fx+2*par-nx+1),sizeof(double));
    double *hogMatrix = (double*) calloc(ny*nx*fz,sizeof(double));
    double mean;
    double * hogVector = (double*) calloc(nx*ny*fz,sizeof(double));
    double std;
    



    for(int k=0; k<(fz-1); k++){
        for(int j=0; j<fx; j++){
            for(int i=0; i<fy; i++){
                features[ [self getPos:i+par ypos:j+par zpos:k xsize: (fy+2*par) ysize: (fx+2*par) zsize: fz ] ] = featuresNoPad[k*fx*fy+fy*j+i];
            }
        }
    }
    


    fx=fx+2*nx;
    fy=fy+2*ny;

    
    for(int i=0; i<(fy-ny+1); i++){
        for(int j=0; j<(fx-nx+1); j++){
            

            //Compute hog vector
            ihog=0;
            jhog=0;
            
            for(int isub=i; isub<(i+ny); isub++){
                for(int jsub=j; jsub<(j+nx); jsub++){
                    for(int ksub = 0; ksub < fz; ksub++){
                        hogMatrix[[self getPos:ihog ypos: jhog zpos:ksub xsize:ny ysize:nx zsize:fz]]=features[[self getPos:isub ypos:jsub zpos:ksub xsize:fy ysize:fx zsize:fz ]];
                    }
                    
                    jhog++;
                }
                
                jhog=0;
                ihog++;
            }
            
            for(int subi=0; subi<fz; subi++){
                for(int subj=0; subj<nx; subj++){
                    for(int subk=0; subk<ny; subk++)
                        hogVector[subi*nx*ny+subj*nx+subk] = hogMatrix[[self getPos:subk ypos:subj zpos:subi xsize:ny ysize:nx zsize:fz]];
                }
            }
            
            mean=0;
            
            for(int r=0; r<ny*nx*fz ;r++) mean=mean+hogVector[r];
            
            mean=mean/(ny*nx*fz);
            
            
            std=0;
            for(int r=0; r<ny*nx*fz ;r++){
                hogVector[r]=hogVector[r]-mean;
                std=std+hogVector[r]*hogVector[r];
            }
            std=sqrt(std);
            
            
            for(int r=0; r<ny*nx*fz ;r++){
                hogVector[r]=hogVector[r]/(std+1E-5);
            }
            
            for(int r=0; r<ny*nx*fz ;r++) windows[[self getPos:0 ypos:r zpos:c xsize:0 ysize:(ny*nx*fz) zsize:(fy-ny+1)*(fx-nx+1)]] = hogVector[r];
            
            c=c+1;
        }
    }
    
    printf("First Part:%f \n",[start timeIntervalSinceNow]);
   double * a = (double*) calloc(1024*(fy-ny+1)*(fx-nx+1),sizeof(double));
    //Solving LASSO Problem

    sol->lassoSolver<double>(ny*nx*fz,(fy-ny+1)*(fx-nx+1),(double *) windows, (double *) a);

    printf("Second Part:%f \n",[start timeIntervalSinceNow]);

    
    
    //Reconstruction matrix
    double * recons = (double*) calloc(3136*(fy-ny+1)*(fx-nx+1),sizeof(double));
    //dgray  -------------> A --------> 3136x1024
    //a      -------------> B --------> 1024x3082
    //recons -------------> C --------> 3136x3082
    
    //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3136, 3082, 1024, 1.0, (double *) dgray, 1024, a, 3082, 0.0, recons, 3082);
    
    //Check of the Solver Output Matrix
    
    
    // (1,2,3,Rows A, Cols B, Cols A, 1, dgray, ,a,colsb, 0,
    double * dgrayT = (double*) calloc(3136*1024,sizeof(double));
    for(int i=0;i <3136; i++) for(int j=0; j<1024; j++) dgrayT[1024*i+j]=dgray[j][i];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3136, (fy-ny+1)*(fx-nx+1), 1024, 1.0, dgrayT, 1024, a, (fy-ny+1)*(fx-nx+1), 0.0, recons, (fy-ny+1)*(fx-nx+1));
    
 
    //Reconstruct
    int im_y =(fy+2)*sbin;
    int im_x =(fx+2)*sbin;
    double * im = (double*) calloc(im_y*im_x,sizeof(double));
    int w_y =(fy+2)*sbin;
    int w_x =(fx+2)*sbin;
    double * weights = (double*) calloc(w_y*w_x,sizeof(double));
    c=0;

    
    
    for(int i=0; i<(fy-ny+1); i++){
        for(int j=0; j<(fx-nx+1); j++){
            
            double * patch = (double*) calloc((ny+2)*sbin*(nx+2)*sbin,sizeof(double));
            
            for(int subj=0 ; subj<(nx+2)*sbin; subj++){
                for(int subi=0 ; subi<(ny+2)*sbin; subi++){
                    patch[[ self getPos:0 ypos:subi zpos:subj xsize:0 ysize:(ny+2)*sbin zsize:(nx+2)*sbin]]= recons[[ self getPos:0 ypos:subj*(ny+2)*sbin+subi zpos:c xsize:0 ysize:3136 zsize:(fy-ny+1)*(fx-nx+1)]];
                }
            }
            
            
            for(int subj=0 ; subj<(nx+2)*sbin; subj++){
                for(int subi=0 ; subi<(ny+2)*sbin; subi++){
                    patch[[ self getPos:0 ypos:subi zpos:subj xsize:0 ysize:(ny+2)*sbin zsize:(nx+2)*sbin]] = patch[[ self getPos:0 ypos:subi zpos:subj xsize:0 ysize:(ny+2)*sbin zsize:(nx+2)*sbin]] * fil[subi][subj];
                }
            }
            
            
            int i_patch=0;
            int j_patch=0;
            for(int subj=j*sbin;subj<(j*sbin+(nx+2)*sbin);subj++,j_patch++){
                for(int subi=i*sbin;subi<(i*sbin+(ny+2)*sbin);subi++,i_patch++){
                    im[[ self getPos:0 ypos:subi zpos:subj xsize:0 ysize:im_y zsize:im_x]] = im[[ self getPos:0 ypos:subi zpos:subj xsize:0 ysize:im_y zsize:im_x]]+patch[[ self getPos:0 ypos:i_patch zpos:j_patch xsize:0 ysize:(ny+2)*sbin zsize:(nx+2)*sbin]];
                    
                    weights[[ self getPos:0 ypos:subi zpos:subj xsize:0 ysize:w_y zsize:w_x]]++;
                }
                i_patch=0;
            }
            
            
            c++;
            
            
        }
    }
    

    
    
    double * imageFinal = (double*) calloc((fy+2)*sbin*(fx+2)*sbin,sizeof(double));
    //Post Processing and clipping
    
    //Computing minimum
    for(int i=0; i<(fy+2)*sbin; i++) for(int j=0; j<(fx+2)*sbin; j++)     imageFinal[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:(fy+2)*sbin zsize:(fx+2)*sbin]] =im[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:im_y zsize:im_x]]/weights[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:w_y zsize:w_x]];
    double min=100000;
    for(int i=0; i<(fy+2)*sbin; i++) for(int j=0; j<(fx+2)*sbin; j++) if(imageFinal[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:(fy+2)*sbin zsize:(fx+2)*sbin]]<min) min=imageFinal[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:(fy+2)*sbin zsize:(fx+2)*sbin]];
    for(int i=0; i<(fy+2)*sbin; i++) for(int j=0; j<(fx+2)*sbin; j++)     imageFinal[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:(fy+2)*sbin zsize:(fx+2)*sbin]] =imageFinal[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:(fy+2)*sbin zsize:(fx+2)*sbin]]-min;
    
    //Computing maximum
    double max=0;
    for(int i=0; i<(fy+2)*sbin; i++) for(int j=0; j<(fx+2)*sbin; j++) if(imageFinal[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:(fy+2)*sbin zsize:(fx+2)*sbin]]>max) max=imageFinal[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:(fy+2)*sbin zsize:(fx+2)*sbin]];
    for(int i=0; i<(fy+2)*sbin; i++) for(int j=0; j<(fx+2)*sbin; j++)     imageFinal[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:(fy+2)*sbin zsize:(fx+2)*sbin]]=imageFinal[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:(fy+2)*sbin zsize:(fx+2)*sbin]]/max;
    
    
    
    double * imageFinalRescale = (double *) calloc(((fy+2)*sbin-2*par*sbin)*((fx+2)*sbin-2*par*sbin),sizeof(double));
    int i_re=0;
    int j_re=0;
    for(int i=(par*sbin-1); i<(fy+2)*sbin-par*sbin-1; i++) {
        for(int j=(par*sbin-1); j<(fx+2)*sbin-par*sbin-1; j++){
            imageFinalRescale[[ self getPos:0 ypos:i_re zpos:j_re xsize:0 ysize:((fy+2)*sbin-2*par*sbin) zsize:((fx+2)*sbin-2*par*sbin)]]=imageFinal[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:(fy+2)*sbin zsize:(fx+2)*sbin]];
            j_re++;
        }
        j_re=0;
        i_re++;
    }
    
    //Load Image
    UInt8 *imageBuffer = (UInt8 *) calloc(((fy+2)*sbin-2*par*sbin)*((fx+2)*sbin-2*par*sbin)*4,sizeof(UInt8));
    for(int i=0; i<((fy+2)*sbin-2*par*sbin); i++)
        for(int j=0; j<((fx+2)*sbin-2*par*sbin); j++){
            imageBuffer[i*((fx+2)*sbin-2*par*sbin)*4+4*j]  =  round(imageFinalRescale[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:((fy+2)*sbin-2*par*sbin) zsize:((fx+2)*sbin-2*par*sbin)]]*255);
            imageBuffer[i*((fx+2)*sbin-2*par*sbin)*4+4*j+1]=  round(imageFinalRescale[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:((fy+2)*sbin-2*par*sbin) zsize:((fx+2)*sbin-2*par*sbin)]]*255);
            imageBuffer[i*((fx+2)*sbin-2*par*sbin)*4+4*j+2]=  round(imageFinalRescale[[ self getPos:0 ypos:i zpos:j xsize:0 ysize:((fy+2)*sbin-2*par*sbin) zsize:((fx+2)*sbin-2*par*sbin)]]*255);
            imageBuffer[i*((fx+2)*sbin-2*par*sbin)*4+4*j+3]=  255;

    }
    CGContextRef context = CGBitmapContextCreate(imageBuffer, //data
                                                 ((fx+2)*sbin-2*par*sbin), //width
                                                 ((fy+2)*sbin-2*par*sbin), //height
                                                 8, //bits per component
                                                 ((fx+2)*sbin-2*par*sbin)*4, //bytes per row
                                                 CGColorSpaceCreateDeviceRGB(),
                                                 kCGImageAlphaPremultipliedLast ); //bitmap info
    
    CGImageRef ima = CGBitmapContextCreateImage(context);
    CGContextRelease(context);
    UIImage *image = [UIImage imageWithCGImage:ima scale:1.0 orientation:UIImageOrientationUp];
    free(imageBuffer);
    CGImageRelease(ima);
    printf("Third Part:%f \n",[start timeIntervalSinceNow]);

    return image;
    
}



@end
