//
//  UIImage+iHOG.h
//  DetectMe
//
//  Created by a on 03/12/13.
//  Copyright (c) 2013 Josep Marc Mingot Hidalgo. All rights reserved.
//

#import <UIKit/UIKit.h>
#import "UIImage+HOG.h"
@interface UIImage (iHOG)

- (int) getPos: (int) x ypos: (int) y zpos: (int) z xsize: (int) sx ysize: (int) sy zsize: (int) sz;
- (UIImage *) convertToIHOG ;

@end
