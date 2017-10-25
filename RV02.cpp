/*
RV02: Median und Sobel

Autor: .....................
HAW-University of Applied Sciences - Hamburg,Germany

*/ 

#include "ltiObject.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <list>
#include <stdio.h>

#include "RV02.h"
#include "ltiTimer.h"
#include "ltiBMPFunctor.h"
#include "ltiViewer.h"
#include "ltiSplitImg.h"
#include "gtk.h"
#include "ltiGtkServer.h"


#define PI 3.14159265

using std::cout;
using std::endl;


namespace lti {

	void RV02::operator()(int argc,char *argv[]) {

		/*********************************************/
		/**** has always to be started (AMei)     ****/
		/**** if program is without gtk-Widgets   ****/
		/*********************************************/

		gtkServer server;
		server.start();

		/******************************************/
		/**** instantiation of used components ****/
		/******************************************/

		/*---------------------*/
		/* loaders and viewers */
		/*---------------------*/
		loadBMP loader;                         // object for loading .bmp-images

		viewer view("Original");                // object for visualizing images
		viewer viewTransformed("Filtered");
		viewer viewGradient("Gradient");
		viewer viewDirection("Direction");

		/*---------------------*/
		/* images & channels   */
		/*---------------------*/
		image img;                              // normalized (color) image

		channel8  src;  // source picture       // 8-bit-image (source)
		channel8  dst;  // destination picture  // 8-bit-image (source) 
		channel8  grad;  // gradient picture  // 8-bit-image (source)
		channel8  direc;  // direction picture  // 8-bit-image (source)


		/*-----------------------------*/
		/* Image processing components */
		/*-----------------------------*/

		// object to split image into hue, saturation and intensity
		// hue        = Farbton
		// saturation = Farbsättigung
		// intensity  = Intensität (Grauwert)
		splitImageToHSI splitter;



		/******************************************/
		/*    the program                         */
		/******************************************/

		// load the source image
		loader.load("circle.bmp",img);

		// extract the intensity channel only
		splitter.getIntensity(img,src);

		// determine image size
		const int rowSize    = src.rows();
		const int columnSize = src.columns();


		// set destination size to source size 
		dst.resize(rowSize,columnSize,0,false,true);

		grad.resize(rowSize,columnSize,0,false,true);
		direc.resize(rowSize,columnSize,0,false,true);



		Median(src,dst,9,9);
		Sobel(src, grad, direc);



		// view pictures
		view.show(src);
	//	viewTransformed.show(dst);
		viewDirection.show(direc);
		viewGradient.show(grad);

		getchar();

	}





	/***************************************************************************/
	/* Function definition: ----- Median-operator----                          */
	/***************************************************************************/
	void RV02::Median(  const	     channel8& sPic, 	// source picture 
		channel8& dPic, 	// destination picture
		const int MaskSizeX,		    // mask size in x-direction
		const int MaskSizeY		 	// mask size in y-direction
		)
	{
		int maskSizeX = MaskSizeX;
		int maskSizeY = MaskSizeY;

		if(maskSizeX%2==0){
			maskSizeX++;
		}
		if(maskSizeY%2==0){
			maskSizeY++;
		}


		const int PicSizeY = sPic.rows();
		const int PicSizeX = sPic.columns();

		int x,y,mx,my,centerMaskX, centerMaskY;


		for(y=0; y < PicSizeY - maskSizeY; y++){
			for(x=0; x<PicSizeX - maskSizeX; x++){
				// create hystogram
				int hisotgram[256];
				for (int k = 0; k<256; k++){
					hisotgram[k] = 0;
				}
				for(int yInMask=y; yInMask < y + maskSizeY; yInMask++){
					for(int xInMask = x; xInMask < x + maskSizeX; xInMask++){
						hisotgram[sPic[yInMask][xInMask]]++;
					}
				}

				int posTillMiddle = (maskSizeY * maskSizeX + 1) /2;

				int i = -1;
				while(posTillMiddle>0){
					i++;
					posTillMiddle = posTillMiddle - hisotgram[i];
				}

				// set new value
				centerMaskX = x + maskSizeX/2;
				centerMaskY = y +maskSizeY/2;
				dPic[centerMaskY][centerMaskX] = i;
			}


		}

	}


	void RV02::Sobel(const channel8& sPic, channel8& GradientPic, channel8& DirectionPic){
		

		int fMaskeX [9] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 }; 
		int fMaskeY [9] = { -1, -2, -1, 0, 0, 0, 1, 2, 1};

		const int PicSizeY = sPic.rows();
		const int PicSizeX = sPic.columns();
		const int maskSize =3;
		int x,y;

		for(y=0; y < PicSizeY - maskSize; y++){
			for(x=0; x<PicSizeX - maskSize; x++){

				float gx = 0;
				float gy = 0;

				for(int yInMask=0; yInMask < maskSize; yInMask++){
					for(int xInMask = 0; xInMask < maskSize; xInMask++){
						int posInMask = yInMask + xInMask * maskSize;
						gx += sPic[y+yInMask][x+xInMask] * fMaskeX[posInMask];
						gy += sPic[y+yInMask][x+xInMask] * fMaskeY[posInMask];
					}
				}

				gx = gx/4;
				gy = gy/4;
				
				int beitragGradient = (int)sqrt(gx*gx+gy*gy);
				GradientPic[y +maskSize/2][x + maskSize/2] = beitragGradient; 

				double richtungGrad = (atan2(gx,gy) * 180 / PI);
				
				if (richtungGrad < -22.5) {
					richtungGrad += 360;
				}

				if (richtungGrad >= -22.5 && richtungGrad <= 22.5)
					richtungGrad = 0;
				else if (richtungGrad >= 22.5 && richtungGrad <= 67.5)
					richtungGrad = 1;
				else if (richtungGrad >= 67.5 && richtungGrad <= 112.5)
					richtungGrad = 2;
				else if (richtungGrad >= 112.5 && richtungGrad <= 157.5)
					richtungGrad = 3;
				else if (richtungGrad >= 157.5 && richtungGrad <= 202.5)
					richtungGrad = 4;
				else if (richtungGrad >= 202.5 && richtungGrad <= 247.5)
					richtungGrad = 5;
				else if (richtungGrad >= 247.5 && richtungGrad <= 292.5)
					richtungGrad = 6;
				else if (richtungGrad >= 292.5 && richtungGrad <= 337.5)
					richtungGrad = 7;
				
				/*
				
				richtungGrad = (richtungGrad+ 22.5) / 45;
				*/
				DirectionPic[y +maskSize/2][x + maskSize/2] = (int) richtungGrad; 

				

				//berechnung des gradients

				//berehcnung der gradientrichtung

			}
		}
	}



};
