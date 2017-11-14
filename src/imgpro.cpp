
// Computer Vision for Digital Post-Production
// Lecturer: Gergely Vass - vassg@vassg.hu
//
// Skeleton Code for programming assigments
// 
// Code originally from Thomas Funkhouser
// main.c
// original by Wagner Correa, 1999
// modified by Robert Osada, 2000
// modified by Renato Werneck, 2003
// modified by Jason Lawrence, 2004
// modified by Jason Lawrence, 2005
// modified by Forrester Cole, 2006
// modified by Tom Funkhouser, 2007
// modified by Chris DeCoro, 2007
//



// Include files
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <assert.h>
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"



// Program arguments
static char options[] =
"  -help\n"
"  -svdTest\n"
"  -sobelX\n"
"  -sobelY\n"
"  -log\n"
"  -harris <real:sigma>\n"
"  -saturation <real:factor>\n"
"  -brightness <real:factor>\n"
"  -blur <real:sigma>\n"
"  -sharpen \n"
"  -matchTranslation <file:other_image>\n"
"  -matchHomography <file:other_image>\n"
"  -trackMarkers <file:referenceImage> <int:numFrames> <int:numMarkers>\n";


static void 
ShowUsage(void)
{
  // Print usage message and exit
  fprintf(stderr, "Usage: imgpro input_image output_image [  -option [arg ...] ...]\n");
  fprintf(stderr, options);
  exit(EXIT_FAILURE);
}



static void 
CheckOption(char *option, int argc, int minargc)
{
  // Check if there are enough remaining arguments for option
  if (argc < minargc)  {
    fprintf(stderr, "Too few arguments for %s\n", option);
    ShowUsage();
    exit(-1);
  }
}



static int 
ReadCorrespondences(char *filename, R2Segment *&source_segments, R2Segment *&target_segments, int& nsegments)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open correspondences file %s\n", filename);
    exit(-1);
  }

  // Read number of segments
  if (fscanf(fp, "%d", &nsegments) != 1) {
    fprintf(stderr, "Unable to read correspondences file %s\n", filename);
    exit(-1);
  }

  // Allocate arrays for segments
  source_segments = new R2Segment [ nsegments ];
  target_segments = new R2Segment [ nsegments ];
  if (!source_segments || !target_segments) {
    fprintf(stderr, "Unable to allocate correspondence segments for %s\n", filename);
    exit(-1);
  }

  // Read segments
  for (int i = 0; i <  nsegments; i++) {

    // Read source segment
    double sx1, sy1, sx2, sy2;
    if (fscanf(fp, "%lf%lf%lf%lf", &sx1, &sy1, &sx2, &sy2) != 4) { 
      fprintf(stderr, "Error reading correspondence %d out of %d\n", i, nsegments);
      exit(-1);
    }

    // Read target segment
    double tx1, ty1, tx2, ty2;
    if (fscanf(fp, "%lf%lf%lf%lf", &tx1, &ty1, &tx2, &ty2) != 4) { 
      fprintf(stderr, "Error reading correspondence %d out of %d\n", i, nsegments);
      exit(-1);
    }

    // Add segments to list
    source_segments[i] = R2Segment(sx1, sy1, sx2, sy2);
    target_segments[i] = R2Segment(tx1, ty1, tx2, ty2);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int 
main(int argc, char **argv)
{
  // Look for help
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-help")) {
      ShowUsage();
    }
	if (!strcmp(argv[i], "-svdTest")) {
      R2Image *image = new R2Image();
	  //image->svdTest();
	  return 0;
    }
	if (!strcmp(argv[i], "-calcHomography")) {
		R2Image *image = new R2Image();
		
		FILE *f = fopen(argv[2], "r");
		
		if (f == NULL) {
			printf("I didn't read this right... \n");
		}

		int input[3], output[3];
		for (int i = 0; i < 8; i += 2) {
			fscanf(f, "%d %d %d %d ", &input[i], &input[i + 1], &output[i], &output[i + 1]);
		}
		image->calcHomographyMatrix(input, output, 8);
		exit(-1);
	}
  }

  // Read input and output image filenames
  if (argc < 3)  ShowUsage();
  argv++, argc--; // First argument is program name

  char *input_image_name = *argv; argv++, argc--; 
  char *output_image_name = *argv; argv++, argc--; 

  // Allocate image
  R2Image *image = new R2Image();

  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    exit(-1);
  }


  // Read input image
  if (!image->Read(input_image_name)) {
    fprintf(stderr, "Unable to read image from %s\n", input_image_name);
    exit(-1);
  }

  // Initialize sampling method
  int sampling_method = R2_IMAGE_POINT_SAMPLING;

  // Parse arguments and perform operations 
  while (argc > 0) {
	  if (!strcmp(*argv, "-brightness")) {
		  CheckOption(*argv, argc, 2);
		  double factor = atof(argv[1]);
		  argv += 2, argc -= 2;
		  image->Brighten(factor);
	  }
	  else if (!strcmp(*argv, "-sobelX")) {
		  argv++, argc--;
		  image->SobelX();
	  }
	  else if (!strcmp(*argv, "-sobelY")) {
		  argv++, argc--;
		  image->SobelY();
	  }
	  else if (!strcmp(*argv, "-log")) {
		  argv++, argc--;
		  image->LoG();
	  }
	  else if (!strcmp(*argv, "-saturation")) {
		  CheckOption(*argv, argc, 2);
		  double factor = atof(argv[1]);
		  argv += 2, argc -= 2;
		  image->ChangeSaturation(factor);
	  }
	  else if (!strcmp(*argv, "-harris")) {
		  CheckOption(*argv, argc, 2);
		  double sigma = atof(argv[1]);
		  argv += 2, argc -= 2;
		  //image->Harris(sigma);
	  }
    else if (!strcmp(*argv, "-blur")) {
      CheckOption(*argv, argc, 2);
      double sigma = atof(argv[1]);
      argv += 2, argc -= 2;
     image->Blur(sigma);
    }
    else if (!strcmp(*argv, "-sharpen")) {
      argv++, argc--;
      image->Sharpen();
    }
	else if (!strcmp(*argv, "-trackMarkers")) {

		//read in the reference image
		CheckOption(*argv, argc, 2);
		R2Image *referenceImageMarker = new R2Image(argv[1]);
		//read in the reference image for the known points
		CheckOption(*argv, argc, 3);
		R2Image *referenceImageKnownPoints = new R2Image(argv[2]);
		//read in the number of frames to track through
		CheckOption(*argv, argc, 4);
		const int numFrames = atoi(argv[3]);
		//read in the number of markers we are working with
		CheckOption(*argv, argc, 5);
		const int numMarkers = atoi(argv[4]);
		//read in how many known locations we have
		CheckOption(*argv, argc, 6);
		const int numKnown = atoi(argv[5]);
		argv += 6, argc -= 6;

		//read in the locations of the reference points
		std::vector<int> knownX;
		std::vector<int> knownY;
		std::vector<int> knownZ;

		for (int i = 0; i < numKnown; i++) {
			printf("Please enter the x y and z measurements for known location %d \n", i);
			scanf("%d %d %d", knownX[i], knownY[i], knownZ[i]);
		}

		image->SetGivenLocs3D(knownX, knownY, knownZ);

		printf("Thank you! Your locations are : \n");
		for (int i = 0; i < numKnown; i++) {
			printf("< %d, %d, %d > \n", knownX[i], knownY[i], knownZ[i]);
		}

		//check to make sure everything is fine
		printf("\nlNUMBER OF FRAMES: %d\n", numFrames);
		printf("NUMBER OF MARKERS: %d\n", numMarkers);
		printf("NUMBER OF KNOWN POINTS: %d\n", numKnown);
		printf("input image name: %s\n", input_image_name);
		printf("output image name: %s\n", output_image_name);

		// extract input and output filepaths
		std::string sInput = input_image_name; //this should be the first frame
		std::string sOutput = output_image_name;

		//find the last instance of "." which should be directly before the file type
		//index is the index of this value
		int index = sInput.find_last_of(".");
		if (index == -1) {
			fprintf(stderr, "Unable to find extension in %s\n", input_image_name);
			exit(-1);
		}
		//splits string between the name and the file type/extension
		std::string extension = sInput.substr(index);

		//search for the frame number (specifically using the first frame
		index = sInput.find("0000001");
		if (index == -1) {
			fprintf(stderr, "Unable to find '0000001' (7-digit padding) in %s\n", input_image_name);
			exit(-1);
		}

		//get the first part of the string until the number, should be the original name thats the same in each input frame
		std::string inputPath = sInput.substr(0, index);

		//search for the frame number in the output file name
		index = sOutput.find("0000001");
		if (index == -1) {
			fprintf(stderr, "Unable to find '0000001' (7-digit padding) in %s\n", output_image_name);
			exit(-1);
		}

		//get the first part of the string till the number, should be original name thats the same in each output frame
		std::string outputPath = sOutput.substr(0, index);

		std::string number; // padding = 7 digits
		const int height = image->Height();
		const int sqRadius = 5;
		int xa, xb, ya, yb;

		// image = first frame
		//TEST
		printf("got here1 \n");
		//TEST
		
		//find the reference points on the screen
		Marker2D* knownPoints2D = image->MarkerDetection(referenceImageKnownPoints);

		//find the markers on the screen
		Marker2D* MarkersA = image->MarkerDetection(referenceImageMarker);

		//set them as the values for the image fields
		std::vector<int> x(numMarkers);
		std::vector<int> y(numMarkers);

		std::vector<int> newX(numKnown);
		std::vector<int> newY(numKnown);

		for (int i = 0; i < numMarkers; i++) {
			printf("found the marker at ( %d, %d )", MarkersA[i].xVal, MarkersA[i].yVal);
			x[i] = MarkersA[i].xVal;
			y[i] = MarkersA[i].yVal;
		}
		for (int i = 0; i < numKnown; i++) {
			printf("found the known spot at ( %d, %d )", knownPoints2D[i].xVal, knownPoints2D[i].yVal);
			newX[i] = knownPoints2D[i].xVal;
			newY[i] = knownPoints2D[i].yVal;
		}
		image->SetMarkerLocs(x, y);
		image->SetGivenLocs2D(newX, newY);

		printf("Found %d markers and %d reference spots in first frame\n", numMarkers, numKnown);

		//draw markers and reference spots
		for (int i = 0; i < numMarkers; i++) {
			image->drawMarkers(i, 1, 0, 1);   
		}
		for (int i = 0; i < numKnown; i++) {
			image->drawReferenceSpots(i, 0, 1, 1);
		}



		// R2Image *imageB = new R2Image(*image);
		// R2Image *tempImage;
		// std::vector<int> markersB;
		// std::vector<double> Hvector;

		// Calculate H using DLTRANSAC between frame(1) and frame(2). reject bad tracks
		// double H[3][3];
		// image->SkyDLTRANSAC(imageB, H);
		// TODO Hvector never set..
		// for (int i = 0; i < 3; i++)
		//   for (int j = 0; j < 3; j++)
		//     Hvector.push_back(H[i][j]);

		// imageB->SetH(Hvector);

		// Translation RANSAC
		// image->RANSAC(imageB);

		// warp and blend sky in frame(1)
		// R2Image *outputOrigImage = new R2Image(*image);
		// outputOrigImage->WarpSkyTranslation(skyImage);

		// Write output image
		if (image->Write(output_image_name)) {
			fprintf(stderr, "Unable to read image from %s\n", output_image_name);
			exit(-1);
		}

		printf("Finished frame 1\n");

		//R2Image *imageA;

		//// iterate through frames
		//// imageA = frame(i-1)
		//// imageB = frame(i)

		//for (int i = 2; i <= numFrames; i++) {
		//	// SETUP
		//	// copy imageB into imageA via operator=
		//	imageA = imageB;

		//	// Calculate new imageB
		//	number = "0000000" + std::to_string(i);
		//	number = number.substr(number.length() - 7);

		//	imageB = new R2Image((inputPath + number + extension).c_str());

		//	// Track features from frame(i-1) to frame(i)
		//	featuresB.clear();
		//	featuresB = imageA->findAFeaturesOnB(imageB, imageA->SkyFeatures(), sqRadius);
		//	imageB->SetSkyFeatures(featuresB);
		//	// Hvector.clear();

		//	// Calculate H between frame(i-1) and frame(i)
		//	// imageA->SkyDLTRANSAC(imageB, H);

		//	// for (int j = 0; j < 3; j++) {
		//	//   for (int k = 0; k < 3; k++) {
		//	//     Hvector.push_back(H[j][k]);
		//	//     printf("%f ", H[j][k]);
		//	//   }
		//	//   printf("\n");
		//	// }

		//	// imageB->SetH(Hvector);

		//	imageA->SkyRANSAC(imageB);

		//	tempImage = new R2Image(*imageB);
		//	tempImage->WarpSkyTranslation(skyImage);

		//	if (!tempImage->Write((outputPath + number + extension).c_str())) {
		//		fprintf(stderr, "Unable to read image from %s\n", (outputPath + number + extension).c_str());
		//		exit(-1);
		//	}

		//	printf("Tracked features from frame%d to frame%d\n", i - 1, i);
		//	delete imageA;
		//	delete tempImage;
		//}
		delete image;
		//delete imageB;
		//image = outputOrigImage;


	}
	
    else {
      // Unrecognized program argument
      fprintf(stderr, "image: invalid option: %s\n", *argv);
      ShowUsage();
    }
  }

  // Write output image
  if (!image->Write(output_image_name)) {
    fprintf(stderr, "Unable to read image from %s\n", output_image_name);
    exit(-1);
  }

  // Delete image
  delete image;

  // Return success
  return EXIT_SUCCESS;
}



