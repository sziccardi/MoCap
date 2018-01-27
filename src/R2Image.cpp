// Source file for image class
#define _USE_MATH_DEFINES

// Include files 
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <cmath>



////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////

//TODO: Are all of these necessary? Need to evaluate.

R2Image::
R2Image(void)
	: pixels(NULL),
	npixels(0),
	width(0),
	height(0),
	markerLocs2D(std::vector<Marker2D>()),
	givenLocs2D(std::vector<Marker2D>()),
	givenLocs3D(std::vector<Marker3D>())


{	
}


R2Image::
R2Image(const char *filename)
	: pixels(NULL),
	npixels(0),
	width(0),
	height(0),
	markerLocs2D(std::vector<Marker2D>()),
	givenLocs2D(std::vector<Marker2D>()),
	givenLocs3D(std::vector<Marker3D>())
{
	// Given the name
	// Read image
	Read(filename);
}



R2Image::
R2Image(const char *filename, int numMarkers, std::vector<Marker3D> known3D)
	: pixels(NULL),
	npixels(0),
	width(0),
	height(0),
	markerLocs2D(std::vector<Marker2D>()),
	givenLocs2D(std::vector<Marker2D>()),
	givenLocs3D(known3D),
	numMarkers(numMarkers)
{
  // Given the name
  // Read image
  Read(filename);
  
}



R2Image::
R2Image(int width, int height, std::vector<Marker3D> known3D)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height),
	markerLocs2D(std::vector<Marker2D>()),
	givenLocs2D(std::vector<Marker2D>()),
	givenLocs3D(known3D)
{
  // Given the size
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p, std::vector<Marker3D> known3D)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height),
	markerLocs2D(std::vector<Marker2D>()),
	givenLocs2D(std::vector<Marker2D>()),
	givenLocs3D(known3D)
{
  // Given the size and the pixels
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height),
	markerLocs2D(image.markerLocs2D),
	givenLocs2D(image.givenLocs2D),
	givenLocs3D(image.givenLocs3D),
	numMarkers(image.numMarkers)
    
{
  // Given an image, make a copy
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}


////////////////////////////////////////////////////////////////////////
// OPERATORS 
////////////////////////////////////////////////////////////////////////

R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Reset marker locations
  markerLocs2D = image.markerLocs2D;

  givenLocs2D = image.givenLocs2D;
  
  givenLocs3D = image.givenLocs3D;
  
  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}

R2Image& R2Image::
operator*(const R2Image image)
{
	// Reset width and height
	width = image.width;
	height = image.height;
	
	R2Image myImg = image;

	// multiply pixels
	for (int c = 0; c < width; c++) {
		for (int r = 0; r < height; r++) {
			
			Pixel(c,r) *= myImg.Pixel(c,r);
			
		}
	}
	
	// Return image
	return *this;
}


////////////////////////////////////////////////////////////////////////
// FUNCTIONS
////////////////////////////////////////////////////////////////////////

//Helpful tools
int* R2Image::matMult(int* vector, int size, double** mat, int rows, int columns) 
{

	int* result = new int[3]; 

	for (int i = 0; i < rows; i++) {
		result[i] = 0;
	}

	if (size != columns) {
		printf("Cannot multiply these. \n");
		return NULL;
	}
	else {		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				result[i] += (int) (vector[j] * mat[j][i]);
			}
		}
		return result; 
	}
}

int* R2Image::Normalize(int* vector, int size) 
{
	int* result = new int[3];
	for (int i = 0; i < size; i++) {
		result[i] = 0;
	}
	for (int i = 0; i < size; i++) {
		result[i] = vector[i] / vector[size - 1];
	}
	return result;
}

bool R2Image::
validPixel(const int x, const int y)
{
	return (x >= 0 && x < width && y >= 0 && y < height);
}

double R2Image::SSD(R2Image* otherImage, int x1, int y1, int x2, int y2, int range) 
{
	//Sum of Squared Difference
	//TODO: write description
	////////////////////////////////////////////////////////////////////////

	double ssd = 0;
	//iterate through a local area defined by range
	for (int i = -range; i < range; i++) {
		for (int j = -range; j < range; j++) {
			//if we dont go off the edge of the image
			if ( validPixel(x1+i, y1+j) && validPixel(x2+i, y1+j)
				&& validPixel(x1+i, y2+j) && validPixel(x2+i, y2+j) ) { 

				//take the difference between the luminances, square and add to sum
				//double diff = Pixel((i + x1), (j + y1)).Luminance() - otherImage->Pixel((i + x2), (j + y2)).Luminance();
				//ssd += diff * diff;

				//take the difference between the r, g, and b values then square and sum those
				double diffR = Pixel((i + x1), (j + y1)).Red() - otherImage->Pixel((i + x2), (j + y2)).Red();
				double diffG = Pixel((i + x1), (j + y1)).Green() - otherImage->Pixel((i + x2), (j + y2)).Green();
				double diffB = Pixel((i + x1), (j + y1)).Blue() - otherImage->Pixel((i + x2), (j + y2)).Blue();
				ssd += ((diffR * diffR) + (diffG * diffG) + (diffB * diffB));
				
			} 
		}
	}

return ssd;
}

//TODO: Combine the box-drawing methods
void R2Image::drawMarkers(int m, float r, float g, float b) 
{
	//draw a box for markers
	//TODO: write description
	////////////////////////////////////////////////////////////////////////

	int x = markerLocs2D[m].xVal;
	int y = markerLocs2D[m].yVal;

	int newx = givenLocs2D[m].xVal;
	int newy = givenLocs2D[m].yVal;

	int i = -20;
	int j = -20;
	while (j < 20) {
		if (i + x >= 0 && j + y >= 0 && i + x < width && j + y < height) {
			Pixel(i + x, j + y).SetBlue(b);
			Pixel(i + x, j + y).SetRed(r);
			Pixel(i + x, j + y).SetGreen(g);

			Pixel(i + newx, j + newy).SetBlue(r);
			Pixel(i + newx, j + newy).SetRed(g);
			Pixel(i + newx, j + newy).SetGreen(b);
		}
		j++;
	}
	while (i < 20) {
		if (i + x >= 0 && j + y >= 0 && i + x < width && j + y < height) {
			Pixel(i + x, j + y).SetBlue(b);
			Pixel(i + x, j + y).SetRed(r);
			Pixel(i + x, j + y).SetGreen(g);
		}
		i++;
	}
	while (j > -20) {
		if (i + x >= 0 && j + y >= 0 && i + x < width && j + y < height) {
			Pixel(i + x, j + y).SetBlue(b);
			Pixel(i + x, j + y).SetRed(r);
			Pixel(i + x, j + y).SetGreen(g);
		}
		j--;
	}
	while (i > -20) {
		if (i + x >= 0 && j + y >= 0 && i + x < width && j + y < height) {
			Pixel(i + x, j + y).SetBlue(b);
			Pixel(i + x, j + y).SetRed(r);
			Pixel(i + x, j + y).SetGreen(g);
		}
		i--;
	}
}


void R2Image::drawReferenceSpots(int m, float r, float g, float b) 
{
	//draw a box for reference spots
	//TODO: write description
	////////////////////////////////////////////////////////////////////////

	int x = givenLocs2D[m].xVal;
	int y = givenLocs2D[m].yVal;

	int i = -20;
	int j = -20;
	while (j < 20) {
		if (i + x >= 0 && j + y >= 0 && i + x < width && j + y < height) {
			Pixel(i + x, j + y).SetBlue(b);
			Pixel(i + x, j + y).SetRed(r);
			Pixel(i + x, j + y).SetGreen(g);
		}
		j++;
	}

	while (i < 20) {
		if (i + x >= 0 && j + y >= 0 && i + x < width && j + y < height) {
			Pixel(i + x, j + y).SetBlue(b);
			Pixel(i + x, j + y).SetRed(r);
			Pixel(i + x, j + y).SetGreen(g);
		}
		i++;
	}
	while (j > -20) {
		if (i + x >= 0 && j + y >= 0 && i + x < width && j + y < height) {
			Pixel(i + x, j + y).SetBlue(b);
			Pixel(i + x, j + y).SetRed(r);
			Pixel(i + x, j + y).SetGreen(g);
		}
		j--;
	}
	while (i > -20) {
		if (i + x >= 0 && j + y >= 0 && i + x < width && j + y < height) {
			Pixel(i + x, j + y).SetBlue(b);
			Pixel(i + x, j + y).SetRed(r);
			Pixel(i + x, j + y).SetGreen(g);
		}
		i--;
	}
}

void R2Image::box(int x, int y, float r, float g, float b)
{
	for (int i = -3; i < 3; i++) {
		for (int j = -3; j < 3; j++) {
			if (i + x >= 0 && j + y >= 0 && i+x < width && j+y < height) {
				Pixel(i + x, j + y).SetBlue(b);
				Pixel(i + x, j + y).SetRed(r);
				Pixel(i + x, j + y).SetGreen(g);
			}
		}
	}
}


void R2Image::line(int x0, int y0, int x1, int y1, float r, float g, float b)
{
	//draw a line for marker tracks with boxes on the features
	//TODO: write description
	////////////////////////////////////////////////////////////////////////

	if (x0>x1) {
		int x = y1;
		y1 = y0;
		y0 = x;

		x = x1;
		x1 = x0;
		x0 = x;		
	}
	int deltax = x1 - x0;
	int deltay = y1 - y0;
	float error = 0;
	float deltaerr = 0.0;

	if (deltax != 0) deltaerr = fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
																		// note that this division needs to be done in a way that preserves the fractional part
	int y = y0;
	for (int x = x0; x <= x1; x++)
	{
		Pixel(x, y).Reset(r, g, b, 1.0);
		error = error + deltaerr;
		if (error >= 0.5)
		{
			if (deltay>0) y = y + 1;
			else y = y - 1;

			error = error - 1.0;
		}
	}
	
	box(x0, y0, r, g, b);
	box(x1, y1, r, g, b);

}


// Tools for Detection ////////////////////////////////////////////////


int compare(const void *p, const void *q)
{
	// qsort struct comparision function (price float field)
	////////////////////////////////////////////////////////////////////////

	Marker2D *a, *b;

	a = (Marker2D*)p;
	b = (Marker2D*)q;
	//compares harris values
	if (((*a).harrisValue) > ((*b).harrisValue)) return -1;
	else if (((*a).harrisValue) < ((*b).harrisValue)) return 1;
	else return 0;
}

// Detection ////////////////////////////////////////////////

Marker2D* R2Image::
MarkerDetection(R2Image* referenceImage)
{
	//implimented for 1 marker right now
	//TODO : iterate through all markers
	
	//radius for search window
	int radius = referenceImage->width/2;
	Marker2D *list;
	list = new Marker2D[numMarkers];

	double minSSD = (double) MAXINT;
	
	//iterate through pixels in main image in chunks of 2*radius 
	for (int j = 0; j < (int)(height / (2*radius)); j++) {
		for (int i = 0; i < (int)(width / (2 * radius)); i++) {

			double tempSSD = SSD(referenceImage, (i * 2 * radius + radius), (j * 2 * radius + radius), radius, radius, radius);

			if (tempSSD < minSSD) {
				minSSD = tempSSD;
				list[0].xVal = (i * 2 * radius + radius);
				list[0].yVal = (j * 2 * radius + radius);
				list[0].harrisValue = minSSD;
			}
		}
	}

	return list;

}


void R2Image::
trackMarkersOntoOtherImage(std::vector<int> originalXLocs, std::vector<int> originalYLocs, R2Image *otherImage, int radius) {

	//TODO: impliment the use of this 
	///////////////////////////////////////////////////////////////////////////////

	std::vector<Marker2D> tempVals(numMarkers);

	//for each marker
	for (int i = 0; i < numMarkers; i++) {

		//search through a window of pixels to find the one that matches
		int searchWindowWidth = (int)width*0.1;
		int searchWindowHeight = (int)height*0.1;

		int currentPixelLocX = originalXLocs[i];
		int currentPixelLocY = originalYLocs[i];

		Marker2D minSSDMarker;
		double minSSD = (double)MAXINT, currentSSD;

		for (int m = -searchWindowWidth; m <= searchWindowWidth; m++) {
			for (int n = -searchWindowHeight; n <= searchWindowHeight; n++) {

				//if the pixel is on frame
				if (validPixel(currentPixelLocX + m, currentPixelLocY + n)) {

					//calcualte the SSD
					currentSSD = SSD(otherImage, currentPixelLocX, currentPixelLocY, m, n, 3);

					//if its the smallest, copy the info into minSSDFeature
					if (currentSSD < minSSD) {
						minSSD = currentSSD;
						minSSDMarker.xVal = currentPixelLocX + m;
						minSSDMarker.yVal = currentPixelLocY + n;
					}

				}

			}
		}

		//copy the min SSD feature into the list
		tempVals[i] = minSSDMarker;
	}

	otherImage->SetMarkerLocs(tempVals);

}


int R2Image::
RANDSAC(int maxIteration, int threshold, Marker2D *originalFeatures, Marker2D *matchedFeatures, double **HMatrix) {

	//for a maxIteration of times,
	//choose 4 random features
	//use these and their calculated matches to produce a Homography matrix
	//iterate through every other feature in the image and calculate the predicted match location using H
	//if the difference between its actual match location and its predicted match location is less than the threshold, then its a good match
	//the H with the highest number of good matches is the accepted H

	///////////////////////////////////////////////////////////////////////////////

	int random[4];
	Match matches[150];

	for (int i = 0; i < maxIteration; i++) {

		int matchCount = 0;

		//choose 4 unique random points
		random[0] = rand() % 150;

		random[1] = rand() % 150;
		while (random[1] == random[0]) { random[1] = rand() % 150; }

		random[2] = rand() % 150;
		while (random[2] == random[0] || random[2] == random[1]) { random[2] = rand() % 150; }

		random[3] = rand() % 150;
		while (random[3] == random[2] || random[3] == random[1] || random[3] == random[0]) { random[3] = rand() % 150; }

		//create input and output vectors
		std::vector<int> input(8);
		input[0] = originalFeatures[random[0]].xVal;
		input[1] = originalFeatures[random[0]].yVal;

		input[2] = originalFeatures[random[1]].xVal;
		input[3] = originalFeatures[random[1]].yVal;

		input[4] = originalFeatures[random[2]].xVal;
		input[5] = originalFeatures[random[2]].yVal;

		input[6] = originalFeatures[random[3]].xVal;
		input[7] = originalFeatures[random[3]].yVal;

		std::vector<int> output(8);
		output[0] = matchedFeatures[random[0]].xVal;
		output[1] = matchedFeatures[random[0]].yVal;

		output[2] = matchedFeatures[random[1]].xVal;
		output[3] = matchedFeatures[random[1]].yVal;

		output[4] = matchedFeatures[random[2]].xVal;
		output[5] = matchedFeatures[random[2]].yVal;

		output[6] = matchedFeatures[random[3]].xVal;
		output[7] = matchedFeatures[random[3]].yVal;

		printf(" mapped < %d %d > to < %d %d > \n mapped < %d %d > to < %d %d > \n mapped < %d %d > to < %d %d > \n mapped < %d %d > to < %d %d > \n \n",
			input[0], input[1], output[0], output[1], input[2], input[3], output[2], output[3], input[4], input[5], output[4], output[5], input[6], input[7], output[6], output[7]);

		HMatrix = calcCameraMatrix(input, output, 4);

		for (int j = 0; j < 150; j++) {
			if (random[0] != j && random[1] != j && random[2] != j && random[3] != j) {

				//calculate the the predicted mapping of the point
				int *vector = new int[3];
				vector[0] = originalFeatures[j].xVal;
				vector[1] = originalFeatures[j].yVal;
				vector[2] = 1;

				int* mappedVector = new int[3];

				printf("before! my vector is < %d , %d , %d > \n", vector[0], vector[1], vector[2]);
				mappedVector = matMult(vector, 3, HMatrix, 3, 3);
				printf("after! my vector is < %d , %d , %d > \n", mappedVector[0], mappedVector[1], mappedVector[2]);

				mappedVector = Normalize(mappedVector, 3);

				double distance = (mappedVector[0] - output[0])*(mappedVector[0] - output[0]) +
					(mappedVector[1] - output[1])*(mappedVector[1] - output[1]) +
					(mappedVector[2] - output[2])*(mappedVector[2] - output[2]);
				distance = sqrt(distance);

				if (distance < threshold) {
					matchCount++;
				}

			}
		}

		//record the number of matches and the map that created them into a list of matches
		matches[i].numMatches = matchCount;
		matches[i].hMatrix = HMatrix;

	}

	Match bestMatch;
	int maxMatches = -1;
	//find the map that has the highest number of matches, stored in good choice
	for (int i = 0; i < maxIteration; i++) {
		if (matches[i].numMatches > maxMatches) {
			maxMatches = matches[i].numMatches;
			bestMatch = matches[i];
		}
	}

	HMatrix = bestMatch.hMatrix;

	return bestMatch.numMatches;
}



double** R2Image::
calcCameraMatrix(std::vector<int> input, std::vector<int> output, int size)
{
	//size = number of matched points
	//input = 3D points
	//output = 2D points

	// build the [size*2 x 12] matrix of equations
	//need the *2 because each output has an x and a y

	///////////////////////////////////////////////////////////////////////////////

	double** linEquations = dmatrix(1, size * 2, 1, 12);
	for (int i = 1; i <= size * 2; i++) {
		
		linEquations[i][1] = (double)((input[i - 1]));
		linEquations[i][2] = (double)(input[i]);
		linEquations[i][3] = (double)(input[i + 1]);
		linEquations[i][4] = 1.0;
		linEquations[i][5] = 0.0;
		linEquations[i][6] = 0.0;
		linEquations[i][7] = 0.0;
		linEquations[i][8] = 0.0;
		linEquations[i][9] = (double)(-output[i] * input[i - 1]);
		linEquations[i][10] = (double)(-output[i] * input[i]);
		linEquations[i][11] = (double)(-output[i] * input[i+1]);
		linEquations[i][12] = (double)(-output[i]);
		i++;
		linEquations[i][1] = 0.0;
		linEquations[i][2] = 0.0;
		linEquations[i][3] = 0.0;
		linEquations[i][4] = 0.0;
		linEquations[i][5] = (double)((input[i - 2]));
		linEquations[i][6] = (double)(input[i - 1]);
		linEquations[i][7] = (double)(input[i]);
		linEquations[i][8] = 1.0;
		linEquations[i][9] = (double)(-output[i] * input[i - 2]);
		linEquations[i][10] = (double)(-output[i] * input[i - 1]);
		linEquations[i][11] = (double)(-output[i] * input[i]);
		linEquations[i][12] = (double)(-output[i]);
		//printf("input : ( %lf , %lf ) output : ( %lf , %lf ) \n", input[i-1], input[i], output[i-1], output[i]);
	}

	// compute the SVD
	double singularValues[13]; // 1..12

	double **nullspaceMatrix;

	nullspaceMatrix = dmatrix(1, 12, 1, 12);
	
	svdcmp(linEquations, (size*2), 12, singularValues, nullspaceMatrix);

	// find the smallest singular value:
	int smallestIndex = 1;
	for (int i = 1; i < 13; i++) { 
		if (singularValues[i] < singularValues[smallestIndex]) { 
			smallestIndex = i; 
		} 
	}
	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	
	double **HomographyMatrix;
	HomographyMatrix = dmatrix(0, 4, 0, 3);

	HomographyMatrix[0][0] = nullspaceMatrix[1][smallestIndex];
	HomographyMatrix[0][1] = nullspaceMatrix[2][smallestIndex];
	HomographyMatrix[0][2] = nullspaceMatrix[3][smallestIndex];
	HomographyMatrix[1][0] = nullspaceMatrix[4][smallestIndex];
	HomographyMatrix[1][1] = nullspaceMatrix[5][smallestIndex];
	HomographyMatrix[1][2] = nullspaceMatrix[6][smallestIndex];
	HomographyMatrix[2][0] = nullspaceMatrix[7][smallestIndex];
	HomographyMatrix[2][1] = nullspaceMatrix[8][smallestIndex];
	HomographyMatrix[2][2] = nullspaceMatrix[9][smallestIndex];
	HomographyMatrix[3][0] = nullspaceMatrix[10][smallestIndex];
	HomographyMatrix[3][1] = nullspaceMatrix[11][smallestIndex];
	HomographyMatrix[3][2] = nullspaceMatrix[12][smallestIndex];

	printf("Camera matrix: \n | %f %f %f |\n | %f %f %f |\n | %f %f %f | \n | %f %f %f | \n",
	HomographyMatrix[0][0], HomographyMatrix[0][1], HomographyMatrix[0][2],
	HomographyMatrix[1][0], HomographyMatrix[1][1], HomographyMatrix[1][2],
	HomographyMatrix[2][0], HomographyMatrix[2][1], HomographyMatrix[2][2],
	HomographyMatrix[3][0], HomographyMatrix[3][1], HomographyMatrix[3][2]);

	return HomographyMatrix;

}

void R2Image::
find3DLocation(std::vector<Marker2D> screenLocs, std::vector<Marker3D> knownLocs, std::vector<Marker2D> unknownLocs, int numKnown, int numUnknown)
{
	//screenLocs = screen locations for the reference spots
	//knownLocs = known locations for the reference spots
	//unknownLocs = screen locs of markers
	//numKnown = number of reference spots
	//numUnknown = number of markers

	////////////////////////////////////////////////////////////////////////

	std::vector<int> output (2 * numKnown);
	std::vector<int> input (3 * numKnown);

	for (int i = 0; i < numKnown; i++) { 

		output[i] = screenLocs[i].xVal;
		i++;
		output[i] = screenLocs[i - 1].yVal;

	}

	for (int i = 0; i < numKnown; i++) {

		input[i] = knownLocs[i].xVal;
		i++;
		input[i] = knownLocs[i-1].yVal;
		i++;
		input[i] = knownLocs[i-2].zVal;

	}
	
	double **HomographyMatrix;
	HomographyMatrix = dmatrix(0, 4, 0, 3);
	HomographyMatrix = calcCameraMatrix(input, output, numKnown);

	FILE * pFile;
	pFile = fopen("3DLocs.txt", "w");

	for (int i = 0; i < numUnknown; i++) {

		//TODO: NEED TO RECHEK THESE SIZES
		int* markerLoc = new int[3];
		markerLoc[0] = unknownLocs[i].xVal;
		markerLoc[1] = unknownLocs[i].yVal;
		markerLoc[2] = 1;
		int* newMarkerLoc = new int[4];
		newMarkerLoc = matMult(markerLoc, 3, HomographyMatrix, 4, 3);

		//output these 3D locations to a file
		fprintf(pFile, "%d %d %d %d\n", newMarkerLoc[0], newMarkerLoc[1], newMarkerLoc[2], newMarkerLoc[3]);

		//might need to normalize?
	}

}


/*	Original Skeleton 
	From here down
*/


////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }  

  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}





////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format 
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format 
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format 
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format 
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format 
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header 
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}


	

int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 95, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
	
  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}






