
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
#include "svd.h"
#include <cmath>

//TODO: evaluate if options[] is even necessary anymore
// Program arguments
static char options[] =
"  -help\n"
"  -svdTest\n"
"  -trackMarkers <file:referenceImage> <int:numFrames> <int:numMarkers>\n";

/* Original Skeleton */
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

int* matMult(int* vector, int size, double** mat, int rows, int columns)
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
				result[i] += (int)(vector[j] * mat[j][i]);
			}
		}
		return result;
	}
}

double** matMult(double** a, double** b, int startA, int rowsA, int columnsA, int startB, int rowsB, int columnsB)
{

	double** result;
	result = dmatrix(0, rowsA, 0, columnsB);

	for (int i = 0; i <= rowsA; i++) {
		for (int j = 0; j <= columnsB; j++) {
			result[i][j] = 0.0;
		}
	}

	if (rowsB != columnsA) {
		printf("CANNOT MULTIPLY THESE \n");
		return NULL;
	}

	for (int i = startA; i < rowsA + startA; i++) {
		for (int j = startB; j < columnsB + startB; j++) {

			for (int k = 0; k < columnsA; k++) {
				result[i - startA][j - startB] += (a[i][startA + k] * b[startB + k][j]);
			}
		}
	}
	
	return result;
}

int* Normalize(int* vector, int size)
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

static int* find3DLocation(double** CameraMatrix, std::vector<Marker2D> unknownLocs, int numUnknown)
{
	//CameraMatrix = transformation matrix for specifc camera
	//unknownLocs = screen locs of markers
	//numUnknown = number of markers

	////////////////////////////////////////////////////////////////////////

	for (int i = 0; i < numUnknown; i++) {

		//TODO: NEED TO RECHEK THESE SIZES
		int* markerLoc = new int[3];
		markerLoc[0] = unknownLocs[i].xVal;
		markerLoc[1] = unknownLocs[i].yVal;
		markerLoc[2] = 1;
		return matMult(markerLoc, 3, CameraMatrix, 4, 3);

		//might need to normalize?
	}
	return NULL;
}

static double**
calcCameraMatrix(std::vector<int> input, std::vector<int> output, int size)
{
	//size = number of matched points
	//input = 3D points
	//output = 2D points

	// build the [size*2 x 12] matrix of equations
	//need the *2 because each output has an x and a y

	///////////////////////////////////////////////////////////////////////////////

	//double** linEquations = dmatrix(1, size * 2, 1, 12);
	//int j = 0, k = 0; //j = input counter, k = output counter
	//for (int i = 1; i <= size * 2; i++) { //i = linEquations counter
	//	
	//	linEquations[i][1] = (double)((input[j]));
	//	linEquations[i][2] = (double)(input[j + 1]);
	//	linEquations[i][3] = (double)(input[j + 2]);
	//	linEquations[i][4] = 1.0;
	//	linEquations[i][5] = 0.0;
	//	linEquations[i][6] = 0.0;
	//	linEquations[i][7] = 0.0;
	//	linEquations[i][8] = 0.0;
	//	linEquations[i][9] = (double)(-output[k] * input[j]);
	//	linEquations[i][10] = (double)(-output[k] * input[j + 1]);
	//	linEquations[i][11] = (double)(-output[k] * input[j + 2]);
	//	linEquations[i][12] = (double)(-output[k]);
	//	k++;
	//	i++;
	//	linEquations[i][1] = 0.0;
	//	linEquations[i][2] = 0.0;
	//	linEquations[i][3] = 0.0;
	//	linEquations[i][4] = 0.0;
	//	linEquations[i][5] = (double)((input[j]));
	//	linEquations[i][6] = (double)(input[j + 1]);
	//	linEquations[i][7] = (double)(input[j + 2]);
	//	linEquations[i][8] = 1.0;
	//	linEquations[i][9] = (double)(-output[k] * input[j]);
	//	linEquations[i][10] = (double)(-output[k] * input[j + 1]);
	//	linEquations[i][11] = (double)(-output[k] * input[j + 2]);
	//	linEquations[i][12] = (double)(-output[k]);
	//	
	//	printf("input : ( %d , %d, %d ) output : ( %d , %d ) \n", input[j], input[j+1], input[j+2], output[k-1], output[k]);
	//	k++;
	//	j += 3;
	//}

	//printf(" %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n ",
	//	linEquations[1][1], linEquations[1][2], linEquations[1][3], linEquations[1][4], linEquations[1][5], linEquations[1][6], linEquations[1][7], linEquations[1][8], linEquations[1][9], linEquations[1][10], linEquations[1][11], linEquations[1][12],
	//	linEquations[2][1], linEquations[2][2], linEquations[2][3], linEquations[2][4], linEquations[2][5], linEquations[2][6], linEquations[2][7], linEquations[2][8], linEquations[2][9], linEquations[2][10], linEquations[2][11], linEquations[2][12],
	//	linEquations[3][1], linEquations[3][2], linEquations[3][3], linEquations[3][4], linEquations[3][5], linEquations[3][6], linEquations[3][7], linEquations[3][8], linEquations[3][9], linEquations[3][10], linEquations[3][11], linEquations[3][12],
	//	linEquations[4][1], linEquations[4][2], linEquations[4][3], linEquations[4][4], linEquations[4][5], linEquations[4][6], linEquations[4][7], linEquations[4][8], linEquations[4][9], linEquations[4][10], linEquations[4][11], linEquations[4][12],
	//	linEquations[5][1], linEquations[5][2], linEquations[5][3], linEquations[5][4], linEquations[5][5], linEquations[5][6], linEquations[5][7], linEquations[5][8], linEquations[5][9], linEquations[5][10], linEquations[5][11], linEquations[5][12],
	//	linEquations[6][1], linEquations[6][2], linEquations[6][3], linEquations[6][4], linEquations[6][5], linEquations[6][6], linEquations[6][7], linEquations[6][8], linEquations[6][9], linEquations[6][10], linEquations[6][11], linEquations[6][12],
	//	linEquations[7][1], linEquations[7][2], linEquations[7][3], linEquations[7][4], linEquations[7][5], linEquations[7][6], linEquations[7][7], linEquations[7][8], linEquations[7][9], linEquations[7][10], linEquations[7][11], linEquations[7][12],
	//	linEquations[8][1], linEquations[8][2], linEquations[8][3], linEquations[8][4], linEquations[8][5], linEquations[8][6], linEquations[8][7], linEquations[8][8], linEquations[8][9], linEquations[8][10], linEquations[8][11], linEquations[8][12],
	//	linEquations[9][1], linEquations[9][2], linEquations[9][3], linEquations[9][4], linEquations[9][5], linEquations[9][6], linEquations[9][7], linEquations[9][8], linEquations[9][9], linEquations[9][10], linEquations[9][11], linEquations[9][12],
	//	linEquations[10][1], linEquations[10][2], linEquations[10][3], linEquations[10][4], linEquations[10][5], linEquations[10][6], linEquations[10][7], linEquations[10][8], linEquations[10][9], linEquations[10][10], linEquations[10][11], linEquations[10][12],
	//	linEquations[11][1], linEquations[11][2], linEquations[11][3], linEquations[11][4], linEquations[11][5], linEquations[11][6], linEquations[11][7], linEquations[11][8], linEquations[11][9], linEquations[11][10], linEquations[11][11], linEquations[11][12],
	//	linEquations[12][1], linEquations[12][2], linEquations[12][3], linEquations[12][4], linEquations[12][5], linEquations[12][6], linEquations[12][7], linEquations[12][8], linEquations[12][9], linEquations[12][10], linEquations[12][11], linEquations[12][12]);


	// compute the SVD
	double singularValues[13]; // 1..12

	double **nullspaceMatrix;

	nullspaceMatrix = dmatrix(1, 12, 1, 12);
	//
	//svdcmp(linEquations, (size*2), 12, singularValues, nullspaceMatrix);
	//


	//// find the smallest singular value:
	//int smallestIndex = 1;
	//printf("singular values : \n");
	//for (int i = 1; i < 13; i++) { 
	//	printf(" %f ", singularValues[i]);
	//	if (singularValues[i] < singularValues[smallestIndex]) { 
	//		smallestIndex = i; 
	//	} 
	//}
	//// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	//printf("\n nullspaceMatrix : \n");
	//for (int i = 1; i < 13; i++) {
	//	printf("  |");
	//	for (int j = 1; j < 13; j++) {
	//		printf(" %lf ", nullspaceMatrix[i][j]);
	//	}
	//	printf("|\n");
	//}

	double **HomographyMatrix;
	HomographyMatrix = dmatrix(0, 3, 0, 4);

	//HomographyMatrix[0][0] = nullspaceMatrix[1][smallestIndex];
	//HomographyMatrix[0][1] = nullspaceMatrix[2][smallestIndex];
	//HomographyMatrix[0][2] = nullspaceMatrix[3][smallestIndex];
	//HomographyMatrix[0][3] = nullspaceMatrix[4][smallestIndex];
	//HomographyMatrix[1][0] = nullspaceMatrix[5][smallestIndex];
	//HomographyMatrix[1][1] = nullspaceMatrix[6][smallestIndex];
	//HomographyMatrix[1][2] = nullspaceMatrix[7][smallestIndex];
	//HomographyMatrix[1][3] = nullspaceMatrix[8][smallestIndex];
	//HomographyMatrix[2][0] = nullspaceMatrix[9][smallestIndex];
	//HomographyMatrix[2][1] = nullspaceMatrix[10][smallestIndex];
	//HomographyMatrix[2][2] = nullspaceMatrix[11][smallestIndex];
	//HomographyMatrix[2][3] = nullspaceMatrix[12][smallestIndex];
	double** rightHandSide;
	rightHandSide = dmatrix(0, 1, 0, 12);
	double** linEquations = dmatrix(1, size * 2, 1, 12);
	int j = 0, k = 0; //j = input counter, k = output counter
	for (int i = 1; i <= size * 2; i++) { //i = linEquations counter

		linEquations[i][1] = (double)((output[j]));
		linEquations[i][2] = (double)(output[j + 1]);
		linEquations[i][3] = 1.0;
		linEquations[i][4] = 0.0;
		linEquations[i][5] = 0.0;
		linEquations[i][6] = 0.0;
		linEquations[i][7] = 0.0;
		linEquations[i][8] = 0.0;
		linEquations[i][9] = 0.0;
		linEquations[i][10] = 0.0;
		linEquations[i][11] = 0.0;
		linEquations[i][12] = 0.0;
		rightHandSide[0][i - 1] = (double)(input[k]);
		i++; k++;
		linEquations[i][1] = 0.0;
		linEquations[i][2] = 0.0;
		linEquations[i][3] = 0.0;
		linEquations[i][4] = (double)((output[j]));
		linEquations[i][5] = (double)(output[j + 1]);
		linEquations[i][6] = 1.0;
		linEquations[i][7] = 0.0;
		linEquations[i][8] = 0.0;
		linEquations[i][9] = 0.0;
		linEquations[i][10] = 0.0;
		linEquations[i][11] = 0.0;
		linEquations[i][12] = 0.0;
		rightHandSide[0][i - 1] = (double)(input[k]);
		i++; k++;
		linEquations[i][1] = 0.0;
		linEquations[i][2] = 0.0;
		linEquations[i][3] = 0.0;
		linEquations[i][4] = 0.0;
		linEquations[i][5] = 0.0;
		linEquations[i][6] = 0.0;
		linEquations[i][7] = (double)(output[j]);
		linEquations[i][8] = (double)(output[j + 1]);
		linEquations[i][9] = 1.0;
		linEquations[i][10] = 0.0;
		linEquations[i][11] = 0.0;
		linEquations[i][12] = 0.0;
		rightHandSide[0][i - 1] = (double)(input[k]);
		i++; k++;
		linEquations[i][1] = 0.0;
		linEquations[i][2] = 0.0;
		linEquations[i][3] = 0.0;
		linEquations[i][4] = 0.0;
		linEquations[i][5] = 0.0;
		linEquations[i][6] = 0.0;
		linEquations[i][7] = 0.0;
		linEquations[i][8] = 0.0;
		linEquations[i][9] = 0.0;
		linEquations[i][10] = (double)(output[j]);
		linEquations[i][11] = (double)(output[j + 1]);
		linEquations[i][12] = 1.0;
		rightHandSide[0][i - 1] = 1.0;

		printf("input : ( %d , %d, %d ) output : ( %d , %d ) \n", input[j], input[j + 1], input[j + 2], output[k - 1], output[k]); // flip flopped them now
		j += 4;
	}


	/******************************************************************************/
	svdcmp(linEquations, 12, 12, singularValues, nullspaceMatrix);
	/*******************************************************************************
	Given a matrix linEquations[1..13][1..13], this routine computes its singular value
	decomposition, linEquations = U.W.VT.  The matrix U replaces linEquations on output.  The diagonal
	matrix of singular values W is output as a vector singularValues[1..13].  The matrix V (not
	the transpose VT) is output as nullspaceMatrix[1..13][1..13].
	*******************************************************************************/

	//pseudo-inverse calculations
	double **pseudoInverse;
	pseudoInverse = dmatrix(1, 12, 1, 12);

	double **U;
	U = dmatrix(1, 12, 1, 12);
	U = linEquations;

	double **W;
	W = dmatrix(1, 12, 1, 12);
	for (int i = 1; i <= 12; i++) {
		W[i][i] = singularValues[i];
	}

	double **V;
	V = dmatrix(1, 12, 1, 12);
	V = nullspaceMatrix;
	//Confirmed normal	

	double **Wi;
	Wi = dmatrix(1, 12, 1, 12);
	for (int i = 1; i <= 12; i++) {
		//printf("my singular value is : %lf \n", singularValues[i]);
		for (int j = 1; j <= 12; j++) {
			if (fabs(singularValues[i])>0.001 && i == j) {
				//printf("I found a non-zero one! \n");
				Wi[i][i] = 1 / singularValues[i];
			}
			else {
				Wi[i][j] = 0.0;
			}
		}
	}

	double **UT;
	UT = dmatrix(1, 12, 1, 12);
	for (int i = 1; i <= 12; i++) {
		for (int j = 1; j <= 12; j++) {
			UT[i][j] = U[j][i];
		}
	}

	pseudoInverse = matMult(V, Wi, 1, 12, 12, 1, 12, 12);
	printf("check a spot %lf \n", pseudoInverse[3][3]);
	pseudoInverse = matMult(pseudoInverse, UT, 1, 12, 12, 1, 12, 12);

	FILE * extraInfo;
	extraInfo = fopen("extraInfo.txt", "w");

	fprintf(extraInfo, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n ",
		pseudoInverse[1][1], pseudoInverse[1][2], pseudoInverse[1][3], pseudoInverse[1][4], pseudoInverse[1][5], pseudoInverse[1][6], pseudoInverse[1][7], pseudoInverse[1][8], pseudoInverse[1][9], pseudoInverse[1][10], pseudoInverse[1][11], pseudoInverse[1][12],
		pseudoInverse[2][1], pseudoInverse[2][2], pseudoInverse[2][3], pseudoInverse[2][4], pseudoInverse[2][5], pseudoInverse[2][6], pseudoInverse[2][7], pseudoInverse[2][8], pseudoInverse[2][9], pseudoInverse[2][10], pseudoInverse[2][11], pseudoInverse[2][12],
		pseudoInverse[3][1], pseudoInverse[3][2], pseudoInverse[3][3], pseudoInverse[3][4], pseudoInverse[3][5], pseudoInverse[3][6], pseudoInverse[3][7], pseudoInverse[3][8], pseudoInverse[3][9], pseudoInverse[3][10], pseudoInverse[3][11], pseudoInverse[3][12],
		pseudoInverse[4][1], pseudoInverse[4][2], pseudoInverse[4][3], pseudoInverse[4][4], pseudoInverse[4][5], pseudoInverse[4][6], pseudoInverse[4][7], pseudoInverse[4][8], pseudoInverse[4][9], pseudoInverse[4][10], pseudoInverse[4][11], pseudoInverse[4][12],
		pseudoInverse[5][1], pseudoInverse[5][2], pseudoInverse[5][3], pseudoInverse[5][4], pseudoInverse[5][5], pseudoInverse[5][6], pseudoInverse[5][7], pseudoInverse[5][8], pseudoInverse[5][9], pseudoInverse[5][10], pseudoInverse[5][11], pseudoInverse[5][12],
		pseudoInverse[6][1], pseudoInverse[6][2], pseudoInverse[6][3], pseudoInverse[6][4], pseudoInverse[6][5], pseudoInverse[6][6], pseudoInverse[6][7], pseudoInverse[6][8], pseudoInverse[6][9], pseudoInverse[6][10], pseudoInverse[6][11], pseudoInverse[6][12],
		pseudoInverse[7][1], pseudoInverse[7][2], pseudoInverse[7][3], pseudoInverse[7][4], pseudoInverse[7][5], pseudoInverse[7][6], pseudoInverse[7][7], pseudoInverse[7][8], pseudoInverse[7][9], pseudoInverse[7][10], pseudoInverse[7][11], pseudoInverse[7][12],
		pseudoInverse[8][1], pseudoInverse[8][2], pseudoInverse[8][3], pseudoInverse[8][4], pseudoInverse[8][5], pseudoInverse[8][6], pseudoInverse[8][7], pseudoInverse[8][8], pseudoInverse[8][9], pseudoInverse[8][10], pseudoInverse[8][11], pseudoInverse[8][12],
		pseudoInverse[9][1], pseudoInverse[9][2], pseudoInverse[9][3], pseudoInverse[9][4], pseudoInverse[9][5], pseudoInverse[9][6], pseudoInverse[9][7], pseudoInverse[9][8], pseudoInverse[9][9], pseudoInverse[9][10], pseudoInverse[9][11], pseudoInverse[9][12],
		pseudoInverse[10][1], pseudoInverse[10][2], pseudoInverse[10][3], pseudoInverse[10][4], pseudoInverse[10][5], pseudoInverse[10][6], pseudoInverse[10][7], pseudoInverse[10][8], pseudoInverse[10][9], pseudoInverse[10][10], pseudoInverse[10][11], pseudoInverse[10][12],
		pseudoInverse[11][1], pseudoInverse[11][2], pseudoInverse[11][3], pseudoInverse[11][4], pseudoInverse[11][5], pseudoInverse[11][6], pseudoInverse[11][7], pseudoInverse[11][8], pseudoInverse[11][9], pseudoInverse[11][10], pseudoInverse[11][11], pseudoInverse[11][12],
		pseudoInverse[12][1], pseudoInverse[12][2], pseudoInverse[12][3], pseudoInverse[12][4], pseudoInverse[12][5], pseudoInverse[12][6], pseudoInverse[12][7], pseudoInverse[12][8], pseudoInverse[12][9], pseudoInverse[12][10], pseudoInverse[12][11], pseudoInverse[12][12]);

	;

	double** resultVector;
	resultVector = dmatrix(0, 0, 0, 11);

	resultVector = matMult(rightHandSide, pseudoInverse, 0, 1, 12, 1, 12, 12);
	printf("result after call : %lf \n", resultVector[0][0]);
	HomographyMatrix[0][0] = resultVector[0][0];
	HomographyMatrix[0][1] = resultVector[0][1];
	HomographyMatrix[0][2] = resultVector[0][2];
	HomographyMatrix[0][3] = resultVector[0][3];
	HomographyMatrix[1][0] = resultVector[0][4];
	HomographyMatrix[1][1] = resultVector[0][5];
	HomographyMatrix[1][2] = resultVector[0][6];
	HomographyMatrix[1][3] = resultVector[0][7];
	HomographyMatrix[2][0] = resultVector[0][8];
	HomographyMatrix[2][1] = resultVector[0][9];
	HomographyMatrix[2][2] = resultVector[0][10];
	HomographyMatrix[2][3] = resultVector[0][11];

	printf("Camera matrix: \n | %f %f %f %f |\n | %f %f %f %f |\n | %f %f %f %f | \n ",
		HomographyMatrix[0][0], HomographyMatrix[0][1], HomographyMatrix[0][2], HomographyMatrix[0][3],
		HomographyMatrix[1][0], HomographyMatrix[1][1], HomographyMatrix[1][2], HomographyMatrix[1][3],
		HomographyMatrix[2][0], HomographyMatrix[2][1], HomographyMatrix[2][2], HomographyMatrix[2][3]);

	return HomographyMatrix;
}

/* */

int 
main(int argc, char **argv)
{
  // Look for help
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-help")) {
      ShowUsage();
    }
	if (!strcmp(argv[i], "-calcHomography")) {
		R2Image *image = new R2Image();
		
		FILE *f = fopen(argv[2], "r");
		
		if (f == NULL) {
			printf("I didn't read this right... \n");
		}

		std::vector<int> input(3), output(3);
		for (int i = 0; i < 8; i += 2) {
			fscanf(f, "%d %d %d %d ", &input[i], &input[i + 1], &output[i], &output[i + 1]);
		}
		calcCameraMatrix(input, output, 8);
		exit(-1);
	}
  }

  // Read input and output image filenames
  if (argc < 3)  ShowUsage();
  argv++, argc--; // First argument is program name

  char *input_image_name_camera1 = *argv; argv++, argc--; 
  char *output_image_name_camera1 = *argv; argv++, argc--;

  char *input_image_name_camera2 = *argv; argv++, argc--;
  char *output_image_name_camera2 = *argv; argv++, argc--;

  // Allocate images
  R2Image *image_camera1 = new R2Image();
  R2Image *image_camera2 = new R2Image();

  if (!image_camera1) {
    fprintf(stderr, "Unable to allocate image\n");
    exit(-1);
  }
  if (!image_camera2) {
	  fprintf(stderr, "Unable to allocate image\n");
	  exit(-1);
  }


  // Read input images
  if (!image_camera1->Read(input_image_name_camera1)) {
    fprintf(stderr, "Unable to read image from %s\n", input_image_name_camera1);
    exit(-1);
  }
  if (!image_camera2->Read(input_image_name_camera2)) {
	  fprintf(stderr, "Unable to read image from %s\n", input_image_name_camera2);
	  exit(-1);
  }

  // Initialize sampling method
  int sampling_method = R2_IMAGE_POINT_SAMPLING;

  // Parse arguments and perform operations 
  while (argc > 0) {
	if (!strcmp(*argv, "-trackMarkers")) {

		//read in the reference image
		CheckOption(*argv, argc, 2);
		R2Image *referenceImageMarker1 = new R2Image(argv[1]);
		CheckOption(*argv, argc, 3);
		R2Image *referenceImageMarker2 = new R2Image(argv[2]);
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

		image_camera1->SetNumMarkers(numMarkers);
		image_camera2->SetNumMarkers(numMarkers);

		//read in the 3D and screen locations of the reference points
		std::vector<Marker3D> known3D(numKnown);
		std::vector<Marker2D> known2D_camera1(numKnown);
		std::vector<Marker2D> known2D_camera2(numKnown);

		FILE * mappingFile3D;
		mappingFile3D = fopen("orig_3DLocs.txt", "r");
		FILE * mappingFile1;
		mappingFile1 = fopen("camera1_mapping.txt", "r");
		FILE * mappingFile2;
		mappingFile2 = fopen("camera2_mapping.txt", "r");

		for (int i = 0; i < numKnown; i++) {
			printf("Please enter the x y and z measurements for known location %d \n", i+1);
			fscanf(mappingFile3D, "%d %d %d", &known3D[i].xVal, &known3D[i].yVal, &known3D[i].zVal);
			printf("Now please enter where they are mapped to on the screen for camera 1\n");
			fscanf(mappingFile1, "%d %d", &known2D_camera1[i].xVal, &known2D_camera1[i].yVal);
			printf("Now please enter where they are mapped to on the screen for camera 2\n");
			fscanf(mappingFile2, "%d %d", &known2D_camera2[i].xVal, &known2D_camera2[i].yVal);

		}

		//TODO: do we need to have these values saved in every image? 
		//image->SetGivenLocs3D(known3D);
		//image->SetGivenLocs2D(known2D);

		printf("Thank you! Your locations are : \n");
		for (int i = 0; i < numKnown; i++) {
			printf("< %d, %d, %d > -> < %d, %d > and < %d, %d >\n", known3D[i].xVal, known3D[i].yVal, known3D[i].zVal, known2D_camera1[i].xVal, known2D_camera1[i].yVal, known2D_camera2[i].xVal, known2D_camera2[i].yVal);
		}

		//check to make sure everything is fine
		printf("\nNUMBER OF FRAMES: %d\n", numFrames);
		printf("NUMBER OF MARKERS: %d\n", numMarkers);
		printf("NUMBER OF KNOWN POINTS: %d\n\n", numKnown);
		printf("input image name for camera 1: %s\n", input_image_name_camera1);
		printf("input image name for camera 2: %s\n", input_image_name_camera2);
		printf("output image name for camera 1: %s\n\n", output_image_name_camera1);
		printf("output image name for camera 2: %s\n\n", output_image_name_camera2);


		// extract input and output filepaths
		std::string sInput1 = input_image_name_camera1; //this should be the first frame
		std::string sOutput1 = output_image_name_camera1;
		
		std::string sInput2 = input_image_name_camera2; //this should be the first frame
		std::string sOutput2 = output_image_name_camera2;

		//find the last instance of "." which should be directly before the file type
		//index is the index of where the "." is
		int index = sInput1.find_last_of(".");
		if (index == -1) {
			fprintf(stderr, "Unable to find extension in %s\n", input_image_name_camera1);
			exit(-1);
		}
		index = sInput2.find_last_of(".");
		if (index == -1) {
			fprintf(stderr, "Unable to find extension in %s\n", input_image_name_camera2);
			exit(-1);
		}

		//splits string between the name and the file type/extension
		std::string extension1 = sInput1.substr(index);
		std::string extension2 = sInput2.substr(index);

		//search for the frame number (specifically using the first frame
		index = sInput1.find("000001");
		if (index == -1) {
			fprintf(stderr, "Unable to find '000001' (6-digit padding) in %s\n", input_image_name_camera1);
			exit(-1);
		}
		index = sInput2.find("000001");
		if (index == -1) {
			fprintf(stderr, "Unable to find '000001' (6-digit padding) in %s\n", input_image_name_camera2);
			exit(-1);
		}

		//get the first part of the string until the number, should be the original name thats the same in each input frame
		std::string inputPath1 = sInput1.substr(0, index);
		std::string inputPath2 = sInput2.substr(0, index);

		//search for the frame number in the output file name
		index = sOutput1.find("000001");
		if (index == -1) {
			fprintf(stderr, "Unable to find '000001' (6-digit padding) in %s\n", output_image_name_camera1);
			exit(-1);
		}
		index = sOutput2.find("000001");
		if (index == -1) {
			fprintf(stderr, "Unable to find '000001' (6-digit padding) in %s\n", output_image_name_camera2);
			exit(-1);
		}

		//get the first part of the string till the number, should be original name thats the same in each output frame
		std::string outputPath1 = sOutput1.substr(0, index);
		std::string outputPath2 = sOutput2.substr(0, index);

		std::string number; // padding = 6 digits

		//TODO: impliment better tracking
		//const int height = image_camera1->Height();
		//const int sqRadius = height*0.1;

		//calculate matrix for camera 1 ///////////////////////////////////
		std::vector<int> output1(2 * numKnown);
		std::vector<int> output2(2 * numKnown);
		std::vector<int> input(3 * numKnown);
		//j = input/output counter
		for (int i = 0, j = 0; i < numKnown; i++, j++) {

			output1[j] = known2D_camera1[i].xVal;
			output2[j] = known2D_camera2[i].xVal;
			j++;
			output1[j] = known2D_camera1[i].yVal;
			output2[j] = known2D_camera2[i].yVal;
		}

		for (int i = 0, j = 0; i < numKnown; i++, j++) {

			input[j] = known3D[i].xVal;
			j++;
			input[j] = known3D[i].yVal;
			j++;
			input[j] = known3D[i].zVal;

		}

		double **CameraMatrix1;
		CameraMatrix1 = dmatrix(0, 4, 0, 3);
		CameraMatrix1 = calcCameraMatrix(input, output1, numKnown);

		double **CameraMatrix2;
		CameraMatrix2 = dmatrix(0, 4, 0, 3);
		CameraMatrix2 = calcCameraMatrix(input, output2, numKnown);
		///////////////////////////////////////////////////////////////////

		// image_camera1 = first frame
		FILE * actualLocs1;
		actualLocs1 = fopen("3DLocs_camera1.txt", "w");
		FILE * actualLocs2;
		actualLocs2 = fopen("3DLocs_camera2.txt", "w");
		
		//MARKER DETECTION////////////////////////////////////////////////////////////
		//find the markers on the screen
		std::vector<Marker2D> MarkersA = image_camera1->MarkerDetection(referenceImageMarker1);
		std::vector<Marker2D> MarkersB = image_camera2->MarkerDetection(referenceImageMarker2);

		//set them as the values for the image fields
		std::vector<Marker2D> points1(numMarkers);
		std::vector<Marker2D> points2(numMarkers);

		for (int i = 0; i < numMarkers; i++) {
			printf("Found marker %d at ( %d, %d ) for camera 1\n", i+1, MarkersA[i].xVal, MarkersA[i].yVal);
			points1[i] = MarkersA[i];
			printf("Found marker %d at ( %d, %d ) for camera 2\n", i+1, MarkersB[i].xVal, MarkersB[i].yVal);
			points2[i] = MarkersB[i];
		}

		image_camera1->SetMarkerLocs(points1);
		image_camera2->SetMarkerLocs(points2);
		printf("Found %d markers in first frames\n", numMarkers);

		//TODO: do I need these loops?
		//draw markers and reference spots
		for (int i = 0; i < numMarkers; i++) {
			image_camera1->drawMarkers(i, 1, 0, 1); 
			image_camera2->drawMarkers(i, 1, 0, 1);
		}
		for (int i = 0; i < numKnown; i++) {
			image_camera1->drawReferenceSpots(known2D_camera1, i, 0, 1, 1);
			image_camera2->drawReferenceSpots(known2D_camera2, i, 0, 1, 1);
		}

		// Write output image
		if (!image_camera1->Write(output_image_name_camera1)) {
			fprintf(stderr, "Unable to read image from %s\n", output_image_name_camera1);
			exit(-1);
		}
		if (!image_camera2->Write(output_image_name_camera2)) {
			fprintf(stderr, "Unable to read image from %s\n", output_image_name_camera2);
			exit(-1);
		}
		//////////////////////////////////////////////////////////////////////////////////////////

		//MAP TO REAL SPACE//////////////////////////////////////////////////////////////////////
		// Write 3D loc for frame 1
		int* newMarkerLoc_camera1 = new int[4];
		newMarkerLoc_camera1 = find3DLocation(CameraMatrix1, image_camera1->MarkerLocs2D(), numMarkers);
		fprintf(actualLocs1, "%d %d %d %d\n", newMarkerLoc_camera1[0], newMarkerLoc_camera1[1], newMarkerLoc_camera1[2], newMarkerLoc_camera1[3]);
		int* newMarkerLoc_camera2 = new int[4];
		newMarkerLoc_camera2 = find3DLocation(CameraMatrix2, image_camera2->MarkerLocs2D(), numMarkers);
		fprintf(actualLocs2, "%d %d %d %d\n", newMarkerLoc_camera2[0], newMarkerLoc_camera2[1], newMarkerLoc_camera2[2], newMarkerLoc_camera2[3]);

		//////////////////////////////////////////////////////////////////////////////////////////

		printf("Finished frame 1!\n");

		// iterate through frames
		// imageA = frame(i-1)
		// imageB = frame(i)

		//initialize imageB as frame 1
		R2Image *imageB1 = new R2Image(*image_camera1);
		R2Image *imageB2 = new R2Image(*image_camera2);
	
		for (int i = 2; i <= numFrames; i++) {
			// SETUP
			R2Image *imageA1;
			R2Image *imageA2;
			// copy imageB into imageA via operator=
			imageA1 = imageB1;
			imageA2 = imageB2;

			// Calculate new imageB
			number = "00000" + std::to_string(i);
			number = number.substr(number.length() - 6);
			imageB1 = new R2Image((inputPath1 + number + extension1).c_str());
			imageB1->SetNumMarkers(numMarkers);
			imageB2 = new R2Image((inputPath2 + number + extension2).c_str());
			imageB2->SetNumMarkers(numMarkers);
			
			// Track features from frame(i-1) to frame(i) and save into imageB
			//imageA->trackMarkersOntoOtherImage(imageA->MarkerLocs2DX(), imageA->MarkerLocs2DY(), imageB, sqRadius);
			//find the markers on the screen
			MarkersA = imageB1->MarkerDetection(referenceImageMarker1);
			MarkersB = imageB2->MarkerDetection(referenceImageMarker2);

			////////
			//set them as the values for the image fields

			for (int i = 0; i < numMarkers; i++) {
				printf("Found marker %d at ( %d, %d ) and ( %d, %d )\n", i+1, MarkersA[i].xVal, MarkersA[i].yVal, MarkersB[i].xVal, MarkersB[i].yVal);
				points1[i] = MarkersA[i];
				points2[i] = MarkersB[i];
			}

			imageB1->SetMarkerLocs(points1);
			imageB2->SetMarkerLocs(points2);
			
			////////
			
			//calculate the 3D location of the marker
			//for (int i = 0; i < numMarkers; i++) {
			//	newLocs[i] = imageB->MarkerLocs2D()[i];
			//}

			newMarkerLoc_camera1 = find3DLocation(CameraMatrix1, imageB1->MarkerLocs2D(), numMarkers);
			fprintf(actualLocs1, "%d %d %d %d\n", newMarkerLoc_camera1[0], newMarkerLoc_camera1[1], newMarkerLoc_camera1[2], newMarkerLoc_camera1[3]);
			newMarkerLoc_camera2 = find3DLocation(CameraMatrix2, imageB2->MarkerLocs2D(), numMarkers);
			fprintf(actualLocs2, "%d %d %d %d\n", newMarkerLoc_camera2[0], newMarkerLoc_camera2[1], newMarkerLoc_camera2[2], newMarkerLoc_camera2[3]);

			//draw markers and reference spots
			for (int i = 0; i < numMarkers; i++) {
				imageB1->drawMarkers(i, 1, 0, 1);
				imageB2->drawMarkers(i, 1, 0, 1);
			}
			for (int i = 0; i < numKnown; i++) {
				imageB1->drawReferenceSpots(known2D_camera1, i, 0, 1, 1);
				imageB2->drawReferenceSpots(known2D_camera2, i, 0, 1, 1);
			}

			if (!imageB1->Write((outputPath1 + number + extension1).c_str())) {
				fprintf(stderr, "Unable to read image from %s\n", (outputPath1 + number + extension1).c_str());
				exit(-1);
			}
			if (!imageB2->Write((outputPath2 + number + extension2).c_str())) {
				fprintf(stderr, "Unable to read image from %s\n", (outputPath2 + number + extension2).c_str());
				exit(-1);
			}

			//printf("Tracked features from frame%d to frame%d\n", i - 1, i);
			delete imageA1;
			delete imageA2;
		}
		delete imageB1;
		delete imageB2;
		//image = outputOrigImage;

	}
	
    else {
      // Unrecognized program argument
      fprintf(stderr, "image: invalid option: %s\n", *argv);
      ShowUsage();
    }
  }

  // Write output image
  if (!image_camera1->Write(output_image_name_camera1)) {
    fprintf(stderr, "Unable to read image from %s\n", output_image_name_camera1);
    exit(-1);
  }
  // Write output image
  if (!image_camera2->Write(output_image_name_camera2)) {
	  fprintf(stderr, "Unable to read image from %s\n", output_image_name_camera2);
	  exit(-1);
  }

  // Delete image
  delete image_camera1;
  delete image_camera2;

  // Return success
  return EXIT_SUCCESS;
}



