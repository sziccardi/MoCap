// Include file for image class
#ifndef R2_IMAGE_INCLUDED
#define R2_IMAGE_INCLUDED

#include <vector>

////////////////////////////////////////////////////////////////////////
// Constant definitions
////////////////////////////////////////////////////////////////////////

/*Originial Skeleton*/
typedef enum {
  R2_IMAGE_RED_CHANNEL,
  R2_IMAGE_GREEN_CHANNEL,
  R2_IMAGE_BLUE_CHANNEL,
  R2_IMAGE_ALPHA_CHANNEL,
  R2_IMAGE_NUM_CHANNELS
} R2ImageChannel;

typedef enum {
  R2_IMAGE_POINT_SAMPLING,
  R2_IMAGE_BILINEAR_SAMPLING,
  R2_IMAGE_GAUSSIAN_SAMPLING,
  R2_IMAGE_NUM_SAMPLING_METHODS
} R2ImageSamplingMethod;

typedef enum {
  R2_IMAGE_OVER_COMPOSITION,
  R2_IMAGE_IN_COMPOSITION,
  R2_IMAGE_OUT_COMPOSITION,
  R2_IMAGE_ATOP_COMPOSITION,
  R2_IMAGE_XOR_COMPOSITION,
} R2ImageCompositeOperation;

/* */

typedef struct Marker2D {
	int xVal;
	int yVal;
	double harrisValue;
} Marker2D;

typedef struct Marker3D {
	int xVal;
	int yVal;
	int zVal;
	double harrisValue;
} Marker3D;

typedef struct Match {
	int numMatches;
	double** hMatrix;
};

////////////////////////////////////////////////////////////////////////
// Class definition
////////////////////////////////////////////////////////////////////////

class R2Image {
 public:

  // Constructors/Destructor
  R2Image(void);
  R2Image(const char *filename); //for reference image
  R2Image(const char *filename, int numMarkers, std::vector<int> givenLocs3DX, std::vector<int> givenLocs3DY, std::vector<int> givenLocs3DZ); //for frames
  R2Image(int width, int height, std::vector<int> givenLocs3DX, std::vector<int> givenLocs3DY, std::vector<int> givenLocs3DZ);
  R2Image(int width, int height, const R2Pixel *pixels, std::vector<int> givenLocs3DX, std::vector<int> givenLocs3DY, std::vector<int> givenLocs3DZ);
  R2Image(const R2Image& image);
  ~R2Image(void);

  // Image properties
  int NumMarkers(void) const;
  int NPixels(void) const;
  int Width(void) const;
  int Height(void) const;

  std::vector<int> MarkerLocs2DX(void) const;
  std::vector<int> MarkerLocs2DY(void) const;

  std::vector<int> GivenLocs2DX(void) const;
  std::vector<int> GivenLocs2DY(void) const;

  std::vector<int> GivenLocs3DX(void) const;
  std::vector<int> GivenLocs3DY(void) const;
  std::vector<int> GivenLocs3DZ(void) const;

  // Pixel access/update
  R2Pixel& Pixel(int x, int y);
  R2Pixel *Pixels(void);
  R2Pixel *Pixels(int row);
  R2Pixel *operator[](int row);
  const R2Pixel *operator[](int row) const;
  void SetPixel(int x, int y,  const R2Pixel& pixel);
  void SetNumMarkers(int num);
  void SetMarkerLocs(std::vector<int> x, std::vector<int> y);
  void SetGivenLocs2D(std::vector<int> x, std::vector<int> y);
  void SetGivenLocs3D(std::vector<int> x, std::vector<int> y, std::vector<int> z);

  // Image processing
  R2Image& operator=(const R2Image& image);
  R2Image& operator*(const R2Image image);

  // Linear filtering operations
  R2Image filterMult(int kernel[3][3]);
  void line(int x0, int x1, int y0, int y1, float r, float g, float b);
  void drawMarkers(int i, float r, float g, float b);
  void drawReferenceSpots(int i, float r, float g, float b);
  void box(int x, int y, float r, float g, float b);
  double SSD(R2Image* otherImage, int x1, int y1, int x2, int y2, int range);
 
  // further operations
  void trackMarkersOntoOtherImage(std::vector<int> originalXLocs, std::vector<int> originalYLocs, R2Image *otherImage, int radius);
  int RANDSAC(int maxIteration, int threshold, Marker2D *originalFeatures, Marker2D *matchedFeatures, double **HMatrix);
  double** calcCameraMatrix(int *input, int * output, int size);
  void find3DLocation(Marker2D *screenLocs, Marker3D *knownLocs, Marker2D *unknownLocs, int numKnown, int numUnknown);
  bool validPixel(const int x, const int y);

  Marker2D* MarkerDetection(R2Image *referenceImage);

  // File reading/writing
  int Read(const char *filename);
  int ReadBMP(const char *filename);
  int ReadPPM(const char *filename);
  int ReadJPEG(const char *filename);
  int Write(const char *filename) const;
  int WriteBMP(const char *filename) const;
  int WritePPM(const char *filename, int ascii = 0) const;
  int WriteJPEG(const char *filename) const;

 private:

  //Fields
  R2Pixel *pixels;
  int numMarkers;
  int npixels;
  int width;
  int height;

  std::vector<int> markerLocs2DX;
  std::vector<int> markerLocs2DY;

  std::vector<int> givenLocs2DX;
  std::vector<int> givenLocs2DY;

  std::vector<int> givenLocs3DX;
  std::vector<int> givenLocs3DY;
  std::vector<int> givenLocs3DZ;

  int* matMult(int* vector, int size, double** mat, int rows, int columns);
  int* Normalize(int* vector, int size);
};

////////////////////////////////////////////////////////////////////////
// GETTERS/SETTERS
////////////////////////////////////////////////////////////////////////

// Setters ////////////////////////////////////////////////////////////////////////
inline void R2Image::
SetPixel(int x, int y, const R2Pixel& pixel)
{
	// Set pixel
	pixels[x*height + y] = pixel;
}

inline void R2Image::SetMarkerLocs(std::vector<int> x, std::vector<int> y)
{
	markerLocs2DX = x;
	markerLocs2DY = y;
}

inline void R2Image::SetGivenLocs2D(std::vector<int> x, std::vector<int> y)
{
	givenLocs2DX = x;
	givenLocs2DY = y;
}

inline void R2Image::SetGivenLocs3D(std::vector<int> x, std::vector<int> y, std::vector<int> z)
{
	givenLocs3DX = x;
	givenLocs3DY = y;
	givenLocs3DZ = z;
}

inline void R2Image::SetNumMarkers(int num) {
	numMarkers = num;
}

// Getters ////////////////////////////////////////////////////////////////////////

inline int R2Image::
NPixels(void) const
{
  // Return total number of pixels in the image
  return npixels;
}

inline int R2Image::
Width(void) const
{
  // Return width
  return width;
}

inline int R2Image::
Height(void) const
{
  // Return height
  return height;
}

inline std::vector<int> R2Image::
MarkerLocs2DX(void) const
{
	return markerLocs2DX;
}

inline std::vector<int> R2Image::
MarkerLocs2DY(void) const
{
	return markerLocs2DY;
}

inline std::vector<int> R2Image::
GivenLocs2DX(void) const
{
	return givenLocs2DX;
}

inline std::vector<int> R2Image::
GivenLocs2DY(void) const
{
	return givenLocs2DY;
}

inline std::vector<int> R2Image::
GivenLocs3DX(void) const
{
	return givenLocs3DX;
}

inline std::vector<int> R2Image::
GivenLocs3DY(void) const
{
	return givenLocs3DY;
}

inline std::vector<int> R2Image::
GivenLocs3DZ(void) const
{
	return givenLocs3DZ;
}

inline int R2Image::
NumMarkers(void) const
{
	return numMarkers;
}

inline R2Pixel& R2Image::
Pixel(int x, int y)
{
  // Return pixel value at (x,y)
  // (pixels start at lower-left and go in row-major order)
  return pixels[x*height + y];
}

inline R2Pixel *R2Image::
Pixels(void)
{
  // Return pointer to pixels for whole image 
  // (pixels start at lower-left and go in row-major order)
  return pixels;
}

inline R2Pixel *R2Image::
Pixels(int x)
{
  // Return pixels pointer for row at x
  // (pixels start at lower-left and go in row-major order)
  return &pixels[x*height];
}

////////////////////////////////////////////////////////////////////////
// OPERATORS 
////////////////////////////////////////////////////////////////////////

inline R2Pixel *R2Image::
operator[](int x) 
{
  // Return pixels pointer for row at x
  return Pixels(x);
}

inline const R2Pixel *R2Image::
operator[](int x) const
{
  // Return pixels pointer for row at x
  // (pixels start at lower-left and go in row-major order)
  return &pixels[x*height];
}


#endif
