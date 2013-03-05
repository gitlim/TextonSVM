#include <string>
//#include <tiffio.h>
#include <vector>
#include "Center.h"
#include "math/matrices/matrix.hh"

using namespace std;

using math::matrices::matrix;

class Image {
private:
  matrix<> *img;
  matrix<short> *mask;
  vector<Center> centers;
  vector<Center> good_centers;
  int resize_factor;
  int X, Y, center_n;
  int couMask;
  string filename;
  
public:
  Image();
  ~Image();

  inline matrix<>* get_img_pointer() {
    return img;
  }
	
  void LoadGoodCenters(string name);
  void LoadCenters(string name);
  void LoadImage(string name);
  void LoadNormImage(string name);
	
  static void WriteNormImage(string name, matrix<>& img, int resize_factor);
	
  void LoadMask(string name);	
  
  void Normalize();

	
	// inline get/set code
  inline int GetResizeFactor() {
    return resize_factor;
  }
  inline void SetResizeFactor(int rf) {
    resize_factor = rf;
  }
  inline vector<Center>& GetCenters() {
    return centers;
  }
  inline vector<Center>& GetGoodCenters() {
    return good_centers;
  }
  inline matrix<>& GetImage() {
    return *img;
  }
  
  inline short GetMask(int x, int y) {
    //    printf("%d %d %f\n", x, y, (*mask)(x,y));
    return (*mask)(x,y);
  }
  inline void SetMask(int x, int y, short val) {
    (*mask)(x,y) = val;
    //printf("s%d %d %f\n", x, y, (*mask)(x,y));
  }
  
  inline double GetImg(int x, int y) {
    // TODO: what should I return if mask is false		
    return (*img)(x, y);
  }
  inline void SetImg(int x, int y, double val) {
    (*img)(x,y) = val;
  }
  
  inline int GetWidth() {
    return X;
  }
  inline int GetHeight() {
    return Y;
  }
};
