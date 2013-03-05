#include <glib.h>
#include <gtk/gtk.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <iostream>

// Library under test.
#include "gtkanimview.h"
#include "gtkimagescrollwin.h"
#include "gtkimagetooldragger.h"

#include "selector.h"
#include "PCAP.h"

//#include "cv.h"
//#include "highgui.h"

#define DEBUG_MODE

//using namespace cv;

short MAXVAL;
short *val_matrix = NULL;
//IplImage *img = NULL;

int nx, ny;
int* dontcare = NULL;

// inline GetVal
inline short GetVal(int x, int y) {
  return val_matrix[x+y*nx];
}

// inline NotDontCare
inline int NotDontCare(int x, int y) {
  return dontcare[x+y*nx] == 0;
}

// inline GetDontCare
inline int GetDontCare(int x, int y) {
  return dontcare[x+y*nx];
}

// inline SetDontCare
inline void SetDontCare(int x, int y, int val) {
  dontcare[x+y*nx] = val;
}

// LoadMRC
IplImage* LoadMRC(char* filename) {
  FILE *fin;
  fin = fopen(filename, "r");
  if (fin == NULL) {
#ifdef DEBUG_MODE
    g_print("File loading error...\n");
#endif
    return NULL;
  }
  
  int nz, mrc_type;
  fread(&nx, sizeof(int), 1, fin);
  fread(&ny, sizeof(int), 1, fin);
  fread(&nz, sizeof(int), 1, fin);
  fread(&mrc_type, sizeof(int), 1, fin);
  
#ifdef DEBUG_MODE
  g_print("%d %d %d %d\n", nx, ny, nz, mrc_type);
#endif
  
  fseek(fin, 1024, SEEK_SET);	

  // TODO: what if nz > 1
  IplImage *tmp_img = cvCreateImage(cvSize(nx, ny), IPL_DEPTH_16S, 1);
  short *tmp_ptr = (short*)tmp_img->imageData;

  for (int i = 0; i < ny; i++)
    for (int j = 0; j < nx; j++) {
      if (mrc_type == 1) {
	short tmp;
	fread(&tmp, sizeof(short), 1, fin);
	tmp_ptr[i*nx+j] = tmp;
      }
      else if (mrc_type == 2) {
	float tmp;
	fread(&tmp, sizeof(float), 1, fin);
	tmp_ptr[i*nx+j] = (short)tmp;
      }
      if (MAXVAL < tmp_ptr[i*nx+j])
	MAXVAL = tmp_ptr[i*nx+j];
    }
  
  fclose(fin);

  return tmp_img;
}

// LoadTIFF
IplImage *LoadTIFF(char *filename) {
  IplImage *tmp_img = cvLoadImage(filename, 2);

  cout << tmp_img->depth << " " << IPL_DEPTH_16S << endl; 
  cout << tmp_img->nChannels << endl;

  CvScalar s;
  for (int i = 0; i < 50; i++) {
    s=cvGet2D(tmp_img,i,0); // get the (i,j) pixel value
    printf("intensity=%f\n",s.val[0]);
  }

  nx = tmp_img->width;
  ny = tmp_img->height; 

#ifdef DEBUG_MODE
  g_print("%d %d\n", nx, ny);
#endif

  uchar *tmp_ptr = (uchar*)tmp_img->imageData;
  IplImage *tmp_img2 = cvCreateImage(cvSize(nx, ny), IPL_DEPTH_16S, 1);
  short *tmp_ptr2 = (short*)tmp_img2->imageData;
  for (int i = 0; i < ny; i++)
    for (int j = 0; j < nx; j++) {
      tmp_ptr2[i*nx+j] = tmp_ptr[i*nx+j];
      if (MAXVAL < tmp_ptr2[i*nx+j])
	MAXVAL = tmp_ptr2[i*nx+j];
    }
  cvReleaseImage(&tmp_img);

  return tmp_img2;
}

// LoadImage
GdkPixbuf* LoadImage(char* filename) {
#ifdef DEBUG_MODE
  g_print("Loading %s\n", filename);
#endif

  IplImage *tmp_img;

  MAXVAL = -1;
  string tmp = filename;
  if (strcmp(tmp.substr(tmp.size()-3, 3).c_str(), "mrc") == 0) {
    tmp_img = LoadMRC(filename);
  }
  else if (strcmp(tmp.substr(tmp.size()-3, 3).c_str(), "img") == 0) {
    // Load image
  } 
  else {
    tmp_img = LoadTIFF(filename);
  }

  nx = nx / RES;
  ny = ny / RES;
  if (img != NULL) {
    cvReleaseImage(&img);
  }
  img = cvCreateImage(cvSize(nx, ny), IPL_DEPTH_16S, 1);
  cvResize(tmp_img, img, CV_INTER_CUBIC);
  val_matrix = (short*)img->imageData;
  cvReleaseImage(&tmp_img);


  GdkPixbuf* anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, nx, ny);
  guchar *pixels = gdk_pixbuf_get_pixels(anim);
  int stride = gdk_pixbuf_get_rowstride(anim);
  int n_chans = gdk_pixbuf_get_n_channels(anim);
  int height = gdk_pixbuf_get_height(anim);
  int width = gdk_pixbuf_get_width(anim);
  // Load (or initiazlie dontcare)
  // TODO: load dontcare
  if (dontcare != NULL) {
    free(dontcare);
  } 
  dontcare = (int*)malloc(sizeof(int)*width*height);

  
#ifdef DEBUG_MODE
  g_print("%d %d %d %d\n", stride, n_chans, height, width);
#endif
  
  int i, j;
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      pixels[j * stride + i * n_chans] = (guchar)(GetVal(i,j) * 255 / MAXVAL);
      pixels[j * stride + i * n_chans+1] = (guchar)(GetVal(i,j) * 255 / MAXVAL);
      pixels[j * stride + i * n_chans+2] = (guchar)(GetVal(i,j) * 255 / MAXVAL);
      SetDontCare(i,j,0);
    }
  }

#ifdef DEBUG_MODE
  g_print("LoadImage is done.\n");
#endif
  
  return anim;
}



// SaveAll
void SaveAll(string path) {
  // TODO: name manipulation
  
  // 1. Save particles
  SaveParticles(path);

  // 2. Save mask and img
  FILE *fout1 = fopen((path+".mask").c_str(), "wt");
  FILE *fout2 = fopen((path+".img").c_str(), "wt");
  fprintf(fout1, "%d %d\n", ny, nx);
  fprintf(fout2, "%d %d\n", ny, nx);
  for (int i = 0; i < ny; i++) {
    for (int j = 0; j < nx; j++) {
      fprintf(fout1, "%d ", GetDontCare(j, i));
      fprintf(fout2, "%d ", GetVal(j, i));
    }
    fprintf(fout1, "\n");
    fprintf(fout2, "\n");
  }  
  fclose(fout1);
  fclose(fout2);
}



// NormalizeValMatrix
GdkPixbuf* NormalizeValMatrix() {
  int i, j, cou = 0;
  
  double mean_val = 0;	
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      if (NotDontCare(i, j)) {
	mean_val = mean_val + GetVal(i, j);
	cou++;
      }	
  mean_val = mean_val / cou;
  g_print("Mean Val: %.3f\n", mean_val);
  
  double std_val = 0;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) 
      if (NotDontCare(i, j)) {
	std_val = std_val + pow(GetVal(i, j) - mean_val, 2);
      }
  std_val = sqrt(std_val / cou);
  g_print("STD Val: %.3f\n", std_val);	
  
  GdkPixbuf* anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, nx, ny);
  guchar *pixels = gdk_pixbuf_get_pixels(anim);
  int stride = gdk_pixbuf_get_rowstride(anim);
  int n_chans = gdk_pixbuf_get_n_channels(anim);
  int height = gdk_pixbuf_get_height(anim);
  int width = gdk_pixbuf_get_width(anim);
  
  guchar val;  
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      if (NotDontCare(i,j)) {
	val = (guchar)(((GetVal(i,j) - mean_val) / std_val / 6.0 + 0.5) * 255);
	if (val < 0) 
	  val = 0;
	else if (val > 0xff)
	  val = 0xff;
	pixels[j * stride + i * n_chans] = val;
	pixels[j * stride + i * n_chans+1] = val;
	pixels[j * stride + i * n_chans+2] = val;
      } else {
	pixels[j * stride + i * n_chans] = 0x00;
	pixels[j * stride + i * n_chans+1] = 0xff;
	pixels[j * stride + i * n_chans+2] = 0x00;
      }
    }
  }
  
  return anim;
}
