#include <glib.h>
#include <gtk/gtk.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

// Library under test.
#include <gtkimageview/gtkanimview.h>
#include <gtkimageview/gtkimagescrollwin.h>
#include <gtkimageview/gtkimagetooldragger.h>

//#include "selector.h"
#include "TEXTONSVM.h"

//#define DEBUG_MODE

float *val_matrix = NULL;
float *det_matrix = NULL;
float *conf_matrix = NULL;
int nx, ny;
int Threshold = 0;
float min_diff, min_x, min_y; 


//int RES = 8;
int RES = -1;

extern GtkWidget *selected_size;
extern GdkPixmap *pixmap;
extern GtkWidget *plot_image_view;

short color_map[64*3];

vector<double> P,R,overallP,overallR;

// inline GetVal
inline short GetVal(float *mat, int x, int y) {
  return mat[x+y*nx];
}

inline void SetVal(float *mat, int x, int y, short val) {
  mat[x+y*nx] = val;
}

/*
// LoadMRC
GdkPixbuf* LoadMRC(char* filename) {
#ifdef DEBUG_MODE
  g_print("Loading %s\n", filename);
#endif
  
  FILE *fin;
  fin = fopen(filename, "r");
  if (fin == NULL) {
#ifdef DEBUG_MODE
    g_print("File loading error...\n");
#endif
    return NULL;
  }
  
#ifdef DEBUG_MODE
  g_print("I am fine here.\n");
#endif
  
  int nz, mrc_type;
  fread(&nx, sizeof(int), 1, fin);
  fread(&ny, sizeof(int), 1, fin);
  fread(&nz, sizeof(int), 1, fin);
  fread(&mrc_type, sizeof(int), 1, fin);
  
#ifdef DEBUG_MODE
  g_print("%d %d %d %d\n", nx, ny, nz, mrc_type);
#endif
  
  fseek(fin, 1024, SEEK_SET);	
  
  if (val_matrix != NULL)
    free(val_matrix);
  val_matrix = (float*)malloc(sizeof(float)*nx*ny*nz);
  
  if (dontcare != NULL)
    free(dontcare);
  dontcare = (int*)malloc(sizeof(int)*nx*ny*nz);
  
  if (mrc_type == 1) {
    fread((short *)val_matrix, sizeof(short), nx*ny*nz, fin);
  }
  else if (mrc_type == 2) {
    fread(val_matrix, sizeof(float), nx*ny*nz, fin);
  }
  
  fclose(fin);
  
  GdkPixbuf* anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, ny, nx);
  guchar *pixels = gdk_pixbuf_get_pixels(anim);
  int stride = gdk_pixbuf_get_rowstride(anim);
  int n_chans = gdk_pixbuf_get_n_channels(anim);
  int height = gdk_pixbuf_get_height(anim);
  int width = gdk_pixbuf_get_width(anim);
  
#ifdef DEBUG_MODE
  g_print("%d %d %d %d\n", stride, n_chans, height, width);
#endif
  
  int i, j;
  for (i = 0; i < width; i++) {
#ifdef DEBUG_MODE
    //		g_print("%d\n", i);
#endif
    for (j = 0; j < height; j++) {
      pixels[j * stride + i * n_chans] = (guchar)(GetVal(val_matrix,j,i) * 255 / 23246);
      pixels[j * stride + i * n_chans+1] = (guchar)(GetVal(val_matrix,j,i) * 255 / 23246);
      pixels[j * stride + i * n_chans+2] = (guchar)(GetVal(val_matrix,j,i) * 255 / 23246);
      SetDontCare(j,i,0);
    }
  }
  
#ifdef DEBUG_MODE
  g_print("LoadMRC is done.\n");
#endif
  
  return anim;
}
*/

GdkPixbuf* LoadImage(string filename) {
#ifdef DEBUG_MODE
  printf("Loading %s\n", filename.c_str());
#endif

  FILE *fin;
  fin = fopen((filename.substr(0, filename.find_last_of(".")) + ".norm").c_str(), "r");
  if (fin == NULL) {
#ifdef DEBUG_MODE
    g_print("File loading error...\n");
#endif
    return NULL;
  }

  int tmp_x, tmp_y;
  // TODO: check dimension mataches with nx and ny
  fscanf(fin, "%d %d %d\n", &tmp_x, &tmp_y, &RES);

  nx = tmp_x;
  ny = tmp_y;

  //  if (val_matrix != NULL)
  //    free(val_matrix);
  //  val_matrix = (float*)malloc(sizeof(float)*nx*ny);
  
  GdkPixbuf* anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, nx, ny);
  guchar *pixels = gdk_pixbuf_get_pixels(anim);
  int stride = gdk_pixbuf_get_rowstride(anim);
  int n_chans = gdk_pixbuf_get_n_channels(anim);
  int height = gdk_pixbuf_get_height(anim);
  int width = gdk_pixbuf_get_width(anim);
  
#ifdef DEBUG_MODE
  g_print("%d %d %d %d\n", stride, n_chans, height, width);
#endif
  
  //  g_print("%d\n", GetVal(val_matrix, 50, 50));

  g_print("%d %d", height, nx);

  int i, j;
  float tmp;
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {      
      fscanf(fin, "%f", &tmp);

      guchar val = (guchar)((tmp+3)*255/6);

      //      val_matrix[i+j*nx] = tmp;
      pixels[j * stride + i * n_chans] = val;
      pixels[j * stride + i * n_chans+1] = val;
      pixels[j * stride + i * n_chans+2] = val;
    }
  }

  fclose(fin);
  
#ifdef DEBUG_MODE
  g_print("LoadDetection is done.\n");
#endif
  
  return anim;
}

// LoadDetection
GdkPixbuf* LoadDetection(string filename) {
#ifdef DEBUG_MODE
  g_print("Loading %s\n", filename.c_str());
#endif
	
  FILE *fin;
  fin = fopen((filename.substr(0,filename.find_last_of(".")) + ".dmap").c_str(), "r");
  if (fin == NULL) {
#ifdef DEBUG_MODE
    g_print("File loading error...\n");
#endif
    return NULL;
  }

  int tmp_x, tmp_y;
  // TODO: check dimension mataches with nx and ny
  fscanf(fin, "%d %d\n", &tmp_x, &tmp_y);  
  
  if (det_matrix != NULL)
    free(det_matrix);
  det_matrix = (float*)malloc(sizeof(float)*nx*ny);
  
  GdkPixbuf* anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, nx, ny);
  guchar *pixels = gdk_pixbuf_get_pixels(anim);
  int stride = gdk_pixbuf_get_rowstride(anim);
  int n_chans = gdk_pixbuf_get_n_channels(anim);
  int height = gdk_pixbuf_get_height(anim);
  int width = gdk_pixbuf_get_width(anim);
  
#ifdef DEBUG_MODE
  g_print("%d %d %d %d\n", stride, n_chans, height, width);
#endif
  
  g_print("%d\n", GetVal(det_matrix, 50, 50));

  g_print("%d %d", height, nx);

  int i, j, tmp2;
  float tmp;
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      fscanf(fin, "%f", &tmp);
      det_matrix[i+j*nx] = tmp;
      tmp2 = (int)tmp / 4;
      // TODO: color_map MAX_LENGTH = 64
      if (tmp2 > 64) {
	tmp2 = 64;
      }
      pixels[j * stride + i * n_chans] = (guchar)color_map[tmp2*3];
      pixels[j * stride + i * n_chans+1] = (guchar)color_map[tmp2*3+1];
      pixels[j * stride + i * n_chans+2] = (guchar)color_map[tmp2*3+2];
    }
  }

  fclose(fin);
  
#ifdef DEBUG_MODE
  g_print("LoadDetection is done.\n");
#endif
  
  return anim;
}

// LoadConfidence
GdkPixbuf* LoadConfidence(string filename) {
#ifdef DEBUG_MODE
  g_print("Loading %s\n", filename.c_str());
#endif
	
  FILE *fin;
  fin = fopen((filename.substr(0, filename.find_last_of(".")) + ".confidence").c_str(), "r");
  if (fin == NULL) {
#ifdef DEBUG_MODE
    g_print("File loading error...\n");
#endif
    return NULL;
  }

  int tmp_x, tmp_y;
  // TODO: check dimension mataches with nx and ny
  fscanf(fin, "%d %d\n", &tmp_x, &tmp_y);  
  
  if (conf_matrix != NULL)
    free(conf_matrix);
  conf_matrix = (float*)malloc(sizeof(float)*nx*ny);
  
  GdkPixbuf* anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, nx, ny);
  guchar *pixels = gdk_pixbuf_get_pixels(anim);
  int stride = gdk_pixbuf_get_rowstride(anim);
  int n_chans = gdk_pixbuf_get_n_channels(anim);
  int height = gdk_pixbuf_get_height(anim);
  int width = gdk_pixbuf_get_width(anim);
  
#ifdef DEBUG_MODE
  g_print("%d %d %d %d\n", stride, n_chans, height, width);
#endif
  
  g_print("%d\n", GetVal(conf_matrix, 50, 50));
  
  int i, j, tmp2;
  float tmp;
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      fscanf(fin, "%f", &tmp);
      conf_matrix[i+j*nx] = tmp;
      tmp2 = (int)tmp / 4;
      // TODO: color_map MAX_LENGTH = 64
      if (tmp2 > 64) {
	tmp2 = 64;
      }
      pixels[j * stride + i * n_chans] = (guchar)color_map[tmp2*3];
      pixels[j * stride + i * n_chans+1] = (guchar)color_map[tmp2*3+1];
      pixels[j * stride + i * n_chans+2] = (guchar)color_map[tmp2*3+2];
    }
  }

  fclose(fin);
  
#ifdef DEBUG_MODE
  g_print("LoadConfidence is done.\n");
#endif
  
  return anim;
}

// ApplyThreshold
int ApplyThreshold() {
  Threshold = gtk_range_get_value((GtkRange*)selected_size) * 255.0 / 100;
  
  GdkPixbuf* anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, nx, ny);
  guchar *pixels = gdk_pixbuf_get_pixels(anim);
  int stride = gdk_pixbuf_get_rowstride(anim);
  int n_chans = gdk_pixbuf_get_n_channels(anim);
  int height = gdk_pixbuf_get_height(anim);
  int width = gdk_pixbuf_get_width(anim);
  
#ifdef DEBUG_MODE
  g_print("%d %d %d %d\n", stride, n_chans, height, width);
#endif
  
  int i, j, tmp;
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      tmp = GetVal(conf_matrix,i,j);	
      tmp = tmp * (tmp > Threshold);
      tmp = tmp / 4;
      pixels[j * stride + i * n_chans] = (guchar)color_map[tmp*3];
      pixels[j * stride + i * n_chans+1] = (guchar)color_map[tmp*3+1];
      pixels[j * stride + i * n_chans+2] = (guchar)color_map[tmp*3+2];
    }
  }
  
#ifdef DEBUG_MODE
  g_print("ApplyThreshold is done.\n");
#endif

  return 0;
}


void helperPlot(guchar *pixels, float x, float y, int height, int width, int stride, int n_chans, guchar R, guchar G, guchar B, int dx = 0, int dy = 0) {
  int i, j; 
  for (i = 0-dx/2; i <= 1+dx/2; i++)
    for (j = 0-dy/2; j <= 1+dy/2; j++) {
      if (height-(y+j) < 0) 
	continue;
      if (x+i < 0)
	continue;
      if (height-(y+j) > height)
	continue;
      if (x+i > width)
	continue;
      pixels[(int)(height-(y+j)) * stride + (int)(x+i) * n_chans] = R;
      pixels[(int)(height-(y+j)) * stride + (int)(x+i) * n_chans+1] = G;
      pixels[(int)(height-(y+j)) * stride + (int)(x+i) * n_chans+2] = B;
    }

  /*
  if (abs(x - Threshold) < min_diff) {
    min_diff = abs(x - Threshold);
    min_x = x;
    min_y = y;
  }
  */
}



//GdkPixbuf* PlotPR(const vector<double>& P, const vector <double>& R, const vector <double>& overallP, const vector<double>& overallR) { 
GdkPixbuf* PlotPR(const vector<double>& P, const vector<double>& R) {
  int nx_p = 160;
  int ny_p = 160;

  int i, j;

  GdkPixbuf* anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, nx_p, ny_p);
  guchar *pixels = gdk_pixbuf_get_pixels(anim);
  int stride = gdk_pixbuf_get_rowstride(anim);
  int n_chans = gdk_pixbuf_get_n_channels(anim);
  int height = gdk_pixbuf_get_height(anim);
  int width = gdk_pixbuf_get_width(anim);

  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      pixels[j * stride + i * n_chans] = (guchar)0x00;
      pixels[j * stride + i * n_chans+1] = (guchar)0x00;
      pixels[j * stride + i * n_chans+2] = (guchar)0x00;
    }
  }

  int n;
  n = P.size();

  if (n == 0) {
    return anim;
  }


  /*
  FILE *fin = fopen((filename.substr(0,filename.find_last_of(".")) + ".perf").c_str(), "r");
  if (!fin) {
    printf("Error: perf file not exist.\n");
    return anim;
  }
    
  fscanf(fin, "%d\n", &n);

  float P[n], R[n];
  for (i = 0; i < n; i++) {
    fscanf(fin, "%f", &(P[i]));
  }
  for (i = 0; i < n; i++) {
    fscanf(fin, "%f", &(R[i]));
  }
  fclose(fin);
  */

  /*
  // sort
  float tmp;
  for (i = 0; i < n; i++)
    for (j = i+1; j < n; j++)
      if (R[i] > R[j]) {
	tmp = P[i];
	P[i] = P[j];
	P[j] = tmp;

	tmp = R[i];
	R[i] = R[j];
	R[j] = tmp;
      }
  */

  /*  
  float x1 = 0, y1 = P[0] * height;
  float x2, y2;
  float step_x = width/n;
  float dx, dy;
  float slope, x, y;
  for (i = 0; i < n-1; i++) {
    x2 = x1 + step_x;
    y2 = P[i+1] * height;

    dx = x2 - x1;
    dy = y2 - y1;
    if (dx > dy) {
      slope = dy/dx;
      y = y1;
      for (x = x1; x < x2; x++) {
	helperPlot(pixels, x, y, height, stride, n_chans, 0, 0, 0xFF);
	y = min(y2, y + slope);	
      }
    } else {
      slope = dx/dy;
      x = x1;
      for (y = y1; y < y2; y++) {
	helperPlot(pixels, x, y, height, stride, n_chans, 0, 0, 0xFF);
	x = min(x2, x + slope);
      }
    }
    x1 = x2;
    y1 = y2;
  }

  x1 = 0, y1 = R[0] * height;
  for (i = 0; i < n-1; i++) {
    x2 = x1 + step_x;
    y2 = R[i+1] * height;

    dx = x2 - x1;
    dy = y2 - y1;
    if (dx > -dy) {
      slope = dy/dx;
      y = y1;
      for (x = x1; x < x2; x++) {
	helperPlot(pixels, x, y, height, stride, n_chans, 0, 0xFF, 0);
	y = max(y2, y + slope);
      }
    } else {
      slope = dx/dy;
      x = x1;
      for (y = y1; y > y2; y--) {
	helperPlot(pixels, x, y, height, stride, n_chans, 0, 0xFF, 0);
	x = min(x2, x + slope);
      }
    }
    x1 = x2;
    y1 = y2;
  }
  */

  float x1, y1;
  float x2, y2;
  float x, y, slope;
  float dx, dy;

  float d_thres;
  //int how_many_away;

  float thres_step = 255.0 / n;

  min_diff = 9999;
  x1 = R[1] * width;
  y1 = P[1] * height;
  min_x = x1;
  min_y = y1;
  for (i = 1; i < n-1; i++) {
    //    printf("%d\n", i);
    x2 = R[i+1] * width;
    y2 = P[i+1] * height;
    if (isnan(y2))
      break;
    if (isnan(x2))
      break;

    dx = x2 - x1;
    dy = y2 - y1;

    if ((std::abs(dx) < 0.000001) && (std::abs(dy) < 0.000001))
      continue;

    d_thres = Threshold - (i * thres_step);
    if ((d_thres >= 0) && (d_thres < thres_step)) {
      min_x = x1;
      min_y = y1;
    }
    //    printf("%.3f %.3f\n", x1, y1);
    //    printf("%.3f %.3f\n", x2, y2);
    if (std::abs(dx) > std::abs(dy)) {
      slope = dy/dx;
      y = y1;
      if (x1 > x2) {
	for (x = x1; x >= x2; x--) {
	  helperPlot(pixels, x, y, height, width, stride, n_chans, 0, 0xFF, 0);
	  y = y - slope; 
	}
      } else {
	for (x = x1; x <= x2; x++) {
	  helperPlot(pixels, x, y, height, width, stride, n_chans, 0, 0xFF, 0);
	  y = y + slope; 
	}
      }
    } else { 
      slope = dx/dy;
      x = x1;
      if (y1 < y2) {
	for (y = y1; y <= y2; y++) {
	  helperPlot(pixels, x, y, height, width, stride, n_chans, 0, 0xFF, 0);
	  x = x + slope;
	}
      } else {
	for (y = y1; y >= y2; y--) {
	  helperPlot(pixels, x, y, height, width, stride, n_chans, 0, 0xFF, 0);
	  x = x - slope;
	}
      }
    }
    x1 = x2;
    y1 = y2;
    //    printf(".\n");
  }

  helperPlot(pixels, min_x, min_y, height, width, stride, n_chans, 0xFF, 0xFF, 0xFF, 8, 8);
  helperPlot(pixels, min_x, min_y, height, width, stride, n_chans, 0xFF, 0, 0);

  return anim;
}

GdkPixbuf* PlotImagePR() {
  return PlotPR(P, R);
}

GdkPixbuf* PlotOverallPR() {
  return PlotPR(overallP, overallR);
}

void compute_perf(string path) {
  // TODO: what should be radius here?
  // TODO: nx/ny should be max of particles' x and y rather?
  calculate_image_PR(nx, ny, P, R);
  calculate_overall_PR(path, nx, ny, overallP, overallR);

  /*
  FILE *f_exist = fopen((path.substr(0, path.find_last_of(".")) + ".box").c_str(), "rt");
  if (!f_exist) {
    int maxF_ind = -1;
    double maxF = -9999999;
    for (int i = 0; i < overallP.size(); i++) {
      double F = 2 * (overallP[i] * overallR[i]) / (overallP[i] + overallR[i]);
      if (F > maxF) {
	maxF = F;
	maxF_ind = i;
      }
    }
    cout << "maxF @ " << maxF_ind << endl;
    update_box_file(path, maxF_ind);
  } else {
    fclose(f_exist);
  }
  */
}
