#include <glib.h>
#include <gtk/gtk.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <imresize.h>

#include <iostream>


#include <selector.h>
#include <TEXTONSVM.h>

#include <tiffio.h>

//#define DEBUG_MODE

static int MAXVAL, MINVAL;
static float *img = NULL;
bool IMG_FLIPPED;

int RES;
static int nx, ny;
static int original_x, original_y;
static int* dontcare = NULL;

// inline GetVal
inline float GetVal(int x, int y) {
  return img[x+y*nx];
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
  if (!IMG_FLIPPED)
    dontcare[x+y*nx] = val;
  else
    dontcare[y+x*nx] = val;
}

float* LoadTIFF(const char *filename) {
  float *tmp_img;

  TIFF *tif = TIFFOpen(filename, "r");
  if (tif) {
    uint32 imagelength;
    uint32 imagewidth;
    tdata_t buf;
    uint32 row;

    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imagewidth);
    buf = _TIFFmalloc(TIFFScanlineSize(tif));
    printf("%d %d %d\n", TIFFScanlineSize(tif), imagelength, imagewidth);

    tmp_img = (float*)malloc(sizeof(float)*imagewidth*imagelength);
    nx = (int)imagewidth;
    ny = (int)imagelength;

    for (row = 0; row < imagelength; row++) {
      TIFFReadScanline(tif, buf, row); 
      for (unsigned int i = 0; i < imagewidth; i++) {
	if ((int)imagewidth < TIFFScanlineSize(tif)) {
	  unsigned char a = ((unsigned char*)buf)[i*2];
	  char b = ((char*)buf)[i*2+1];
	  //((float*)(tmp_img2->imageData + row*tmp_img2->widthStep))[i*tmp_img2->nChannels+0]=((float)b*256+a + 32768)/65536;
	  //	  ((float*)(tmp_img2->imageData + row*tmp_img2->widthStep))[i*tmp_img2->nChannels+0]=((float)b*256+a);
	  tmp_img[row*imagewidth+i] = (float)b*256+a;
	  if (b*256+a>MAXVAL) {
	    MAXVAL = b*256+a;
	  }
	  if (b*256+a<MINVAL) {
	    MINVAL = b*256+a;
	  } 
	}
      }
    }
    _TIFFfree(buf);
    TIFFClose(tif);
  }

  return tmp_img;
}

// TODO: convert everything of TEXTONSVM to MRC? or not
float* LoadMRCorTEXTONSVM(const char *filename) {
  //  IplImage *tmp_img2;
  float *tmp_img;

  FILE *fin = fopen(filename, "r");
  if (fin == NULL) {
#ifdef DEBUG_MODE
    g_print("File loading error...\n");
#endif
    return NULL;
  }

  int tmp_nx, tmp_ny, tmp_nz;
  fread(&tmp_nx, 1,sizeof(int), fin);
  fread(&tmp_ny, 1,sizeof(int), fin);
  fread(&tmp_nz, 1,sizeof(int), fin);

  // TODO: nz
  nx = (int)tmp_ny;
  ny = (int)tmp_nx;

  g_print("%d %d\n", nx, ny);

  int type;
  fread(&type, 1,sizeof(int), fin);

  g_print("Type: %d\n", type);

  fseek(fin, 1024, SEEK_SET);

  if (type == 1) {
    int tmp_var = nx;
    nx = ny;
    ny = tmp_var;

    int row;
    short *tmp = (short*)malloc(sizeof(short)*nx*ny);
    fread(tmp, nx*ny,sizeof(short), fin);
    int pointer=0;

    /*
    tmp_img = (float*)malloc(sizeof(float)*nx*ny);
    for (int col = 0; col < nx; col++) {
      for (row = 0; row < ny; row++) {
	// TODO: normalize images?
	//	((float*)(tmp_img2->imageData + row*tmp_img2->widthStep))[i*tmp_img2->nChannels+0]=tmp[pointer];
	tmp_img[row*nx+col] = tmp[pointer];
	if (tmp[pointer] > MAXVAL) {
	  MAXVAL = tmp[pointer];
	}
	if (tmp[pointer] < MINVAL) {
	  MINVAL = tmp[pointer];
	} 

	pointer++;
      }
    }
    */
    tmp_img = (float*)malloc(sizeof(float)*nx*ny);
    for (row = 0; row < ny; row++) {
      for (int i = 0; i < nx; i++) {
        // TODO: normalize images?                                                                                                                                                                       
        //      ((float*)(tmp_img2->imageData + row*tmp_img2->widthStep))[i*tmp_img2->nChannels+0]=tmp[pointer];                                                                                         
        tmp_img[row*nx+i] = tmp[pointer];
        if (tmp[pointer] > MAXVAL) {
          MAXVAL = tmp[pointer];
        }
        if (tmp[pointer] < MINVAL) {
          MINVAL = tmp[pointer];
        }

        pointer++;
      }
    }
    free(tmp);
  }
  else if (type == 2) {
    int tmp_var = nx;
    nx = ny;
    ny = tmp_var;

    int row;
    float *tmp = (float*)malloc(sizeof(float)*nx*ny);
    fread(tmp, nx*ny,sizeof(float), fin);
    int pointer=0;

    //    tmp_img2 = cvCreateImage(cvSize(nx, ny), IPL_DEPTH_32F, 1);
    tmp_img = (float*)malloc(sizeof(float)*nx*ny);
    for (row = 0; row < ny; row++) {
      for (int i = 0; i < nx; i++) {
        // TODO: normalize images?	
	//      ((float*)(tmp_img2->imageData + row*tmp_img2->widthStep))[i*tmp_img2->nChannels+0]=tmp[pointer]; 
	tmp_img[row*nx+i] = tmp[pointer];
        if (tmp[pointer] > MAXVAL) {
          MAXVAL = tmp[pointer];
        }
        if (tmp[pointer] < MINVAL) {
          MINVAL = tmp[pointer];
        }

        pointer++;
      }
    }
    /*
    for (int col = 0; col < nx; col++) {
      for (row = 0; row < ny; row++) {
	// TODO: normalize images?
	//	((float*)(tmp_img2->imageData + row*tmp_img2->widthStep))[i*tmp_img2->nChannels+0]=tmp[pointer];
	tmp_img[row*nx+col] = tmp[pointer];
	if (tmp[pointer] > MAXVAL) {
	  MAXVAL = tmp[pointer];
	}
	if (tmp[pointer] < MINVAL) {
	  MINVAL = tmp[pointer];
	} 

	pointer++;
      }
    }
    */
    free(tmp);
  }
  // TEXTONSVM type
  else if (type == 99) {
    //    int tmp_var = nx;
    //    nx = ny;
    //    ny = tmp_var;

    int row;
    float *tmp = (float*)malloc(sizeof(float)*nx*ny);
    fread(tmp, nx*ny,sizeof(float), fin);
    int pointer=0;

    //    tmp_img2 = cvCreateImage(cvSize(nx, ny), IPL_DEPTH_32F, 1);
    tmp_img = (float*)malloc(sizeof(float)*nx*ny);
    for (row = 0; row < ny; row++) {
      for (int i = 0; i < nx; i++) {
	// TODO: normalize images?
	//	((float*)(tmp_img2->imageData + row*tmp_img2->widthStep))[i*tmp_img2->nChannels+0]=tmp[pointer];
	tmp_img[row*nx+i] = tmp[pointer];
	if (tmp[pointer] > MAXVAL) {
	  MAXVAL = tmp[pointer];
	}
	if (tmp[pointer] < MINVAL) {
	  MINVAL = tmp[pointer];
	} 

	pointer++;
      }
    }
    free(tmp);
  }

  fclose(fin);

  return tmp_img;
}

// TODO: convert everything of TEXTONSVM to MRC? or not
float* LoadMRC_new(const char *filename) {
  int startSlice = 1; //default input
  int numSlices = 1; //default input
  int test = 0; //default input

  float *tmp_img;

  FILE *fin = fopen(filename, "r");
  if (fin == NULL) {
#ifdef DEBUG_MODE
    g_print("File loading error...\n");
#endif
    return NULL;
  }

  int a[10];
  fread(a, sizeof(int), 10, fin);
  /* debug output
  for (int i = 0; i < 10; i++)
    cout << a[i] << endl;
  */

  // TODO: ??
  if (abs(a[0]) > 1e5) {
    cout << "NO!!!!!!!" << endl;
    fclose(fin);
    fin=fopen(filename, "r");
  }

  int mode = a[3];

  // Get the next 12
  float b[12];  
  fread(b, sizeof(float), 12, fin);
  /* debug output 
  for (int i = 0; i < 12; i++)
    cout << b[i] << endl;
  */

  // b(4,5,6) angles
  float mi = b[9];
  float ma = b[10];
  float av = b[11];

  double s_rez = (double)b[0];

  //cout << mi << " " << ma << " " << av << " " << s_rez << endl;

  int c[30];
  fread(c, sizeof(int), 30, fin);
  /* debug output
  for (int i = 0; i < 30; i++) {
    cout << c[i] << endl;
  }
  */

  char s_chars[8];
  fread(s_chars, sizeof(char), 8, fin);
  /* debug output
  for (int i = 0; i < 8; i++)
    cout << d[i] << endl;
  */

  int e[2];
  fread(e, sizeof(int), 2, fin);
  //cout << e[0] << " " << e[1] << endl;

  int ns = 10;
  if (ns > e[1])
    ns = e[1];

  char s_header[ns][80];
  char garbage[80];
  for (int i = 0; i < 10; i++) {
    if (i < ns)
      fread(s_header[i], sizeof(char), 80, fin);
    else
	fread(garbage, sizeof(char), 80, fin);
    //cout << str[i] << endl;
  }

  int s_nx = a[0];
  int s_ny = a[1];
  int s_nz = a[2];
  cout << s_nx << " " << s_ny << " " << s_nz << endl;

  char ex_header[c[1]];
  if (c[1] > 0) {
    fread(ex_header, sizeof(char), c[1], fin);
    cout << ex_header << endl;
  }

  int pixbytes;
  switch (mode) {
  case 0:
    pixbytes = 1;
    break;
  case 1:
    pixbytes = 2;
    break;
  case 2:
    pixbytes = 4;
    break;
  case 6:
    break;
  default:
    pixbytes = 0;
  }

  int skipbytes = 0;
  int nz = s_nz;
  if (startSlice>1) {
    skipbytes = (startSlice-1)*s_nx*s_ny*pixbytes;
    fseek(fin, skipbytes, SEEK_CUR);
    nz = min(numSlices, s_nz-(startSlice-1));
  }
  int ndata = s_nx * s_ny * nz;
  // cout << ndata << " " << pixbytes << endl;

  nx = s_nx;
  ny = s_ny;

  if (mode == 0) {
    // int8
    int row;
    char *tmp = (char*)malloc(sizeof(char)*nx*ny);
    fread(tmp, nx*ny,sizeof(char), fin);
    int pointer=0;

    tmp_img = (float*)malloc(sizeof(float)*nx*ny);
    for (row = 0; row < ny; row++) {
      for (int i = 0; i < nx; i++) {
        tmp_img[row*nx+i] = tmp[pointer];
        if (tmp[pointer] > MAXVAL) {
          MAXVAL = tmp[pointer];
        }
        if (tmp[pointer] < MINVAL) {
          MINVAL = tmp[pointer];
        }

        pointer++;
      }
    }
    free(tmp);    
  } else if (mode == 1) {
    // int16
    int row;
    short *tmp = (short*)malloc(sizeof(short)*nx*ny);
    fread(tmp, nx*ny,sizeof(short), fin);
    int pointer=0;

    tmp_img = (float*)malloc(sizeof(float)*nx*ny);
    for (row = 0; row < ny; row++) {
      for (int i = 0; i < nx; i++) {
        tmp_img[row*nx+i] = tmp[pointer];
        if (tmp[pointer] > MAXVAL) {
          MAXVAL = tmp[pointer];
        }
        if (tmp[pointer] < MINVAL) {
          MINVAL = tmp[pointer];
        }

        pointer++;
      }
    }
    free(tmp);    
  } else if (mode == 2) {
    // float32
    int row;
    float *tmp = (float*)malloc(sizeof(float)*nx*ny);
    fread(tmp, nx*ny,sizeof(float), fin);
    int pointer=0;

    tmp_img = (float*)malloc(sizeof(float)*nx*ny);
    for (row = 0; row < ny; row++) {
      for (int i = 0; i < nx; i++) {
        tmp_img[row*nx+i] = tmp[pointer];
        if (tmp[pointer] > MAXVAL) {
          MAXVAL = tmp[pointer];
        }
        if (tmp[pointer] < MINVAL) {
          MINVAL = tmp[pointer];
        }

        pointer++;
      }
    }
    free(tmp);    
  } else if (mode == 6) {
    // error
  } else {
    // error
  }    

  fclose(fin);

  return tmp_img;
}

void convert_resize_image(string src_file, string target_file, string src_file2, string target_file2, int RES2) {
  cout << "Loading " + src_file + "   " + src_file2 + "\n";


  float *tmp_img;

  // TODO: handle RES2 properly
  MAXVAL = -1;
  if (strcmp(src_file.substr(src_file.size()-4, 4).c_str(), ".mrc") == 0) {
    //tmp_img = LoadMRCorTEXTONSVM(src_file.c_str());
    tmp_img = LoadMRC_new(src_file.c_str());
  }
  else {
    tmp_img = LoadTIFF(src_file.c_str());
  }

  if (img != NULL) {
    free(img);
  }
  if (RES2 != 1) {
    img = imresize(tmp_img, nx, ny, 1.0/RES2);
    free(tmp_img);
  } else {
    img = tmp_img;
  }

  cout << nx << " " << ny << endl;

  // Actual Image
  if (1) {
    float *tmp = (float*)malloc(sizeof(float)*nx*ny);
    int pointer=0;
    for (int i = 0; i < ny; i++) {
      for (int j = 0; j < nx; j++) {
	tmp[pointer++]= img[j+i*nx];
      }
    }  
    
    FILE *fout2 = fopen((target_file).c_str(), "wt");
    fwrite(&ny, 1,sizeof(int), fout2);
    fwrite(&nx, 1,sizeof(int), fout2);
    int tmp_nz = 1;
    fwrite(&tmp_nz, 1,sizeof(int), fout2);
    int tmp_type = 99; //pcap
    fwrite(&tmp_type, 1,sizeof(int), fout2);
    fwrite(&RES2, 1, sizeof(int), fout2);

    fseek(fout2, 1024, SEEK_SET);

    fwrite(tmp,nx*ny,sizeof(float), fout2);
    fclose(fout2);
    free(tmp);
  }

  // Mask
  if (1) {
    FILE *fin = fopen(src_file2.c_str(), "rt");
    if (fin) {
      int original_x, original_y;
      fread(&original_y, 1, sizeof(int), fin);
      fread(&original_x, 1, sizeof(int), fin);
      short *tmp= (short*)malloc(sizeof(short)*original_x*original_y);
      
      fseek(fin, 1024, SEEK_SET);
      fread(tmp, sizeof(short), original_x*original_y, fin);
      fclose(fin);
      
      short *tmp_out = (short*)malloc(sizeof(short)*nx*ny);
      int pointer=0;
      for (int i = 0; i < ny; i++) {
	for (int j = 0; j < nx; j++) {
	  tmp_out[pointer++] = tmp[(j*RES2)+(i*RES2)*original_x];
	}
      }
      
      FILE *fout = fopen(target_file2.c_str(), "wt");
      fwrite(&ny, 1, sizeof(int), fout);
      fwrite(&nx, 1, sizeof(int), fout);
      int tmp_nz = 1;
      fwrite(&tmp_nz, 1, sizeof(int), fout);
      int tmp_type = 98; //pcap mask
      fwrite(&tmp_type, 1, sizeof(int), fout);
      
      fseek(fout, 1024, SEEK_SET);
      fwrite(tmp_out, nx*ny, sizeof(short), fout);
      fclose(fout);
      
      free(tmp_out);
      free(tmp);
    }
  }
}


// LoadImage
GdkPixbuf* LoadImage(char* filename) {
#ifdef DEBUG_MODE
  g_print("Loading %s\n", filename);
#endif

  float *tmp_img;

  // TODO: handle RES2 properly
  int RES2 = 8;
  MAXVAL = -1;
  string tmp = filename;
  if (strcmp(tmp.substr(tmp.size()-4, 4).c_str(), ".mrc") == 0) {
    //tmp_img = LoadMRCorTEXTONSVM(filename);
    tmp_img = LoadMRC_new(filename);
    RES2 = 4;
  }
  else if (strcmp(tmp.substr(tmp.size()-5, 5).c_str(), ".tsvm") == 0) {
    tmp_img = LoadMRCorTEXTONSVM(filename);
    RES2 = 1;
  }
  else {
    tmp_img = LoadTIFF(filename);
    RES2 = 8;
  }
  RES = RES2;
  original_x = nx;
  original_y = ny;

  // TODO: handle RES properly

  if (img != NULL) {
    free(img);
  }
  if (RES2 != 1) {
    img = imresize(tmp_img, nx, ny, 1.0/RES2);
    free(tmp_img);
  } else {
    img = tmp_img;
  }
  g_print("%d %d %d %d\n", nx, ny, original_x, original_y);

  // Load (or initiazlie dontcare)
  // TODO: load dontcare
  if (dontcare != NULL) {
    free(dontcare);
  } 
  dontcare = (int*)malloc(sizeof(int)*nx*ny);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      SetDontCare(i,j,0);
    }
  }

  IMG_FLIPPED = false;

  GdkPixbuf* anim;
  if (!IMG_FLIPPED) {
    anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, nx, ny);
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
      for (j = 0; j < height; j++) {
	pixels[j * stride + i * n_chans] = (guchar)((GetVal(i,j)-MINVAL) * 255 / (MAXVAL-MINVAL));
	pixels[j * stride + i * n_chans+1] = pixels[j*stride+i*n_chans];
	pixels[j * stride + i * n_chans+2] = pixels[j*stride+i*n_chans];
      }
    }
  }
  else {
    anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, ny, nx);
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
      for (j = 0; j < height; j++) {
	pixels[j * stride + i * n_chans] = (guchar)((GetVal(j,i)-MINVAL) * 255 / (MAXVAL-MINVAL));
	pixels[j * stride + i * n_chans+1] = pixels[j*stride+i*n_chans];
	pixels[j * stride + i * n_chans+2] = pixels[j*stride+i*n_chans];
      }
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
  if (path.size() > 0) {
    int dot_p = path.find_last_of(".");
    int slash_p = path.find_last_of("/");
    if (dot_p > slash_p) {
      path = path.substr(0, dot_p);
    }
  }
  
  // 1. Save particles
  SaveParticles(path);

  // 2. Save img (this seems to be useless at this point)
  if (0) {
    float *tmp = (float*)malloc(sizeof(float)*nx*ny);
    int pointer=0;
    for (int i = 0; i < ny; i++) {
      for (int j = 0; j < nx; j++) {
	tmp[pointer++]= GetVal(j,i);
      }
    }  
    
    FILE *fout2 = fopen((path+".tsvm").c_str(), "wt");
    fwrite(&ny, 1,sizeof(int), fout2);
    fwrite(&nx, 1,sizeof(int), fout2);
    int tmp_nz = 1;
    fwrite(&tmp_nz, 1,sizeof(int), fout2);
    int tmp_type = 99; //pcap
    fwrite(&tmp_type, 1,sizeof(int), fout2);

    fseek(fout2, 1024, SEEK_SET);

    fwrite(tmp,nx*ny,sizeof(float), fout2);
    fclose(fout2);
    free(tmp);
  }

  // 3. Save Mask
  if (1) {
    short *tmp = (short*)malloc(sizeof(short)*original_x*original_y);
    int pointer=0;
    for (int i = 0; i < original_y; i++) {
      for (int j = 0; j < original_x; j++) {
	tmp[pointer]= (short)GetDontCare(j/RES,i/RES);
	tmp[pointer] = !tmp[pointer];
	pointer++;
      }
    }  
    
    FILE *fout2 = fopen((path+".mask").c_str(), "wt");
    fwrite(&original_y, 1,sizeof(int), fout2);
    fwrite(&original_x, 1,sizeof(int), fout2);
    int tmp_nz = 1;
    fwrite(&tmp_nz, 1,sizeof(int), fout2);
    // TODO: set type for mask.. maybe 98?
    int tmp_type = 98; //pcap
    fwrite(&tmp_type, 1,sizeof(int), fout2);    

    fseek(fout2, 1024, SEEK_SET);

    fwrite(tmp,original_x*original_y,sizeof(short), fout2);
    fclose(fout2);
    free(tmp);
  }
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


  GdkPixbuf* anim;
  if (!IMG_FLIPPED) {
    anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, nx, ny);
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
  }
  else {
    anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, ny, nx);
    guchar *pixels = gdk_pixbuf_get_pixels(anim);
    int stride = gdk_pixbuf_get_rowstride(anim);
    int n_chans = gdk_pixbuf_get_n_channels(anim);
    int height = gdk_pixbuf_get_height(anim);
    int width = gdk_pixbuf_get_width(anim);
    
    guchar val;
    for (i = 0; i < width; i++) {
      for (j = 0; j < height; j++) {
	if (NotDontCare(j,i)) {
	  val = (guchar)(((GetVal(j,i) - mean_val) / std_val / 6.0 + 0.5) * 255);
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
  }
    
  return anim;
}




// Flip image
GdkPixbuf* FlipImage(GdkPixbuf* cur_buf) {
  IMG_FLIPPED = !IMG_FLIPPED;

  GdkPixbuf* anim = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, gdk_pixbuf_get_height(cur_buf), gdk_pixbuf_get_width(cur_buf));
  guchar *pixels = gdk_pixbuf_get_pixels(anim);
  int stride = gdk_pixbuf_get_rowstride(anim);
  int n_chans = gdk_pixbuf_get_n_channels(anim);
  int height = gdk_pixbuf_get_height(anim);
  int width = gdk_pixbuf_get_width(anim);

  guchar *pixels_src = gdk_pixbuf_get_pixels(cur_buf);
  int stride_src = gdk_pixbuf_get_rowstride(cur_buf);
  int n_chans_src = gdk_pixbuf_get_n_channels(cur_buf);
  
  guchar val;
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      pixels[j * stride + i * n_chans] = pixels_src[i*stride_src + j*n_chans_src];
      pixels[j * stride + i * n_chans+1] = pixels_src[i*stride_src + j*n_chans_src+1];
      pixels[j * stride + i * n_chans+2] = pixels_src[i*stride_src + j*n_chans_src+2];
    }
  }
  return anim;
}


