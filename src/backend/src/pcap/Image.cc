#include "pcap/Image.hh"
#include <iostream>

#include <cstdlib>
#include <cstdio>

using namespace std;

#define SQR(x) (x * x)
#define DEBUG_MODE 0

Image::Image() {
  img = NULL;
  mask = NULL;
  X = Y = couMask = -1;
  filename = "";
  resize_factor = -1;
}

Image::~Image() {
  if (img)
    delete img;
  if (mask)
    delete mask;
  
#if DEBUG_MODE	
  cout << "Image deconstructor" << endl;
#endif
}

void Image::LoadGoodCenters(string name) {
  FILE* fin = fopen((name+".gbox").c_str(), "rt");
  if (fin == NULL) {
    printf("Fild loading error.\n");
    exit(1);
  }
  if (resize_factor == -1) {
    printf("Resize Factor not set.\n");
    exit(1);
  }
  printf("%d\n", resize_factor);
  //	fscanf(fin, "%d\n", &center_n);
  //	for (int i = 0; i < center_n; i++) {
  int x,y,w,h,garbage;
  while (fscanf(fin, "%d %d %d %d %d\n", &x, &y, &w, &h, &garbage)==5) {
    x = (x+(w/2)) / resize_factor;
    y = (y+(h/2)) / resize_factor;
    //TODO: coordinate
    good_centers.push_back(Center(y-1, x-1));
  }
  /* old format
  while (1) {
    int x, y;
    fscanf(fin, "%d %d\n", &x, &y);
    x = x / resize_factor;
    y = y / resize_factor;
    if ((x == 0) && (y == 0))
      break;
    //TODO: coordinate
    good_centers.push_back(Center(x-1, y-1));
  }
  */
  fclose(fin);
}

void Image::LoadCenters(string name) {
  FILE* fin = fopen((name+".box").c_str(), "rt");
  if (fin == NULL) {
    printf("Fild loading error.\n");
    exit(1);
  }
  if (resize_factor == -1) {
    printf("Resize Factor not set.\n");
    exit(1);
  }
  //	fscanf(fin, "%d\n", &center_n);
  //	for (int i = 0; i < center_n; i++) {
  int x,y,w,h,garbage;
  while (fscanf(fin, "%d %d %d %d %d\n", &x, &y, &w, &h, &garbage)==5) {
    x = (x+(w/2)) / resize_factor;
    y = (y+(h/2)) / resize_factor;
    //TODO: coordinate
    centers.push_back(Center(y-1, x-1));
  }
  /* old format
  while (1) {
    int x, y;
    fscanf(fin, "%d %d\n", &x, &y);
    x = x / resize_factor;
    y = y / resize_factor;
    if ((x == 0) && (y == 0))
      break;
    //TODO: coordinate
    centers.push_back(Center(x-1, y-1));
  }
  */
  fclose(fin);
}

void Image::LoadImage(string name) {
  if (img) {
    free(img);
    img = NULL;
  }	

  FILE* fin = fopen((name+".tsvm").c_str(), "rt");
  if (fin == NULL) {
    printf("File loading error: %s.\n", (name+".tsvm").c_str());
    exit(1);
  }
  int tmp_z, type;
  fread(&X, 1,sizeof(int), fin);
  fread(&Y, 1,sizeof(int), fin);
  fread(&tmp_z, 1,sizeof(int), fin);
  fread(&type, 1,sizeof(int), fin);
  fread(&resize_factor, 1, sizeof(int), fin);
  if (type != 99) {
    printf("File type does not match.\n");
    exit(1);
  }

  float *tmp = (float*)malloc(sizeof(float)*X*Y);
  fseek(fin, 1024, SEEK_SET);
  fread(tmp, X*Y,sizeof(float), fin);

  img = new matrix<>(X,Y);

  int pointer = 0;
  for (int i = 0; i < X; i++) {
    for (int j = 0; j < Y; j++) {
      SetImg(i, j, tmp[pointer]);
      pointer++;
    }
  }

  fclose(fin);
  free(tmp);
}

void Image::LoadMask(string name) {
  if (mask) {
    free(mask);
    mask = NULL;
  }

  FILE* fin = fopen((name+".mask").c_str(), "rt");
  if (!fin) {
    printf("Warning: Mask file does not exist.\n");
    printf("   %s\n", (name+".mask").c_str()); 

    if (X == -1) {
      LoadImage(name);
    }
    mask = new matrix<short>(X,Y);

    for (int i = 0; i < X; i++) {
      for (int j = 0; j < Y; j++) {
	SetMask(i, j, 1);
      }
    }
  } else {
    fread(&X, 1, sizeof(int), fin);
    fread(&Y, 1, sizeof(int), fin);
    
    short *tmp = (short*)malloc(sizeof(short)*X*Y);
    fseek(fin, 1024, SEEK_SET);
    fread(tmp, X*Y, sizeof(short), fin);
    
    mask = new matrix<short>(X,Y);
    
    int pointer = 0;
    for (int i = 0; i < X; i++) {
      for (int j = 0; j < Y; j++) {
	SetMask(i, j, tmp[pointer]);
	pointer++;
      }
    }
    
    fclose(fin);
    free(tmp);
  }
}	

void Image::LoadNormImage(string name) {
  if (img) {
    free(img);
    img = NULL;
  }	
  
  FILE* fin = fopen((name+".norm").c_str(), "rt");
  if (fin == NULL) {
    printf("Fild loading error.\n");
    exit(1);
  }
  int tmp_resize, tmp_X, tmp_Y;
  fread(&tmp_X, 1, sizeof(unsigned long), fin); X = static_cast<int>tmp_X;
  fread(&tmp_Y, 1, sizeof(unsigned long), fin); Y = static_cast<int>tmp_Y;
  fread(&tmp_resize, 1, sizeof(int), fin);

  fseek(fin, 1024, SEEK_SET);
  if (resize_factor == -1) {
    resize_factor = tmp_resize;
  }
  img = new matrix<>(X,Y);
  for (int i = 0; i < X; i++) {
    for (int j = 0; j < Y; j++) {
      double tmp;
      fread(&tmp, 1, sizeof(double), fin);
      SetImg(i, j, tmp);
    }
  }
  /* old format
  int tmp_resize;
  fscanf(fin, "%d %d %d\n", &X, &Y, &tmp_resize);
  if (resize_factor == -1) {
    resize_factor = tmp_resize;
  }
  img = new matrix<>(X, Y);
  for (int i = 0; i < X; i++) {
    for (int j = 0; j < Y; j++) {
      double tmp;
      fscanf(fin, "%lf", &tmp);
      SetImg(i, j, tmp);
    }
  }
  */
  fclose(fin);

  LoadMask(name);
}

void Image::WriteNormImage(string name, matrix<>& image, int resize_factor) {
  FILE* fout = fopen((name+".norm").c_str(), "wt");
  if (fout == NULL) {
    printf("Fild loading error.\n");
    exit(1);
  }
  unsigned long s0 = static_cast<unsigned long>image.size(0), s1 = static_cast<unsigned long>image.size(1);
  fwrite(&s0, 1, sizeof(unsigned long), fout);
  fwrite(&s1, 1, sizeof(unsigned long), fout);
  fwrite(&resize_factor, 1, sizeof(int), fout);
  fseek(fout, 1024, SEEK_SET);
  for (int i = 0; i < (long)image.size(0); i++) {
    for (int j = 0; j < (long)image.size(1); j++) {
      fwrite(&(image(i,j)), 1, sizeof(double), fout);
    }
  }
  /* old format
  fprintf(fin, "%ld %ld %d\n", image.size(0), image.size(1), resize_factor);
  
  for (int i = 0; i < (long)image.size(0); i++) {
    for (int j = 0; j < (long)image.size(1); j++) {
      fprintf(fin, "%.9lf ", image(i,j));
    }
    fprintf(fin, "\n");
  }
  */
  fclose(fout);
}

void Image::Normalize() {
  float avg = 0;
  for (int x = 0; x < X; x++)
    for (int y = 0; y < Y; y++)
      if (GetMask(x, y) == 1)
	avg += GetImg(x, y);		
  avg /= couMask;
  
  float std = 0;
  for (int x = 0; x < X; x++)
    for (int y = 0; y < Y; y++)
      std += SQR(GetImg(x,y) - avg);
  std /= couMask;
  
  for (int x = 0; x < X; x++)
    for (int y = 0; y < Y; y++)
      if (GetMask(x, y) == 1)
	SetImg(x, y, (GetImg(x, y) - avg) / std);
}
