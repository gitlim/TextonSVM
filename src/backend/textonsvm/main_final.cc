// TEXTONSVM project (v1.0)
// UC Berkeley - Computer Vision group

// Author: Joseph J. Lim (lim@csail.mit.edu)
// Date: Mar 7th, 2013

// TODO
// 1. variable globalization (or localization) including SVM_PARAM, DICT
// 2. check if passing SVM_PARAM globally rather than "loading from parameters.mat" is fine. ( i think it's ok )
// 3. mask effect
// 12. Change all temporary .img to binary.. it will speed up some


#include <dirent.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
#include <queue>
#include <ctime>
#include <stdlib.h>

#include <libgen.h>

#include "pcap/libsvm/svm.hh"
#include "pcap/Image.hh"
#include "pcap/Examplar.h"
#include "pcap/safe_vector.h"
#include "lang/string.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/array.hh"
#include "collections/list.hh"
#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "math/matrices/matrix.hh"
#include "mlearning/clustering/clusterers/kmeans/matrix_clusterer.hh"
#include "math/matrices/functors/matrix_distance_functors.hh"
#include "io/formats/image/jpeg.hh"
#include "io/formats/image/png.hh"

#define DETAIL (2)
#define BRIEF (1)
#define OFF (0)
//#define DEBUG_MODE DETAIL
//#define DEBUG_MODE BRIEF
#define DEBUG_MODE OFF

//#define WINDOW_RADIUS (30)
//#define WINDOW_RADIUS (36)
//#define MAX_SAMPLES (10000)
#define MAX_SAMPLES (300000)

// k-means
//#define MAX_ITER (3)
#define MAX_ITER (30)
//#define MAX_ITER (100)
//#define MAX_ITER (5000)

// train_004
//#define TEXTON_NBINS (64)
//#define GRAYHIST_NBINS (64)
//#define NBCLS (64)
#define WEIGHT_FEATURE (0.5)

const int dx[] = {-1, 0, 1, -1, 1, -1, 0, 1};
const int dy[] = {-1, -1, -1, 0, 0, 1, 1, 1};

string svm_model_name;
string RESULT_DIR;
string IMG_DIR;
string OUTPUT_DIR;

bool TEST_ONLY;
bool SEARCH_OPTIMAL;
int WINDOW_RADIUS = 30;
int TEXTON_RADIUS = 4;
int K = 1000;

int TEXTON_NBINS = 64;
int GRAYHIST_NBINS = 64;
int NBCLS = 64;
//double WEIGHT_FEATURE = 0.5;

int rescale_factor = -1;
int default_particle_radius = -1;
int default_window_radius = -1;
int default_thresh_Hmax = 50;

char* DICT_NAME;

FILE *random_seed_F;

using namespace std;

using lang::array;
using collections::list;
using collections::array_list;
using collections::pointers::auto_collection;
using math::matrices::matrix;
using mlearning::clustering::clusterers::kmeans::matrix_clusterer;
using math::matrices::functors::matrix_L2_distance;
using io::formats::image::jpeg;
using io::formats::image::png;
using lang::pointers::auto_ptr;

time_t t0, t1;

inline int c3to1(const int a, const int b, const int c, const int mb, const int mc) {
  if (b >=  mb) {
#if (DEBUG_MODE == DETAIL)
      cout << "ERROR" << endl;
#endif
  }
  if (c >= mc) {
#if (DEBUG_MODE == DETAIL)
      cout << "ERROR" << endl;
#endif
  }

  return a*mb*mc + b*mc + c; 
}

void tic() {
  t0 = time(NULL);
}

void toc() {
  t1 = time(NULL);
  printf("[[TIME: %ld seconds ellapsed]]\n", (t1 - t0));
}

bool exist_file(string name) {
  ifstream inp;
  inp.open(name.c_str(), ifstream::in);
  inp.close();
  return !inp.fail();
}

void touchFile(string dir, string name) {
  FILE* fout;
#if (DEBUG_MODE == DETAIL)
    cout << dir+name << endl;
#endif
  fout = fopen((dir+name).c_str(), "wt");
  fclose(fout);
}

int min_f(int a, int b) {
  if (a > b) {
    return b;
  }
  return a;
}
double min_f(double a, double b) {
  if (a > b)
    return b;
  return a;
}

void normalize_img(const vector<string>& image_names,
		   const int window_radius,
		   const int center_offset = 0) {
  for (unsigned int i = 0; i < image_names.size(); i++) {
#if (DEBUG_MODE)
      cout << "Image " << (i+1) << ": " << image_names[i] << endl;
#endif

    if (exist_file(OUTPUT_DIR + image_names[i] + ".norm"))
      continue;

    Image image;
    image.LoadImage((IMG_DIR+image_names[i]));

    const int w = image.GetWidth(), h = image.GetHeight();
    const double *img_ptr = &((*image.get_img_pointer())(0,0));
    //    cout << image.GetImg(100, 500) << " " << *(img_ptr+ (100*h + 500)) << " " << *(img_ptr + (100*w + 500)) << endl;
    
    matrix<> patch((2*window_radius+1)*(2*window_radius+1), 1);
    double *patch_ptr = &(patch(0,0));

    matrix<> image_norm(w, h);

#if (DEBUG_MODE)
    cout << "Normalizing Img: [";
#endif
    for (int ii = 0+window_radius; ii < w - window_radius; ii++) {
      if (ii % 100 == 0) {
#if (DEBUG_MODE == DETAIL)
	cout << ii * 100.0 / w << endl;
#endif
      }
      for (int jj = 0+window_radius; jj < h - window_radius; jj++) {
	int cou = 0;
	for (int x = ii-window_radius; x <= ii+window_radius; x++) {
	  for (int y = jj-window_radius; y <= jj+window_radius; y++) {
	    *(patch_ptr + (cou++)) = *(img_ptr + x*h + y);
	  }
	}

	const double patch_mean = mean(patch);
	double patch_std = 0;
	for (int kk = 0; kk < cou; kk++) {
	  const double tmp = *(patch_ptr+kk) - patch_mean;
	  patch_std += tmp*tmp; 
	}
	patch_std = sqrt(patch_std / (cou - 1));

	if (patch_std < 0.00001) {
	  image_norm(ii,jj) = *(img_ptr + ii*h + jj) - patch_mean;
	}
	else {
	  image_norm(ii,jj) = (*(img_ptr + ii*h + jj) - patch_mean) / patch_std;
	}
      }
#if (DEBUG_MODE) 
      cout << ".";
#endif
    }
#if (DEBUG_MODE)
    cout << "]" << endl;
#endif

    Image::WriteNormImage(OUTPUT_DIR+image_names[i], image_norm, image.GetResizeFactor());

#if (DEBUG_MODE)
    toc();
#endif
  }
}

void collect_exemplars(const int window_radius,
		       const vector<string>& image_names,
		       safe_vector<matrix<>*>& examplars) {
  if (exist_file(RESULT_DIR+"examplars.denoise"))
    return ;

  for (unsigned int j = 0; j < image_names.size(); j++) {
#if (DEBUG_MODE)
      cout << "Image " << (j+1) << ": " << image_names[j] << endl;
#endif

    Image image;
    image.LoadNormImage(OUTPUT_DIR+image_names[j]);
    image.LoadGoodCenters(IMG_DIR+image_names[j]);

    // eroding
    for (int x = 0; x < image.GetWidth(); x++) {
      for (int y = 0; y < image.GetHeight(); y++) {
	//	if ((image.GetMask(x, y) == 0) || (x == 0) || (y == 0) || (x == image.GetWidth()-1) || (y == image.GetHeight()-1)) {
	if ((image.GetMask(x, y) == 0) || (x <= window_radius) || (y <= window_radius) || (x >= image.GetWidth()-1-window_radius) || (y >= image.GetHeight()-1-window_radius)) {
	  for (int nx = x - window_radius; nx <= x + window_radius; nx++) {
	    for (int ny = y - window_radius; ny <= y + window_radius; ny++) {
	      if ((nx < 0) || (nx >= image.GetWidth()) || (ny < 0) || (ny >= image.GetHeight()))
		continue;
	      
	      image.SetMask(nx, ny, -1);
	    }
	  }
	}
      }
    }
    for (int x = 0; x < image.GetWidth(); x++)
      for (int y = 0; y < image.GetHeight(); y++) {
	image.SetMask(x, y, image.GetMask(x, y) > 0);
      }


    const vector<Center>& centers = image.GetGoodCenters();
#if (DEBUG_MODE == DETAIL)
      cout << centers.size() << endl;
#endif

#if DEBUG_MODE
    printf("image %d is loaded\n", j+1);
#endif
    for (unsigned int i = 0; i < centers.size(); i++) {
      int xx = centers[i].GetX();
      int yy = centers[i].GetY();

      if (!image.GetMask(xx, yy))
	continue;

#if (DEBUG_MODE == DETAIL)
      cout << xx << " " << yy << endl;
#endif


      //      if ((xx - (window_radius*2) < 0+198) || (yy - (window_radius*2) < 0) ||
      //	  (xx + (window_radius*2+1) >= image.GetWidth()) || (yy + (window_radius*2+1) >= image.GetHeight()))
      //	continue;
      
      examplars.push_back(new matrix<>(window_radius*2+1, window_radius*2+1));
      for (int x = -window_radius; x <= window_radius; x++) {
	for (int y = -window_radius; y <= window_radius ; y++) {
	  (*examplars.back())(x+window_radius, y+window_radius) = image.GetImg(x+xx, y+yy);
	}
      }
    }
  }
}

int partition(double a[], int b[], int p, int r) {
  double x = a[r];
  int y = b[r];

  int j = p - 1;

  for (int i = p; i < r; i++) {
    if (x > a[i]) {
      j = j + 1;
      double temp = a[j];
      a[j] = a[i];
      a[i] = temp;

      int temp2 = b[j];
      b[j] = b[i];
      b[i] = temp2;
    }
  }

  a[r] = a[j + 1];
  a[j+1] = x;

  b[r] = b[j + 1];
  b[j+1] = y;

  return j + 1;
}

void quicksort(double a[], int b[], int left, int right) {
  if (left < right) {
    int pivot = partition(a, b, left, right);
    quicksort(a, b, left, pivot-1);
    quicksort(a, b, pivot+1, right);
  }
}

void randperm(int n, int* rp) {
  double tmp[n];

  for (int i = 0; i < n; i++) {
    if ((!random_seed_F) || (feof(random_seed_F))) {
      tmp[i] = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
    } else {
      fread(&(tmp[i]), sizeof(double), 1, random_seed_F);
    }
    rp[i] = i;

    //    cout << tmp[i] << " ";
  }
  //  cout << endl<<endl<<endl;

  quicksort(tmp, rp, 0, n-1);
#if (DEBUG_MODE==DETAIL)
  for (int i = 0; i < n; i ++) {
    cout << rp[i] << " ";
  }
  cout << endl << endl;
#endif
  //  exit(1);

  /*
  for (int i = 0; i < n; i++)
    rp[i] = i;  

  for (int i = n-1; i>0; i--) {
    int j = rand() % (i+1); 
    int tmp = rp[i];
    rp[i] = rp[j];
    rp[j] = tmp;
  }
  */

  /*
  int tmp[n];
  for (int i = 0; i < n; i++) {
    tmp[i] = rand();
    rp[i] = i;
  }
  
  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      if (tmp[i] > tmp[j]) {
	int swap = tmp[i];
	tmp[i] = tmp[j];
	tmp[j] = swap;
	
	swap = rp[i];
	rp[i] = rp[j];
	rp[j] = swap;
      }
    }
  }
  */
}

void computeMembership(const array_list<matrix<> > data, const array_list<matrix<> > means, int *membership, double& rms) {
  const int n = data.size();
  const int d = data[0].size();

  rms = 0;

  const int k = means.size();

  double *mean_ptr[k];
  for (int i = 0; i < k; i++)
    mean_ptr[i] = &(means[i](0,0));

  for (int i = 0; i < n; i++) {
    //    const matrix<>& datum = data[i];
    const double *datum = &(data[i](0,0));

    double min_dist = 99999999;
    int min_ind = -1;
    
    for (int j = 0; j < k; j++) {
      const double *mean = mean_ptr[j];
      //      const double *mean = &(means[j](0,0));
      //	const matrix<>& mean = means[j];

      double dist = 0;
      for (int l = 0; l < d; l++) {
	const double tmp = *(datum+l) - *(mean+l);
	dist += tmp * tmp;
	//	if (dist > min_dist) {
	//	  break;
	//	}
      }      
      if (min_dist > dist) {
	min_dist = dist;
	min_ind = j;
      }
    }
	   
    membership[i] = min_ind;
    rms = rms + min_dist;

  }
  rms = sqrt(rms / n);
}
  
void computeMeans(const int k, const array_list<matrix<> > data, const int* membership, auto_collection<matrix<>, array_list<matrix<> > >& means) {
  int n = data.size();
  int d = data[0].size();

  for (int i = 0; i < k; i++)
    for (int j = 0; j < d; j++)
      (*means)[i](j,0) = 0;

  int cou[k];
  for (int i = 0; i < k; i++)
    cou[i] = 0;
  for (int j = 0; j < n; j++) {
    for (int l = 0; l < d; l++) {
      (*means)[membership[j]](l,0) += data[j](l,0);
    }
    cou[membership[j]]++;
  }

  for (int i = 0; i < k; i++) {
    if (cou[i] == 0) {
      continue;
    }
    for (int l = 0; l < d; l++) {      
      (*means)[i](l,0) /= cou[i];
    }
  }
}

void kmeansInternal(const int k, const array_list<matrix<> > data, const int maxIter, const double dtol, const double etol, const bool ml, const bool verbose, const bool retry, auto_collection<matrix<>, array_list<matrix<> > >& means) {
  int n = data.size();
  int d = data[0].size();
  if (verbose > 0) {
#if (DEBUG_MODE==DETAIL)
    cout << n << "  " << d << endl;
#endif
  }
  
  int rp[n];
  randperm(n, rp);

  // compute initial means (line 57)
  int rate = 3;
  int minN = 50;
  int coarseN = n / rate;

  if ((!ml) || (coarseN < k) || (coarseN < minN)) {    
    // pick random points as means
    means.reset(new array_list<matrix<> >());
    for (int i = 0; i < k; i++) {
      means->add(*(new matrix<>(d, 1)));
      for (int j = 0; j < d; j++) {
	(*means).tail()(j,0) = data[rp[i]](j,0);
      }
    }
  }
  else {
    auto_collection< matrix<>, array_list<matrix<> > > coarse_data(new array_list<matrix<> >());
    for (int i = 0; i < coarseN; i++) {
      coarse_data->add(*(new matrix<>(data[rp[i]])));
    }
    kmeansInternal(k, *coarse_data, maxIter, dtol, etol, ml, verbose, 0, means);
  }

  // Iterate (line 71)
  int iter = 0;
  double rms = 99999999;
  double rmsPrev;
  double rmsPctChange;
  double maxMoved;
  int membership[n];

  if (verbose > 0) {
#if (DEBUG_MODE)
    printf("kmeansML: n=%d d=%d k=%d [", n, d, k);
#endif
  }
  auto_collection<matrix<>, array_list<matrix<> > > prev_means(new array_list<matrix<> >()); 
  for (int i = 0; i < k; i++) {
    prev_means->add(*(new matrix<>(d, 1)));    
  }
  while (iter < maxIter) {
    if (verbose > 0) {
#if (DEBUG_MODE)
      cout << ".";
#endif
    }
    iter++;

    rmsPrev = rms;
    computeMembership(data, *means, membership, rms);
    for (int i = 0; i < k; i++) {
      //double* prev_mean_ptr = &((*prev_means)[i](0,0));
      //      const double* mean_ptr = &((*means)[i](0,0));
      for (int j = 0; j < d; j++) 
	(*prev_means)[i](j,0) = (*means)[i](j,0);
	//*(prev_mean_ptr+j) = *(mean_ptr+j);
    }
    computeMeans(k, data, membership, means);

    if (rms > rmsPrev) {
#if (DEBUG_MODE==DETAIL)
      cout << "ERROR" << endl;
#endif
      break;
    }

    rmsPctChange = 2*(rmsPrev - rms)/(rmsPrev +rms+ 0.000001);
    double max_sum = 0;
    for (int i = 0; i < k; i++) {
      double sum_up = 0;
      for (int j = 0; j < d; j++) {
	const double tmp = (*means)[i](j,0) - (*prev_means)[i](j,0);
	sum_up += tmp * tmp;
      }
      if (sum_up > max_sum)
	max_sum = sum_up;
    }
    maxMoved = sqrt(max_sum);

#if (DEBUG_MODE==DETAIL)
    cout << rmsPctChange << "   " << maxMoved << endl;
#endif
    if ((rmsPctChange <= etol) && (maxMoved <= dtol))
      break;
  }

  if (verbose > 0) {
    printf("] rms=%.5g\n", rms);
  }
}

void kmeansML(const int k, const int max_iter, const array_list<matrix<> > feature, auto_collection<matrix<>, array_list<matrix<> > >& textons) {
  double dtol = 0;
  double etol = 0;
  bool ml = true;
  
  kmeansInternal(k, feature, max_iter, dtol, etol, ml, 1, 1, textons);
}

void kmeans_denoise_textons(const safe_vector<matrix<>*>& examplars,
			    const int texton_radius,
			    const int k,
			    const int max_sample,
			    auto_collection<matrix<>, array_list<matrix<> > >& textons) {
  int tx, ty, n = examplars.size();
  tx = ty = examplars.at(0)->size(0);

  int max_sample_ex = max_sample/n;
  if (max_sample_ex > tx*ty)
    max_sample_ex = tx*ty;
  
  int rp[tx*ty];
  randperm(tx*ty, rp);

  bool mask[tx][ty];
  for (int i = 0; i < tx; i++)
    for (int j = 0; j < ty; j++)
      mask[i][j] = false;
  for (int j = 0; j < max_sample_ex; j++) {
    mask[rp[j]%tx][rp[j]/tx] = true;
  }

  //  int x, y;
  int vector_n = (2*texton_radius+1)*(2*texton_radius+1);

#if (DEBUG_MODE == DETAIL)
  cout << examplars.size() << "   " << max_sample_ex << endl;
#endif

  auto_collection< matrix<>, array_list<matrix<> > >
    feature(new array_list<matrix<> >());
  for (int i = 0; i < examplars.size(); i++) {
    for (int x = 0; x < tx; x++) {
      for (int y = 0; y < ty; y++) {
	if (!mask[x][y])
	  continue;

      //    for (int j = 0; j < max_sample_ex; j++) {
      //      x = rp[j] % tx;
      //      y = rp[j] / tx;
      
        //TODO: boundary condition
	if ((x - texton_radius < 0) || (y - texton_radius < 0) ||
	    (x + texton_radius >= tx) || (y + texton_radius >= ty))
	  continue;
	
	feature->add(*(new matrix<>(vector_n, 1)));
	
	int cou = 0;
	for (int ii = -texton_radius; ii <= texton_radius; ii++)
	  for (int jj = -texton_radius; jj <= texton_radius; jj++) {
	    (*feature).tail()(cou++, 0) = (*examplars.at(i))(x+ii, y+jj);
	  }
      }
    }
    //    cout <<(*feature).tail()(0,0) << endl;
  }

#if (DEBUG_MODE)
  cout << "Starting Kmeans..." << endl;
#endif
  kmeansML(k, MAX_ITER, *feature, textons);

  /*
  matrix_clusterer<double> clusterer(k, MAX_ITER);
  array<unsigned long> assign = clusterer.cluster(*feature, textons);
  cout << assign.size() << endl;
  */

  // Store textons
  FILE *fout = fopen((RESULT_DIR+"textons.dat").c_str(), "wt");
  fprintf(fout, "%lu\n", textons->size());
  fprintf(fout, "%d\n", vector_n);
  for (unsigned int i = 0; i < textons->size(); i++) {
    for (int j = 0; j < vector_n; j++) {
      fprintf(fout, "%.5f ", (*textons)[i](j,0));
    }
    fprintf(fout, "\n"); 
  }
  fclose(fout);
}

bool load_textons(const string& output_dir,
		  const int texton_radius,
		  const int k,
		  auto_collection<matrix<>, array_list<matrix<> > >& textons) {
  if (exist_file(output_dir+"textons.dat")) {
    FILE *fin = fopen((output_dir+"textons.dat").c_str(), "rt");

    int cluster_k, vector_n;
    fscanf(fin, "%d\n", &cluster_k);
    fscanf(fin, "%d\n", &vector_n);
    if ((cluster_k != k) || ((texton_radius*2+1)*(texton_radius*2+1) != vector_n)) {
      fclose(fin);
      return false;
    }
    textons.reset(new array_list<matrix<> >());
    double tmp;
    for (int i = 0; i < cluster_k; i++) {
      textons->add(*(new matrix<>(vector_n, 1)));
      for (int j = 0; j < vector_n; j++) {
	fscanf(fin, "%lf", &tmp);
	(*textons).tail()(j, 0) = tmp;
      }
    } 
      
    fclose(fin);

    printf("Texton loaded.\n");
    return true;
  }
  return false;
}

void kmeans_denoise_rec(const matrix<>& examplar, 
			const auto_collection<matrix<>, array_list<matrix<> > >& textons,
			const int texton_radius,
			matrix<>& reconstruction) {
  int tx = examplar.size(0), ty = examplar.size(1);
  int vector_n = (2*texton_radius+1)*(2*texton_radius+1);
  int x, y;

  matrix_L2_distance<double> mat_dist;
  matrix<> patch(vector_n, 1);

  const int l_size = (texton_radius*2+1)*(texton_radius*2+1);
  const double *examplar_ptr = &(examplar(0,0));
  double *patch_ptr = &(patch(0,0));

#if (DEBUG_MODE)
  cout << "[";
#endif  
  for (unsigned int i = 0; i < examplar.size(); i++) {
    x = i % tx;
    y = i / tx;
#if (DEBUG_MODE==DETAIL)
    if (i % 50000 == 49999) {
      toc();
      cout << i*100.0/examplar.size() << endl;
    }
#endif

    //TODO: MASK
    //TODO: boundary condition
    if ((x - texton_radius < 0) || (y - texton_radius < 0) ||
	(x + texton_radius >= tx) || (y + texton_radius >= ty))
      continue;
    
    int cou = 0;
    for (int ii = -texton_radius; ii <= texton_radius; ii++)
      for (int jj = -texton_radius; jj <= texton_radius; jj++) {
	//	if (*(examplar_ptr + (x+ii)*ty + y+jj) - examplar(x+ii, y+jj) != 0) {
	//	  exit(1);
	//	}
	  
	*(patch_ptr + (cou++)) = *(examplar_ptr + (x+ii)*ty + y+jj);
	//	patch(cou++, 0) = examplar(x+ii, y+jj);
      }

    int min_id = 0;
    double min_D = 999999;
    for (unsigned int j = 0; j < textons->size(); j++) {
      const double* texton_ptr = &((*textons)[j](0,0));
      double tmp = 0;
      for (int l = 0; l < l_size; l++) {
	const double ddd = *(patch_ptr + l) - *(texton_ptr+l);
	tmp += ddd * ddd;
	if (tmp > min_D)
	  break;
      }

      //      double tmp = mat_dist(patch, (*textons)[j]);
      //      cout << sqrt(tmp) << " " << mat_dist(patch, (*textons)[j]) << endl;
      if (tmp < min_D) {
	min_D = tmp;
	min_id = j;
      }
    }
    
    reconstruction(x,y) = (*textons)[min_id](vector_n/2, 0);
  }
#if (DEBUG_MODE)
  cout << "]";
#endif
}

bool load_examplars(const string& output_dir,
		    const string& filename,
		    const int texton_radius,
		    const int k,
		    const int window_radius,
		    safe_vector<matrix<>*>& examplars_denoise) {
  if (exist_file(output_dir + filename)) {
    FILE *fin = fopen((output_dir + filename).c_str(), "rt");

    int cluster_k, stored_t_r, stored_w_r;
    fscanf(fin, "%d\n", &cluster_k);
    fscanf(fin, "%d\n", &stored_t_r);
    fscanf(fin, "%d\n", &stored_w_r);
    if ((cluster_k != k) || 
	(stored_t_r != texton_radius) ||
	(window_radius*2+1 != stored_w_r)) {
      fclose(fin);
      return false;
    }
    int n;
    fscanf(fin, "%d\n", &n);
    double tmp;
    for (int i = 0; i < n; i++) {
      examplars_denoise.push_back(new matrix<>(window_radius*2+1,window_radius*2+1));
      for (int j = 0; j < window_radius*2+1; j++)
	for (int k = 0; k < window_radius*2+1; k++) {
	  fscanf(fin, "%lf", &tmp); 
	  (*examplars_denoise.at(i))(j,k) = tmp;
	}
    }
    fclose(fin);

    printf("Denoised Examplars loaded.\n");    
    return true;
  }
  return false;
}

void save_examplars(const string& output_dir,
		    const string& filename,
		    const int texton_radius,
		    const int k,
		    safe_vector<matrix<>*>& examplars_denoise) {
  FILE *fout = fopen((output_dir + filename).c_str(), "wt");
  // TODO: remove
  //"denoised_examplars.dat").c_str(), "wt");

  fprintf(fout, "%d\n", k);
  fprintf(fout, "%d\n", texton_radius);
  fprintf(fout, "%lu\n", examplars_denoise.at(0)->size(0));
  fprintf(fout, "%d\n", examplars_denoise.size());
  for (int i = 0; i < examplars_denoise.size(); i++) {
    matrix<>* tmp = examplars_denoise.at(i);
    for (unsigned int j = 0; j < tmp->size(0); j++) {
      for (unsigned int k = 0; k < tmp->size(1); k++) {
	fprintf(fout, "%.5f ", (*tmp)(j, k));
      }
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
}

/*
void save_denoised_examplars(const string& output_dir,
			     const int texton_radius,
			     const int k,
			     safe_vector<matrix<>*>& examplars_denoise) {
  FILE *fout = fopen((output_dir + "denoised_examplars.dat").c_str(), "wt");

  fprintf(fout, "%d\n", k);
  fprintf(fout, "%d\n", texton_radius);
  fprintf(fout, "%lu\n", examplars_denoise.at(0)->size(0));
  fprintf(fout, "%d\n", examplars_denoise.size());
  for (int i = 0; i < examplars_denoise.size(); i++) {
    matrix<>* tmp = examplars_denoise.at(i);
    for (unsigned int j = 0; j < tmp->size(0); j++) {
      for (unsigned int k = 0; k < tmp->size(1); k++) {
	fprintf(fout, "%.5f ", (*tmp)(j, k));
      }
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
}
*/

void save_image(const matrix<>& img, const string& work_dir, const string& filename, const array<int>& extra_info = array<int>()) {
#if (DEBUG_MODE == DETAIL)
  cout << "save_image_double 1" << endl;
#endif
  FILE *fout = fopen((work_dir+filename).c_str(), "wt");

#if (DEBUG_MODE == DETAIL)
  cout << img.size(0) << " " << img.size(1) << endl;
#endif

  fprintf(fout, "%lu\n", img.size(0));
  fprintf(fout, "%lu\n", img.size(1));
  for (unsigned int ii = 0; ii < img.size(0); ii++) {
    for (unsigned int jj = 0; jj < img.size(1); jj++)
      fprintf(fout, "%.5f ", img(ii,jj));
    fprintf(fout, "\n");
  }
  fprintf(fout, "%lu\n", extra_info.size());
  for (unsigned int i = 0; i < extra_info.size(); i++)
    fprintf(fout, "%d\n", extra_info[i]);

  fclose(fout);
#if (DEBUG_MODE==DETAIL)
  cout << "save_image 0" << endl;
#endif
}

void save_image(const matrix<int>& img, const string& work_dir, const string& filename, const array<int>& extra_info = array<int>()) {
#if (DEBUG_MODE==DETAIL)
  cout << "save_image 1" << endl;
#endif
  FILE *fout = fopen((work_dir+filename).c_str(), "wt");

  fprintf(fout, "%lu\n", img.size(0));
  fprintf(fout, "%lu\n", img.size(1));
  for (unsigned int ii = 0; ii < img.size(0); ii++) {
    for (unsigned int jj = 0; jj < img.size(1); jj++)
      fprintf(fout, "%d ", img(ii,jj));
    fprintf(fout, "\n");
  }
  fprintf(fout, "%lu\n", extra_info.size());
  for (unsigned int i = 0; i < extra_info.size(); i++)
    fprintf(fout, "%d\n", extra_info[i]);

  fclose(fout);
#if (DEBUG_MODE==DETAIL)
  cout << "save_image 0" << endl;
#endif
}

void denoise_image_assign(const vector<string>& image_names,
			  const auto_collection<matrix<>, array_list<matrix<> > >& textons,
			  const int texton_radius) {
  for (unsigned int i = 0; i < image_names.size(); i++) {
#if DEBUG_MODE
    cout << "Image " << (i+1) << ": " << image_names[i] << endl;
#endif

    if (exist_file(RESULT_DIR + image_names[i] + "_denoise.img"))
	continue;

    Image image;
    image.LoadNormImage(OUTPUT_DIR+image_names[i]);

    matrix<> image_denoise(image.GetWidth(), image.GetHeight());
    kmeans_denoise_rec(image.GetImage(), textons, texton_radius,
		       image_denoise);

    save_image(image_denoise, RESULT_DIR, image_names[i]+"_denoise.img");

    toc();
  }
}

array<int> load_image(const string& image_name,
		      const string& work_dir,
		      matrix<>& image_denoise,
		      bool cut_range) {
  FILE *fin = fopen((work_dir+image_name).c_str(), "rt");

  int width, height;
  double tmp;
  fscanf(fin, "%d\n%d\n", &width, &height);
#if (DEBUG_MODE == DETAIL)
  cout << width << " " << height << endl;
#endif
  image_denoise.resize(static_cast<long unsigned>(width), static_cast<long unsigned>(height));
  for (unsigned int ii = 0; ii < image_denoise.size(0); ii++) {
    for (unsigned int jj = 0; jj < image_denoise.size(1); jj++) {
      fscanf(fin, "%lf ", &tmp);
      if (cut_range) {
	if (tmp <= -2.0)
	  tmp = -2.0;
	else if (tmp >= 2.0)
	  tmp = 2.0;
      }
      image_denoise(ii,jj) = tmp;
    }
    fscanf(fin, "\n");
  }

  int n;
  fscanf(fin, "%d\n", &n);
  array<int> extra_info(n);
  for (int i = 0; i < n; i++) {
    fscanf(fin, "%d\n", &(extra_info[i]));
  }
  fclose(fin);

#if (DEBUG_MODE == DETAIL)
  cout << "load_image: done." << endl;
#endif

  return extra_info;
}

array<int> load_denoise_image(const string& image_name,
			      matrix<>& image_denoise) {
  return load_image(image_name+"_denoise.img", "", image_denoise, true /*cut_range*/);
}

array<int> load_denoise_image(const string& image_name,
			      const string& work_dir,
			      matrix<>& image_denoise) {
  return load_image(image_name+"_denoise.img", work_dir, image_denoise, true /*cut_range*/);
}

class SVM_PARAM_ABSTRACT {
public:
  int texton_radius; 
  int texton_nbins;
  int param_K;
  int nbcls;
  int clusterMethod;

  int grayhist_nbins;
  int particle_radius;
  int window_radius;

  double weight_feature;
  int disType;
  int winType;
  int SVMopts;

  int thresh_Hmax;

  string imgName;
  string dataFile;
  string modelsFile;

  int tx_out;
  int ty_out;
  int overlap;

  int x_init;
  int y_init;

  string batchDir;
  string piecesDir;
  string pieces_detDir;
  string pieces_origDir;
  string pieces_clsDir;

  void print() {
#if (DEBUG_MODE==DETAIL)
    cout << texton_radius << " " << texton_nbins << " " << param_K << " " << nbcls << " " << clusterMethod << " " << grayhist_nbins << " " << particle_radius << " "
	 << window_radius << " " << weight_feature << " " << disType << " " << winType << " " << SVMopts << " " << thresh_Hmax << " " << dataFile << " " << modelsFile << " "
	 << tx_out << " " << ty_out << " " << overlap << " " << x_init << " " << y_init << " " << batchDir << " " << piecesDir << endl; 
#endif
  }
};

class SVM_PARAM : public SVM_PARAM_ABSTRACT {
public:
  SVM_PARAM(int particles_radius,
	    int window_radius,
	    int thresh_Hmax) {
    texton_radius = TEXTON_RADIUS;
    texton_nbins = TEXTON_NBINS;
    param_K = K;
    nbcls = NBCLS;
    // clusterMethod = 'average'

    grayhist_nbins = GRAYHIST_NBINS;

    weight_feature = WEIGHT_FEATURE;

    this->particle_radius = particles_radius;
    this->window_radius = window_radius;

    // disk = 0, square = 1
    winType = 1;

    this->thresh_Hmax = thresh_Hmax;

    //    dataFile = "006_Fer";
    dataFile = DICT_NAME;
    modelsFile = "models.dat";
    imgName = DICT_NAME; 
    //    imgName = "006_Fer";

    tx_out = 500;
    ty_out = tx_out;
    overlap = 2*(window_radius+texton_radius)+1;
    
    x_init = 1;
    y_init = 1;

    batchDir = RESULT_DIR + "batch/";
    mkdir(batchDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

    piecesDir = RESULT_DIR + "pieces/";
    mkdir(piecesDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

    pieces_detDir = piecesDir + "det/";
    mkdir(pieces_detDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

    pieces_origDir = piecesDir + "orig/";
    mkdir(pieces_origDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

    pieces_clsDir = piecesDir + "cls/";
    mkdir(pieces_clsDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
  }
};

class SVM_PARAM_ALL : public SVM_PARAM_ABSTRACT {
public:
  SVM_PARAM_ALL(int texton_radius,
		int texton_nbins,
		int grayhist_nbins,
		double weight_feature,
		int particle_radius,
		int window_radius,
		int param_K,
		int nbcls,
		int thresh_Hmax,
		string imgName) {
    this->texton_radius = texton_radius;
    this->texton_nbins = texton_nbins;
    this->grayhist_nbins = grayhist_nbins;
    this->weight_feature = weight_feature;
    this->particle_radius = particle_radius;
    this->window_radius = window_radius;
    this->param_K = param_K;
    this->nbcls = nbcls;
    this->thresh_Hmax = thresh_Hmax;
    this->dataFile = RESULT_DIR + imgName;
    this->imgName = imgName;

    // clusterMethod = 'average'

    // disk = 0, square = 1
    winType = 1;

    modelsFile = "models.dat";

    tx_out = 500;
    ty_out = tx_out;
    overlap = 2*(window_radius+texton_radius)+1;
    
    x_init = 1;
    y_init = 1;

    //batchDir = RESULT_DIR + imgName + "/" + "batch/";
    batchDir = RESULT_DIR + "results/" + "batch/";
    mkdir(batchDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

    //piecesDir = RESULT_DIR + imgName + "/" + "pieces/";
    piecesDir = RESULT_DIR + "results/" + "pieces/";
    mkdir(piecesDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

    pieces_detDir = piecesDir + "det/";
    mkdir(pieces_detDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

    pieces_origDir = piecesDir + "orig/";
    mkdir(pieces_origDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

    pieces_clsDir = piecesDir + "cls/";
    mkdir(pieces_clsDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
  }
};

class DICT {
public:
  double *thresh;
  int thresh_n;
  auto_collection<matrix<>, array_list<matrix<> > > textons;

  DICT() {
    thresh = NULL;
  }

  ~DICT() {
    if (thresh) {
      delete[] thresh;
    }
  }

  void save(string path, string filename) {
#if (DEBUG_MODE)
    cout << "DICT: save starts" << endl;
#endif

    FILE *fout = fopen((path+filename).c_str(), "wt");
    fprintf(fout, "%lu\n", textons->size());
    for (unsigned int i = 0; i < textons->size(); i++) {
      matrix<>& mat = (*textons)[i];
      fprintf(fout, "%lu %lu\n", mat.size(0), mat.size(1));
      for (unsigned int ii = 0; ii < mat.size(0); ii++) {
	for (unsigned int jj = 0; jj < mat.size(1); jj++) {
	  fprintf(fout, "%.5f ", mat(ii, jj));
	} 
	fprintf(fout, "\n");
      }
    } 
    fprintf(fout, "%d\n", thresh_n);
    for (int i = 0; i < thresh_n; i++)
      fprintf(fout, "%.5f ", thresh[i]);

    fclose(fout);

#if (DEBUG_MODE)
    cout << "DICT: save ends" << endl;
#endif
  }

  void load(string path, string filename) {
#if (DEBUG_MODE)
    cout << "DICT: load starts" << endl;
#endif

    FILE *fin = fopen((path+filename).c_str(), "rt");

    int n;
    fscanf(fin, "%d\n", &n);
#if (DEBUG_MODE==DETAIL)
    cout << n << endl;
#endif

    textons.reset(new array_list<matrix<> >());
    for (int i = 0; i < n; i++) {
      int tx, ty;
      fscanf(fin, "%d %d\n", &tx, &ty);

      textons->add(*(new matrix<>(tx, ty)));
      matrix<>& mat = (*textons)[i];

      double tmp;
      for (int ii = 0; ii < tx; ii++)
	for (int jj = 0; jj < ty; jj++) {
	  fscanf(fin, "%lf ", &tmp);
	  if (tmp >= 2.0)
	    tmp = 2.0;
	  if (tmp <= -2.0)
	    tmp = -2.0;
	  mat(ii, jj) = tmp;
	}
    }
    
    fscanf(fin, "%d\n", &thresh_n);
#if (DEBUG_MODE == DETAIL)
    cout << thresh_n << endl;
#endif
    if (thresh)
      delete[] thresh;
    thresh = new double[thresh_n];
    for (int i = 0; i < thresh_n; i++) {
      double tmp;
      fscanf(fin, "%lf", &tmp);
      thresh[i] = tmp;
    }

    fclose(fin);

#if (DEBUG_MODE)
    cout << "DICT: load ends" << endl;
#endif
  }
};

void compute_textons(const auto_collection<matrix<>, array_list<matrix<> > >& textons, const SVM_PARAM& param, DICT& dictionary) {
  float min_textons = 99999;
  float max_textons = -99999;
  for (int i = 0; i < param.param_K; i++) {
    const matrix<>& texton = (*textons)[i];
    for (unsigned int j = 0; j < texton.size(); j++) {
      if (min_textons > texton(j, 0))
	min_textons = texton(j, 0);
      if (max_textons < texton(j, 0))
	max_textons = texton(j, 0);
    }
  }
  /*
  // TODO: is this right? also, should I just cut textons range?
  if (min_textons <= -2.0)
    min_textons = -2.0;
  if (max_textons >= 2.0)
    max_textons = 2.0;
  */

  double step = (max_textons - min_textons) / (param.texton_nbins - 1);
  dictionary.thresh = new double[param.texton_nbins+1];
  dictionary.thresh_n = param.texton_nbins;
  for (int k = 0; k < param.texton_nbins; k++)
    dictionary.thresh[k] = min_textons + step*k;
  dictionary.thresh[param.texton_nbins] = dictionary.thresh[param.texton_nbins-1] + step/100; 

  // collect histograms
  matrix<> hists(param.texton_nbins, param.param_K);
  for (int i = 0; i < param.param_K; i++) {
    const matrix<>& texton = (*textons)[i];
    for (unsigned int j = 0; j < texton.size(); j++) {
      for (int k = 0; k < param.texton_nbins; k++) {
	if ((texton(j,0) >= dictionary.thresh[k]) && (texton(j,0) < dictionary.thresh[k+1])) {
	  hists(k, i)++;
	  break;
	}	  
      }
    }
  }

  // dissimilarities matrix
  double chi2;
  matrix<> diss(param.param_K*2, param.param_K*2);
  for (int i = 0; i < param.param_K; i++) {
    diss(i,i) = 0;
    for (int j = i + 1; j < param.param_K; j++) {
      chi2 = 0;
      for (unsigned int k = 0; k < hists.size(0); k++) {
	if (hists(k, i) + hists(k, j) == 0)
	  chi2 += (hists(k,i) - hists(k,j)) * (hists(k,i) - hists(k,j));
	else
	  chi2 += (hists(k,i) - hists(k,j)) * (hists(k,i) - hists(k,j)) / (hists(k,i) + hists(k,j));
      }
      diss(i,j) = 0.5 * chi2;
      diss(j,i) = 0.5 * chi2;
    }
  }

#if (DEBUG_MODE==DETAIL)
  for (int i = 0; i < 10; i++)
    cout << diss(0,i) << " ";
  cout << endl;
#endif

  // linkage
  bool flag[param.param_K*2];
  int cou[param.param_K*2];
  int linkage[param.param_K*2][2];
  int height[param.param_K*2];
  double linkage_d[param.param_K*2];
  for (int i = 0; i < param.param_K; i++) {
    flag[i] = true;
    cou[i] = 1;
    linkage[i][0] = i;
    linkage[i][1] = i;
    linkage_d[i] = 0;
  }
  for (int i = param.param_K; i < param.param_K*2; i++) {
    flag[i] = false;
    cou[i] = 0;
  }

  int n = param.param_K;
  int min_i1 = -1, min_i2 = -1;
  double min;
  while (true) {
    min = 9999999;
    for (int i = 0; i < n; i++) {
      if (!flag[i])
	continue;
      for (int j = i+1; j < n; j++) {
	if (!flag[j])
	  continue;
	if (diss(i,j) < min) {
	  min = diss(i,j);
	  min_i1 = i;
	  min_i2 = j;
	}
      }
    }

    cou[n] = cou[min_i1] + cou[min_i2];
    flag[min_i1] = false;
    flag[min_i2] = false;
    flag[n] = true;
    linkage[n][0] = min_i1;
    linkage[n][1] = min_i2;
    linkage_d[n] = min;
#if (DEBUG_MODE==DETAIL)
    cout << min_i1 << " " << min_i2 << " " << min << endl;
#endif
    for (int i = 0; i < n; i++) {
      // average linkage method
      diss(n,i) = (diss(min_i1,i) * cou[min_i1] + diss(min_i2,i) * cou[min_i2]) / cou[n];
      diss(i,n) = diss(n,i);
    }
    n++;
    if (cou[n-1] == param.param_K)
      break;
  }


  bool tmp[n];
  for (int i = 0; i < n; i++) {
    tmp[i] = false;
  }
  int vt_cou = 0, vt_cls[param.nbcls];
  for (int i = n-param.nbcls+1; i < n; i++) {
    if (!tmp[linkage[i][0]])
      vt_cls[vt_cou++] = linkage[i][0];
    if (!tmp[linkage[i][1]])
      vt_cls[vt_cou++] = linkage[i][1];
    tmp[i] = true;
  }

  // PABLO VERISON
  double tot_dist[param.param_K];
  for (int i = 0; i < param.param_K; i++) {
    tot_dist[i] = 0;
    for (int j = 0; j < param.param_K; j++)
      tot_dist[i] += diss(i, j);
#if (DEBUG_MODE==DETAIL)
    cout << tot_dist[i] << " ";
#endif
  }

  dictionary.textons.reset(new array_list<matrix<> >());
  int min_ind = -1;
  for (int i = 0; i < vt_cou; i++) {
    min = 99999;

    for (int j = 0; j < n; j++)
      tmp[j] = false;
    tmp[linkage[vt_cls[i]][0]] = true;
    tmp[linkage[vt_cls[i]][1]] = true;

    for (int j = n-1;  j >= 0; j--) {
      if (tmp[j]) {
	tmp[linkage[j][0]] = true;
	tmp[linkage[j][1]] = true;
	
	if (j < param.param_K) {
	  //	  if (diss(vt_cls[i], j) < min) {
	  if (tot_dist[j] < min) {
	    //	    min = diss(vt_cls[i], j);
	    min = tot_dist[j];
	    min_ind = j;
	  }
	}
      }
    }    

#if (DEBUG_MODE==DETAIL)
    cout << vt_cls[i] << " " << min_ind << " " << endl;
#endif
    dictionary.textons->add(*(new matrix<>((*textons)[min_ind])));
  }

#if (DEBUG_MODE)
  cout << "compute_textons: done." << endl;
#endif
}

void compute_texton_map(const matrix<>& img, const DICT& dictionary, const SVM_PARAM_ABSTRACT& param, matrix<>& texton_map) {
  int tx, ty;
  tx = img.size(0);
  ty = img.size(1);

  int u[param.texton_radius*2+1][param.texton_radius*2+1], v[param.texton_radius*2+1][param.texton_radius*2+1];
  for (int i = 0; i < param.texton_radius*2+1; i++)
    for (int j = 0; j < param.texton_radius*2+1; j++) {
      u[i][j] = j - param.texton_radius;
      v[i][j] = i - param.texton_radius;
    }

  bool NHOOD[param.texton_radius*2+1][param.texton_radius*2+1]; 

  if (param.winType == 0) {
    for (int i = 0; i < param.texton_radius*2+1; i++)
      for (int j = 0; j < param.texton_radius*2+1; j++)
	NHOOD[i][j] = (u[i][j]*u[i][j] + v[i][j]*v[i][j] <= param.texton_radius*param.texton_radius); 
  } else {
    for (int i = 0; i < param.texton_radius*2+1; i++)
      for (int j = 0; j < param.texton_radius*2+1; j++)
	NHOOD[i][j] = (abs(u[i][j]) <= param.texton_radius) && (abs(v[i][j]) <= param.texton_radius);
  }

  // histc
  matrix<> hists_ref(param.texton_nbins, dictionary.textons->size());
  for (unsigned int i = 0; i < dictionary.textons->size(); i++) {
    const matrix<>& texton = (*dictionary.textons)[i];
    for (unsigned int j = 0; j < texton.size(); j++) {
      for (int k = 0; k < param.texton_nbins; k++) {
	if ((texton(j,0) >= dictionary.thresh[k]) && (texton(j,0) < dictionary.thresh[k+1])) {
	  hists_ref(k, i)++;
	  break;
	}	  
      }
    }
  }

  double patch[param.texton_radius*2+1][param.texton_radius*2+1];
  for (int x = param.texton_radius; x < tx - param.texton_radius; x++) {
    for (int y = param.texton_radius; y < ty - param.texton_radius; y++) {
      // patch = part of img
      for (int i = x - param.texton_radius; i <= x + param.texton_radius; i++)
	for (int j = y - param.texton_radius; j <= y + param.texton_radius; j++) {
	  patch[i - x + param.texton_radius][j - y + param.texton_radius] = img(i, j); 
	}
      
      // histc
      int hist[param.texton_nbins];
      for (int k = 0; k < param.texton_nbins; k++)
	hist[k] = 0;
      for (int i = 0; i < param.texton_radius*2+1; i++)
	for (int j = 0; j < param.texton_radius*2+1; j++) {
	  if (!NHOOD[i][j])
	    continue;

	  for (int k = 0; k < param.texton_nbins; k++) {
	    if ((patch[i][j] >= dictionary.thresh[k]) && (patch[i][j] < dictionary.thresh[k+1])) {
	      hist[k]++;
	      break;
	    }
	  }
	}

      // chi2
      double min = 99999;
      int min_ind = -1;
      for (unsigned int i = 0; i < hists_ref.size(1); i++) {
	double chi2 = 0;
	for (int k = 0; k < param.texton_nbins; k++) {
	  if (hist[k] + hists_ref(k,i) != 0)
	    chi2 += (hist[k] - hists_ref(k,i)) * (hist[k] - hists_ref(k,i)) / (hist[k] + hists_ref(k,i));
	}
	if (chi2 < min) {
	  min = chi2;
	  min_ind = i;
	}
      }
      texton_map(x, y) = min_ind; 
    }
  }

  //  cout << "compute_texton_map: done." << endl;
}

void orig2pieces(const matrix<>& orig,
		 const string& piecesDir,
		 const int tx_out,
		 const int ty_out,
		 const int overlap,
		 const int x_init,
		 const int y_init) {
  int tx, ty;
  tx = orig.size(0);
  ty = orig.size(1);

  //  cout << tx << " " << ty << endl;
  //  cout << tx_out << " " << ty_out << endl;
  //  cout << overlap << endl;

  int k = 1;
  int x0 = x_init;
  while (x0 < tx) {
    int xf = min_f(x0+tx_out-1, tx);
    int l = 1;
    int y0 = y_init;
    while (y0 < ty) {
      int yf = min_f(y0+ty_out-1, ty);
      array<unsigned long> start(2);
      start[0] = x0-1;
      start[1] = y0-1;
      array<unsigned long> end(2);
      end[0] = xf-1;
      end[1] = yf-1;

      matrix<> img = orig.submatrix(start, end);

      int tx_piece, ty_piece;
      tx_piece = img.size(0);
      ty_piece = img.size(1);

      char buffer[50];
      sprintf(buffer, "%02d-%02d.img", k, l);
      
      array<int> extra_info(6);
      extra_info[0] = x0;
      extra_info[1] = y0;
      extra_info[2] = tx;
      extra_info[3] = ty;
      extra_info[4] = tx_piece;
      extra_info[5] = ty_piece;

      save_image(img, piecesDir, string(buffer), extra_info);

      y0 = y0+ty_out-overlap;
      l=l+1;
      if (y0+overlap >= ty) {
	break;
      }
    }
    x0 = x0+tx_out-overlap;
    k=k+1;
    if (x0+overlap >= tx) {
      break;
    }
  }

#if (DEBUG_MODE)
  cout << "orig2pieces: done." << endl;
#endif
}

void TEXTONSVM_batch_dict(const SVM_PARAM& param) {
  // load param file (probably I will just pass this from TEXTONSVM_batch_svm_opt)
  
  // mkdirs.. (pieces_det, orig, clsDir)

#if (DEBUG_MODE)
  cout << "Computing Dictionary..." << endl;
#endif

  // load denoised image
#if (DEBUG_MODE==DETAIL)
  cout << param.dataFile << endl;
#endif
  matrix<> img;
  load_denoise_image(param.dataFile, RESULT_DIR, img);

  // load denoised examplars
  safe_vector<matrix<>*> examplars_denoise;
  load_examplars(RESULT_DIR, string("denoised_examplars.dat"), TEXTON_RADIUS, K, WINDOW_RADIUS, examplars_denoise);

  // load textons
  auto_collection<matrix<>, array_list<matrix<> > > textons;
  load_textons(RESULT_DIR, TEXTON_RADIUS, K, textons);

  // compute textons
  DICT dictionary;
  compute_textons(textons, param, dictionary);

  // store dict (model.mat)
  dictionary.save(RESULT_DIR, param.modelsFile);

  // compute texton map
  safe_vector<matrix<>*> examplars_tmap;
  int tx, ty;
  tx = examplars_denoise.at(0)->size(0);
  ty = examplars_denoise.at(0)->size(1);
#if (DEBUG_MODE)
  cout << "[";
#endif
  for (int i = 0; i < examplars_denoise.size(); i++) {
    examplars_tmap.push_back(new matrix<>(tx, ty));
    compute_texton_map(*examplars_denoise.at(i), dictionary, param, *examplars_tmap.at(i));

#if (DEBUG_MODE)
    cout << ".";
#endif
  }
#if (DEBUG_MODE)
  cout << "]";
#endif
  cout << endl;

  // store examplars_tmap (examplars_denoise.mat append)
  save_examplars(RESULT_DIR, string("examplars_tmap_denoise.dat"), TEXTON_RADIUS, K, examplars_tmap);

  orig2pieces(img, param.pieces_origDir, param.tx_out, param.ty_out, param.overlap, param.x_init, param.y_init);

  touchFile(param.batchDir, "dict.touch");

  cout << "TEXTONSVM_batch_dict: done." << endl;
}

void create_piece_texton_map(const DICT& dictionary,
			     const SVM_PARAM_ABSTRACT& param,
			     const string& pieceName,
			     const string& origDir,
			     const string& pieces_clsDir) {
  cout << "tmap: " << origDir << pieceName << endl;

  matrix<> img;
  array<int> extra_info;

  extra_info = load_image(pieceName, origDir, img, true /*cut_range*/);
  
  //  cout << img.size(0) << " " << img.size(1) << endl;
  //  cout << extra_info(2) << " " << extra_info(3) << endl;
  matrix<> texton_map(img.size(0), img.size(1));
  compute_texton_map(img, dictionary, param, texton_map);

  array<int> extra_info2(4);
  extra_info2(0) = extra_info(0); // x0
  extra_info2(1) = extra_info(1); // y0
  extra_info2(2) = extra_info(2); // tx
  extra_info2(3) = extra_info(3); // ty
  
  save_image(texton_map, pieces_clsDir, pieceName, extra_info2);
}

void pieces2texton_map(const string& texton_mapDir, const int& intwin_rad, matrix<>& texton_map_global) {
  cout << "pieces2texton_map started" << endl;

  DIR *dirp = opendir(texton_mapDir.c_str());
  if (dirp) {
    int iter = 0;
    struct dirent *entry;
    while ((entry = readdir(dirp))) {
      if (entry->d_type == 8) {
	string pieceName = string(entry->d_name);
	if (strcmp(pieceName.substr(pieceName.size()-3, 3).c_str(), "img") != 0)
	  continue;

	matrix<> img;
	array<int> extra_info = load_image(pieceName, texton_mapDir, img, false /*cut_range*/);
	if (iter == 0)
	  texton_map_global = matrix<>(extra_info[2], extra_info[3]);

	int tx_piece, ty_piece;
	tx_piece = img.size(0);
	ty_piece = img.size(1);

	int x0 = extra_info[0];
	int y0 = extra_info[1]; 

	// TODO: check coordinate..
	array<int> start(2);
	start[0] = x0+intwin_rad;
	start[1] = y0+intwin_rad;
	array<int> end(2);
	end[0] = x0+tx_piece-1-intwin_rad;
	end[1] = y0+ty_piece-1-intwin_rad;

	array<int> start2(2);
	start2[0] = intwin_rad;
	start2[1] = intwin_rad;
	array<int> end2(2);
	end2[0] = tx_piece-intwin_rad-1;
	end2[1] = ty_piece-intwin_rad-1;

#if DEBUG_MODE
	printf("global: %d %d\n", extra_info[2], extra_info[3]);
	printf("start end: %d %d %d %d\n", start[0], start[1], end[0], end[1]);
	printf("piece tx ty: %d %d\n", tx_piece, ty_piece);
	printf("start2 end2: %d %d %d %d\n", start2[0], start2[1], end2[0], end2[1]);
#endif

	texton_map_global.subassign(start, end, img.submatrix(start2, end2));
	iter++;

#if DEBUG_MODE
	printf("Iter %d done.\n", iter);
#endif
      }
    }
    closedir(dirp);
  }

  cout << "pieces2texton_map ended" << endl;
}

void texton_map2pieces(const matrix<>& texton_map_global,
		       const string& piecesDir,
		       const int tx_out,
		       const int ty_out,
		       const int overlap,
		       const int x_init,
		       const int y_init) {
  int tx, ty;
  tx = texton_map_global.size(0);
  ty = texton_map_global.size(1);
  int k = 1;
  int x0 = x_init;
  while (x0 < tx) {
    int xf = min_f(x0+tx_out-1, tx);
    int l = 1;
    int y0 = y_init;
    while (y0 < ty) {
      int yf = min(y0+ty_out-1, ty);
      
      array<int> start(2);
      start[0] = x0-1;
      start[1] = y0-1;
      array<int> end(2);
      end[0] = xf-1;
      end[1] = yf-1;
      matrix<> texton_map = texton_map_global.submatrix(start, end);
      
      array<int> extra_info(4);
      extra_info[0] = x0;
      extra_info[1] = y0;
      extra_info[2] = tx;
      extra_info[3] = ty;

      char buffer[50];
      sprintf(buffer, "%02d-%02d.img", k, l);
      
      save_image(texton_map, piecesDir, string(buffer), extra_info);
      y0 += ty_out - overlap;
      l++;
      if (y0 + overlap >= ty) {
	break;
      }
    }
    x0 += tx_out - overlap;
    k++;
    if (x0 + overlap >= tx) {
      break;
    }
  }

  cout << "texton_map2pieces: done." << endl;
}

void compute_svm_model(const safe_vector<matrix<>*>& examplars,
		       const SVM_PARAM& param,
		       const safe_vector<matrix<>*>& examplars_tmap,
		       svm_model*& model,
		       vector<double>& thresh) {
  // TODO: if nargin < 3

  // meshgrid
  int u[param.window_radius*2+1][param.window_radius*2+1], v[param.window_radius*2+1][param.window_radius*2+1];
  for (int i = 0; i < param.window_radius*2+1; i++)
    for (int j = 0; j < param.window_radius*2+1; j++) {
      u[i][j] = j - param.window_radius;
      v[i][j] = i - param.window_radius;
    }

  bool NHOOD[param.window_radius*2+1][param.window_radius*2+1];
  bool CROWN[param.window_radius*2+1][param.window_radius*2+1];

  if (param.winType == 0) {
    for (int i = 0; i < param.window_radius*2+1; i++)
      for (int j = 0; j < param.window_radius*2+1; j++) {
	NHOOD[i][j] = (u[i][j]*u[i][j] + v[i][j]*v[i][j] <= param.texton_radius*param.texton_radius); 
	// TODO: crown
      }
  } else {
    for (int i = 0; i < param.window_radius*2+1; i++)
      for (int j = 0; j < param.window_radius*2+1; j++) {
	NHOOD[i][j] = (abs(u[i][j]) <= param.particle_radius) && (abs(v[i][j]) <= param.particle_radius);
	CROWN[i][j] = (abs(u[i][j]) > param.particle_radius) || (abs(v[i][j]) > param.particle_radius);
      }
  }

  int nb_examplars = examplars.size();
  int tx = (*examplars.at(0)).size(0), ty = (*examplars.at(0)).size(1);

  if ((tx < param.window_radius) || (ty < param.window_radius)) {
    //TODO: error
  }

  int cx = (tx-1)/2, cy = (ty - 1)/2;

  array<double> Mxs(nb_examplars), Mns(nb_examplars);
  double min_exmp = 0, max_exmp = 0;
  for (int e = 0; e < nb_examplars; e++) {
    Mxs[e] = max(*examplars.at(e));
    max_exmp += Mxs[e];

    Mns[e] = min(*examplars.at(e));
    min_exmp += Mns[e];
  }
  max_exmp /= nb_examplars;
  min_exmp /= nb_examplars;

  svm_problem data;
  data.l = nb_examplars*4;
  data.y = static_cast<double*>(malloc(sizeof(double)*data.l));
  data.x = static_cast<struct svm_node**>(malloc(sizeof(struct svm_node*)*data.l));
  for (int i = 0; i < data.l; i++) {
    data.x[i] = static_cast<struct svm_node*>(malloc(sizeof(struct svm_node)*(2*param.grayhist_nbins+2*param.nbcls+1)));
  }

  //  double trd[nb_examplars*4][(2*param.grayhist_nbins)+(2*param.nbcls)];
  //  double trl[nb_examplars*4];

  double step = (max_exmp - min_exmp) / (param.grayhist_nbins-1);
  //  double thresh[param.grayhist_nbins+1];
  thresh.resize(static_cast<unsigned long>(param.grayhist_nbins)+1);
  for (int k = 0; k < param.grayhist_nbins; k++) 
    thresh[k] = min_exmp + step*k;
  thresh[param.grayhist_nbins] = thresh[param.grayhist_nbins-1] + step/100;

  for (int e = 0; e < nb_examplars; e++) {
    double grayhist_part[param.grayhist_nbins], grayhist_bgd[param.grayhist_nbins];
    double texthist_part[param.nbcls], texthist_bgd[param.nbcls];
    int grayhist_part_cou, grayhist_bgd_cou, texthist_part_cou, texthist_bgd_cou;
    for (int k = 0; k < param.grayhist_nbins; k++)
      grayhist_part[k] = grayhist_bgd[k] = 0;
    for (int k = 0; k < param.nbcls; k++)
      texthist_part[k] = texthist_bgd[k] = 0;
    grayhist_part_cou = grayhist_bgd_cou = texthist_part_cou = texthist_bgd_cou = 0;

    array<int> start(2);
    start[0] = cx - param.window_radius;
    start[1] = cy - param.window_radius;
    array<int> end(2);
    end[0] = cx + param.window_radius;
    end[1] = cy + param.window_radius;

    // collect positive examples.. which will be also used for negative examples
    matrix<> patch = (*examplars.at(e)).submatrix(start,end);
    matrix<> patch2 = (*examplars_tmap.at(e)).submatrix(start,end);
    for (int x = 0 ; x < param.window_radius*2+1; x++)
      for (int y = 0; y < param.window_radius*2+1; y++) {
	if (NHOOD[x][y]) {
	  for (int k = 0; k < param.grayhist_nbins; k++)
	    if ((patch(x,y) >= thresh[k]) && (patch(x,y) < thresh[k+1])) {
	      grayhist_part[k]++;
	      grayhist_part_cou++; 
	      break;
	    }
	  for (int k = 0; k < param.nbcls; k++)
	    if ((patch2(x,y) >= k - 0.001) && (patch2(x,y) <= k + 0.001)) {
	      texthist_part[k]++;
	      texthist_part_cou++; 
	      break;
	    }
	}
	if (CROWN[x][y]) {
	  for (int k = 0; k < param.grayhist_nbins; k++)
	    if ((patch(x,y) >= thresh[k]) && (patch(x,y) < thresh[k+1])) {
	      grayhist_bgd[k]++;
	      grayhist_bgd_cou++; 
	      break;
	    }
	  for (int k = 0; k < param.nbcls; k++)
	    if ((patch2(x,y) >= k - 0.001) && (patch2(x,y) <= k + 0.001)) {
	      texthist_bgd[k]++;
	      texthist_bgd_cou++;
	      break;
	    }
	}
      }

    // normalize hist with their cou
    for (int k = 0; k < param.grayhist_nbins; k++) {
      grayhist_part[k] /= grayhist_part_cou;
      grayhist_bgd[k] /= grayhist_bgd_cou;
    }
    for (int k = 0; k < param.nbcls; k++) { 
      texthist_part[k] /= texthist_part_cou;
      texthist_bgd[k] /= texthist_bgd_cou;
    }

    // create a real training dataset. YAY!!
    data.y[e*4] = 1;
    data.y[e*4+1] = 0;
    data.y[e*4+2] = 0;
    data.y[e*4+3] = 0;
    for (int k = 0; k < param.grayhist_nbins; k++) {
      data.x[e*4][k].value = grayhist_part[k] * param.weight_feature / 2;
      data.x[e*4][param.grayhist_nbins+k].value = grayhist_bgd[k] * param.weight_feature / 2;
      data.x[e*4][k].index = k+1;
      data.x[e*4][param.grayhist_nbins+k].index = param.grayhist_nbins+k+1;

      data.x[e*4+1][k].value = grayhist_bgd[k] * param.weight_feature / 2;
      data.x[e*4+1][param.grayhist_nbins+k].value = grayhist_part[k] * param.weight_feature / 2;
      data.x[e*4+1][k].index = k+1;
      data.x[e*4+1][param.grayhist_nbins+k].index = param.grayhist_nbins+k+1;

      data.x[e*4+2][k].value = grayhist_part[k] * param.weight_feature / 2;
      data.x[e*4+2][param.grayhist_nbins+k].value = grayhist_part[k] * param.weight_feature / 2;
      data.x[e*4+2][k].index = k+1;
      data.x[e*4+2][param.grayhist_nbins+k].index = param.grayhist_nbins+k+1;

      data.x[e*4+3][k].value = grayhist_bgd[k] * param.weight_feature / 2;
      data.x[e*4+3][param.grayhist_nbins+k].value = grayhist_bgd[k] * param.weight_feature / 2;
      data.x[e*4+3][k].index = k+1;
      data.x[e*4+3][param.grayhist_nbins+k].index = param.grayhist_nbins+k+1;
    }
    for (int k = 0; k < param.nbcls; k++) {
      data.x[e*4][param.grayhist_nbins*2+k].value = texthist_part[k] * (1-param.weight_feature) / 2;
      data.x[e*4][param.grayhist_nbins*2+param.nbcls+k].value = texthist_bgd[k] * (1-param.weight_feature) / 2;
      data.x[e*4][param.grayhist_nbins*2+k].index = param.grayhist_nbins*2+k+1;
      data.x[e*4][param.grayhist_nbins*2+param.nbcls+k].index = param.grayhist_nbins*2+param.nbcls+k+1;

      data.x[e*4+1][param.grayhist_nbins*2+k].value = texthist_bgd[k] * (1-param.weight_feature) / 2;
      data.x[e*4+1][param.grayhist_nbins*2+param.nbcls+k].value = texthist_part[k] * (1-param.weight_feature) / 2;
      data.x[e*4+1][param.grayhist_nbins*2+k].index = param.grayhist_nbins*2+k+1;
      data.x[e*4+1][param.grayhist_nbins*2+param.nbcls+k].index = param.grayhist_nbins*2+param.nbcls+k+1;

      data.x[e*4+2][param.grayhist_nbins*2+k].value = texthist_part[k] * (1-param.weight_feature) / 2;
      data.x[e*4+2][param.grayhist_nbins*2+param.nbcls+k].value = texthist_part[k] * (1-param.weight_feature) / 2;
      data.x[e*4+2][param.grayhist_nbins*2+k].index = param.grayhist_nbins*2+k+1;
      data.x[e*4+2][param.grayhist_nbins*2+param.nbcls+k].index = param.grayhist_nbins*2+param.nbcls+k+1;

      data.x[e*4+3][param.grayhist_nbins*2+k].value = texthist_bgd[k] * (1-param.weight_feature) / 2;
      data.x[e*4+3][param.grayhist_nbins*2+param.nbcls+k].value = texthist_bgd[k] * (1-param.weight_feature) / 2;
      data.x[e*4+3][param.grayhist_nbins*2+k].index = param.grayhist_nbins*2+k+1;
      data.x[e*4+3][param.grayhist_nbins*2+param.nbcls+k].index = param.grayhist_nbins*2+param.nbcls+k+1;
    }
    data.x[e*4][param.grayhist_nbins*2+param.nbcls*2].index = -1;
    data.x[e*4+1][param.grayhist_nbins*2+param.nbcls*2].index = -1;
    data.x[e*4+2][param.grayhist_nbins*2+param.nbcls*2].index = -1;
    data.x[e*4+3][param.grayhist_nbins*2+param.nbcls*2].index = -1;
  }

  // now let's train FINALLY!!!!!
  svm_parameter svm_param;
  svm_param.svm_type = C_SVC;
  svm_param.kernel_type = SUMMIN;
  svm_param.degree = 3;
  svm_param.gamma = 0;	// 1/k
  //  svm_param.gamma = 1/1000;
  svm_param.coef0 = 0;
  svm_param.nu = 0.5;
  svm_param.cache_size = 100;
  svm_param.C = 1;
  svm_param.eps = 1e-3;
  svm_param.p = 0.1;
  svm_param.shrinking = 1;
  svm_param.probability = 0;
  svm_param.nr_weight = 2;
  svm_param.weight_label = static_cast<int*>(malloc(sizeof(int)*2));
  svm_param.weight_label[0] = 1;
  svm_param.weight_label[1] = 0;
  svm_param.weight = static_cast<double*>(malloc(sizeof(double)*2));
  svm_param.weight[0] = 3;
  svm_param.weight[1] = 1;
  
  cout << "TEXTONSVM: svm_train" << endl;
  model = svm_train(&data, &svm_param);

  svm_destroy_param(&svm_param);

  svm_save_model((RESULT_DIR + "tmp.svm").c_str(), model);
  svm_destroy_model(model);

  free(data.y);
  for (int i = 0; i < data.l; i++)
    free(data.x[i]);
  free(data.x);

  // TODO: later you should load this model somewhere else more appropriate
  model = svm_load_model((RESULT_DIR + "tmp.svm").c_str());
}

void TEXTONSVM_batch_update(const SVM_PARAM& param, svm_model*& model, vector<double>& thresh) {
  cout << "TEXTONSVM_batch_upodate: start" << endl;

  // orig2pieces
  cout << param.dataFile << endl;
  matrix<> img;
  load_denoise_image(param.dataFile, RESULT_DIR, img);
  orig2pieces(img, param.pieces_origDir, param.tx_out, param.ty_out, param.overlap, param.x_init, param.y_init);

  // texton_map2pieces
  // TODO: I think I need to uncomment and change other part
  //  matrix<> texton_map;
  //  load_image(string("texton_map.img"), resultsDir, texton_map, false /*cur_range*/);
  //  texton_map2pieces(texton_map, param.pieces_clsDir, param.tx_out, param.ty_out, param.overlap, param.x_init, param.y_init);

  // load examplars_denoise
  safe_vector<matrix<>*> examplars_denoise;
  safe_vector<matrix<>*> examplars_tmap;
  load_examplars(RESULT_DIR, string("denoised_examplars.dat"), TEXTON_RADIUS, K, WINDOW_RADIUS, examplars_denoise);
  load_examplars(RESULT_DIR, string("examplars_tmap_denoise.dat"), TEXTON_RADIUS, K, WINDOW_RADIUS, examplars_tmap);

  // TODO: [model, thresh] = compute_svm_model
  compute_svm_model(examplars_denoise, param, examplars_tmap, model, thresh);

#if (DEBUG_MODE==DETAIL)
  cout << model->nr_class << endl;
  cout << model->l << endl;
#endif

  // TODO: save -append model, thresh

  // touchFile
  touchFile(param.batchDir, "update.touch");

  cout << "TEXTONSVM_batch_update: done." << endl;
}

void compute_svm_distance(const matrix<>& img, const svm_model* model, const SVM_PARAM_ABSTRACT& param, const vector<double>& thresh, const matrix<>& texton_map, matrix<>& detections) {
  cout << "compute_svm_distance: started." << endl;

  int tx = img.size(0), ty = img.size(1);

  bool NHOOD[param.window_radius*2+1][param.window_radius*2+1];
  bool CROWN[param.window_radius*2+1][param.window_radius*2+1];

  double* IH_gray = static_cast<double*>(malloc(param.grayhist_nbins*tx*ty*sizeof(double)*2));
  double* CTD = static_cast<double*>(malloc(param.nbcls*tx*ty*sizeof(double)*2));

#if (DEBUG_MODE==DETAIL)
  cout << thresh.size() << endl;
#endif

  if (param.winType == 0) {
    // meshgrid
    int u[param.window_radius*2+1][param.window_radius*2+1], v[param.window_radius*2+1][param.window_radius*2+1];
    for (int i = 0; i < param.window_radius*2+1; i++)
      for (int j = 0; j < param.window_radius*2+1; j++) {
	u[i][j] = j - param.window_radius;
	v[i][j] = i - param.window_radius;
      }

    for (int i = 0; i < param.window_radius*2+1; i++)
      for (int j = 0; j < param.window_radius*2+1; j++) {
	NHOOD[i][j] = (u[i][j]*u[i][j] + v[i][j]*v[i][j] <= param.texton_radius*param.texton_radius); 
	// TODO: crown
      }
  } else {
    int c;

    for (c = 1; c < param.grayhist_nbins-1; c++) {
      IH_gray[c3to1(c,0,0,tx,ty)] = (thresh[c] <= img(0,0)) && (img(0,0) < thresh[c+1]);
      for (int y = 1; y < ty; y++)
	IH_gray[c3to1(c,0,y,tx,ty)] = IH_gray[c3to1(c,0,y-1,tx,ty)] + ((thresh[c] <= img(0,y)) && (img(0,y) < thresh[c+1])) ;
      for (int x = 1; x < tx; x++) {
	IH_gray[c3to1(c,x,0,tx,ty)] = IH_gray[c3to1(c,x-1,0,tx,ty)] + ((thresh[c] <= img(x,0)) && (img(x,0) < thresh[c+1]));
	for (int y = 1; y < ty; y++) {
	  IH_gray[c3to1(c,x,y,tx,ty)] = IH_gray[c3to1(c,x,y-1,tx,ty)] + IH_gray[c3to1(c,x-1,y,tx,ty)] - IH_gray[c3to1(c,x-1,y-1,tx,ty)] + ((thresh[c] <= img(x,y)) && (img(x,y) < thresh[c+1]));
	}
      }
    }

    // TODO: bug
    c = 0;
    IH_gray[c3to1(c,0,0,tx,ty)] = (img(0,0) < thresh[0]);
    for (int y = 1; y < ty; y++)
      IH_gray[c3to1(c,0,y,tx,ty)] = IH_gray[c3to1(c,0,y-1,tx,ty)] + (img(0,y) < thresh[c]);
    for (int x = 1; x < tx; x++) {
      IH_gray[c3to1(c,x,0,tx,ty)] = IH_gray[c3to1(c,x-1,0,tx,ty)] + (img(x,0) < thresh[0]);
      for (int y = 1; y < ty; y++)
	IH_gray[c3to1(c,x,y,tx,ty)] = IH_gray[c3to1(c,x,y-1,tx,ty)] + IH_gray[c3to1(c,x-1,y,tx,ty)] - IH_gray[c3to1(c,x-1,y-1,tx,ty)] + (img(x,y) < thresh[0]);
    }

    c = param.grayhist_nbins-1;
    IH_gray[c3to1(c,0,0,tx,ty)] = (img(0,0) >= thresh[c]);
    for (int y = 1; y < ty; y++)
      IH_gray[c3to1(c,0,y,tx,ty)] = IH_gray[c3to1(c,0,y-1,tx,ty)] + (img(0,y) >= thresh[c]);
    for (int x = 1; x < tx; x++) {
      IH_gray[c3to1(c,x,0,tx,ty)] = IH_gray[c3to1(c,x-1,0,tx,ty)] + (img(x,0) >= thresh[c]);
      for (int y = 1; y < ty; y++)
	IH_gray[c3to1(c,x,y,tx,ty)] = IH_gray[c3to1(c,x,y-1,tx,ty)] + IH_gray[c3to1(c,x-1,y,tx,ty)] - IH_gray[c3to1(c,x-1,y-1,tx,ty)] + (img(x,y) >= thresh[c]);
    }

    for (c = 0; c < param.nbcls; c++) {
      CTD[c3to1(c,0,0,tx,ty)] = texton_map(0,0) == c;
	for (int y = 1; y < ty; y++)
	  CTD[c3to1(c,0,y,tx,ty)] = CTD[c3to1(c,0,y-1,tx,ty)] + (texton_map(0,y) == c);
	for (int x = 1; x < tx; x++) {
	  CTD[c3to1(c,x,0,tx,ty)] = CTD[c3to1(c,x-1,0,tx,ty)] + (texton_map(x,0) == c);
	  for (int y = 1; y < ty; y++) {
	    CTD[c3to1(c,x,y,tx,ty)] = CTD[c3to1(c,x,y-1,tx,ty)] + CTD[c3to1(c,x-1,y,tx,ty)] - CTD[c3to1(c,x-1,y-1,tx,ty)] + (texton_map(x,y) == c);
	  }
	}
    }
  }

  detections.resize(static_cast<long unsigned>(tx), static_cast<long unsigned>(ty));
  int radius = param.window_radius + param.texton_radius;

  // TODO
  double min = 99999, max = -99999;

  struct svm_node* data = static_cast<struct svm_node*>(malloc(sizeof(struct svm_node)*(2*param.grayhist_nbins+2*param.nbcls+1)));
  for (int i = 0; i < 2*param.grayhist_nbins+2*param.nbcls+1; i++)
    data[i].value = 0;

  for (int x = radius; x < tx - radius; x++) {
    for (int y = radius; y < ty - radius; y++) {
      double h1[param.grayhist_nbins], h2[param.grayhist_nbins], h3[param.nbcls], h4[param.nbcls];
      if (param.winType == 0) {
      } else {
	double sum_h1,sum_h2,sum_h3,sum_h4;
	sum_h1=sum_h2=sum_h3=sum_h4=0;

	for (int c = 0; c < param.grayhist_nbins; c++) {
	  // internal square
	  h1[c] = IH_gray[c3to1(c, x+param.particle_radius, y+param.particle_radius, tx,ty)] + IH_gray[c3to1(c, x-param.particle_radius-1, y-param.particle_radius-1, tx,ty)] -
	    IH_gray[c3to1(c, x+param.particle_radius, y-param.particle_radius-1, tx,ty)] - IH_gray[c3to1(c, x-param.particle_radius-1, y+param.particle_radius, tx,ty)];
	  // external square
	  h2[c] = IH_gray[c3to1(c, x+param.window_radius, y+param.window_radius, tx,ty)] + IH_gray[c3to1(c, x-param.window_radius-1, y-param.window_radius-1, tx,ty)] -
	    IH_gray[c3to1(c, x+param.window_radius, y-param.window_radius-1, tx,ty)] - IH_gray[c3to1(c, x-param.window_radius-1, y+param.window_radius, tx,ty)];
	  // crown
	  h2[c] = h2[c] - h1[c];

	  sum_h1 += h1[c];
	  sum_h2 += h2[c];
	}

	for (int c = 0; c < param.nbcls; c++) {
	  // internal square
	  h3[c] = CTD[c3to1(c, x+param.particle_radius, y+param.particle_radius, tx,ty)] + CTD[c3to1(c, x-param.particle_radius-1, y-param.particle_radius-1, tx,ty)] -
	    CTD[c3to1(c, x+param.particle_radius, y-param.particle_radius-1, tx,ty)] - CTD[c3to1(c, x-param.particle_radius-1, y+param.particle_radius, tx,ty)];
	  // external square
	  h4[c] = CTD[c3to1(c, x+param.window_radius, y+param.window_radius, tx,ty)] + CTD[c3to1(c, x-param.window_radius-1, y-param.window_radius-1, tx,ty)] -
	    CTD[c3to1(c, x+param.window_radius, y-param.window_radius-1, tx,ty)] - CTD[c3to1(c, x-param.window_radius-1, y+param.window_radius, tx,ty)];
	  // crown
	  h4[c] = h4[c] - h3[c];

	  sum_h3 += h3[c];
	  sum_h4 += h4[c];
	}

	if (sum_h1 == 0) sum_h1 = 1;
	if (sum_h2 == 0) sum_h2 = 1;
	if (sum_h3 == 0) sum_h3 = 1;
	if (sum_h4 == 0) sum_h4 = 1;

	for (int c = 0; c < param.grayhist_nbins; c++) {
	  data[c].value = h1[c] / sum_h1 * param.weight_feature / 2;
	  data[c+param.grayhist_nbins].value = h2[c] / sum_h2 * param.weight_feature / 2;
	  data[c].index = c+1;
	  data[c+param.grayhist_nbins].index = c+param.grayhist_nbins+1;
	}
	for (int c = 0; c < param.nbcls; c++) {
	  data[c+param.grayhist_nbins*2].value = h3[c] / sum_h3 * (1-param.weight_feature) / 2;
	  data[c+param.grayhist_nbins*2+param.nbcls].value = h4[c] / sum_h4 * (1-param.weight_feature) / 2;
	  data[c+param.grayhist_nbins*2].index = c+param.grayhist_nbins*2+1;
	  data[c+param.grayhist_nbins*2+param.nbcls].index = c+param.grayhist_nbins*2+param.nbcls+1;
	}

	data[2*param.grayhist_nbins+2*param.nbcls].index = -1;

	double dec = 0;
	//TODO: fast predict? 
	svm_predict_values(model, data, &dec);
	detections(x, y) = dec;

	if (dec < min)
	  min = dec;
	if (dec > max)
	  max = dec;
      }
    }
  }

  free(IH_gray);
  free(CTD);
  free(data);

#if (DEBUG_MODE==DETAIL)
  cout << min << " " << max << endl;
#endif

  for (int y = 0; y < ty; y++) {
    for (int x = 0; x < radius; x++)
      detections(x, y) = min;
    for (int x = tx - radius; x < tx; x++) 
      detections(x, y) = min;
  }
  for (int x = 0; x < tx; x++) {
    for (int y = 0; y < radius; y++)
      detections(x, y) = min;
    for (int y = ty - radius; y < ty; y++)
      detections(x, y) = min;
  }

  cout << "compute_svm_distance: done." << endl;
}

void create_piece_det(const string& pieceName, const SVM_PARAM_ABSTRACT& param, const svm_model *model, const vector<double>& thresh) {
  matrix<> img, texton_map;
  array<int> extra_info;

  cout << "create_piece_det started." << endl;
#if (DEBUG_MODE==DETAIL)
  cout << "original: " << param.pieces_origDir << endl;
  cout << "cls: " << param.pieces_clsDir << pieceName << endl;
  cout << "number of SVs: " << model->l << endl;
  cout << "number of thresh: " << thresh.size() << endl;
#endif

  extra_info = load_image(pieceName, param.pieces_origDir, img, true /*cut_range*/);
  load_image(pieceName, param.pieces_clsDir, texton_map, false /*cut_range*/);

  matrix<> detections;
  compute_svm_distance(img, model, param, thresh, texton_map, detections);

  //TODO: saving extra_info; in Pablo's code, he takes it from the original image again
  save_image(detections, param.pieces_detDir, pieceName, extra_info);

  cout << "create_piece_det done." << endl << endl; 
}

void pieces2det(const string& resultsDir, const int& radius, matrix<>& det_global) {
  cout << "pieces2det started" << endl;

  matrix<> img;
  array<int> extra_info;
  bool first = true;

  DIR *dirp = opendir(resultsDir.c_str());
  if (dirp) {
    struct dirent *entry;
    while ((entry = readdir(dirp))) {
      if (entry->d_type == 8) {
	string pieceName = string(entry->d_name);
	if (strcmp(pieceName.substr(pieceName.size()-3, 3).c_str(), "img") != 0)
	  continue;
	cout << pieceName << endl;

	extra_info = load_image(pieceName, resultsDir, img, false);
	if (first) {
	  det_global.resize(static_cast<long unsigned>(extra_info[2]), static_cast<long unsigned>(extra_info[3])); 
	  first = false;
	}

	int x0 = extra_info[0]-1, y0 = extra_info[1]-1;
	int tx = extra_info[4], ty = extra_info[5];
	if ((tx != static_cast<int>(img.size(0))) || (ty != static_cast<int>(img.size(1)))) {
	  cout << "ERROR!!";
	  exit(1);
	}
	//	cout << x0 << " " << y0 << " " << tx << " " << ty << endl; 
	for (int mx = radius; mx < tx-radius; mx++)
	  for (int my = radius; my < ty-radius; my++) {
	    det_global(x0+mx, y0+my) = img(mx, my);
	  }
      }
    }
  }

  cout << "pieces2det finished" << endl;
}

void im_normalize(matrix<>& detections, const double& range_min, const double& range_max) {
  for (unsigned int x = 0; x < detections.size(0); x++) {
    for (unsigned int y = 0; y < detections.size(1); y++) {
      detections(x, y) = round((detections(x, y) - range_min) / (range_max - range_min) * 255);
    }
  }
} 

void morphological_op(const matrix<>& detections, const SVM_PARAM_ABSTRACT& param, matrix<int>& labels, int& maxima_n) {
  cout << "morphological operation started" << endl;

  int tx = detections.size(0), ty = detections.size(1);

  labels.resize(static_cast<long unsigned>(tx), static_cast<long unsigned>(ty));

  // TODO:
  //  double thresh_Hmax = param.thresh_Hmax;
  double thresh_Hmax = param.thresh_Hmax;

  int next_x, next_y;
  matrix<> I(static_cast<long unsigned>(tx), static_cast<long unsigned>(ty)), J(static_cast<long unsigned>(tx), static_cast<long unsigned>(ty));
  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      I(x,y) = detections(x, y);
      J(x,y) = I(x,y) - thresh_Hmax;
    }

  double max_val;
  for (int y = 0; y < ty; y++)
    for (int x = 0; x < tx; x++) {
      max_val = J(x,y);
      for (int k = 0; k < 4; k++) {
	next_x = x + dx[k];
	next_y = y + dy[k];
	if ((next_x >= 0) && (next_x < tx) && (next_y >= 0) && (next_y < ty)) 
	  if (J(next_x,next_y) > max_val)
	    max_val = J(next_x,next_y);
      }

      if (max_val < I(x,y))
	J(x,y) = max_val;
      else
	J(x,y) = I(x,y);
    }

  queue<int> queue_x;
  queue<int> queue_y;
  for (int y = ty-1; y >= 0; y--)
    for (int x = tx-1; x >= 0; x--) {
      max_val = J(x,y);
      for (int k = 4; k < 8; k++) {
	next_x = x + dx[k];
	next_y = y + dy[k];
	if ((next_x >= 0) && (next_x < tx) && (next_y >= 0) && (next_y < ty))
	  if (J(next_x, next_y) > max_val)
	    max_val = J(next_x, next_y);
      }

      if (max_val < I(x,y))
	J(x,y) = max_val;
      else
	J(x,y) = I(x,y);

      for (int k = 4; k < 8; k++) {
	next_x = x + dx[k];
	next_y = y + dy[k];
	if ((next_x >= 0) && (next_x < tx) && (next_y >= 0) && (next_y < ty))
	  if ((J(next_x, next_y) < J(x, y)) && (J(next_x, next_y) < I(x,y))) {
	    queue_x.push(x);
	    queue_y.push(y);
	    break;
	  }
      }
    }

  int cur_x, cur_y;
  while (!queue_x.empty()) {
    cur_x = queue_x.front();
    queue_x.pop();
    cur_y = queue_y.front();
    queue_y.pop();

    for (int k = 0; k < 8; k++) {
      next_x = cur_x + dx[k];
      next_y = cur_y + dy[k];
      if ((next_x >= 0) && (next_x < tx) && (next_y >= 0) && (next_y < ty)) {
	if ((J(next_x, next_y) < J(cur_x, cur_y)) && (I(next_x, next_y) != J(next_x, next_y))) {
	  if (J(cur_x, cur_y) < I(next_x, next_y))
	    J(next_x, next_y) = J(cur_x, cur_y);
	  else
	    J(next_x, next_y) = I(next_x, next_y);
	  queue_x.push(next_x);
	  queue_y.push(next_y);
	}
      }      
    }
  }

  // sanity check
  while (!queue_x.empty())
    queue_x.pop();
  while (!queue_y.empty())
    queue_y.pop();

  int group = 0;
  matrix<int> group_n(tx, ty);
  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++)
      group_n(x, y) = -1;

  bool maxima;
  double cur_value;

  vector<bool> maxima_group(tx*ty);
  int maxima_group_n = 0;

  for (int x = 0; x < tx; x++) {
    for (int y = 0; y < ty; y++) {
      if (group_n(x, y) == -1) {
	queue_x.push(x);
	queue_y.push(y);
	maxima = true; cur_value = J(x, y);
	group_n(x, y) = group;
	maxima_group[group] = false;

	while (!queue_x.empty()) {
	  cur_x = queue_x.front();
	  queue_x.pop();
	  cur_y = queue_y.front();
	  queue_y.pop();

	  for (int k = 0; k < 8; k++) {
	    next_x = cur_x + dx[k];
	    next_y = cur_y + dy[k];
	    if ((next_x >= 0) && (next_x < tx) && (next_y >= 0) && (next_y < ty)) {
	      if (J(next_x, next_y) == cur_value) {
		if (group_n(next_x, next_y) == -1) {
		  group_n(next_x, next_y) = group;
		  queue_x.push(next_x);
		  queue_y.push(next_y);
		}
	      } else if (J(next_x, next_y) > cur_value) {
		maxima = false;
	      }
	    }
	  }
	}

	maxima_group[group] = maxima;
	if (maxima)
	  maxima_group_n++;
	group++;
      }	
    }
  }

  vector<int> maxima_correspondence(group);
  maxima_n = 1;
  for (int i = 0; i < group; i++)
    if (maxima_group[i]) {
      maxima_correspondence[i] = maxima_n++;
    } else {
      maxima_correspondence[i] = 0;
    }

  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      labels(x, y) = maxima_correspondence[group_n(x, y)];
    }

#if (DEBUG_MODE==DETAIL)
  cout << "Groups detected after reconstruction: " << group << endl;
  cout << "Local maxima groups: " << maxima_group_n << endl;
#endif

  cout << "morphological operation ended" << endl; 
}

class EvalRes {
public:
  vector<int> vgt;
  vector<int> TP, FP, FN;
  vector<double> P, R, F;
  double Area_PR;
};

class EX05_result {
public:
  EvalRes evalRes;
  matrix<int> confidence;
  matrix<int> particle_centers;
  matrix<> detections;
  double rangeDetMin, rangeDetMax;
  int particle_radius, window_radius, thresh_Hmax;

  svm_model* model;
  vector<double> thresh;

  EX05_result() : confidence(1,1), particle_centers(2,2), detections(1,1) {
  }

  void save(string path, string filename, int particle_r, int window_r, int thresh_H) {
    FILE *fout;
    fout = fopen((path + filename + ".dat").c_str(), "wt");

    //TODO: save everything? (e.g. evalRes)
    fprintf(fout, "%d %d %d\n", particle_r, window_r, thresh_H);
    fprintf(fout, "%.5f\n", evalRes.Area_PR);
    fprintf(fout, "%.5f %.5f\n", rangeDetMin, rangeDetMax);
    fprintf(fout, "%lu\n", thresh.size());
    for (unsigned int i = 0; i < thresh.size(); i++) {
      fprintf(fout, "%.5f ", thresh[i]);
    }
    fprintf(fout, "\n");

    fclose(fout);

    svm_save_model((path + filename + ".svm").c_str(), model);
  }

  void load(string path, string filename) {
    FILE *fin;
    fin = fopen((path + filename + ".dat").c_str(), "rt");

    //TODO: load everything?
    fscanf(fin, "%d %d %d\n", &particle_radius, &window_radius, &thresh_Hmax);
    fscanf(fin, "%lf\n", &(evalRes.Area_PR));
    fscanf(fin, "%lf %lf\n", &rangeDetMin, &rangeDetMax);

    int n;
    fscanf(fin, "%d\n", &n);
    double tmp;
    thresh.resize(n);
    for (int i = 0; i < n; i++) {
      fscanf(fin, "%lf", &tmp);
      thresh[i] = tmp;
    }

    fclose(fin);

    model = svm_load_model((path + filename + ".svm").c_str());
  }
}; 

class EX07_result {
public:
  EvalRes evalRes;
  matrix<int> confidence;
  matrix<int> particle_centers;
  matrix<> detections;
  double rangeDetMin, rangeDetMax;

  EX07_result() : confidence(1,1), particle_centers(2,2), detections(1,1) {
  }

  void save_perf(string path, string filename) {
    FILE *fout;
    fout = fopen((path+ filename + ".perf").c_str(), "wt");

    fprintf(fout, "%d\n", 255);
    for (int i = 1; i <= 255; i++) {
      fprintf(fout, "%f\n", evalRes.P[i]);
    }
    for (int i = 1; i <= 255; i++) {
      fprintf(fout, "%f\n", evalRes.R[i]);
    }

    fclose(fout);
  }

  void save(string path, string filename, const double min, const double max, const int particle_radius) { 
    FILE *fout;

    fout = fopen((path + filename).c_str(), "wt");

    //TODO: sort?
    fprintf(fout, "%lu %d\n", particle_centers.size(0), particle_radius*2);

    for (unsigned int i = 0; i < particle_centers.size(0); i++) {
      //TODO: coordinate
      //fprintf(fout, "%d %d %.2f\n", particle_centers(i, 0)+1, particle_centers(i, 1)+1, ((particle_centers(i, 2) / 255.0 * (max - min) + min) - rangeDetMin) / (rangeDetMax - rangeDetMin)*100);
      fprintf(fout, "%d %d %.2f\n", particle_centers(i, 1)+1, particle_centers(i, 0)+1, ((particle_centers(i, 2) / 255.0 * (max - min) + min) - rangeDetMin) / (rangeDetMax - rangeDetMin)*100);
      //fprintf(fout, "%d %d %.2f\n",
      //(particle_centers(i, 1)+1)*rescale_factor, 
      //(particle_centers(i, 0)+1)*rescale_factor,
      //((particle_centers(i, 2) / 255.0 * (max - min) + min) - rangeDetMin) / (rangeDetMax - rangeDetMin)*100);
    }

    fclose(fout);
  }
}; 


void calculate_PR_pcap(const matrix<int>& det_centers, const vector<Center>& GT_centers, const int& radius, const int& tx, const int& ty, EvalRes& evalRes) {
  cout << "calculate_PR_pcap started" << endl;

  // Line 10 - 19
  matrix<int> lab2indcent(tx, ty);
  matrix<bool> gt_balls(tx, ty);

  //#if (DEBUG_MODE==DETAIL)
#if (DEBUG_MODE)
  cout << tx << " " << ty << endl;
  cout << det_centers.size(0) << " " << GT_centers.size() << endl;
  cout << radius << endl;
#endif

  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      lab2indcent(x, y) = -1;
      gt_balls(x, y) = false;
    }

#if DEBUG_MODE
  cout << "step a" << endl;
#endif

  for (int mx = -radius; mx < radius; mx++) {
    for (int my = -radius; my < radius; my++) {
      if (round(sqrt(mx*mx + my*my)) > radius)
	continue;
      for (unsigned int i = 0; i < det_centers.size(0); i++) {
	int x = det_centers(i, 0);
	int y = det_centers(i, 1);
	if ((x + mx < 0) || (x + mx >= tx) || (y + my < 0) || (y + my >= ty))
	  continue;
	if ((lab2indcent(x+mx, y+my) > -1) && (det_centers(lab2indcent(x+mx, y+my), 2) > det_centers(i, 2)))
	  continue;
	else
	  lab2indcent(x+mx, y+my) = i;
      }
    }
  }

#if DEBUG_MODE
  cout << "step a" << endl;
#endif

  // Line 22-43  
  vector<int> &vgt = evalRes.vgt;
  vgt.resize(GT_centers.size());
  for (unsigned int p = 0; p < GT_centers.size(); p++) {
#if DEBUG_MODE
    cout << "step a-1" << endl;
#endif

    int x = GT_centers[p].GetX();
    int y = GT_centers[p].GetY();

#if DEBUG_MODE
    cout << x << " " << y << endl;
#endif

    if (lab2indcent(x, y) > -1)
      vgt[p] = det_centers(lab2indcent(x,y), 2);
    else
      vgt[p] = 0;

#if DEBUG_MODE
    cout << "step a-2" << endl;
#endif

    for (int mx = -radius; mx < radius; mx++) {
      for (int my = -radius; my < radius; my++) {
	if (round(sqrt(mx*mx + my*my)) > radius)
	  continue;
	if ((x + mx < 0) || (x + mx >= tx) || (y + my < 0) || (y + my >= ty))
	  continue;
	gt_balls(x+mx,y+my) = true;
      }
    }

#if DEBUG_MODE
    cout << "step a-3" << endl;
#endif
  }

#if DEBUG_MODE
  cout << "step b" << endl;
#endif

  vector<int> vfp(det_centers.size(0));
  for (unsigned int l = 0; l < vfp.size(); l++) {
    if (!gt_balls(det_centers(l,0), det_centers(l, 1)))
      vfp[l] = det_centers(l, 2);
    else
      vfp[l] = 0;
  }

#if DEBUG_MODE
  cout << "step c" << endl;
#endif

 // Line 48-
  vector<int> &TP = evalRes.TP, &FP = evalRes.FP, &FN = evalRes.FN, temp, temp2;
  temp.resize(256);
  temp2.resize(256);
  TP.resize(256);
  FP.resize(256);
  FN.resize(256);

  for (int l = 0; l <= 255; l++) {
    temp[l] = 0;
    temp2[l] = 0;
  }

  int totally_missed = 0;
  for (unsigned int l = 0; l < vgt.size(); l++) {
    if (vgt[l] <= 0) {
      totally_missed++;
      continue;
    }
    if (vgt[l] > 255)
      temp[255]++;
    else
      temp[vgt[l]]++;
  }
  for (unsigned int l = 0; l < vfp.size(); l++) {
    if (vfp[l] <= 0)
      continue;
    if (vfp[l] > 250)
      temp2[255]++; 
    else
      temp2[vfp[l]]++;
  }

  FN[1] = totally_missed;
  TP[255] = temp[255];
  FP[255] = temp2[255];
  for (int l = 2; l <= 255; l++) {
    FN[l] = FN[l-1] + temp[l-1];
    TP[256-l] = TP[257-l] + temp[256-l];
    FP[256-l] = FP[257-l] + temp2[256-l];
  }

#if DEBUG_MODE
  cout << "step d" << endl;
#endif


  double P2[257], R2[257];
  vector<double> &P = evalRes.P, &R = evalRes.R, &F = evalRes.F;
  P.resize(256);
  R.resize(256);
  F.resize(256);
  for (int l = 1; l <= 255; l++) {
    P[l] = TP[l] * 1.0 / (TP[l] + FP[l]);
    R[l] = TP[l] * 1.0 / (TP[l] + FN[l]); 
    F[l] = 2 * P[l] * R[l] / (P[l] + R[l]);
    printf("%.3f %.3f\n", P[l], R[l]);

    P2[l] = P[l];
    R2[l] = R[l];
  }

  // TODO: More stuff (e.g. bestF, Area_PR, etc)
  for (int i = 1; i <= 255; i++)
    for (int j = i+1; j <= 255; j++) {
      if (R2[i] > R2[j]) {
	double tmp_swap = R2[i];
	R2[i] = R2[j];
	R2[j] = tmp_swap;

	tmp_swap = P2[i];
	P2[i] = P2[j];
	P2[j] = tmp_swap;
      }
    }

  P2[0] = 0;
  P2[256] = 0;
  R2[0] = 0;
  R2[256] = 1;
  for (int i = 255; i >= 0; i--) {
    if (P2[i+1]>P2[i])
      P2[i] = P2[i+1];
  }
  double &Area_PR = evalRes.Area_PR;
  Area_PR = 0;
  for (int i = 1; i <= 256; i++) {
    Area_PR += (R2[i] - R2[i-1])*P2[i];
  }

  cout << "calculate_PR_pcap ended" << endl;
}

void extract_particles(const string name, const int tx, const int ty, const int group_n, const matrix<int>& labels, const matrix<>& detections, matrix<int>& particle_centers ) {
  cout << "extract_particle started" << endl;

  matrix<int> tmp_centers;

  tmp_centers.resize(static_cast<unsigned long>(group_n), static_cast<unsigned long>(4));
  for (unsigned int i = 0; i < tmp_centers.size(0); i++)
    tmp_centers(i, 2) = -1;

  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      if (labels(x,y) == 0)
	continue;

      if (detections(x, y) > tmp_centers(labels(x, y), 2)) {
	tmp_centers(labels(x,y), 0) = x;
	tmp_centers(labels(x,y), 1) = y;
	tmp_centers(labels(x,y), 2) = detections(x,y);
	tmp_centers(labels(x,y), 3) = labels(x,y);   // TODO: necessary?
      }	
    }

  // TODO: rather than -5, you should use particle radius probably...

  Image img;
  img.LoadMask(name);

  // prune them if they are out-of-boundary
  int cou = 0;
  for (unsigned int i = 0; i < tmp_centers.size(0); i++) {
    if ((tmp_centers(i,0) < 5) || (tmp_centers(i,0) > tx - 5) ||
	(tmp_centers(i,1) < 5) || (tmp_centers(i,1) > ty - 5) ||
	(!img.GetMask(tmp_centers(i,0),tmp_centers(i,1)))) {
    } else {
      cou++;
    }
  }

  particle_centers.resize(cou, static_cast<unsigned long>(4));
  cou = 0;
  for (unsigned int i = 0; i < tmp_centers.size(0); i++) {
    if ((!img.GetMask(tmp_centers(i,0),tmp_centers(i,1))) ||
	(tmp_centers(i,0) < 5) || (tmp_centers(i,0) > tx - 5) ||
	(tmp_centers(i,1) < 5) || (tmp_centers(i,1) > ty - 5)) {
    } else {
      for (int k = 0; k < 4; k++)
	particle_centers(cou, k) = tmp_centers(i, k);
      cou++;
    }
  }

  cout << "extract_particle ended" << endl;
}

void TEXTONSVM_batch_eval(const SVM_PARAM& param, EX05_result& ex05_result) {
  EvalRes& evalRes = ex05_result.evalRes;
  matrix<>& detections = ex05_result.detections;
 
  pieces2det(param.pieces_detDir, param.window_radius+param.texton_radius, detections);
  double &range_min = ex05_result.rangeDetMin, &range_max = ex05_result.rangeDetMax;
  range_min = range_max = detections(0, 0);
  int tx = detections.size(0), ty = detections.size(1);
  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      if (range_min > detections(x, y))
	range_min = detections(x, y);
      if (range_max < detections(x, y))
	range_max = detections(x, y);
    }
  im_normalize(detections, range_min, range_max);

  // Load ground truth
  Image tmp;
  tmp.SetResizeFactor(rescale_factor);
  // TODO: resultsDir name misleading
  // TODO: have to compare with mask.. probably reimplementation will be better here
  // (also, loadcenters is better to be reimplemented)
  tmp.LoadCenters(IMG_DIR + param.dataFile);
  vector<Center> centers = tmp.GetCenters();
  vector<Center> GT_centers;    
  for (unsigned int i = 0; i < centers.size(); i++) {
    // TODO: MASK
    if (centers[i].GetX() <= 233)
      continue;
    GT_centers.push_back(centers[i]);
  }
  
  matrix<int> labels;
  int group_n;
  morphological_op(detections, param, labels, group_n);

  matrix<int> &particle_centers = ex05_result.particle_centers;
  extract_particles(IMG_DIR + param.dataFile, tx, ty, group_n, labels, detections, particle_centers);

  /*
  particle_centers.resize(static_cast<unsigned long>(group_n), static_cast<unsigned long>(4));
  for (unsigned int i = 0; i < particle_centers.size(0); i++)
    particle_centers(i, 2) = -1;

  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      if (labels(x,y) == 0)
	continue;

      if (detections(x, y) > particle_centers(labels(x, y), 2)) {
	particle_centers(labels(x,y), 0) = x;
	particle_centers(labels(x,y), 1) = y;
	particle_centers(labels(x,y), 2) = static_cast<int>(detections(x,y));
	particle_centers(labels(x,y), 3) = labels(x,y);   // TODO: necessary?
      }	
    }
  */

  matrix<int>& confidence = ex05_result.confidence;
  confidence.resize(static_cast<unsigned long>(tx), static_cast<unsigned long>(ty));
  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      if ((labels(x,y) > 0) && (labels(x,y) < particle_centers.size(0)))
	confidence(x, y) = particle_centers(labels(x,y), 2);
      else
	confidence(x, y) = 0;
    }

  // TOWORK
  // maybe I should just get group_n rather than labels
  // read TEXTONSVM_batch_eval before you proceed further

  calculate_PR_pcap(particle_centers, GT_centers, param.particle_radius, tx, ty, evalRes);
}

void TEXTONSVM_batch_svm_opt(const SVM_PARAM& param, EX05_result& ex05_result) {
  string touchFileName = param.batchDir + "dict.touch";
  if (!exist_file(touchFileName)) {
    TEXTONSVM_batch_dict(param);
    toc();
  }

  // load dictionary
  DICT dictionary;
  dictionary.load(RESULT_DIR, param.modelsFile);

  {
    matrix<> img;
    load_denoise_image(param.dataFile, RESULT_DIR, img);
    orig2pieces(img, param.pieces_origDir, param.tx_out, param.ty_out, param.overlap, param.x_init, param.y_init);
  }
  
  touchFileName = param.batchDir + "tmap.touch";
  if (!exist_file(touchFileName)) {
    DIR *dirp = opendir(param.pieces_origDir.c_str());
    if (dirp) {
      struct dirent *entry;
      while ((entry = readdir(dirp))) {
	if (entry->d_type == 8) {
	  string pieceName = string(entry->d_name);
	  if (strcmp(pieceName.substr(pieceName.size()-3, 3).c_str(), "img") != 0)
	    continue;

	  create_piece_texton_map(dictionary, param, pieceName, param.pieces_origDir, param.pieces_clsDir);
	}
      }
      closedir(dirp);
    } else {
      cout << "Error to open directory:"
	   << param.pieces_origDir << endl;
      exit(1);
    }

    matrix<> texton_map;
    pieces2texton_map(param.pieces_clsDir, param.texton_radius, texton_map);

    /*
    {
      // texton_map2pieces
      texton_map2pieces(texton_map, param.pieces_clsDir, param.tx_out, param.ty_out, param.overlap, param.x_init, param.y_init);
    }
    */

    save_image(texton_map, RESULT_DIR, "texton_map.img");

    touchFile(param.batchDir, "tmap.touch");
  }



  // update pieces and model
  touchFileName = param.batchDir + "update.touch";
  // TODO: should remove
  //  if (!exist_file(touchFileName)) {
    TEXTONSVM_batch_update(param, ex05_result.model, ex05_result.thresh);
    //    touchFile(param.batchDir, "update.touch");
    //  } 


  // detection
  touchFileName = param.batchDir + "det.touch";
  if (!exist_file(touchFileName)) {
    DIR *dirp = opendir(param.pieces_origDir.c_str());
    if (dirp) {
      struct dirent *entry;
      while ((entry = readdir(dirp))) {
	if (entry->d_type == 8) {
	  string pieceName = string(entry->d_name);
	  if (strcmp(pieceName.substr(pieceName.size()-3, 3).c_str(), "img") != 0)
	    continue;
	  create_piece_det(pieceName, param, ex05_result.model, ex05_result.thresh);
	}
      }
    }
    touchFile(param.batchDir, "det.touch");
  }

  // PR
  touchFileName = param.batchDir + "eval.touch";
  if (!exist_file(touchFileName)) {
    TEXTONSVM_batch_eval(param, ex05_result);
    touchFile(param.batchDir, "eval.touch");
  }

  // TODO: clean up batch
  system(("rm -rRf " + param.batchDir + "det.touch").c_str());
  system(("rm -rRf " + param.batchDir + "update.touch").c_str());
  system(("rm -rRf " + param.batchDir + "tmap.touch").c_str());
  system(("rm -rRf " + param.batchDir + "eval.touch").c_str());
  system(("rm -rRf " + param.piecesDir).c_str());
}

void ex05() {
  double bestArea = 0.0;

  if (SEARCH_OPTIMAL) {
    // TODO: what should I do if search_optimal is true and you have some default particel_radius
    for (int particle_radius = default_particle_radius - 3; particle_radius <= default_particle_radius + 3; particle_radius++) {
      //    for (int particle_radius = 4; particle_radius <= 10; particle_radius++) {
      for (int window_radius = particle_radius + 1; window_radius <= particle_radius + 10; window_radius++) {
      //      for (int window_radius = particle_radius + 8; window_radius <= particle_radius + 10; window_radius++) {
	for (int thresh_Hmax = 50; thresh_Hmax <= 50; thresh_Hmax += 10) {
	  //	particle_radius = 7;
	  //	window_radius = 11;
	  //	particle_radius = 8;
	  //	window_radius = 18;
#if DEBUG_MODE
	  printf("\n-------------Setup %d %d %d----\n\n", particle_radius, window_radius, thresh_Hmax);
#endif
	  
	  EX05_result result;
	  SVM_PARAM param(particle_radius, window_radius, thresh_Hmax); 
	  
	  TEXTONSVM_batch_svm_opt(param, result);
	  
	  cout << "Area_PR: " << result.evalRes.Area_PR << endl;
	  
	  char buffer[50];
	  sprintf(buffer,"%d_%d_%d_train", particle_radius, window_radius, thresh_Hmax);
	  result.save(RESULT_DIR, buffer, particle_radius, window_radius, thresh_Hmax);
	  
	  if (bestArea < result.evalRes.Area_PR) {
	    bestArea = result.evalRes.Area_PR;
	    result.save(RESULT_DIR, "FINAL_TRAIN_BEST", particle_radius, window_radius, thresh_Hmax);
	  }
	  
	  /*
	    {
	    result.save(RESULT_DIR, "FINAL_TRAIN_BEST");
	    return ;
	    }
	  */
	  
	  svm_destroy_model(result.model);
	  toc();
	  
	  //	return ;
	}
      }
    }
  } else {
    // TODO: is this all we need to change?
    // TODO: how about for generating texton?
    EX05_result result;
    SVM_PARAM param(default_particle_radius, default_window_radius, default_thresh_Hmax); 
#if DEBUG_MODE
    printf("\n-------------Setup %d %d %d----\n\n", default_particle_radius, default_window_radius, default_thresh_Hmax);
#endif
    
    TEXTONSVM_batch_svm_opt(param, result);
    
    cout << "Area_PR: " << result.evalRes.Area_PR << endl;
    
    char buffer[50];
    sprintf(buffer,"%d_%d_%d_train", default_particle_radius, default_window_radius, default_thresh_Hmax);
    result.save(RESULT_DIR, buffer, default_particle_radius, default_window_radius, default_thresh_Hmax);
    
    bestArea = result.evalRes.Area_PR;
    result.save(RESULT_DIR, "FINAL_TRAIN_BEST", default_particle_radius, default_window_radius, default_thresh_Hmax);
	  
    svm_destroy_model(result.model);
    toc();
  }

  cout << "ex05: done." << endl; 
}

void TEXTONSVM_batch_eval_all(const string& SVM_ALL_DIR, const SVM_PARAM_ALL& param, const EX05_result& ex05_result, EX07_result& ex07_result) {
  EvalRes& evalRes = ex07_result.evalRes;
  //  matrix<>& detections = ex05_result.detections;
  
  matrix<> detections;
  pieces2det(param.pieces_detDir, param.window_radius+param.texton_radius, detections);
  //  save_image(detections, resultsDir, "detection_original.img");

  double &range_min = ex07_result.rangeDetMin, &range_max = ex07_result.rangeDetMax;
  range_min = range_max = detections(0, 0);
  int tx = detections.size(0), ty = detections.size(1);
  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      if (range_min > detections(x, y))
	range_min = detections(x, y);
      if (range_max < detections(x, y))
	range_max = detections(x, y);
    }
  im_normalize(detections, ex05_result.rangeDetMin, ex05_result.rangeDetMax);

  //  save_image(detections, SVM_ALL_DIR, "pieces2det_result.img");
  save_image(detections, SVM_ALL_DIR, (param.imgName+".dmap").c_str()); 

  //TODO:
  // Load ground truth
  Image tmp; 
  tmp.SetResizeFactor(rescale_factor);
  // TODO: resultsDir name misleading
  // TODO: have to compare with mask.. probably reimplementation will be better here
  // (also, loadcenters is better to be reimplemented)
  tmp.LoadCenters(IMG_DIR + param.imgName);
  vector<Center> centers = tmp.GetCenters();
  vector<Center> GT_centers;    
  for (unsigned int i = 0; i < centers.size(); i++) {
    // TODO: MASK
    if (centers[i].GetX() <= 233)
      continue;
    GT_centers.push_back(centers[i]);
  }
  
  matrix<int> labels;
  int group_n;
  morphological_op(detections, param, labels, group_n);

  //  save_image(labels, SVM_ALL_DIR, "morphological_op.img");
  matrix<int> &particle_centers = ex07_result.particle_centers;
  extract_particles(IMG_DIR + param.imgName, tx, ty, group_n, labels, detections, particle_centers);

  /*
  // TODO:
  particle_centers.resize(static_cast<unsigned long>(group_n), static_cast<unsigned long>(4));
  for (unsigned int i = 0; i < particle_centers.size(0); i++)
    particle_centers(i, 2) = -1;

  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      if (labels(x,y) == 0)
	continue;

      if (detections(x, y) > particle_centers(labels(x, y), 2)) {
	particle_centers(labels(x,y), 0) = x;
	particle_centers(labels(x,y), 1) = y;
	particle_centers(labels(x,y), 2) = detections(x,y);
	particle_centers(labels(x,y), 3) = labels(x,y);   // TODO: necessary?
      }	
    }
  */
  
  matrix<int>& confidence = ex07_result.confidence;
  confidence.resize(static_cast<unsigned long>(tx), static_cast<unsigned long>(ty));
  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      if ((labels(x,y) > 0) && (labels(x,y) < particle_centers.size(0)))
	confidence(x, y) = particle_centers(labels(x,y), 2);
      else
	confidence(x, y) = 0;
    }
  save_image(confidence, SVM_ALL_DIR, (param.imgName+".confidence").c_str());

  ex07_result.save(SVM_ALL_DIR, (param.imgName+".det_cen").c_str(), ex05_result.rangeDetMin, ex05_result.rangeDetMax, param.particle_radius);

  // TOWORK
  // maybe I should just get group_n rather than labels
  // read TEXTONSVM_batch_eval before you proceed further


  // TODO: is PR necessary?
  calculate_PR_pcap(particle_centers, GT_centers, param.particle_radius, tx, ty, evalRes);
  ex07_result.save_perf(SVM_ALL_DIR, param.imgName);
}

void TEXTONSVM_batch_svm_all(const string& SVM_ALL_DIR, const SVM_PARAM_ALL& param, EX07_result& ex07_result) {
  matrix<> img;
  load_denoise_image(param.dataFile, img);
  orig2pieces(img, param.pieces_origDir, param.tx_out, param.ty_out, param.overlap, param.x_init, param.y_init);

  EX05_result ex05_result;
  // TODO: what if FINAL_TRAIN_BEST is not there? In fact, that means something is wrong.
  //  ex05_result.load(RESULT_DIR, "CUR_TRAIN_BEST");
  ex05_result.load(RESULT_DIR, "FINAL_TRAIN_BEST");


  // load dictionary
  DICT dictionary;
  dictionary.load(RESULT_DIR, param.modelsFile);

  string touchFileName = param.batchDir + "tmap.touch";
  if (!exist_file(touchFileName)) {
    DIR *dirp = opendir(param.pieces_origDir.c_str());
    if (dirp) {
      struct dirent *entry;
      while ((entry = readdir(dirp))) {
	if (entry->d_type == 8) {
	  string pieceName = string(entry->d_name);
	  if (strcmp(pieceName.substr(pieceName.size()-3, 3).c_str(), "img") != 0)
	    continue;

	  create_piece_texton_map(dictionary, param, pieceName, param.pieces_origDir, param.pieces_clsDir);
	}
      }
      closedir(dirp);
    } else {
      cout << "Error to open directory:"
	   << param.pieces_origDir << endl;
      exit(1);
    }

    matrix<> texton_map;
    pieces2texton_map(param.pieces_clsDir, param.texton_radius, texton_map);

    save_image(texton_map, SVM_ALL_DIR, "texton_map.img");

    /*
    {
      // texton_map2pieces
      texton_map2pieces(texton_map, param.pieces_clsDir, param.tx_out, param.ty_out, param.overlap, param.x_init, param.y_init);
    }
    */

    touchFile(param.batchDir, "tmap.touch");
  }

  // detection
  touchFileName = param.batchDir + "det.touch";
  if (!exist_file(touchFileName)) {
    DIR *dirp = opendir(param.pieces_origDir.c_str());
    if (dirp) {
      struct dirent *entry;
      while ((entry = readdir(dirp))) {
	if (entry->d_type == 8) {
	  string pieceName = string(entry->d_name);
	  if (strcmp(pieceName.substr(pieceName.size()-3, 3).c_str(), "img") != 0)
	    continue;
	  create_piece_det(pieceName, param, ex05_result.model, ex05_result.thresh);
	}
      }
    }
    touchFile(param.batchDir, "det.touch");
  }


  // PR
  touchFileName = param.batchDir + "eval.touch";
  if (!exist_file(touchFileName)) {
    TEXTONSVM_batch_eval_all(SVM_ALL_DIR, param, ex05_result, ex07_result);
    //    touchFile(param.batchDir, "eval.touch");
  }

  //TODO: remove files
  system(("rm -rRf " + param.piecesDir).c_str());
  system(("rm -rRf " + param.batchDir).c_str());
}

void ex07(vector<string> filenames) {
  int texton_radius = TEXTON_RADIUS;
  int param_K = K;
  int texton_nbins = TEXTON_NBINS;
  int grayhist_nbins = GRAYHIST_NBINS;
  int nbcls = NBCLS;
  double weight_feature = WEIGHT_FEATURE;

  // TODO: should I store default param on ex05 and copy to ex07??

  EX05_result ex05_result;
  // TODO: what if FINAL_TRAIN_BEST is not there? In fact, that means something is wrong.
  //  ex05_result.load(RESULT_DIR, "CUR_TRAIN_BEST");
  ex05_result.load(RESULT_DIR, "FINAL_TRAIN_BEST");
  int particle_radius = ex05_result.particle_radius;
  int window_radius = ex05_result.window_radius;
  int thresh_Hmax = ex05_result.thresh_Hmax;
#if (DEBUG_MODE==DETAIL)
  cout << "p,w,t: " << particle_radius << " " << window_radius << " " << thresh_Hmax << endl;
#endif
  
  for (unsigned int i = 0; i < filenames.size(); i++) {
    cout << filenames[i] << endl;

    string imgName = filenames[i];
    //string resultsDir = RESULT_DIR + imgName + "/";
    string resultsDir = RESULT_DIR + "results/";
    mkdir(resultsDir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

    SVM_PARAM_ALL param(texton_radius, texton_nbins, grayhist_nbins, weight_feature, particle_radius, window_radius, param_K, nbcls, thresh_Hmax, imgName);

    EX07_result ex07_result;
    TEXTONSVM_batch_svm_all(resultsDir, param, ex07_result);
    // TODO:

    // TODO: should I move or copy?
    string cmd;
    cmd = "cp " + RESULT_DIR + "/" + imgName + ".tsvm " + resultsDir + ".";
    system(cmd.c_str());
    cmd = "cp " + RESULT_DIR + "/" + imgName + ".box " + resultsDir + ".";
    system(cmd.c_str());
    cmd = "cp " + RESULT_DIR + "/" + imgName + ".norm " + resultsDir + ".";
    system(cmd.c_str());

    cout << endl;
    toc();
    cout << endl << endl << endl;
  }

  // post processing (.box and overall.perf)
  cout << "Post processing started.." << endl;

  int maxX = -1, maxY = -1;
  vector<Center> overall_gt_center;
  matrix<int> overall_det_center;
  unsigned int det_center_size[filenames.size()];
  unsigned int overall_det_cou = 0;
#if DEBUG_MODE
  cout << "step a" << endl;;
#endif
  for (unsigned int i = 0; i < filenames.size(); i++) {
    string resultsDir = RESULT_DIR + "results/";

    FILE *fin;
    int tmp_size;

    fin = fopen((resultsDir+filenames[i]+".det_cen").c_str(), "rt");
    fscanf(fin, "%d", &tmp_size);
    fclose(fin);
    det_center_size[i] = tmp_size;
    overall_det_cou += tmp_size;
  }
#if DEBUG_MODE
  cout << "step b" << endl;;
#endif
  overall_det_center.resize(static_cast<unsigned long>(overall_det_cou), static_cast<unsigned long>(3));
  overall_det_cou = 0;
  for (unsigned int i = 0; i < filenames.size(); i++) {
    string resultsDir = RESULT_DIR + "results/";

    FILE *fin;

#if DEBUG_MODE
    cout << "step c" << endl;
#endif
    {
      int tmp_size, tmp_radius;
      fin = fopen((resultsDir+filenames[i]+".det_cen").c_str(), "rt");
      fscanf(fin, "%d %d", &tmp_size, &tmp_radius);
      for (int j = 0; j < tmp_size; j++) {
#if DEBUG_MODE
	cout << j << " / " << tmp_size << endl;
#endif
	int tmp1, tmp2;
	double tmp3;
	fscanf(fin, "%d %d %lf", &tmp2, &tmp1, &tmp3);
	overall_det_center(overall_det_cou, 0) = tmp1;
	overall_det_center(overall_det_cou, 1) = tmp2;
	if (maxX < tmp1)
	  maxX = tmp1;
	if (maxY < tmp2)
	  maxY = tmp2;
	overall_det_center(overall_det_cou, 2) = static_cast<int>(tmp3/100*255);
	overall_det_cou++;
      }
      fclose(fin);
    }

#if DEBUG_MODE
    cout << "step d" << endl;
#endif
    {
      fin = fopen((resultsDir+filenames[i]+".box").c_str(), "rt");
      int tmp1, tmp2, tmp3, tmp4, garbage;
      while (fscanf(fin, "%d %d %d %d %d\n", &tmp1, &tmp2, &tmp3, &tmp4, &garbage)==5) {
	Center tmp_center(tmp2+tmp4/2, tmp1+tmp3/2);
	if (maxX < tmp2+tmp4/2)
	  maxX = tmp2+tmp4/2;
	if (maxY < tmp1+tmp3/2)
	  maxY = tmp1+tmp3/2;
	overall_gt_center.push_back(tmp_center);
      }
      fclose(fin);
    }   
  }
  EvalRes overall_evalRes;
  calculate_PR_pcap(overall_det_center, overall_gt_center, particle_radius, maxX+1, maxY+1, overall_evalRes);

  int maxF_ind = -1;
  double maxF = -999999;
  for (int i = 0; i < 255; i++) {
#if DEBUG_MODE
    cout << i << endl;
#endif
    double F = 2*(overall_evalRes.P[i]*overall_evalRes.R[i])/(overall_evalRes.P[i] + overall_evalRes.R[i]);
    if (F > maxF) {
      maxF = F;
      maxF_ind = i;
    }
  }
  cout << "maxF @ " << maxF_ind << endl;

  for (unsigned int i = 0; i < filenames.size(); i++) {
    string resultsDir = RESULT_DIR + "results/";

    FILE *fin;

    matrix<int> det_center;
    det_center.resize(static_cast<unsigned long>(det_center_size[i]), static_cast<unsigned long>(3));
    int det_cou = 0;
    {
      int tmp_size, tmp_radius;
      fin = fopen((resultsDir+filenames[i]+".det_cen").c_str(), "rt");
      fscanf(fin, "%d %d", &tmp_size, &tmp_radius);
      for (int j = 0; j < tmp_size; j++) {
	int tmp1, tmp2;
	double tmp3;
	fscanf(fin, "%d %d %lf", &tmp2, &tmp1, &tmp3);
	det_center(det_cou, 0) = tmp2;
	det_center(det_cou, 1) = tmp1;
	if (maxX < tmp1)
	  maxX = tmp1;
	if (maxY < tmp2)
	  maxY = tmp2;
	det_center(det_cou, 2) = static_cast<int>(tmp3/100*255);
	det_cou++;
      }
      fclose(fin);

      for (int s1 = 0; s1 < tmp_size; s1++) {
	for (int s2 = s1+1; s2 < tmp_size; s2++) {
	  if (det_center(s1, 2) < det_center(s2, 2)) {
	    for (int s3 = 0; s3 < 3; s3++) {
	      int swap = det_center(s1, s3);
	      det_center(s1, s3) = det_center(s2, s3);
	      det_center(s2, s3) = swap;
	    }
	  }
	}
      }
    }

    vector<Center> gt_center;
    {
      fin = fopen((resultsDir+filenames[i]+".box").c_str(), "rt");
      int tmp1, tmp2, tmp3, tmp4, garbage;
      while (fscanf(fin, "%d %d %d %d %d\n", &tmp1, &tmp2, &tmp3, &tmp4, &garbage)==5) {
	int a = (tmp1+tmp3/2)/rescale_factor;
	int b = (tmp2+tmp4/2)/rescale_factor;
	Center tmp_center(a, b);
	if (maxX < a)
	  maxX = a;
	if (maxY < b)
	  maxY = b;
	gt_center.push_back(tmp_center);
      }
      fclose(fin);

      /*
      fin = fopen((resultsDir+filenames[i]+".centers").c_str(), "rt");
      while (fin) {
	int tmp1, tmp2;
	fscanf(fin, "%d %d", &tmp1, &tmp2);
	if ((tmp1 == 0) && (tmp2 == 0)) {
	  break;
	}
	Center tmp_center(tmp1, tmp2);
	if (maxX < tmp1)
	  maxX = tmp1;
	if (maxY < tmp2)
	  maxY = tmp2;
	gt_center.push_back(tmp_center);
      }
      fclose(fin);
      */
    }   

    EvalRes tmp_evalRes;
    calculate_PR_pcap(det_center, gt_center, particle_radius, maxX+1, maxY+1, tmp_evalRes);

    FILE *fout = fopen((resultsDir+filenames[i]+".det").c_str(), "w");

    int SOMERADIUS = particle_radius*4;
    for (unsigned int i = 0; i < det_center.size(0); i++) {
      int a,b,c,d;
      a = det_center(i,1)-SOMERADIUS/2;
      b = det_center(i,0)-SOMERADIUS/2;
      c = SOMERADIUS;
      d = SOMERADIUS;
      //fprintf(fout, "%d %d %d %d -3\n", a*rescale_factor, b*rescale_factor, c*rescale_factor, d*rescale_factor);
      //fprintf(fout, "%d %d %d %d -3\n", b, a, d, c);
      fprintf(fout, "%d %d %d %d -3\n", b*rescale_factor, a*rescale_factor, d*rescale_factor, c*rescale_factor);
    }
    for (unsigned int i = 0; i < tmp_evalRes.vgt.size(); i++) {
      //if (tmp_evalRes.vgt[i] < maxF_ind) {
      if (tmp_evalRes.vgt[i] < 1) {
	int a,b,c,d;
	a = gt_center[i].GetX()-SOMERADIUS/2;
	b = gt_center[i].GetY()-SOMERADIUS/2;
	c = SOMERADIUS;
	d = SOMERADIUS;
	//fprintf(fout, "%d %d %d %d -3\n", a*rescale_factor, b*rescale_factor, c*rescale_factor, d*rescale_factor);
	//fprintf(fout, "%d %d %d %d -3\n", b, a, d, c);
	//fprintf(fout, "%d %d %d %d -3\n", a, b, c, d);
	fprintf(fout, "%d %d %d %d %d\n", a*rescale_factor, b*rescale_factor, c*rescale_factor, d*rescale_factor, -3); //tmp_evalRes.vgt[i]);
      }
    }

    fclose(fout);
  }

  cout << "ex07: done." << endl;
}

int main(int argc, const char* argv[]) {
  string work_dir;

  safe_vector<matrix<>*> examplars_denoise;
  vector<string> image_names;
  vector<string> train_image_names;
  vector<string> test_image_names;
  
  // TODO: have to somehow load the list of files and centers

  random_seed_F = fopen("random.gen", "rt");

  //  image_names.push_back("001_Fer.tif");

  work_dir = "./input/";
  TEST_ONLY = false;

  if (argc < 2) {
    cout << "Usage: main [options]" << endl;
    cout << "Options:" << endl;
    cout << "  -d                         Default option" << endl;
    cout << "  -w                         Work directory" << endl;
    cout << "  -i N filename1 ... filenameN" << endl
	 << "                             Train images. Filenames should be without extensions" << endl
	 << "                               (either .tiff or ETC). Also, it requires ground" << endl
	 << "                               thruth particles (.centers) and mask files (.mask)" << endl;
    cout << "  -t N filename1 ... filenameN" << endl
	 << "                             Test images. Filenames should be without extensions" << endl
	 << "                                (either .tiff or ETC)." << endl;
    cout << "  -s filename                SVM filename without any extension" << endl;
    cout << "  -p #1 #2 #3                #1: Particle size, #2: Window radius," << endl
         << "                               #3: Flag to search for optimal [0: off, 1: on]" << endl;
    cout << "  -r [0-1]" << endl
	 << "                             Test only option. At least one training should have been" << endl
	 << "                               done with SVM_FILENAME. If this option is off, all" << endl
	 << "                               computations will be done again." << endl
	 << "                               [0: off (default), 1: on]" << endl;
    cout << endl;
    cout << "Report bugs to <lim@csail.mit.edu>" << endl;
    return 0;
  }

  // process input
  for (int i = 1;  i < argc; i++) {
#if (DEBUG_MODE==DETAIL)
    cout << argv[i] << endl;
#endif

    if (strcmp(argv[i], "-d") == 0) { // default option
      cout << "Default option" << endl;

      /*
      image_names.push_back("./input/001_Fer");
      image_names.push_back("./input/002_Fer");
      image_names.push_back("./input/003_Fer");
      image_names.push_back("./input/004_Fer");
      image_names.push_back("./input/005_Fer");
      */
      image_names.push_back("./input/006_Fer");
      /*
      image_names.push_back("./input/007_Fer");
      image_names.push_back("./input/008_Fer");
      image_names.push_back("./input/009_Fer");
      image_names.push_back("./input/010_Fer");
      image_names.push_back("./input/011_Fer");
      image_names.push_back("./input/012_Fer");
      image_names.push_back("./input/013_Fer");
      image_names.push_back("./input/014_Fer");
      image_names.push_back("./input/015_Fer");
      image_names.push_back("./input/016_Fer");
      image_names.push_back("./input/017_Fer");
      image_names.push_back("./input/018_Fer");
      image_names.push_back("./input/019_Fer");
      */
      train_image_names.push_back("./input/006_Fer");
      /*
      test_image_names.push_back("./input/001_Fer");
      test_image_names.push_back("./input/002_Fer");
      test_image_names.push_back("./input/003_Fer");
      test_image_names.push_back("./input/004_Fer");
      test_image_names.push_back("./input/005_Fer");
      test_image_names.push_back("./input/006_Fer");
      test_image_names.push_back("./input/007_Fer");
      test_image_names.push_back("./input/008_Fer");
      test_image_names.push_back("./input/009_Fer");
      test_image_names.push_back("./input/010_Fer");
      test_image_names.push_back("./input/011_Fer");
      test_image_names.push_back("./input/012_Fer");
      test_image_names.push_back("./input/013_Fer");
      test_image_names.push_back("./input/014_Fer");
      test_image_names.push_back("./input/015_Fer");
      test_image_names.push_back("./input/016_Fer");
      test_image_names.push_back("./input/017_Fer");
      test_image_names.push_back("./input/018_Fer");
      test_image_names.push_back("./input/019_Fer");
      */

      svm_model_name = "006_Fer_only";

      TEST_ONLY = 0;
    } else if (strcmp(argv[i], "-w") == 0) { // work_dir
      work_dir = argv[i+1];
      if (argv[i+1][strlen(argv[i+1])-1] != '/')
	work_dir = work_dir + '/';
      i++;
    } else if (strcmp(argv[i], "-i") == 0) { // train_images
      int n = atoi(argv[i+1]);
      for (int j = i+2; j < i+2+n; j++) {
	//	image_names.push_back(argv[j]);
	train_image_names.push_back(argv[j]); 
      }
      i+=(n+1);
    } else if (strcmp(argv[i], "-t") == 0) { // test images
      int n = atoi(argv[i+1]);
      cout << "train_n: " << n << endl;
      for (int j = i+2; j < i+2+n; j++) {
	//	image_names.push_back(argv[j]);
	test_image_names.push_back(argv[j]);
      }
      i+=(n+1);
    } else if (strcmp(argv[i], "-s") == 0) { // svm model filename
      svm_model_name = argv[i+1];
      // TODO: check extension of svm_model_name (it should be w/o).
      i++;
    } else if (strcmp(argv[i], "-r") == 0) {
      TEST_ONLY = argv[i+1][0] - '0';
      i++;
    } else if (strcmp(argv[i], "-p") == 0) {
      // TODO: rename this and make it global for ex05 to use
      // TODO: how about for generating texton?
      rescale_factor = atoi(argv[i+1]);
      default_particle_radius = atoi(argv[i+2]);
      default_window_radius = atoi(argv[i+3]);
      SEARCH_OPTIMAL = (argv[i+4][0] != '0');
      WINDOW_RADIUS = atoi(argv[i+5]);
      TEXTON_RADIUS = atoi(argv[i+6]);
      K = atoi(argv[i+7]);
      TEXTON_NBINS = atoi(argv[i+8]);
      GRAYHIST_NBINS = atoi(argv[i+9]);
      NBCLS = atoi(argv[i+10]);
      //WEIGHT_FEATURE = atof(argv[i+11]);
      i+=10;
    }
  }

  //  IMG_DIR = work_dir;
  //  OUTPUT_DIR = work_dir;
  RESULT_DIR = work_dir + svm_model_name + "/";
  IMG_DIR = RESULT_DIR;
  OUTPUT_DIR = RESULT_DIR;

  //  mkdir(work_dir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
  //  mkdir(OUTPUT_DIR.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
  mkdir(RESULT_DIR.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

#if (DEBUG_MODE==DETAIL)
  cout << "IMG_DIR: " << IMG_DIR << endl;
  cout << "OUTPUT_DIR: " << OUTPUT_DIR << endl;
  cout << "RESULT_DIR: " << RESULT_DIR << endl;
#endif

  string cmd;
  for (unsigned int i = 0; i < train_image_names.size(); i++) {
    //    cout << "Copying train files to the work dir [" << i+1 << "/" << train_image_names.size() << "]" << endl;
 
    char* copy_str = static_cast<char*>(malloc(train_image_names[i].size()+1));
    strcpy(copy_str, train_image_names[i].c_str());
    char* filename = basename(copy_str);
    //    cmd = "cp " + train_image_names[i] + ".pcap " + IMG_DIR + filename + ".pcap";
    //    system(cmd.c_str());
    //    cmd = "cp " + train_image_names[i] + ".good_centers " + IMG_DIR + filename + ".good_centers";
    //    system(cmd.c_str());
    //    cmd = "cp " + train_image_names[i] + ".centers " + IMG_DIR + filename + ".centers";
    //    system(cmd.c_str());

    train_image_names[i] = filename;
    image_names.push_back(filename);
    free(copy_str);
  }
  for (unsigned int i = 0; i < test_image_names.size(); i++) {
    //    cout << "Copying test files to the work dir [" << i+1 << "/" << test_image_names.size() << "]" << endl;
 
    char* copy_str = static_cast<char*>(malloc(test_image_names[i].size()+1));
    strcpy(copy_str, test_image_names[i].c_str());
    char* filename = basename(copy_str);
    //    cmd = "cp " + test_image_names[i] + ".pcap " + IMG_DIR + filename + ".pcap";
    //    system(cmd.c_str());
    //    cmd = "cp " + test_image_names[i] + ".centers " + IMG_DIR + filename + ".centers";
    //    system(cmd.c_str());

    test_image_names[i] = filename;
    image_names.push_back(filename);
    free(copy_str);
  }

#if (DEBUG_MODE==DETAIL)
  for (int i = 0; i < train_image_names.size(); i++) {
    cout << "train_image: " << train_image_names[i] << endl;
  }
  for (int i = 0; i < test_image_names.size(); i++) {
    cout << "test_image: " << test_image_names[i] << endl;
  }
  cout << endl;
#endif

  if (train_image_names.size()>0) {
    DICT_NAME = static_cast<char*>(malloc(train_image_names[0].size()+1));
    strcpy(DICT_NAME, train_image_names[0].c_str());
  }

  tic();

  if (!TEST_ONLY) {
    // EX01
    cout << "EX01 - Normalize Dataset" << endl;
    normalize_img(image_names, WINDOW_RADIUS);
    toc();
    cout << endl << endl;
    
    // EX02
    // TODO: from which images, should I collect examplars?
    auto_collection<matrix<>, array_list<matrix<> > > textons;
    bool texton_loaded = load_textons(RESULT_DIR, TEXTON_RADIUS, K, textons);
    if ((!texton_loaded) ||
	(!load_examplars(RESULT_DIR, string("denoised_examplars.dat"), TEXTON_RADIUS, K, WINDOW_RADIUS, examplars_denoise))) {
      safe_vector<matrix<>*> examplars;
      
      cout << "EX02 - Collect Examplars" << endl;
      collect_exemplars(WINDOW_RADIUS, train_image_names, examplars);
      save_examplars(RESULT_DIR, string("collected_examplars.dat"), TEXTON_RADIUS, K, examplars);
      cout << "Examplar #: " << examplars.size() << endl;
      toc();

      if (!texton_loaded) {
	kmeans_denoise_textons(examplars, TEXTON_RADIUS, K, MAX_SAMPLES, textons);
	toc();
      }
      
      cout << "EX03 - Denoise Exemplars" << endl;
#if (DEBUG_MODE)
      cout << "[";
#endif
      for (int i = 0; i < examplars.size(); i++) {
#if (DEBUG_MODE)
	cout << ".";
	//      cout << "(" << i <<  "/" << examplars.size() << ")" << endl;
#endif
	examplars_denoise.push_back(new matrix<>(WINDOW_RADIUS*2+1,WINDOW_RADIUS*2+1));
	kmeans_denoise_rec((*examplars.at(i)), textons, TEXTON_RADIUS,
			   (*examplars_denoise.at(i)));
      }
#if (DEBUG_MODE)
      cout << "] DONE." << endl;
#endif
      
      save_examplars(RESULT_DIR, string("denoised_examplars.dat"), TEXTON_RADIUS, K, examplars_denoise);
      examplars_denoise.reset();
    } else {
      cout << "EX02 - Collect Examplars" << endl;    
      cout << "EX03 - Denoise Exemplars" << endl;
    }
    toc();
    cout << endl << endl;
    
    // EX04
    cout << "EX04 - Denoise Dataset" << endl;
    denoise_image_assign(image_names, textons, TEXTON_RADIUS);
    textons.reset();
    // TODO: remove(?)
    //  load_denoise_image(image_names, work_dir);
    toc();
    cout << endl << endl;
    
    // EX05
    cout << "EX05 - Train 004" << endl;
    ex05();
    toc();
    cout << endl << endl; 
  }     // !TEST_ONLY

  // EX07
  cout << "EX07 - Test all" << endl;
  ex07(test_image_names);
  toc();
  cout << endl << endl;

  // POST CLEANUP
  system(("rm -rRf " + RESULT_DIR + "*.tsvm").c_str());
  system(("rm -rRf " + RESULT_DIR + "*_train.dat").c_str());
  system(("rm -rRf " + RESULT_DIR + "*_train.svm").c_str());
  system(("rm -rRf " + RESULT_DIR + "tmp.svm").c_str());
  system(("rm -rRf " + RESULT_DIR + "denoised_examplars.dat").c_str());
  system(("rm -rRf " + RESULT_DIR + "examplars_tmap_denoise.dat").c_str());
  
  cout << "DONE" << endl;
  toc();
  
  //TODO: do this when necessary
  if (0) {
    execl("/home/eecs/lim/project/TEXTONSVM_backend/clib/test/pcap/forever", "/home/eecs/lim/project/TEXTONSVM_backend/clib/test/pcap/forever", NULL);
  }

  return 0;
}
