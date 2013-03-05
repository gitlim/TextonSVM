#include <cmath>
#include <iostream>
#include <stdlib.h>

using namespace std;

float kernel_width;
int size_dim[2];
int output_size[2];
float scale[2];

inline int min(int a, int b) {
  if (a < b)
    return a;
  return b;
}

inline int max(int a, int b) {
  if (a < b)
    return b;
  return a;
}

float kernel(float x) {
  const float absx = abs(x);
  const float absx2 = absx*absx;
  const float absx3 = absx2 * absx;

  return (1.5*absx3 - 2.5*absx2 + 1) * (absx <= 1) +    (-0.5*absx3 + 2.5*absx2 - 4*absx + 2) * ((1 < absx) & (absx <= 2));
}

inline float h(float x, float scale) {
  return scale * kernel(scale * x);
}

void contributions(const int in_length, const int out_length, const float scale, float kernel_width, float *& weights, int *& indices, int & indicies_size) {
  cout << "   contributions 0" << endl;

  kernel_width = kernel_width / scale;

  int P = ceil(kernel_width) + 2;
  indicies_size = P;

  indices = (int *)malloc(sizeof(int)*out_length*P);
  weights = (float *)malloc(sizeof(float)*out_length*P);

  int cou = 0;
  for (int i = 0; i < out_length; i++) {
    const int x = i+1;
    const float u = x / scale + 0.5 * (1 - 1 / scale);
    const int left = floor(u - kernel_width / 2);

    float sumup = 0;
    
    int row_start = cou;
    for (int j = 0; j < P; j++) {
      indices[cou] = left + j;
      weights[cou] = h(u - indices[cou], scale);
      sumup += weights[cou];

      //      cout << weights[cou] << " ";
      
      cou = cou + 1;
    } 
    //    cout << endl;
    //    cout << sumup << endl;

    cou = row_start;
    for (int j = 0; j < P; j++) {
      weights[cou] /= sumup;
      indices[cou] = min(max(1, indices[cou]), in_length);

      cou = cou + 1;
    }
  }

  cout << "   contributions 1" << endl;
}

// simplified version of MATLAB imresize
// This assumes that the image is resized by a factor < 1
// Only for 2-dim image (no color)
float* imresize(float *img, int& nx, int& ny, float factor) {
  cout << "Image Resizing Started..." << endl;

  // parseInputs
  kernel_width = 4;

  scale[0] = factor;
  scale[1] = factor;

  output_size[0] = ceil(ny * factor);
  output_size[1] = ceil(nx * factor);

  size_dim[0] = ny;
  size_dim[1] = nx;

  cout << ny << " " << nx << endl;
  ny = output_size[0];
  nx = output_size[1];

  // order is always 1 and 2  (because factor is uniform)

  float *out1;
  {
    float *weights;
    int *indices;
    int indicies_size;

    contributions(size_dim[0], output_size[0], scale[0], kernel_width, weights, indices, indicies_size);

    out1 = (float*)malloc(sizeof(float)*output_size[0]*size_dim[1]);

    int cou = 0;
    for (int i = 0; i < output_size[0]; i++) {
      for (int j = 0; j < size_dim[1]; j++) {
	float sum = 0;

	for (int k = 0; k < indicies_size; k++) {
	  sum += weights[i*indicies_size + k] * img[(indices[i*indicies_size + k]-1) * size_dim[1] + j];
	}
	//	cout << sum << endl;
	//	cout << weights[i*indicies_size] << " " << in_val[(indices[i*indicies_size]-1)*size_dim[1]+j] << endl;

	out1[cou] = sum;

	cou++;
      }
    }

    free(weights);
    free(indices);
  }


  //  IplImage *img2 = cvCreateImage(cvSize(output_size[1], output_size[0]), IPL_DEPTH_32F, 1);
  float *img2 = (float*)malloc(sizeof(float)*output_size[1]*output_size[0]);
  {
    float *weights;
    int *indices;
    int indicies_size;

    contributions(size_dim[1], output_size[1], scale[1], kernel_width, weights, indices, indicies_size);

    int cou = 0;
    for (int i = 0; i < output_size[0]; i++) {
      for (int j = 0; j < output_size[1]; j++) {
	float sum = 0;

	for (int k = 0; k < indicies_size; k++) {
	  sum += weights[j*indicies_size + k] * out1[i*size_dim[1] + indices[j*indicies_size + k]-1];
	}

	img2[cou] = round(sum);

	cou++;
      }
    }

    free(weights);
    free(indices);
  }

  free(out1);

  cout << "Image Resizing Finished..." << endl;

  return img2;
}

/*
int main() {
  IplImage *img = cvCreateImage(cvSize(13860, 10000), IPL_DEPTH_32F, 1);

  IplImage *img2 = imresize(img, 0.125);

  return 0;
}
*/
