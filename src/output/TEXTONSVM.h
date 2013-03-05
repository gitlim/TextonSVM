#include <glib.h>
#include <gtk/gtk.h>
#include <assert.h>
#include <stdio.h>
#include <string>
#include <vector>

// Library under test.
#include <gtkimageview/gtkanimview.h>
#include <gtkimageview/gtkimagescrollwin.h>
#include <gtkimageview/gtkimagetooldragger.h>

#define MAX_W (110)
#define MAX_H (110)
#define MARGIN (5)
#define SELECTED_PER_ROW (2)
#define SELECTED_ZOOM (1.5)

using namespace std;

struct MyPoint {
  int x;
  int y;
  int w;
  int h;
  double conf;
};
void MarkDetectedPoint(int);
void MarkGTPoint(int);
void HighlightSelected(int old_selected, int new_selected, vector<struct MyPoint>& center, vector<GtkImageView*> centers);
void ShowDetectedCenter();
void ShowGTCenter();

void UpdateConfidence();

// Load Image
GdkPixbuf* LoadImage(string);
GdkPixbuf* LoadDetection(string);
GdkPixbuf* LoadConfidence(string);
int ApplyThreshold();
int TotalDetected();
int TotalShown();
void LoadCenters(string);

// Tools
static GtkIImageTool *dragger = NULL;
static GtkIImageTool *selector = NULL;
static GtkIImageTool *painter = NULL;


inline short GetVal(float *, int, int);

void PlotGTCenters(GdkPixbuf*, bool);
void PlotDetCenters(GdkPixbuf*);

//GdkPixbuf* PlotPR(string);
//GdkPixbuf* PlotPR(const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&);
//GdkPixbuf* PlotPR();
GdkPixbuf* PlotImagePR();
GdkPixbuf* PlotOverallPR();
void update_box_file(string, int);
//void calculate_PR_pcap(const vector<struct MyPoint> &, const vector<struct MyPoint> &, const int&, const int&, const int&, vector<double>&, vector<double>&);
void calculate_image_PR(const int&, const int&, vector<double>&, vector<double>&);
void calculate_overall_PR(const string&, const int&, const int&, vector<double>&, vector<double>&);
void compute_perf(string);
