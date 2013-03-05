#include <glib.h>
#include <gtk/gtk.h>
#include <assert.h>
#include <stdio.h>
#include <string>

// Library under test.
#include <gtkimageview/gtkanimview.h>
#include <gtkimageview/gtkimagescrollwin.h>
#include <gtkimageview/gtkimagetooldragger.h>

using namespace std;

#define MAX_W (110)
#define MAX_H (110)
#define MARGIN (5)
#define SELECTED_PER_ROW (2)
#define SELECTED_ZOOM (1.5)

struct MyPoint {
  int x;
  int y;
  int w;
  int h;
  bool good_center;
};
void MarkGoodParticlePoint();
void MarkSelectedPoint(int x, int y);
void RemoveSelectedPoint();
void AddSelectedPoint(int x, int y, int w, int h, bool good, bool flipped);
struct MyPoint GetSelectedPoint(int i);
void LoadParticles(string name);
void SaveParticles(string name);
void ResetSelectedPoint();

void MoveSelected(int dx, int dy);
void HighlightSelected(int , int);
void ShowSelected();
void RedrawSelected(GdkPixbuf* clone_buf);

// MRC Load
GdkPixbuf* LoadImage(char*);
GdkPixbuf* NormalizeValMatrix();
GdkPixbuf* FlipImage(GdkPixbuf*);
void SaveAll(string);

void convert_resize_image(string, string, string, string, int);


// Tools
static GtkIImageTool *dragger = NULL;
static GtkIImageTool *selector = NULL;
static GtkIImageTool *painter = NULL;
//static GtkWidget *sub_scroll_win = NULL;


extern bool IMG_FLIPPED;

// Dontcare
inline int NotDontCare(int x, int y);
inline void SetDontCare(int x, int y, int val);

// Constant values
//const int RES = 8;
//int RES;
//const int RES = 1;
