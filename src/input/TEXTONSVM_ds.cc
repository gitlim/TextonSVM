#include <vector>
#include <glib.h>
#include <gtk/gtk.h>
#include <assert.h>
#include <stdio.h>
#include <iostream>

// Library under test.
//#include "gtkanimview.h"
//#include "gtkimagescrollwin.h"
//#include "gtkimagetooldragger.h"

#include "selector.h"
#include <TEXTONSVM.h>

#define DEBUG_MODE

extern GtkWidget *view;
extern GtkWidget *selected_view, *sub_label;
GdkPixbuf *cur_buf = NULL;

static vector<struct MyPoint>  point_selected;
//static struct MyPoint point_selected[10000];
//static int selected_cou = 0;
static int selected_ind = 0;

extern int RES;

void MarkGoodParticlePoint() {
  if (selected_ind >= (int)point_selected.size()) {
    return ;
  }

  point_selected[selected_ind].good_center = !point_selected[selected_ind].good_center;

  HighlightSelected(selected_ind, selected_ind);
}

void MarkSelectedPoint(int x, int y) {
  int old_ind = selected_ind;

  selected_ind = (x / MAX_W) + (y / MAX_H) * SELECTED_PER_ROW;

#ifdef DEBUG_MODE
  g_print("%d %d %d \n", MAX_W, MAX_H, selected_ind);
#endif

  HighlightSelected(old_ind, selected_ind);
}

void RemoveSelectedPoint() {
  if (selected_ind >= (int)point_selected.size()) {
    return ;
  }

  point_selected.erase(point_selected.begin()+selected_ind);
  selected_ind = point_selected.size() - 1;

  ShowSelected();

  // TODO: if user clicks to add and zoom, only one activates.. 
  // TODO: should I change this to gtk_text?
  char tmp[50];
  sprintf(tmp, "Number of selected particles: %4lu", point_selected.size());
  gtk_label_set_text((GtkLabel*)sub_label, tmp);
}

void AddSelectedPoint(int x, int y, int w, int h, bool good, bool flipped) {
  int rx = x, ry = y;
  int ny, nx;
  if (IMG_FLIPPED) {
    ny = gdk_pixbuf_get_width(cur_buf);
    nx = gdk_pixbuf_get_height(cur_buf);
  }
  else {
    nx = gdk_pixbuf_get_width(cur_buf);
    ny = gdk_pixbuf_get_height(cur_buf);
  }
  if (flipped) {
    rx = y;
    ry = x;
  }

  if ((x < 0) || (x + w >= nx) ||
      (y < 0) || (y + h >= ny)) {
#ifdef DEBUG_MODE
    g_print("AddSelectedPoint: out of range\n");
    g_print("     %d %d %d %d\n", x, y, w, h);
#endif
    return ;
  }

  for (unsigned int i = 0; i < point_selected.size(); i++) {
    struct MyPoint &point = point_selected[i];
    if ((point.x == rx) && (point.y == ry) && (point.w == w) && (point.h == h)) {
      if (!point.good_center)
	point.good_center = good;
      return ;
    }
  }

  struct MyPoint new_point;
  new_point.x = rx;
  new_point.y = ry;
  new_point.w = w;
  new_point.h = h;
  new_point.good_center = good;
  point_selected.push_back(new_point);

#ifdef DEBUG_MODE
  g_print("AddSelectedPoint: Pointed Added: %d.\n", point_selected.size());
  g_print("   %d %d %d %d\n", x, y, w, h);
#endif

  // TODO: if user clicks to add and zoom, only one activates.. 
  // TODO: should I change this to gtk_text?
  char tmp[50];
  sprintf(tmp, "Number of selected particles: %4lu", point_selected.size());
  gtk_label_set_text((GtkLabel*)sub_label, tmp);
}

struct MyPoint GetSelectedPoint(int i) {
	return point_selected[i];
}

void LoadParticles(string name) {
  string ext = name.substr(name.find_last_of(".") + 1);

  cout << "Resolution: " << RES << endl;

  // what size should I add? (i.e. should I resize?)
  if (ext == "spi") {
    // load spi
    FILE *fin = fopen(name.c_str(), "rt");
    int tmp, tmp2;
    double y, x;
    while (fscanf(fin, "%d %d %lf %lf\n", &tmp, &tmp2, &x, &y) == 4) {
#ifdef DEBUG_MODE
      cout << y << " " << x << endl;
#endif
      AddSelectedPoint((int)(x/RES), (int)(y/RES), 30, 30, false, false);
    }
    fclose(fin);
  }
  else if (ext == "box") {
    // load centers
    FILE *fin = fopen(name.c_str(), "rt");
    int y, x, w, h, garbage;
    //while (fscanf(fin, "%d %d %d %d %d\n", &y, &x, &h, &w, &garbage) == 5) {
    while (fscanf(fin, "%d %d %d %d %d\n", &x, &y, &w, &h, &garbage) == 5) {
      AddSelectedPoint((int)(x+w/2)/RES, (int)(y+h/2)/RES, (int)w/RES, (int)h/RES, false, false);
      //      AddSelectedPoint(x-10, y-10, 20, 20);
    }
    fclose(fin);
  }
  /* old format*/
  else if (ext == "centers") {
    // load centers
    FILE *fin = fopen(name.c_str(), "rt");
    int y, x;
    while (fscanf(fin, "%d %d\n", &y, &x) == 2) {
      if ((x==0) && (y == 0)) {
	break;
      }
      // TODO: particle size
      AddSelectedPoint((int)x/RES, (int)y/RES, 30, 30, false, false);
      //      AddSelectedPoint(x-10, y-10, 20, 20);
    }
    fclose(fin);
  }
  else {
    cout << "File extension didn't match." << endl;
  }

  {
    string good_centers = name.substr(0, name.find_last_of(".")) + ".gbox";
    FILE *fin = fopen(good_centers.c_str(), "rt");
    if (fin) {
      int y, x, w, h, garbage;
      //while (fscanf(fin, "%d %d %d %d %d\n", &y, &x, &h, &w, &garbage) == 5) {
      while (fscanf(fin, "%d %d %d %d %d\n", &x, &y, &w, &h, &garbage) == 5) {
	AddSelectedPoint((int)(x+w/2)/RES, (int)(y+h/2)/RES, (int)w/RES, (int)h/RES, true, false);
	//      AddSelectedPoint(x-10, y-10, 20, 20);
      }
    }
    fclose(fin);
  }
  /* old format
  {
    string good_centers = name.substr(0, name.find_last_of(".")) + ".good_centers";
    FILE *fin = fopen(good_centers.c_str(), "rt");
    if (fin) {
      int y, x;
      while (fscanf(fin, "%d %d\n", &y, &x) == 2) {
	if ((x==0) && (y == 0)) {
	  break;
	}
	// TODO: particle size
	AddSelectedPoint((int)x/RES, (int)y/RES, 30, 30, true, false);
	//      AddSelectedPoint(x-10, y-10, 20, 20);
      }
    }
  }
  */
}

void SaveParticles(string name) {
  FILE *fout = fopen((name+".box").c_str(), "wt");  
  FILE *fout2 = fopen((name+".gbox").c_str(), "wt");

  for (unsigned int i = 0; i < point_selected.size(); i++) {
    fprintf(fout, "%d %d %d %d -3\n", (point_selected[i].x-point_selected[i].w/2)*RES, (point_selected[i].y-point_selected[i].h/2)*RES, point_selected[i].w*RES, point_selected[i].h*RES);
    if (point_selected[i].good_center)
      fprintf(fout2, "%d %d %d %d -3\n", (point_selected[i].x-point_selected[i].w/2)*RES, (point_selected[i].y-point_selected[i].h/2)*RES, point_selected[i].w*RES, point_selected[i].h*RES);
  }

  fclose(fout2);
  fclose(fout);  


  /* old format
  FILE *fout = fopen((name+".centers").c_str(), "wt");  
  FILE *fout2 = fopen((name+".good_centers").c_str(), "wt");

  // TODO: should I stroe resolution?
  for (unsigned int i = 0; i < point_selected.size(); i++) {
    //    fprintf(fout, "%d %d %d %d\n", point_selected[i].x, point_selected[i].y, point_selected[i].w, point_selected[i].h); 
    //    fprintf(fout, "%d %d\n", point_selected[i].y, point_selected[i].x);

    // TODO: shouldn't I store point size?
    //    fprintf(fout, "%d %d\n", point_selected[i].y+point_selected[i].h/2, point_selected[i].x+point_selected[i].w/2);
    fprintf(fout, "%d %d\n", point_selected[i].y*RES, point_selected[i].x*RES);
    if (point_selected[i].good_center)
      fprintf(fout2, "%d %d\n", point_selected[i].y*RES, point_selected[i].x*RES);
  }

  fprintf(fout, "0 0\n");
  fprintf(fout2, "0 0\n");

  fclose(fout2);
  fclose(fout);  
  */
}

void ResetSelectedPoint() {
  point_selected.clear();
  char tmp[50];
  sprintf(tmp, "Number of selected particles: %4d", 0);
  gtk_label_set_text((GtkLabel*)sub_label, tmp);
  selected_ind = 0;
}

// TODO: should I highlight which one is currently shown?
// If clicked, should I choose that selected tile?
// (if so, it might be a bit complicated, as mouse event is already included)
void RedrawSelected(GdkPixbuf* clone_buf) {
  guchar *pixels = gdk_pixbuf_get_pixels(clone_buf);
  int n_chans = gdk_pixbuf_get_n_channels(clone_buf);
  int stride = gdk_pixbuf_get_rowstride(clone_buf);
  int x, y;
  int ofs;

  for (unsigned int k = 0; k < point_selected.size(); k++) {
    short R = 0xFF, G = 0x00, B = 0x00;
    if (point_selected[k].good_center) {
      R = 0x00;
      G = 0xFF;
    }
    if ((int)k == selected_ind) {
      R = 0xFF;
      G = 0xFF;
    }

    for (y = point_selected[k].y - point_selected[k].h/2; y < point_selected[k].y + point_selected[k].h/2; y++) {
      for (x = point_selected[k].x - point_selected[k].w/2; x < point_selected[k].x - point_selected[k].w/2 + 5; x++) {
	int rx = x, ry = y;
	if (IMG_FLIPPED) {
	  rx = y;
	  ry = x;
	}
	if ((rx < 0) || (ry < 0))
	  continue;

	ofs = ry * stride + rx * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }
      
      for (x = point_selected[k].x + point_selected[k].w/2 - 5; x <point_selected[k].x + point_selected[k].w/2; x++) {
	int rx = x, ry = y;
	if (IMG_FLIPPED) {
	  rx = y;
	  ry = x;
	}
	if ((rx < 0) || (ry < 0))
	  continue;

	ofs = ry * stride + rx * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }
    }
    
    for (x = point_selected[k].x - point_selected[k].w/2; x < point_selected[k].x + point_selected[k].w/2; x++) {
      for (y = point_selected[k].y - point_selected[k].h/2; y < point_selected[k].y - point_selected[k].h/2 + 5; y++) {
	int rx = x, ry = y;
	if (IMG_FLIPPED) {
	  rx = y;
	  ry = x;
	}
	if ((rx < 0) || (ry < 0))
	  continue;

	ofs = ry * stride + rx * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }
      
      for (y = point_selected[k].y + point_selected[k].h/2 - 5; y < point_selected[k].y + point_selected[k].h/2; y++) {
	int rx = x, ry = y;
	if (IMG_FLIPPED) {
	  rx = y;
	  ry = x;
	}
	if ((rx < 0) || (ry < 0))
	  continue;

	ofs = ry * stride + rx * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }		
    }
  }
}

void MoveSelected(int dx, int dy) {
  if (point_selected.size() == 0) {
#ifdef DEBUG_MODE
		g_print("MoveSelected: no selected pt\n");
#endif
		return ;
	}
	
  if ((point_selected[selected_ind].x + dx < 0) ||
      (point_selected[selected_ind].x + point_selected[selected_ind].w + dx >= gdk_pixbuf_get_width(cur_buf)) ||
      (point_selected[selected_ind].y + dy < 0) ||
      (point_selected[selected_ind].y + point_selected[selected_ind].h + dy >= gdk_pixbuf_get_height(cur_buf))) {
#ifdef DEBUG_MODE
    g_print("MoveSelected: out of range\n");
#endif
    return ;
  }
	
  point_selected[selected_ind].x += dx;
  point_selected[selected_ind].y += dy;
		
  GdkPixbuf* clone_buf = gdk_pixbuf_copy(cur_buf);
  RedrawSelected(clone_buf);
  gtk_image_view_set_pixbuf((GtkImageView*)view, clone_buf, FALSE);
  g_object_unref(clone_buf);
}

// Highlight in the selected view
void HighlightSelected(int old_selected, int new_selected) {
#ifdef DEBUG_MODE
  g_print("Highlight %d %d\n", old_selected, new_selected);
#endif
  //int how_many_lines = ((int)point_selected.size() - 1) / SELECTED_PER_ROW + 1;
  if (point_selected.size() <= 0)
    return ; 

  //  GdkPixbuf* tmp = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, SELECTED_PER_ROW * MAX_W, how_many_lines*MAX_H);
  //  gtk_image_view_set_pixbuf((GtkImageView*)selected_view, tmp, FALSE);
  GdkPixbuf* tmp = gtk_image_view_get_pixbuf((GtkImageView*)selected_view);
  guchar *pixels = gdk_pixbuf_get_pixels(tmp);
  int stride = gdk_pixbuf_get_rowstride(tmp);
  int n_chans = gdk_pixbuf_get_n_channels(tmp);

  int k = old_selected;

  int px = (k % SELECTED_PER_ROW) * MAX_W;
  int py = (k / SELECTED_PER_ROW) * MAX_H;

  int x, y;

  short R, G, B;
  if (point_selected[k].good_center) {
    R = 0x00;
    G = 0xFF;
    B = 0x00;
  } else {
    R = 0xDC;
    G = 0xDC;
    B = 0xDC;
  }
  for (x = 0; x < MAX_W; x++) {
    for (int margin = 0; margin < MARGIN; margin++) {
      pixels[(py+margin)*stride+(px+x)*n_chans+0] = R;
      pixels[(py+margin)*stride+(px+x)*n_chans+1] = G;
      pixels[(py+margin)*stride+(px+x)*n_chans+2] = B;
      pixels[(py+MAX_H-1-margin)*stride+(px+x)*n_chans+0] = R;
      pixels[(py+MAX_H-1-margin)*stride+(px+x)*n_chans+1] = G;
      pixels[(py+MAX_H-1-margin)*stride+(px+x)*n_chans+2] = B;
    }
  }
  for (y = 0; y < MAX_H; y++) {
    for (int margin = 0; margin < MARGIN; margin++) {
      pixels[(py+y)*stride+(px+margin)*n_chans+0] = R;
      pixels[(py+y)*stride+(px+margin)*n_chans+1] = G;
      pixels[(py+y)*stride+(px+margin)*n_chans+2] = B;
      pixels[(py+y)*stride+(px+MAX_W-1-margin)*n_chans+0] = R;
      pixels[(py+y)*stride+(px+MAX_W-1-margin)*n_chans+1] = G;
      pixels[(py+y)*stride+(px+MAX_W-1-margin)*n_chans+2] = B;
    }
  }
  

  k = new_selected;

  px = (k % SELECTED_PER_ROW) * MAX_W;
  py = (k / SELECTED_PER_ROW) * MAX_H;

  R = 0xFF;
  G = 0xFF;
  B = 0x00;
  for (x = 0; x < MAX_W; x++) {
    for (int margin = 1; margin < MARGIN; margin++) {
      pixels[(py+margin)*stride+(px+x)*n_chans+0] = R;
      pixels[(py+margin)*stride+(px+x)*n_chans+1] = G;
      pixels[(py+margin)*stride+(px+x)*n_chans+2] = B;
      pixels[(py+MAX_H-1-margin)*stride+(px+x)*n_chans+0] = R;
      pixels[(py+MAX_H-1-margin)*stride+(px+x)*n_chans+1] = G;
      pixels[(py+MAX_H-1-margin)*stride+(px+x)*n_chans+2] = B;
    }
  }
  for (y = 0; y < MAX_H; y++) {
    for (int margin = 1; margin < MARGIN; margin++) {
      pixels[(py+y)*stride+(px+margin)*n_chans+0] = R;
      pixels[(py+y)*stride+(px+margin)*n_chans+1] = G;
      pixels[(py+y)*stride+(px+margin)*n_chans+2] = B;
      pixels[(py+y)*stride+(px+MAX_W-1-margin)*n_chans+0] = R;
      pixels[(py+y)*stride+(px+MAX_W-1-margin)*n_chans+1] = G;
      pixels[(py+y)*stride+(px+MAX_W-1-margin)*n_chans+2] = B;
    }
  }

  gtk_image_view_set_pixbuf((GtkImageView*)selected_view, tmp, FALSE);
  //  g_object_unref(tmp);
}

void ShowSelected() {
  int how_many_lines = (point_selected.size() - 1) / SELECTED_PER_ROW + 1;
  if (point_selected.size() <= 0)
    how_many_lines = 1;
	
  GdkPixbuf* tmp = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, SELECTED_PER_ROW * MAX_W, how_many_lines*MAX_H);
  guchar *pixels = gdk_pixbuf_get_pixels(tmp);
  int stride = gdk_pixbuf_get_rowstride(tmp);
  int n_chans = gdk_pixbuf_get_n_channels(tmp);
  
  int y, x, n;
  
  if (cur_buf == NULL)
    return ;
  guchar *im_pixels = gdk_pixbuf_get_pixels(cur_buf);
  //int im_height = gdk_pixbuf_get_height(cur_buf);
  //int im_width = gdk_pixbuf_get_width(cur_buf);
  int im_n_chans = gdk_pixbuf_get_n_channels(cur_buf);
  int im_stride = gdk_pixbuf_get_rowstride(cur_buf);
  
#ifdef DEBUG_MODE
  g_print("TEXTONSVM_ds (ShowSelected): %d %d %d\n", im_pixels, im_n_chans, im_stride);
#endif

  int dx, dy;
  int sx, sy;
  int px, py;
  short R, G, B;
  for (unsigned int k = 0; k < point_selected.size(); k++) {
#ifdef DEBUG_MODE
    g_print("TEXTONSVM_ds: %d %d\n", point_selected[k].w, point_selected[k].h);
    g_print("   %d %d\n", point_selected[k].x, point_selected[k].y);
#endif
    px = (k % SELECTED_PER_ROW) * MAX_W;
    py = (k / SELECTED_PER_ROW) * MAX_H;

    int dupx = MAX_W/point_selected[k].w, dupy = MAX_H/point_selected[k].h;

    int ox, oy;
    if (IMG_FLIPPED) {
      ox = point_selected[k].y;
      oy = point_selected[k].x;
    }
    else {
      ox = point_selected[k].x;
      oy = point_selected[k].y;
    }
    for (x = 0; x < point_selected[k].w; x++) {
      for (y = 0; y < point_selected[k].h; y++)
	for (n = 0; n < n_chans; n++) {
	  sx = ox+x-point_selected[k].w/2;
	  sy = oy+y-point_selected[k].h/2;
	  for (dx = 0; dx < dupx; dx++) {
	    for (dy = 0; dy < dupy; dy++) {
	      pixels[(py + y*dupy+dy + MARGIN) * stride + (px + x*dupx+dx + MARGIN)*n_chans+n] = 
		im_pixels[sy*im_stride + sx*im_n_chans+n];
	    }
	  }
	}
    }
    for (x = 0; x < point_selected[k].w*dupx; x++) {
      for (y = 0; y < MARGIN; y++)
	for (n = 0; n < n_chans; n++) {
	  pixels[(py + y) * stride + (px + x + MARGIN)*n_chans+n] = 0x00;
	}
      for (y = point_selected[k].h*dupy + MARGIN; y < MAX_H; y++)
	for (n = 0; n < n_chans; n++) {
	  //		sx = point_selected[k].x+x;
	  //		sy = point_selected[k].y+y;	
	  pixels[(py + y) * stride + (px + x + MARGIN)*n_chans+n] = 0x00;
	}
    }
    for (x = 0; x < MARGIN; x++) {
      for (y = 0; y < MAX_H; y++)
	for (n = 0; n < n_chans; n++) {
	  //		sx = point_selected[k].x+x;
	  //		sy = point_selected[k].y+y;	
	  pixels[(py + y) * stride + (px + x)*n_chans+n] = 0x00;
	}
    }
    for (x = point_selected[k].w*dupx + MARGIN; x < MAX_W; x++) {
      for (y = 0; y < MAX_H; y++)
	for (n = 0; n < n_chans; n++) {
	  //		sx = point_selected[k].x+x;
	  //		sy = point_selected[k].y+y;	
	  pixels[(py + y) * stride + (px + x)*n_chans+n] = 0x00;
	}
    }

    if (point_selected[k].good_center) {
      R = 0x00;
      G = 0xFF;
      B = 0x00;
    } else {
      R = 0xDC;
      G = 0xDC;
      B = 0xDC;
    }
    for (x = 0; x < MAX_W; x++) {
      for (int margin = 0; margin < MARGIN; margin++) {
	pixels[(py+margin)*stride+(px+x)*n_chans+0] = R;
	pixels[(py+margin)*stride+(px+x)*n_chans+1] = G;
	pixels[(py+margin)*stride+(px+x)*n_chans+2] = B;
	pixels[(py+MAX_H-1-margin)*stride+(px+x)*n_chans+0] = R;
	pixels[(py+MAX_H-1-margin)*stride+(px+x)*n_chans+1] = G;
	pixels[(py+MAX_H-1-margin)*stride+(px+x)*n_chans+2] = B;
      }
    }
    for (y = 0; y < MAX_H; y++) {
      for (int margin = 0; margin < MARGIN; margin++) {
	pixels[(py+y)*stride+(px+margin)*n_chans+0] = R;
	pixels[(py+y)*stride+(px+margin)*n_chans+1] = G;
	pixels[(py+y)*stride+(px+margin)*n_chans+2] = B;
	pixels[(py+y)*stride+(px+MAX_W-1-margin)*n_chans+0] = R;
	pixels[(py+y)*stride+(px+MAX_W-1-margin)*n_chans+1] = G;
	pixels[(py+y)*stride+(px+MAX_W-1-margin)*n_chans+2] = B;
      }
    }
  }
  
  int additional = (SELECTED_PER_ROW - point_selected.size()) % SELECTED_PER_ROW;
  if ((additional < 0) || (point_selected.size() <= 0))
    additional += SELECTED_PER_ROW;
  for (unsigned int k = point_selected.size(); k < point_selected.size() + additional; k++) {
    px = (k % SELECTED_PER_ROW) * MAX_W;
    py = (k / SELECTED_PER_ROW) * MAX_H;
    for (x = 0; x < MAX_W; x++) 
      for (y = 0; y < MAX_H; y++)
	for (n = 0; n < n_chans; n++) {
	  //		sx = point_selected[k].x+x;
	  //		sy = point_selected[k].y+y;	
	  pixels[(py + y) * stride + (px + x)*n_chans+n] = 0xDC;
	//	pixels[(py + y) * stride + (px + x)*n_chans+0] = 0xDC;
	//	pixels[(py + y) * stride + (px + x)*n_chans+1] = 0xDA;
	//	pixels[(py + y) * stride + (px + x)*n_chans+2] = 0xD6;
	}
  }

  gtk_image_view_set_pixbuf((GtkImageView*)selected_view, tmp, FALSE);
  g_object_unref(tmp);

  HighlightSelected(selected_ind, selected_ind);
}
