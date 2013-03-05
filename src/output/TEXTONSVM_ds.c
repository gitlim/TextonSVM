#include <glib.h>
#include <gtk/gtk.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
//#include <vector>
#include <iostream>
#include <cmath>
#include <dirent.h>

// Library under test.
#include <gtkimageview/gtkanimview.h>
#include <gtkimageview/gtkimagescrollwin.h>
#include <gtkimageview/gtkimagetooldragger.h>

//#include "selector.h"
#include "TEXTONSVM.h"

//#define DEBUG_MODE

const int thickness = 4;

//struct MyPoint *gt_center = NULL;
//int* det_center = NULL;
vector<struct MyPoint> gt_center, det_center;
vector<int> vgt;

int WINDOW_RADIUS = -1;

extern float *conf_matrix;

vector<GtkImageView*> centers;
extern GtkWidget *view;
extern GtkWidget *det_center_box;
extern GtkWidget *det_scroll_win;
extern GtkWidget *sub_label;
extern int currentTool;
GdkPixbuf *cur_buf = NULL;

extern int RES;
extern int Threshold;

static int selected_det = 0;
//static int selected_gt = 0;
static int selected_gt = -99;

int TotalDetected() {
  return det_center.size();
}

int TotalShown() {
  int cou = 0;
  for (unsigned int i = 0; i < det_center.size(); i++) {
    if (det_center[i].conf >= Threshold) {
      cou++;
    }
  }
  return cou;
}

void MarkDetectedPoint(int new_selected) {
  int old_ind = selected_det;
  selected_det = new_selected;

  HighlightSelected(old_ind, selected_det, det_center, centers);
}

/*
void MarkGTPoint(int x, int y) {
  int old_ind = selected_gt;
  selected_gt = (x / MAX_W) + (y / MAX_H) * SELECTED_PER_ROW;

  HighlightSelected(old_ind, selected_gt, gt_center, gt_center_view);
}
*/


///////
// mouse callback
///////
static void det_button_press_callback(GtkWidget *widget,
                                  GdkEventButton *event,
                                  gpointer dat)
{
  int new_selected=(long)dat;
  MarkDetectedPoint(new_selected);

    GdkPixbuf * clone_buf = gdk_pixbuf_copy(cur_buf);
    if (currentTool == 10) {
        PlotGTCenters(clone_buf, false);
        PlotDetCenters(clone_buf);
        gtk_image_view_set_pixbuf((GtkImageView*)view, clone_buf, FALSE);
        g_object_unref(clone_buf);
    }
    else if (currentTool == 30) {
        PlotGTCenters(clone_buf, true);
        PlotDetCenters(clone_buf);
        gtk_image_view_set_pixbuf((GtkImageView*)view, clone_buf, FALSE);
        g_object_unref(clone_buf);
    }
}


void HighlightSelected(int old_selected, int new_selected, vector<struct MyPoint>& center, vector<GtkImageView*> centers) {
#ifdef DEBUG_MODE
  g_print("Highlight %d %d\n", old_selected, new_selected);
#endif
  if (center.size() <= 0)
    return ; 

  GdkPixbuf* tmp = gtk_image_view_get_pixbuf((GtkImageView*)centers[old_selected]);
  guchar *pixels = gdk_pixbuf_get_pixels(tmp);
  int stride = gdk_pixbuf_get_rowstride(tmp);
  int n_chans = gdk_pixbuf_get_n_channels(tmp);

  int x, y;

  short R, G, B;
  R = 0xDC;
  G = 0xDC;
  B = 0xDC;
  for (x = 0; x < MAX_W; x++) {
    for (int margin = 0; margin < MARGIN; margin++) {
      pixels[(margin)*stride+(x)*n_chans+0] = R;
      pixels[(margin)*stride+(x)*n_chans+1] = G;
      pixels[(margin)*stride+(x)*n_chans+2] = B;
      pixels[(MAX_H-1-margin)*stride+(x)*n_chans+0] = R;
      pixels[(MAX_H-1-margin)*stride+(x)*n_chans+1] = G;
      pixels[(MAX_H-1-margin)*stride+(x)*n_chans+2] = B;
    }
  }
  for (y = 0; y < MAX_H; y++) {
    for (int margin = 0; margin < MARGIN; margin++) {
      pixels[(y)*stride+(margin)*n_chans+0] = R;
      pixels[(y)*stride+(margin)*n_chans+1] = G;
      pixels[(y)*stride+(margin)*n_chans+2] = B;
      pixels[(y)*stride+(MAX_W-1-margin)*n_chans+0] = R;
      pixels[(y)*stride+(MAX_W-1-margin)*n_chans+1] = G;
      pixels[(y)*stride+(MAX_W-1-margin)*n_chans+2] = B;
    }
  }
  gtk_image_view_set_pixbuf(centers[old_selected], tmp, FALSE);  


  tmp = gtk_image_view_get_pixbuf((GtkImageView*)centers[new_selected]);
  pixels = gdk_pixbuf_get_pixels(tmp);
  stride = gdk_pixbuf_get_rowstride(tmp);
  n_chans = gdk_pixbuf_get_n_channels(tmp);

  R = 0xFF;
  G = 0xFF;
  B = 0x00;
  for (x = 0; x < MAX_W; x++) {
    for (int margin = 1; margin < MARGIN; margin++) {
      pixels[(margin)*stride+(x)*n_chans+0] = R;
      pixels[(margin)*stride+(x)*n_chans+1] = G;
      pixels[(margin)*stride+(x)*n_chans+2] = B;
      pixels[(MAX_H-1-margin)*stride+(x)*n_chans+0] = R;
      pixels[(MAX_H-1-margin)*stride+(x)*n_chans+1] = G;
      pixels[(MAX_H-1-margin)*stride+(x)*n_chans+2] = B;
    }
  }
  for (y = 0; y < MAX_H; y++) {
    for (int margin = 1; margin < MARGIN; margin++) {
      pixels[(y)*stride+(margin)*n_chans+0] = R;
      pixels[(y)*stride+(margin)*n_chans+1] = G;
      pixels[(y)*stride+(margin)*n_chans+2] = B;
      pixels[(y)*stride+(MAX_W-1-margin)*n_chans+0] = R;
      pixels[(y)*stride+(MAX_W-1-margin)*n_chans+1] = G;
      pixels[(y)*stride+(MAX_W-1-margin)*n_chans+2] = B;
    }
  }

  gtk_image_view_set_pixbuf(centers[new_selected], tmp, FALSE);
  //  g_object_unref(tmp);
}

GtkWidget *entire_box = NULL;
void ShowSelected(vector<struct MyPoint>& center, GtkWidget* center_box) {
  g_print("ShowSelected...\n");

  int how_many_lines = (center.size() - 1) / SELECTED_PER_ROW + 1;

  int y, x, n; //, ofs;
  
  if (cur_buf == NULL)
    return ;
  guchar *im_pixels = gdk_pixbuf_get_pixels(cur_buf);
  int im_height = gdk_pixbuf_get_height(cur_buf);
  int im_width = gdk_pixbuf_get_width(cur_buf);
  int im_n_chans = gdk_pixbuf_get_n_channels(cur_buf);
  int im_stride = gdk_pixbuf_get_rowstride(cur_buf);
  
#ifdef DEBUG_MODE
  g_print("TEXTONSVM_ds (ShowSelected): %d %d %d %d %d\n", im_pixels, im_height, im_width, im_n_chans, im_stride);
#endif

  unsigned int k = 0;

  centers.clear();
  if (entire_box) {
    gtk_container_remove((GtkContainer*)center_box, entire_box);
  }
  gtk_widget_set_usize(center_box, -1, how_many_lines*(MAX_H*SELECTED_ZOOM+9));
  entire_box = gtk_vbox_new(FALSE, 0);
  gtk_widget_set_usize(entire_box, -1, how_many_lines*(MAX_H*SELECTED_ZOOM+9));
  for (int i = 0; i < how_many_lines; i++) {
    GtkWidget *tmp3 = gtk_hbox_new(FALSE, 0);
    for (int j = 0; j < SELECTED_PER_ROW; j++) {
      GtkWidget *tmp_img_conf = gtk_vbox_new(FALSE, 0);

      GtkImageView* tmp2 = GTK_IMAGE_VIEW(gtk_image_view_new());
      gtk_image_view_set_fitting((GtkImageView*)tmp2, FALSE);
      gtk_image_view_set_zoom (GTK_IMAGE_VIEW (tmp2), SELECTED_ZOOM);
      g_signal_connect(G_OBJECT(tmp2), "button_press_event",
		       G_CALLBACK (det_button_press_callback), (gpointer)k);
      

      GdkPixbuf* tmp = gdk_pixbuf_new((GdkColorspace)0, FALSE, 8, MAX_W, MAX_H);
      guchar *pixels = gdk_pixbuf_get_pixels(tmp);
      int stride = gdk_pixbuf_get_rowstride(tmp);
      int n_chans = gdk_pixbuf_get_n_channels(tmp);
      //int height = gdk_pixbuf_get_height(tmp);
      //int width = gdk_pixbuf_get_width(tmp);

      int dx, dy;
      int sx, sy;
#ifdef DEBUG_MODE
      g_print("TEXTONSVM_ds: %d %d\n", center[k].w, center[k].h);
      g_print("%d %d\n", center[k].x, center[k].y);
#endif
      int dupx = (MAX_W-MARGIN*2)/center[k].w, dupy = (MAX_H-MARGIN*2)/center[k].h;

      for (x = 0; x < center[k].w; x++) {
	for (y = 0; y < center[k].h; y++)
	  for (n = 0; n < n_chans; n++) {
	    sx = center[k].x+x-center[k].w/2;
	    sy = center[k].y+y-center[k].h/2;
	    for (dx = 0; dx < dupx; dx++)
	      for (dy = 0; dy < dupy; dy++) {
		pixels[(y*dupy + dy + MARGIN) * stride + (x*dupx + dx + MARGIN)*n_chans+n] = im_pixels[sy*im_stride + sx*im_n_chans+n];
	      }
	  }
      }
      for (x = 0; x < center[k].w*dupx; x++) {
	for (y = 0; y < MARGIN; y++)
	  for (n = 0; n < n_chans; n++) {
	    pixels[(y) * stride + (x + MARGIN)*n_chans+n] = 0x00;
	  }
	for (y = center[k].h*dupy + MARGIN; y < MAX_H; y++)
	  for (n = 0; n < n_chans; n++) {
	    pixels[(y) * stride + (x + MARGIN)*n_chans+n] = 0x00;
	  }
      }
      for (x = 0; x < MARGIN; x++) {
	for (y = 0; y < MAX_H; y++)
	  for (n = 0; n < n_chans; n++) {
	    pixels[(y) * stride + (x)*n_chans+n] = 0x00;
	  }
      }
      for (x = center[k].w*dupx + MARGIN; x < MAX_W; x++) {
	for (y = 0; y < MAX_H; y++)
	  for (n = 0; n < n_chans; n++) {
	    pixels[(y) * stride + (x)*n_chans+n] = 0x00;
	  }
      }

      short R = 0xDC;
      short G = 0xDC;
      short B = 0xDC;

      for (x = 0; x < MAX_W; x++) {
	for (int margin = 0; margin < MARGIN; margin++) {
	  pixels[(margin)*stride+(x)*n_chans+0] = R;
	  pixels[(margin)*stride+(x)*n_chans+1] = G;
	  pixels[(margin)*stride+(x)*n_chans+2] = B;
	  pixels[(MAX_H-1-margin)*stride+(x)*n_chans+0] = R;
	  pixels[(MAX_H-1-margin)*stride+(x)*n_chans+1] = G;
	  pixels[(MAX_H-1-margin)*stride+(x)*n_chans+2] = B;
	}
      }
      for (y = 0; y < MAX_H; y++) {
	for (int margin = 0; margin < MARGIN; margin++) {
	  pixels[(y)*stride+(margin)*n_chans+0] = R;
	  pixels[(y)*stride+(margin)*n_chans+1] = G;
	  pixels[(y)*stride+(margin)*n_chans+2] = B;
	  pixels[(y)*stride+(MAX_W-1-margin)*n_chans+0] = R;
	  pixels[(y)*stride+(MAX_W-1-margin)*n_chans+1] = G;
	  pixels[(y)*stride+(MAX_W-1-margin)*n_chans+2] = B;
	}
      }

      gtk_image_view_set_pixbuf(tmp2, tmp, FALSE);
      gtk_widget_set_usize((GtkWidget*)tmp2, MAX_W*SELECTED_ZOOM, MAX_H*SELECTED_ZOOM);
      gtk_box_pack_start((GtkBox*)tmp_img_conf, (GtkWidget*)tmp2, FALSE, FALSE, 0);
      gtk_widget_show((GtkWidget*)tmp2);
      centers.push_back(tmp2);

      char buffer[10];
      sprintf(buffer, "%d", int(center[k].conf/2.55));
      GtkWidget *label = gtk_label_new(buffer);
      gtk_widget_set_usize(label, -1, 9);
      gtk_box_pack_start((GtkBox*)tmp_img_conf, label, FALSE, FALSE, 0);
      gtk_widget_show(label);
      
      gtk_box_pack_start((GtkBox*)tmp3, tmp_img_conf, FALSE, FALSE, 0);
      gtk_widget_show(tmp_img_conf);

      k = k + 1;
      if (k >= center.size()) {
	break;
      }
    }
    gtk_box_pack_start((GtkBox*)entire_box, (GtkWidget*)tmp3, FALSE, FALSE, 0);
    gtk_widget_show((GtkWidget*)tmp3);
    if (k >= center.size()) {
      break;
    }
  }
  gtk_box_pack_start((GtkBox*)center_box, entire_box, FALSE, FALSE, 0);
  gtk_widget_show(entire_box);
}

/*
void ShowGTCenter() {
  ShowSelected(gt_center, gt_center_view);
  HighlightSelected(selected_gt, selected_gt, gt_center, gt_center_view);
}
*/

void ShowDetectedCenter() {
  ShowSelected(det_center, det_center_box);
  HighlightSelected(selected_det, selected_det, det_center, centers);
}

void LoadCenters(string filename) {
#ifdef DEBUG_MODE
  g_print("LoadCenters\n");
#endif

  FILE *fin;

  fin = fopen((filename.substr(0, filename.find_last_of(".")) + ".det_cen").c_str(), "r");
  det_center.clear();

  if (fin) {
    struct MyPoint point;
    int tmp;
    fscanf(fin, "%d %d\n", &tmp, &WINDOW_RADIUS);
    WINDOW_RADIUS *= 2;
    while (fscanf(fin, "%d %d %lf\n", &point.y, &point.x, &point.conf) == 3) {
      //      fscanf(fin, "%d", &point.x);
      //      fscanf(fin, "%d", &point.y);
      //      fscanf(fin, "%f", &tmp); // confidence  TODO: should I use this?
      //      fscanf(fin, "%lf", &point.conf);
      point.conf *= 2.55;
      point.w = WINDOW_RADIUS;
      point.h = WINDOW_RADIUS;

      det_center.push_back(point);
    }
#ifdef DEBUG_MODE
    printf("Detected Centers: %d\n", det_center.size());
#endif

    for (unsigned int i = 0; i < det_center.size(); i++) {
      for (unsigned int j = i+1; j < det_center.size(); j++) {
	if (det_center[i].conf < det_center[j].conf) {
	  MyPoint tmp = det_center[i];
	  det_center[i] = det_center[j];
	  det_center[j] = tmp;
	  //	  cout << det_center[i].conf << "  " << det_center[j].conf << endl;
	}
      }
    }

    fclose(fin);
  }

#ifdef DEBUG_MODE
  for (unsigned int i = 0; i < det_center.size(); i++) {
    cout << det_center[i].x << " " << det_center[i].y << " " << det_center[i].conf << endl;
  }
#endif
  
  fin = fopen((filename.substr(0, filename.find_last_of(".")) + ".box").c_str(), "r");
  gt_center.clear();

  if (fin) {    
    struct MyPoint point;
    int x, y, w, h, garbage;
    while (fscanf(fin, "%d %d %d %d %d\n", &y, &x, &h, &w, &garbage) == 5) {
      x = x + w / 2;
      y = y + h / 2;
      
      point.w = WINDOW_RADIUS;
      point.h = WINDOW_RADIUS;
      point.x = x / RES;
      point.y = y / RES;
      if ((point.x - point.w < 0) || (point.y - point.h < 0)) {
	continue;
      }
      
      gt_center.push_back(point);
    }
#ifdef DEBUG_MODE
    printf("Ground Truth: %d\n", gt_center.size());
#endif
    
    fclose(fin);
  }


#ifdef DEBUG_MODE
  g_print("LoadCenters is done.\n");
#endif
}

void calculate_PR_pcap(const vector<struct MyPoint> &gt_center, const vector<struct MyPoint> &det_center, const int& tx, const int& ty, const int& radius, vector<double> &P, vector<double> &R, vector<int> &vgt) {
  cout << "calculate_PR_pcap started" << endl;

  // Line 10 - 19
  //  matrix<int> lab2indcent(tx, ty);
  // matrix<bool> gt_balls(tx, ty);

  vector< vector<int> > lab2indcent, gt_balls;
  lab2indcent.resize(tx);
  gt_balls.resize(tx);
  for (int i = 0; i < tx; i++) {
    lab2indcent[i].resize(ty);
    gt_balls[i].resize(ty);
  }

#ifdef DEBUG_MODE
  cout << tx << " " << ty << endl;
  cout << det_center.size() << " " << gt_center.size() << endl;
  cout << radius << endl;
#endif

  for (int x = 0; x < tx; x++)
    for (int y = 0; y < ty; y++) {
      lab2indcent[x][y] = -1;
      gt_balls[x][y] = false;
    }

  for (int mx = -radius; mx < radius; mx++) {
    for (int my = -radius; my < radius; my++) {
      if (round(sqrt(mx*mx + my*my)) > radius)
        continue;
      for (unsigned int i = 0; i < det_center.size(); i++) {
        int x = det_center[i].x;
        int y = det_center[i].y;
        if ((x + mx < 0) || (x + mx >= tx) || (y + my < 0) || (y + my >= ty))
          continue;
        if ((lab2indcent[x+mx][y+my] > -1) && (det_center[lab2indcent[x+mx][y+my]].conf > det_center[i].conf))
          continue;
        else
          lab2indcent[x+mx][y+my] = i;
      }
    }
  }

  // Line 22-43  
  vgt.resize(gt_center.size());
  //  vector<int> vgt(gt_center.size());
  for (unsigned int p = 0; p < gt_center.size(); p++) {
    int x = gt_center[p].x;
    int y = gt_center[p].y;
    if (lab2indcent[x][y] > -1)
      vgt[p] = det_center[lab2indcent[x][y]].conf;
    else
      vgt[p] = 0;

    for (int mx = -radius; mx < radius; mx++) {
      for (int my = -radius; my < radius; my++) {
        if (round(sqrt(mx*mx + my*my)) > radius)
          continue;
        if ((x + mx < 0) || (x + mx >= tx) || (y + my < 0) || (y + my >= ty))
          continue;
        gt_balls[x+mx][y+my] = true;
      }
    }
  }

  vector<int> vfp(det_center.size());
  for (unsigned int l = 0; l < vfp.size(); l++) {
    if (!gt_balls[det_center[l].x][det_center[l].y])
      vfp[l] = det_center[l].conf;
    else
      vfp[l] = 0;
  }

  // Line 48-
  vector<int> TP, FP, FN, temp, temp2;
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


  double P2[256], R2[256];
  vector<double> F;
  P.resize(256);
  R.resize(256);
  F.resize(256);
  for (int l = 1; l <= 255; l++) {
    P[l] = TP[l] * 1.0 / (TP[l] + FP[l]);
    R[l] = TP[l] * 1.0 / (TP[l] + FN[l]); 
    F[l] = 2 * P[l] * R[l] / (P[l] + R[l]);
    //printf("%.3f %.3f\n", P[l], R[l]);

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

  cout << "calculate_PR_pcap ended" << endl;
}

void calculate_image_PR(const int &nx, const int &ny, vector<double> &P, vector<double> &R) {
  const int radius = WINDOW_RADIUS / 4;
#ifdef DEBUG_MODE
  cout << "radius = " << radius << endl;
#endif
  calculate_PR_pcap(gt_center, det_center, nx, ny, radius, P, R, vgt);
}

void calculate_overall_PR(const string &path, const int &nx, const int &ny, vector<double> &overallP, vector<double> &overallR) {
  const int radius = WINDOW_RADIUS / 4;

  string dir = path.substr(0, path.find_last_of("/")+1);

  DIR *dirp;
  struct dirent *entry;

  FILE *fin;

  vector<struct MyPoint> overall_gt_center, overall_det_center;

  // cout << dir << endl;
  if ((dirp = opendir(dir.c_str()))) {
    while ((entry = readdir(dirp))) {
      if (entry->d_type == 8) {
	//printf("%s\n", entry->d_name);
	string filename = entry->d_name;

	if (filename.substr(filename.find_last_of(".")+1) == "box") {
	  fin = fopen((dir + filename).c_str(), "r");
	  
	  if (fin) {    
	    struct MyPoint point;
	    int y, x, h, w, garbage;
	    while (fscanf(fin, "%d %d %d %d %d\n", &y, &x, &h, &w, &garbage) == 5) {
	      y = y + h/2;
	      x = x + w/2;

	      //point.w = w;
	      //point.h = h;
	      point.w = WINDOW_RADIUS;
	      point.h = WINDOW_RADIUS;
	      point.x = x / RES;
	      point.y = y / RES;      
	      
	      if ((point.x - point.w < 0) || (point.y - point.h < 0)) {
		continue;
	      }
      
	      overall_gt_center.push_back(point);
	    }
    
	    fclose(fin);
	  }
	} else if (filename.substr(filename.find_last_of(".")+1) == "det_cen") {
	  fin = fopen((dir+filename).c_str(), "r");

	  if (fin) {
	    struct MyPoint point;
	    int tmp;
	    fscanf(fin, "%d %d\n", &tmp, &tmp);
	    while (fscanf(fin, "%d %d %lf\n", &point.y, &point.x, &point.conf) == 3) {
	      point.conf *= 2.55;
	      point.w = WINDOW_RADIUS;
	      point.h = WINDOW_RADIUS;

	      overall_det_center.push_back(point);
	    }
	    fclose(fin);
	  }
	}
      }
    }
    closedir(dirp);
  }

  /*
  fin = fopen(path.substr(0, path.find_last_of("/")).c_str(), "r");
  fclose(fin);
  */

  vector<int> overall_vgt;
  calculate_PR_pcap(overall_gt_center, overall_det_center, nx, ny, radius, overallP, overallR, overall_vgt);

  FILE *fout = fopen((dir+"overall.perf").c_str(), "wt");
  for (unsigned int i = 0; i < overallP.size(); i++) {
    fprintf(fout, "%.3f %.3f\n", overallP[i], overallR[i]);
  }
  fclose(fout);
}

void update_box_file(string path, int threshold) {
  FILE *fout = fopen((path.substr(0, path.find_last_of(".")) + ".det").c_str(), "w");

  for (unsigned int i = 0; i < det_center.size(); i++) {
    if (det_center[i].conf >= threshold) {
      int a,b,c,d;
      a = det_center[i].y-det_center[i].h/2;
      b = det_center[i].x-det_center[i].w/2;
      c = det_center[i].h;
      d = det_center[i].w;
      fprintf(fout, "%d %d %d %d -3\n", a*RES, b*RES, c*RES, d*RES);
    }
  }
  for (unsigned int i = 0; i < vgt.size(); i++) {
    if (vgt[i] < threshold) {
      int a,b,c,d;
      //a = gt_center[i].y-gt_center[i].h/2;
      //b = gt_center[i].x-gt_center[i].w/2;
      //c = gt_center[i].h;
      //d = gt_center[i].w;
      //fprintf(fout, "%d %d %d %d -3\n", b*RES, a*RES, d*RES, c*RES);
      a = gt_center[i].y*RES-gt_center[i].h*RES/2;
      b = gt_center[i].x*RES-gt_center[i].w*RES/2;
      c = gt_center[i].h*RES;
      d = gt_center[i].w*RES;
      fprintf(fout, "%d %d %d %d -3\n", a, b, c, d);
    }
  }

  fclose(fout);

  fout = fopen((path.substr(0, path.find_last_of(".")) + ".th").c_str(), "w");
  fprintf(fout, "%d\n", (int)(threshold/255.0*100));
  fclose(fout);
}

void PlotGTCenters(GdkPixbuf* img, bool darkgreen) {
#ifdef DEBUG_MODE
  g_print("PlotGTCenters...\n");
#endif
  
  guchar *pixels = gdk_pixbuf_get_pixels(img);
  int stride = gdk_pixbuf_get_rowstride(img);
  int n_chans = gdk_pixbuf_get_n_channels(img);
  //int height = gdk_pixbuf_get_height(img);
  //int width = gdk_pixbuf_get_width(img);
  
  int x, y, ofs;
  //  int wx = 60, wy = 50;
  // TODO: what size of wx, wy?
  int wx = 30, wy = 20;
  
  for (unsigned int i = 0; i < gt_center.size(); i++) {
    GdkRectangle rect = {gt_center[i].x -wx/2, gt_center[i].y -wy/2, wx, wy};

    short R, G, B;
    if ((int)i == selected_gt) {
      R = 0xFF;
      G = 0x00;
      B = 0x00;
    } else {
      if (darkgreen) {
	R = 0x00;
	G = 0x64;
	B = 0x00;
      } else {
	R = 0x7f;
	G = 0xff;
	B = 0x00;
      }
    }

    for (y = rect.y; y < rect.y + rect.height; y++) {
      for (x = rect.x; x < rect.x + thickness; x++) {
	ofs = y * stride + x * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }
      
      for (x = rect.x + rect.width - thickness; x < rect.x + rect.width; x++) {
	ofs = y * stride + x * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }
    }
    
    for (x = rect.x; x < rect.x + rect.width; x++) {
      for (y = rect.y; y < rect.y + thickness; y++) {
	ofs = y * stride + x * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }
      
      for (y = rect.y + rect.height - thickness; y < rect.y+rect.height; y++) {
	ofs = y * stride + x * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }		
    }
  }
}

void PlotDetCenters(GdkPixbuf* img) {
#ifdef DEBUG_MODE
  g_print("PlotDetCenters...\n");
#endif

  guchar *pixels = gdk_pixbuf_get_pixels(img);
  int stride = gdk_pixbuf_get_rowstride(img);
  int n_chans = gdk_pixbuf_get_n_channels(img);
  //int height = gdk_pixbuf_get_height(img);
  //int width = gdk_pixbuf_get_width(img);
  
  int x, y, ofs;
  //  int wx = 50, wy = 60;
  // TODO: what size of wx, wy?
  int wx = 20, wy = 25;

  GdkRectangle pb_rect = {0, 0, gdk_pixbuf_get_width (img), gdk_pixbuf_get_height (img)};
  
  for (unsigned int i = 0; i < det_center.size(); i++) {
    if (det_center[i].conf < Threshold)
	continue;
    //if (GetVal(conf_matrix, det_center[i].x, det_center[i].y) < Threshold)
    //continue;
    
    GdkRectangle rect = {det_center[i].x -wx/2, det_center[i].y -wy/2, wx, wy};		
    gdk_rectangle_intersect (&pb_rect, &rect, &rect);

    short R, G, B;
    if ((int)i == selected_det) {
      R = 0xFF;
      G = 0x00;
      B = 0x00;
    } else {
      R = 0xFF;
      G = 0x00;
      B = 0xFF;
    }
    
    for (y = rect.y; y < rect.y + rect.height; y++) {
      for (x = rect.x; x < rect.x + thickness; x++) {
	ofs = y * stride + x * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }
      
      for (x = rect.x + rect.width - thickness; x < rect.x + rect.width; x++) {
	ofs = y * stride + x * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }
    }
    
    for (x = rect.x; x < rect.x + rect.width; x++) {
      for (y = rect.y; y < rect.y + thickness; y++) {
	ofs = y * stride + x * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }
      
      for (y = rect.y + rect.height - thickness; y < rect.y+rect.height; y++) {
	ofs = y * stride + x * n_chans;
	pixels[ofs] = R;
	pixels[ofs+1] = G;
	pixels[ofs+2] = B;
      }		
    }
  }
}

