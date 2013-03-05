/* -*- Mode: C; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; coding: utf-8 -*- */

#include <glib.h>
#include <gtk/gtk.h>
#include <assert.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

// Library under test.
//#include <gtkanimview.h>
//#include <gtkimagescrollwin.h"
//#include <gtkimagetooldragger.h"

#include "painter.h"
#include "selector.h"
#include "TEXTONSVM.h"

#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdlib>

// Defines for backwards compatibility with GTK+ 2.6
#if !GTK_CHECK_VERSION(2,8,0)
#define GTK_STOCK_FULLSCREEN ""
#endif


// TODO:
// (1) Load/Save function (data + dontcare) : for saving, shoudl I copy original images to the target places too?
// (2) Should ask desired resize (default = 1/4)
// (3) Resolution
// (4) Should convert_resize support *.spi (non .centers file)?

using namespace std;

extern GdkPixbuf* cur_buf;

//////////////////////////////////////////////////////////////////////
///// Global data ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static GtkFileChooserDialog *open_dialog = NULL, *load_dialog = NULL, *save_dialog = NULL, *particle_dialog = NULL, *train_dialog = NULL, *test_dialog = NULL, *model_dialog = NULL;
static GtkWidget *image_param_dialog = NULL, *backend_param_dialog = NULL, *start_backend_dialog = NULL;
static GtkWidget *res_textbox = NULL, *radius_textbox = NULL;
static GtkWidget *input_resize_factor = NULL, *input_particle_size = NULL, *input_window_radius = NULL, *input_particle_window_radius = NULL, *input_texton_radius = NULL, *input_K = NULL;
static GtkWidget *input_texton_nbins = NULL, *input_grayhist_nbins = NULL, *input_nbcls = NULL; //, *input_weight_feature = NULL;
static GtkWidget *input_search_optimal = NULL, *delete_particle_button = NULL, *good_particle_button = NULL;
static GtkWindow *main_window = NULL; //, *progress_window = NULL;
GtkWidget* sub_scroll_win = NULL;

static GtkActionGroup *default_group = NULL;
static GtkActionGroup *image_group = NULL;
static GtkActionGroup *transform_group = NULL;
static gboolean is_fullscreen = FALSE;
static GtkWidget *statusbar = NULL;
GtkWidget *view;
GtkWidget *selected_size, *selected_label, *sub_label;
GtkWidget *selected_view;

// Label that displays the active selection.
static GtkWidget *sel_info_label = NULL;

// Context ID:s for the Statusbar
int help_msg_cid = -1;
int image_info_cid = -1;
int currentTool = -1;

string CUR_APP_PATH = "";

static void error_message_window(char* err_msg) {
    GtkWidget* dialog = gtk_message_dialog_new (main_window,
                                                GTK_DIALOG_DESTROY_WITH_PARENT,
                                                GTK_MESSAGE_ERROR,
                                                GTK_BUTTONS_CLOSE,
                                                err_msg);
    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
}

//////////////////////////////////////////////////////////////////////
///// Opener dialog //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void init_open_dialog ()
{
	open_dialog = (GtkFileChooserDialog *)
		gtk_file_chooser_dialog_new ("Open Image",
									 main_window,
                                     GTK_FILE_CHOOSER_ACTION_OPEN,
									 GTK_STOCK_CANCEL,
									 GTK_RESPONSE_CANCEL,
									 GTK_STOCK_OPEN,
									 GTK_RESPONSE_ACCEPT,
									 NULL);

    GtkFileFilter *imgFilter = gtk_file_filter_new();
    gtk_file_filter_set_name(imgFilter, "IMG");
    gtk_file_filter_add_pattern(imgFilter, "*.tiff");
    gtk_file_filter_add_pattern(imgFilter, "*.TIFF");
    gtk_file_filter_add_pattern(imgFilter, "*.tif");
    gtk_file_filter_add_pattern(imgFilter, "*.TIF");
    gtk_file_filter_add_pattern(imgFilter, "*.MRC");
    gtk_file_filter_add_pattern(imgFilter, "*.mrc");
    gtk_file_filter_add_pattern(imgFilter, "*.IMG");
    gtk_file_filter_add_pattern(imgFilter, "*.img");
    //    gtk_file_filter_add_pattern(imgFilter, "*.TEXTONSVM");
    //    gtk_file_filter_add_pattern(imgFilter, "*.pcap");

    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(open_dialog), imgFilter);
}

//////////////////////////////////////////////////////////////////////
///// Loader dialog //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void init_load_dialog ()
{
	load_dialog = (GtkFileChooserDialog *)
	gtk_file_chooser_dialog_new ("Load Particles",
								 main_window,
								 GTK_FILE_CHOOSER_ACTION_OPEN,
								 GTK_STOCK_CANCEL,
								 GTK_RESPONSE_CANCEL,
								 GTK_STOCK_OPEN,
								 GTK_RESPONSE_ACCEPT,
								 NULL);

    GtkFileFilter *imgFilter = gtk_file_filter_new();
    gtk_file_filter_set_name(imgFilter, "Centers");
    //gtk_file_filter_add_pattern(imgFilter, "*.centers");
    gtk_file_filter_add_pattern(imgFilter, "*.box");
    gtk_file_filter_add_pattern(imgFilter, "*.spi");

    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(load_dialog), imgFilter);

}

//////////////////////////////////////////////////////////////////////
///// Particle dialog //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void init_particle_dialog ()
{
	particle_dialog = (GtkFileChooserDialog *)
        gtk_file_chooser_dialog_new ("Choose file(s) to get good particles",
								 main_window,
								 GTK_FILE_CHOOSER_ACTION_OPEN,
								 GTK_STOCK_CANCEL,
								 GTK_RESPONSE_CANCEL,
								 GTK_STOCK_OPEN,
								 GTK_RESPONSE_ACCEPT,
								 NULL);
    GtkFileFilter *imgFilter = gtk_file_filter_new();
    gtk_file_filter_set_name(imgFilter, "IMG");
    gtk_file_filter_add_pattern(imgFilter, "*.tiff");
    gtk_file_filter_add_pattern(imgFilter, "*.TIFF");
    gtk_file_filter_add_pattern(imgFilter, "*.tif");
    gtk_file_filter_add_pattern(imgFilter, "*.TIF");
    gtk_file_filter_add_pattern(imgFilter, "*.MRC");
    gtk_file_filter_add_pattern(imgFilter, "*.mrc");

    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(particle_dialog), imgFilter);

    gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(particle_dialog), true);
}

//////////////////////////////////////////////////////////////////////
///// Train dialog //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void init_train_dialog ()
{
	train_dialog = (GtkFileChooserDialog *)
        gtk_file_chooser_dialog_new ("Choose file(s) to train",
								 main_window,
								 GTK_FILE_CHOOSER_ACTION_OPEN,
								 GTK_STOCK_CANCEL,
								 GTK_RESPONSE_CANCEL,
								 GTK_STOCK_OPEN,
								 GTK_RESPONSE_ACCEPT,
								 NULL);
    GtkFileFilter *imgFilter = gtk_file_filter_new();
    gtk_file_filter_set_name(imgFilter, "IMG");
    gtk_file_filter_add_pattern(imgFilter, "*.tiff");
    gtk_file_filter_add_pattern(imgFilter, "*.TIFF");
    gtk_file_filter_add_pattern(imgFilter, "*.tif");
    gtk_file_filter_add_pattern(imgFilter, "*.TIF");
    gtk_file_filter_add_pattern(imgFilter, "*.MRC");
    gtk_file_filter_add_pattern(imgFilter, "*.mrc");

    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(train_dialog), imgFilter);

    gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(train_dialog), false);
}

//////////////////////////////////////////////////////////////////////
///// Test dialog //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void init_test_dialog ()
{
	test_dialog = (GtkFileChooserDialog *)
        gtk_file_chooser_dialog_new ("Choose file(s) to test",
								 main_window,
								 GTK_FILE_CHOOSER_ACTION_OPEN,
								 GTK_STOCK_CANCEL,
								 GTK_RESPONSE_CANCEL,
								 GTK_STOCK_OPEN,
								 GTK_RESPONSE_ACCEPT,
								 NULL);




    GtkFileFilter *imgFilter = gtk_file_filter_new();
    gtk_file_filter_set_name(imgFilter, "IMG");
    gtk_file_filter_add_pattern(imgFilter, "*.tiff");
    gtk_file_filter_add_pattern(imgFilter, "*.TIFF");
    gtk_file_filter_add_pattern(imgFilter, "*.tif");
    gtk_file_filter_add_pattern(imgFilter, "*.TIF");
    gtk_file_filter_add_pattern(imgFilter, "*.MRC");
    gtk_file_filter_add_pattern(imgFilter, "*.mrc");

    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(test_dialog), imgFilter);

    gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(test_dialog), true);
}

//////////////////////////////////////////////////////////////////////
///// Model dialog //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void init_model_dialog ()
{
	model_dialog = (GtkFileChooserDialog *)
        gtk_file_chooser_dialog_new ("Choose directory to store results",
                                 main_window,
                                     GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER,
                                 GTK_STOCK_CANCEL,
								 GTK_RESPONSE_CANCEL,
                                     GTK_STOCK_OK,
                                     GTK_RESPONSE_OK,
								 NULL);

    gtk_file_chooser_set_create_folders(GTK_FILE_CHOOSER(model_dialog), true);
}

//////////////////////////////////////////////////////////////////////
///// Image Param dialog /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void init_image_param_dialog ()
{
    /* Create the widgets */ 
    image_param_dialog = gtk_dialog_new_with_buttons("Scale Options", NULL, 
                                               GTK_DIALOG_DESTROY_WITH_PARENT,
                                               GTK_STOCK_OK, GTK_RESPONSE_ACCEPT, GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT, NULL);
    GtkWidget *res_label, *radius_label;
    res_label = gtk_label_new ("Rescaling factor (1/n)");
    radius_label = gtk_label_new("Particle radius");

    res_textbox = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(res_textbox), "8");
    radius_textbox = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(radius_textbox), "30");
   
    // Add the label, and show everything we've added to the dialog.
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(image_param_dialog))),
                       res_label);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(image_param_dialog))),
                       res_textbox);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(image_param_dialog))),
                       radius_label);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(image_param_dialog))),
                       radius_textbox);

    gtk_widget_show_all(gtk_dialog_get_content_area(GTK_DIALOG(image_param_dialog)));
    //    gtk_widget_show_all (image_param_dialog);
    //    gtk_dialog_run(GTK_DIALOG(image_param_dialog));
}

//////////////////////////////////////////////////////////////////////
///// Backend Param dialog ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void init_backend_param_dialog ()
{
    GtkWidget *label1, *label2, *label3, *label4, *label5, *label6, *label7, *label8, *label9; 
   
    /* Create the widgets */ 
    backend_param_dialog = gtk_dialog_new_with_buttons("Advanced Options", main_window,
                                               GTK_DIALOG_DESTROY_WITH_PARENT,
                                               GTK_STOCK_OK, GTK_RESPONSE_ACCEPT, GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT, NULL);

    label1 = gtk_label_new ("Image resize factor");
    input_resize_factor = gtk_entry_new();
    gtk_entry_set_text((GtkEntry*)input_resize_factor, "8");    

    label2 = gtk_label_new ("Particle size");
    input_particle_size = gtk_entry_new();
    gtk_entry_set_text((GtkEntry*)input_particle_size, "7");

    label3 = gtk_label_new ("Window radius");
    input_window_radius = gtk_entry_new();
    gtk_entry_set_text((GtkEntry*)input_window_radius, "11");

    label4 = gtk_label_new ("Particle Window radius");
    input_particle_window_radius = gtk_entry_new();
    gtk_entry_set_text((GtkEntry*)input_particle_window_radius, "30");

    label5 = gtk_label_new ("Texton radius");
    input_texton_radius = gtk_entry_new();
    gtk_entry_set_text((GtkEntry*)input_texton_radius, "4");

    label6 = gtk_label_new ("K");
    input_K = gtk_entry_new();
    gtk_entry_set_text((GtkEntry*)input_K, "1000");

    label7 = gtk_label_new("Texton Number of Bins");
    input_texton_nbins = gtk_entry_new();
    gtk_entry_set_text((GtkEntry*)input_texton_nbins, "64");

    label8 = gtk_label_new("Grayhist Number of Bins");
    input_grayhist_nbins = gtk_entry_new();
    gtk_entry_set_text((GtkEntry*)input_grayhist_nbins, "64");

    label9 = gtk_label_new("NBCLS");
    input_nbcls = gtk_entry_new();
    gtk_entry_set_text((GtkEntry*)input_nbcls, "64");

    /*
    label10 = gtk_label_new("Weight Feature");
    input_weight_feature = gtk_entry_new();
    gtk_entry_set_text((GtkEntry*)input_weight_feature, "0.5");
    */

    input_search_optimal = gtk_check_button_new_with_label("Check to search for the optimal particle size/window radius");
    gtk_toggle_button_set_active((GtkToggleButton*)input_search_optimal, true);
   
    // Ensure that the dialog box is destroyed when the user clicks ok.
    //    g_signal_connect_swapped (param_dialog,
    //                              "response", 
    //                              G_CALLBACK (gtk_widget_destroy),
    //                              param_dialog);

    // Add the label, and show everything we've added to the dialog.
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       label1);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_resize_factor);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       label2);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_particle_size);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       label3);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_window_radius);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       label4);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_particle_window_radius);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       label5);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_texton_radius);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       label6);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_K);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       label7);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_texton_nbins);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       label8);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_grayhist_nbins);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       label9);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_nbcls);
    /*
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       label10);
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_weight_feature);
    */
    gtk_container_add (GTK_CONTAINER (gtk_dialog_get_content_area(GTK_DIALOG(backend_param_dialog))),
                       input_search_optimal);

    gtk_widget_show(label1);
    gtk_widget_show(input_resize_factor);
    gtk_widget_show(label2);
    gtk_widget_show(input_particle_size);
    gtk_widget_show(label3);
    gtk_widget_show(input_window_radius);
    gtk_widget_show(label4);
    gtk_widget_show(input_particle_window_radius);
    gtk_widget_show(label5);
    gtk_widget_show(input_texton_radius);
    gtk_widget_show(label6);
    gtk_widget_show(input_K);
    gtk_widget_show(label7);
    gtk_widget_show(input_texton_nbins);
    gtk_widget_show(label8);
    gtk_widget_show(input_grayhist_nbins);
    gtk_widget_show(label9);
    gtk_widget_show(input_nbcls);
    //    gtk_widget_show(label10);
    //gtk_widget_show(input_weight_feature);
    gtk_widget_show(input_search_optimal);

    //    gtk_widget_show_all (backend_param_dialog);
    //    gtk_dialog_run(GTK_DIALOG(backend_param_dialog));
}

static void init_start_backend_dialog() {
    GtkWidget *label = gtk_label_new ("Press OK to start pre-processing");

    start_backend_dialog = gtk_dialog_new_with_buttons("Starting Backend...", main_window,
                                               GTK_DIALOG_DESTROY_WITH_PARENT,
                                                       GTK_STOCK_OK, GTK_RESPONSE_ACCEPT, NULL);
    gtk_container_add(GTK_CONTAINER(gtk_dialog_get_content_area(GTK_DIALOG(start_backend_dialog))), label);
    gtk_widget_show(label);
}



//////////////////////////////////////////////////////////////////////
///// Saver dialog //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void init_save_dialog ()
{
	save_dialog = (GtkFileChooserDialog *)
	gtk_file_chooser_dialog_new ("Save Data",
								 main_window,
								 GTK_FILE_CHOOSER_ACTION_SAVE,
								 GTK_STOCK_CANCEL,
								 GTK_RESPONSE_CANCEL,
								 GTK_STOCK_SAVE,
								 GTK_RESPONSE_ACCEPT,
								 NULL);

    gtk_file_chooser_set_create_folders(GTK_FILE_CHOOSER(save_dialog), true);
}

///////
// mouse callback
///////
static void button_press_callback(GtkWidget *widget,
                                  GdkEventButton *event,
                                  gpointer dat)
{
    GdkRectangle image_rect;

    GdkRectangle viewport, draw_rect;
    if (!gtk_image_view_get_viewport((GtkImageView*)selected_view, &viewport))
        return ;
    if (!gtk_image_view_get_draw_rect((GtkImageView*)selected_view, &draw_rect))
        return ;

    image_rect.x = (int)((viewport.x - draw_rect.x + event->x)/SELECTED_ZOOM);
    image_rect.y = (int)((viewport.y - draw_rect.y + event->y)/SELECTED_ZOOM);

    g_print("%d %d\n", image_rect.x, image_rect.y);

    MarkSelectedPoint(image_rect.x, image_rect.y);
    MoveSelected(0, 0);
}

////////
// delete particle callback
////////
static void delete_particle_callback(GtkWidget *widget, GdkEventButton *event, gpointer dat)
{
    RemoveSelectedPoint();
    MoveSelected(0, 0);
}

///////
// good particle callback
///////
static void good_particle_callback(GtkWidget *widget, GdkEventButton *event, gpointer dat)
{
    MarkGoodParticlePoint();
}

//////////////////////////////////////////////////////////////////////
///// ImageViewerApp /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void
push_image_info (char               *basename,
                 GdkPixbuf *anim)
{
	int width = gdk_pixbuf_get_width (anim);
	int height = gdk_pixbuf_get_height (anim);
	char *msg = g_strdup_printf ("%s, %d x %d pixels",
								 basename, width, height);
	gtk_statusbar_push (GTK_STATUSBAR (statusbar), image_info_cid, msg);
	g_free (msg);
}


static void load_filename (char *path)
{
	GdkPixbuf *anim = LoadImage(path);
//	anim = NormalizeValMatrix();

//    LoadParticles(path);


//    string filename(path);
    //    filename  = filename.substr(filename.find_last_of("/")+1);
    //    cout << filename << endl;

    if (!anim)
    {
        printf ("No anim!\n");
        return;
    }

	//g_print("%d\n", gdk_pixbuf_get_pixels(anim));
	gtk_image_view_set_pixbuf((GtkImageView*)view, anim, TRUE);

    if (cur_buf) {
        g_object_unref(cur_buf);
    }
	cur_buf = gdk_pixbuf_copy(anim);

	g_object_unref(anim);


	
    char *basename = g_path_get_basename (path);
    gtk_window_set_title (main_window, basename);
    push_image_info (basename, anim);
    g_free (basename);

    gtk_action_group_set_sensitive (image_group, TRUE);

	ResetSelectedPoint();
    ShowSelected();
	
    /* Only active the transform_group if the loaded object is a single
       image -- transformations cannot be applied to animations. */
    gboolean is_image = TRUE;
    gtk_action_group_set_sensitive (transform_group, is_image);	
}

//////////////////////////////////////////////////////////////////////
///// Callbacks //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void change_image_tool_cb (GtkAction      *action,
                      GtkRadioAction *current) {
    // TODO: somewhere has to redraw selected red boxes + selected particles on the right panel

	if (currentTool == 30) {
		GdkPixbuf *anim = NormalizeValMatrix();		
		if (!anim)
		{
			printf ("No anim!\n");
			return;
		}
		
        if (cur_buf) {
            g_object_unref(cur_buf);
        }
		cur_buf = gdk_pixbuf_copy(anim);
		RedrawSelected(anim);
        ShowSelected();
		gtk_image_view_set_pixbuf((GtkImageView*)view, anim, FALSE);
        g_object_unref(anim);
	}
	
    int value = gtk_radio_action_get_current_value (current);
    GtkIImageTool *tool = dragger;
    if (value == 10)
        tool = dragger;
	else if (value == 20)
		tool = selector;
    else if (value == 30)
        tool = painter;
	currentTool = value;
    gtk_image_view_set_tool (GTK_IMAGE_VIEW (view), tool);
}

static void zoom_in_cb () {
	gtk_image_view_zoom_in (GTK_IMAGE_VIEW (view));
}

static void zoom_out_cb () {
	gtk_image_view_zoom_out (GTK_IMAGE_VIEW (view));
}

static void zoom_100_cb () {
	gtk_image_view_set_zoom (GTK_IMAGE_VIEW (view), 1.0);
}

static void zoom_to_fit_cb () {
	gtk_image_view_set_fitting (GTK_IMAGE_VIEW (view), TRUE);
}

static void open_image_cb (GtkAction *action) {
	if (!open_dialog)
		init_open_dialog ();
    if (!image_param_dialog)
        init_image_param_dialog();

	if (gtk_dialog_run (GTK_DIALOG (open_dialog)) == GTK_RESPONSE_ACCEPT)
	{
        if (gtk_dialog_run(GTK_DIALOG(image_param_dialog)) == GTK_RESPONSE_ACCEPT) {
            char *fname;
            fname = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (open_dialog));
            load_filename (fname);
            g_free (fname);
        }
        gtk_widget_destroy(GTK_WIDGET(image_param_dialog));
        image_param_dialog = NULL;
    }
	gtk_widget_destroy (GTK_WIDGET (open_dialog));
    open_dialog = NULL;
}

static void load_data_cb(GtkAction *action) {
    // if no image is opened yet
    if (!cur_buf) {
        error_message_window("No image is opened yet");
        return ;
    }

	if (!load_dialog)
		init_load_dialog ();

	if (gtk_dialog_run (GTK_DIALOG (load_dialog)) == GTK_RESPONSE_ACCEPT)
	{
        char *fname;
        fname = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (load_dialog));
            
        g_print("%s\n", fname);
        LoadParticles(fname); 
        ShowSelected();
        GdkPixbuf* clone_buf = gdk_pixbuf_copy(cur_buf);
        RedrawSelected(clone_buf);
        gtk_image_view_set_pixbuf((GtkImageView*)view, clone_buf, FALSE);
        g_object_unref(clone_buf);
        g_free (fname);
    }
	gtk_widget_destroy (GTK_WIDGET (load_dialog));	
    load_dialog = NULL;
}

static void launch_backend_cb(GtkAction *action) {
    if (!particle_dialog) {
        init_particle_dialog();
    }
    if (!train_dialog) {
        init_train_dialog();
    }
    if (!test_dialog) {
        init_test_dialog();
    }
    if (!model_dialog) {
        init_model_dialog();
    }
    if (!backend_param_dialog) {
        init_backend_param_dialog();
    }
    if (!start_backend_dialog) {
        init_start_backend_dialog();
    }

    string option = "";
    ostringstream oss;
    int n_index;
    if (gtk_dialog_run(GTK_DIALOG(particle_dialog)) == GTK_RESPONSE_ACCEPT) {
        if (gtk_dialog_run(GTK_DIALOG(train_dialog)) == GTK_RESPONSE_ACCEPT) {
            if (gtk_dialog_run(GTK_DIALOG(test_dialog)) == GTK_RESPONSE_ACCEPT) {
                if (gtk_dialog_run(GTK_DIALOG(model_dialog)) == GTK_RESPONSE_OK) {
                    if (gtk_dialog_run(GTK_DIALOG(backend_param_dialog)) == GTK_RESPONSE_ACCEPT) {
                        gtk_widget_hide(GTK_WIDGET(backend_param_dialog));
                        gtk_widget_hide(GTK_WIDGET(model_dialog));
                        gtk_widget_hide(GTK_WIDGET(test_dialog));
                        gtk_widget_hide(GTK_WIDGET(train_dialog));
                        gtk_widget_hide(GTK_WIDGET(particle_dialog));

                        gtk_dialog_run(GTK_DIALOG(start_backend_dialog));
                        
                        gtk_widget_destroy(GTK_WIDGET(start_backend_dialog));
                        start_backend_dialog = NULL;
                        
                        gtk_widget_hide(GTK_WIDGET(main_window));
                        
                        /*
                          progress_window = (GtkWindow*)gtk_window_new(GTK_WINDOW_TOPLEVEL);
                          gtk_window_set_default_size(progress_window, 300, 300);
                          gtk_widget_show_all (GTK_WIDGET (progress_window));
                          gtk_widget_show (GTK_WIDGET (progress_window));
                          gtk_window_present(progress_window);
                        */
                        
                        while (gtk_events_pending()) {
                            gtk_main_iteration();
                        }
                        
                        int ui_resize_factor = atoi(gtk_entry_get_text((GtkEntry*)input_resize_factor));
                        int ui_particle_size = atoi(gtk_entry_get_text((GtkEntry*)input_particle_size));
                        int ui_window_radius = atoi(gtk_entry_get_text((GtkEntry*)input_window_radius));
                        int ui_particle_window_radius = atoi(gtk_entry_get_text((GtkEntry*)input_particle_window_radius));
                        int ui_texton_radius = atoi(gtk_entry_get_text((GtkEntry*)input_texton_radius));
                        int ui_K = atoi(gtk_entry_get_text((GtkEntry*)input_K));
                        int ui_texton_nbins = atoi(gtk_entry_get_text((GtkEntry*)input_texton_nbins));
                        int ui_grayhist_nbins = atoi(gtk_entry_get_text((GtkEntry*)input_grayhist_nbins));
                        int ui_nbcls = atoi(gtk_entry_get_text((GtkEntry*)input_nbcls));
                        //double ui_weight_feature = atof(gtk_entry_get_text((GtkEntry*)input_weight_feature));
                        bool ui_search_optimal = gtk_toggle_button_get_active((GtkToggleButton*)input_search_optimal);                    
                        
                        /*
                          cout << ui_resize_factor << " "
                          << ui_particle_size << " "
                          << ui_window_radius << " "
                          << ui_particle_window_radius << " "
                          << ui_search_optimal << " "
                          << ui_texton_radius << " "
                          << ui_K << endl;
                        */
                        if ((ui_resize_factor < 1) || (ui_particle_size < 1) || (ui_window_radius < 1) || (ui_particle_window_radius < 1) || (ui_search_optimal < 0) || (ui_texton_radius < 1) || (ui_K < 1) ||
                            (ui_texton_nbins < 1) || (ui_grayhist_nbins < 1) || (ui_nbcls < 1)) {
                            error_message_window("Parameters are set wrong. Either they are non-positive or containing non-numerical character");
                        } else {
                            char* param[10000];
                            int cou = 1;
                            
                            // copy
                            // resize to pcap
                            // run backend
                            GFile* model_file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(model_dialog));
                            
                            //string work_dir = g_file_get_path(g_file_get_parent(model_file));
                            string work_dir = g_file_get_path(model_file) + (string)("/");
                            
                            
                            // training param
                            param[cou++] = "-i";
                            n_index = cou++;
                            string train_file;
                            GSList* train_files = gtk_file_chooser_get_files(GTK_FILE_CHOOSER(train_dialog));
                            while (train_files) {
                                // copy from get_path to work_dir+get_basename
                                string src_file = g_file_get_path((GFile*)train_files->data);
                                string target_file = work_dir + "/" + g_file_get_basename((GFile*)train_files->data);
                                target_file = target_file.substr(0,target_file.find_last_of(".")) + ".tsvm";
                                train_file = src_file;
                                
                                string src_file2 = src_file.substr(0,src_file.find_last_of(".")) + ".mask";
                                string target_file2 = target_file.substr(0,target_file.find_last_of(".")) + ".mask";
                                
                                string cmd = (string)("cp ") + (src_file.substr(0,src_file.find_last_of(".")) + ".*box") + " " + work_dir;
                                cout << cmd << endl;
                                system(cmd.c_str());
                                
                                // convert+resize target_file
                                convert_resize_image(src_file, target_file, src_file2, target_file2, ui_resize_factor);
                                
                                string tmp = target_file.substr(0,target_file.find_last_of("."));
                                param[cou] = (char*)malloc(tmp.size()+1);
                                strcpy(param[cou], tmp.c_str());
                                cou++; 
                                
                                cout << param[cou-1] << endl;
                                
                                train_files = train_files->next;
                            }
                            GSList* particle_files = gtk_file_chooser_get_files(GTK_FILE_CHOOSER(particle_dialog));
                            while (particle_files) {
                                // copy from get_path to work_dir+get_basename
                                string src_file = g_file_get_path((GFile*)particle_files->data);
                                string target_file = work_dir + "/" + g_file_get_basename((GFile*)particle_files->data);
                                target_file = target_file.substr(0,target_file.find_last_of(".")) + ".tsvm";
                                if (train_file.compare(src_file) == 0) {
                                    printf("duplicated train / particle file\n");
                                    printf("skippign this file.......\n");
                                    particle_files = particle_files->next;
                                    continue;
                                }
                                
                                string src_file2 = src_file.substr(0,src_file.find_last_of(".")) + ".mask";
                                string target_file2 = target_file.substr(0,target_file.find_last_of(".")) + ".mask";
                                
                                string cmd = (string)("cp ") + (src_file.substr(0,src_file.find_last_of(".")) + ".*box") + " " + work_dir;
                                cout << cmd << endl;
                                system(cmd.c_str());
                                
                                // convert+resize target_file
                                convert_resize_image(src_file, target_file, src_file2, target_file2, ui_resize_factor);
                                
                                string tmp = target_file.substr(0,target_file.find_last_of("."));
                                param[cou] = (char*)malloc(tmp.size()+1);
                                strcpy(param[cou], tmp.c_str());
                                cou++; 
                                
                                cout << param[cou-1] << endl;
                                
                                particle_files = particle_files->next;
                            }
                            oss << cou - n_index - 1;
                            param[n_index] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[n_index], oss.str().c_str());
                            oss.str("");
                            g_slist_free(train_files);
                            g_slist_free(particle_files);
                            
                            // testing param
                            param[cou++] = "-t";
                            n_index = cou++;
                            GSList* test_files = gtk_file_chooser_get_files(GTK_FILE_CHOOSER(test_dialog));
                            while (test_files) {
                                // copy from get_path to work_dir+get_basename
                                string src_file = g_file_get_path((GFile*)test_files->data);
                                string target_file = work_dir + "/" + g_file_get_basename((GFile*)test_files->data);
                                target_file = target_file.substr(0,target_file.find_last_of(".")) + ".tsvm";
                                
                                string src_file2 = src_file.substr(0,src_file.find_last_of(".")) + ".mask";
                                string target_file2 = target_file.substr(0,target_file.find_last_of(".")) + ".mask";
                                
                                string cmd = (string)("cp ") + (src_file.substr(0,src_file.find_last_of(".")) + ".*box") + " " + work_dir;
                                cout << cmd << endl;
                                system(cmd.c_str());
                                
                                // convert+resize target_file
                                convert_resize_image(src_file, target_file, src_file2, target_file2, ui_resize_factor);
                                
                                string tmp = target_file.substr(0,target_file.find_last_of("."));
                                param[cou] = (char*)malloc(tmp.size()+1);
                                strcpy(param[cou], tmp.c_str());
                                cou++; 
                                
                                cout << param[cou-1] << endl;
                                
                                test_files = test_files->next;
                            }
                            oss << cou - n_index - 1;
                            param[n_index] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[n_index], oss.str().c_str());
                            oss.str("");
                            g_slist_free(test_files);
                            
                            // model param                
                            param[cou++] = "-s";
                            //                oss << "-s ";
                            param[cou++] = g_file_get_basename(model_file);
                            //                oss << file << " ";
                            
                            param[cou++] = "-w";
                            //                oss << "-w ";
                            param[cou++] = g_file_get_path(g_file_get_parent(model_file));
                            //                oss << file << " ";
                            
                            param[cou++] = "-p";
                            oss << ui_resize_factor;
                            param[cou] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[cou], oss.str().c_str());
                            cou++;
                            oss.str("");
                            oss << ui_particle_size;
                            param[cou] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[cou], oss.str().c_str());
                            cou++;
                            oss.str("");
                            oss << ui_window_radius;
                            param[cou] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[cou], oss.str().c_str());
                            cou++;
                            oss.str("");
                            if (ui_search_optimal) {
                                param[cou++] = "1";
                            }
                            else {
                                param[cou++] = "0";
                            }
                            oss.str("");
                            oss << ui_particle_window_radius;
                            param[cou] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[cou], oss.str().c_str());
                            cou++;
                            oss.str("");
                            oss << ui_texton_radius;
                            param[cou] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[cou], oss.str().c_str());
                            cou++;
                            oss.str("");
                            oss << ui_K;
                            param[cou] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[cou], oss.str().c_str());
                            cou++;
                            
                            oss.str("");
                            oss << ui_texton_nbins;
                            param[cou] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[cou], oss.str().c_str());
                            cou++;
                            oss.str("");
                            oss << ui_grayhist_nbins;
                            param[cou] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[cou], oss.str().c_str());
                            cou++;
                            oss.str("");
                            oss << ui_nbcls;
                            param[cou] = (char*)malloc(oss.str().size()+1);
                            strcpy(param[cou], oss.str().c_str());
                            cou++;
                            //                            oss.str("");
                            //                            oss << ui_weight_feature;
                            //                            param[cou] = (char*)malloc(oss.str().size()+1);
                            //                            strcpy(param[cou], oss.str().c_str());
                            //                            cou++;
                            
                            
                            //                    param[cou++] = "-f";
                            //                    oss << resize_factor;
                            //                    param[cou] = (char*)malloc(oss.str().size()+1);
                            //                    strcpy(param[cou], oss.str().c_str());
                            //                    oss.str("");
                            //                    cou++;
                            
                            //                option = oss.str();
                            //cout << option << endl;
                            //                g_print("%s\n", option);
                            
                            // TODO: relative path
                            //                string backend_path = CUR_APP_PATH.substr(0,CUR_APP_PATH.find_last_of("/")) + "/../TEXTONSVM_backend/clib/test/pcap/main";
                            string backend_path = CUR_APP_PATH.substr(0,CUR_APP_PATH.find_last_of("/")) + "/main_final";
                            char* app_name = (char*)malloc(backend_path.size()+1);
                            strcpy(app_name, backend_path.c_str());
                            
                            param[0] = app_name;
                            param[cou] = NULL;
                            
                            for (int i = 0; i < cou; i++)
                                g_print("%s\n", param[i]);
                            
                            execv(app_name, param);
                            free(app_name);
                            
                            exit(1);
                            
                            g_print("quit\n");
                        }
                    }
                    gtk_widget_destroy(GTK_WIDGET(backend_param_dialog));
                    backend_param_dialog = NULL;
                }
                gtk_widget_destroy(GTK_WIDGET(model_dialog));
                model_dialog = NULL;
            }
            gtk_widget_destroy(GTK_WIDGET(test_dialog));
            test_dialog = NULL;
        }
        gtk_widget_destroy(GTK_WIDGET(train_dialog));
        train_dialog = NULL;
    }
    gtk_widget_destroy(GTK_WIDGET(particle_dialog));
    particle_dialog = NULL;
}

static void save_data_cb(GtkAction *action) {
	if (!save_dialog)
		init_save_dialog ();
	if (gtk_dialog_run (GTK_DIALOG (save_dialog)) == GTK_RESPONSE_ACCEPT)
	{
		char *fname;
		fname = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (save_dialog));
		g_print("%s\n", fname);

        SaveAll(fname);

		g_free (fname);
	}
	gtk_widget_destroy (GTK_WIDGET (save_dialog));	
    save_dialog = NULL;
}

static void flipimg_cb ()
{
    GdkPixbuf* anim = FlipImage(cur_buf);
    if (!anim)
		{
			printf ("No anim!\n");
			return;
		}
		
    if (cur_buf) {
        g_object_unref(cur_buf);
    }
    cur_buf = gdk_pixbuf_copy(anim);
    RedrawSelected(anim);
    ShowSelected();
    gtk_image_view_set_pixbuf((GtkImageView*)view, anim, FALSE);
    g_object_unref(anim);
}   



static void menu_item_select_cb (GtkMenuItem *proxy)
{
	GtkAction *action = (GtkAction*) g_object_get_data (G_OBJECT (proxy), "gtk-action");

	char *msg;
	g_object_get (G_OBJECT (action), "tooltip", &msg, NULL);

	if (msg)
	{
		gtk_statusbar_push (GTK_STATUSBAR (statusbar), help_msg_cid, msg);
		g_free (msg);
	}
}

static void menu_item_deselect_cb (GtkMenuItem *item)
{
	gtk_statusbar_pop (GTK_STATUSBAR (statusbar), help_msg_cid);
}

static void
connect_proxy_cb (GtkUIManager *ui,
				  GtkAction    *action,
				  GtkWidget    *proxy)
{
	if (!GTK_IS_MENU_ITEM (proxy))
        return;
    g_signal_connect (proxy, "select", G_CALLBACK (menu_item_select_cb), NULL);
    g_signal_connect (proxy, "deselect",
                      G_CALLBACK (menu_item_deselect_cb), NULL);
}

static void disconnect_proxy_cb (GtkUIManager *ui,
					 GtkAction    *action,
					 GtkWidget    *proxy)
{
	if (!GTK_IS_MENU_ITEM (proxy))
        return;
    g_signal_handlers_disconnect_by_func (proxy,
                                          (void*)G_CALLBACK (menu_item_select_cb),
                                          NULL);
    g_signal_handlers_disconnect_by_func (proxy,
                                          (void*)G_CALLBACK (menu_item_deselect_cb),
                                          NULL);
}

static void zoom_changed_cb (GtkImageView *view,  GtkLabel     *label)
{
	gdouble zoom = gtk_image_view_get_zoom (view);
	char *text = g_strdup_printf ("%d%%", (int)(zoom * 100.0));
	gtk_label_set_text (label, text);
	g_free (text);
}

static void kill_app_cb (void) {    
	/* Kill the widgets. */
	gtk_widget_destroy (GTK_WIDGET (main_window));
	if (open_dialog)
		gtk_widget_destroy (GTK_WIDGET (open_dialog));
    if (save_dialog)
        gtk_widget_destroy (GTK_WIDGET (save_dialog));
    if (particle_dialog)
        gtk_widget_destroy (GTK_WIDGET (particle_dialog));
    if (train_dialog)
        gtk_widget_destroy (GTK_WIDGET (train_dialog));
    if (test_dialog)
        gtk_widget_destroy (GTK_WIDGET (test_dialog));
    if (model_dialog)
        gtk_widget_destroy (GTK_WIDGET (model_dialog));
    if (load_dialog)
        gtk_widget_destroy (GTK_WIDGET (load_dialog));

	gtk_main_quit ();
}

//////////////////////////////////////////////////////////////////////
///// MainWindow /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static GtkWindow *main_window_new (GtkWidget *widget,
                 int        width,
                 int        height)
{
    GtkWindow *window = (GtkWindow *)gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_default_size (window, width, height);
    gtk_container_add (GTK_CONTAINER (window), widget);
    g_signal_connect (G_OBJECT (window), "delete_event",
                      G_CALLBACK (kill_app_cb), NULL);
    return window;
}

//////////////////////////////////////////////////////////////////////
///// UI Setup ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static const GtkActionEntry default_actions[] = {
    {"FileMenu", NULL, "_File"},
    {
        "Open",
        GTK_STOCK_OPEN,
        "_Open image",
        NULL,
        "Open an image",
        G_CALLBACK (open_image_cb)
    },
	{
		"Load",
		GTK_STOCK_OPEN,
		"_Load particles",
		NULL,
		"Load particle data",
		G_CALLBACK (load_data_cb)
	},
	{
		"Save",
		GTK_STOCK_SAVE,
		"_Save data",
		NULL,
		"Save data",
		G_CALLBACK (save_data_cb)
	},
    {
        "LaunchBackend",
        GTK_STOCK_QUIT,
        "_Launch Backend",
        NULL,
        "Launch Backend",
        G_CALLBACK (launch_backend_cb)
    },
    {
        "Quit",
        GTK_STOCK_QUIT,
        "_Quit",
        NULL,
        "Quit the program",
        G_CALLBACK (kill_app_cb)
    },
    {"EditMenu", NULL, "_Edit"},
    {"ViewMenu", NULL, "_View"},
    //    {"TranspMenu", NULL, "_Transparency"}
};

static const GtkActionEntry image_actions[] = {
    {
        "ZoomIn",
        GTK_STOCK_ZOOM_IN,
        "Zoom _In",
		"<control>plus",
        "Zoom in one step",
		G_CALLBACK (zoom_in_cb)
    },
    {
        "ZoomOut",
        GTK_STOCK_ZOOM_OUT,
        "Zoom _Out",
		"<control>minus",
        "Zoom out one step",
		G_CALLBACK (zoom_out_cb)
    },
	{
		"ZoomNormal",
		GTK_STOCK_ZOOM_100,
		"_Normal Size",
		"<control>0",
		"Set zoom to natural size of the image",
		G_CALLBACK (zoom_100_cb)
	},
	{
		"ZoomFit",
		GTK_STOCK_ZOOM_FIT,
		"Best _Fit",
		NULL,
		"Adapt zoom to fit image",
		G_CALLBACK (zoom_to_fit_cb)
	},
	{
		"FlipImg",
		GTK_STOCK_REFRESH,
		"_Flip Image",
		NULL,
		"Flip Image",
		G_CALLBACK (flipimg_cb)
	}
};

static const GtkRadioActionEntry image_tools[] = {
    {
        "DraggerTool",
        GTK_STOCK_REFRESH,
        "_Move",
        NULL,
        "Use the hand tool",
        10
    },
    {
        "PainterTool",
        GTK_STOCK_MEDIA_PLAY,
        "_Paint Unused Regions",
        NULL,
        "Use the painter tool",
        30
    },
    {
        "SelectorTool",
        GTK_STOCK_MEDIA_PAUSE,
        "_Mark Particles",
        NULL,
        "Use the rectangular selection tool",
        20
    }
};

gchar *ui_info =
    "<ui>"
    "  <menubar name = 'MenuBar'>"
    "    <menu action = 'FileMenu'>"
    "      <menuitem action = 'Open'/>"
	"      <menuitem action = 'Load'/>"
    "      <separator/>" 
	"	   <menuitem action = 'Save'/>"
    "      <separator/>" 
    "      <menuitem action = 'LaunchBackend'/>"
    "      <menuitem action = 'Quit'/>"
    "    </menu>"
    "    <menu action = 'EditMenu'>"
    //    "      <menuitem action = 'Transform'/>"
    //    "      <separator/>" 
    "      <menuitem action = 'DraggerTool'/>"
    "      <menuitem action = 'PainterTool'/>"
    "      <menuitem action = 'SelectorTool'/>"
    "    </menu>"
    "    <menu action = 'ViewMenu'>"
	"      <menuitem action = 'Fullscreen'/>"
	"      <separator/>"
    "      <menuitem action = 'ZoomIn'/>"
    "      <menuitem action = 'ZoomOut'/>"
	"      <menuitem action = 'ZoomNormal'/>"
	"      <menuitem action = 'ZoomFit'/>"
    "      <separator/>"
    "    </menu>"
    "  </menubar>"
    "  <toolbar name = 'ToolBar'>"
    //    "    <toolitem action='Quit'/>"
    "    <toolitem action='Open'/>"
    "    <toolitem action='Load'/>" 
    "    <separator/>"
    "    <toolitem action='DraggerTool'/>"
    "    <toolitem action='PainterTool'/>"
    "    <toolitem action='SelectorTool'/>"
    "    <separator/>"
    "    <toolitem action='ZoomIn'/>"
    "    <toolitem action='ZoomOut'/>"
	"    <toolitem action='ZoomNormal'/>"
	"    <toolitem action='ZoomFit'/>"
    "    <separator/>"
    "    <toolitem action='FlipImg'/>"
    "  </toolbar>"
    "</ui>";


static void parse_ui (GtkUIManager *uimanager) {
    GError *err;
    if (!gtk_ui_manager_add_ui_from_string (uimanager, ui_info, -1, &err))
    {
        g_warning ("Unable to create menus: %s", err->message);
        g_free (err);
    }
}

static void add_action_groups (GtkUIManager *uimanager) {
    // Setup the default group.
    default_group = gtk_action_group_new ("default");
    gtk_action_group_add_actions (default_group,
                                  default_actions,
                                  G_N_ELEMENTS (default_actions),
                                  NULL);
    gtk_action_group_add_radio_actions (default_group,
                                        image_tools,
                                        G_N_ELEMENTS (image_tools),
                                        10,
                                        G_CALLBACK (change_image_tool_cb),
                                        NULL);
    gtk_ui_manager_insert_action_group (uimanager, default_group, 0);

    // Setup the image group.
    image_group = gtk_action_group_new ("image");
    gtk_action_group_add_actions (image_group,
                                  image_actions,
                                  G_N_ELEMENTS (image_actions),
                                  NULL);
	gtk_action_group_set_sensitive (image_group, FALSE);
    gtk_ui_manager_insert_action_group (uimanager, image_group, 0);
}

static GtkWidget * setup_layout (GtkUIManager *uimanager) {
    GtkWidget *box = gtk_vbox_new (FALSE, 0);
    
    GtkWidget *menu = gtk_ui_manager_get_widget (uimanager, "/MenuBar");
    gtk_box_pack_start (GTK_BOX (box), menu, FALSE, FALSE, 0);

    GtkWidget *toolbar = gtk_ui_manager_get_widget (uimanager, "/ToolBar");
    gtk_box_pack_start (GTK_BOX (box), toolbar, FALSE, FALSE, 0);

	GtkWidget *scroll_win = gtk_image_scroll_win_new (GTK_IMAGE_VIEW (view));

	gtk_box_pack_start (GTK_BOX (box), scroll_win, TRUE, TRUE, 0); 

	statusbar = gtk_statusbar_new ();

    // A label in the statusbar that displays the current selection if
    // there is one.
    GtkWidget *sel_info_frame = gtk_frame_new (NULL);
    gtk_frame_set_shadow_type (GTK_FRAME (sel_info_frame), GTK_SHADOW_IN);

    sel_info_label = gtk_label_new ("");
    gtk_container_add (GTK_CONTAINER (sel_info_frame), sel_info_label);

    gtk_box_pack_start (GTK_BOX (statusbar), sel_info_frame, FALSE, FALSE, 0);

    // A label in the statusbar that displays the current zoom. It
    // updates its text when the zoom-changed signal is fired from the
    // view.
	GtkWidget *zoom_info_frame = gtk_frame_new (NULL);
	gtk_frame_set_shadow_type (GTK_FRAME (zoom_info_frame), GTK_SHADOW_IN);

	GtkWidget *zoom_info_label = gtk_label_new ("100%");
	gtk_container_add (GTK_CONTAINER (zoom_info_frame), zoom_info_label);

	g_signal_connect (G_OBJECT (view), "zoom_changed",
					  G_CALLBACK (zoom_changed_cb), zoom_info_label);

	gtk_box_pack_start (GTK_BOX (statusbar), zoom_info_frame, FALSE, FALSE, 0);
	
	gtk_box_pack_end (GTK_BOX (box), statusbar, FALSE, FALSE, 0);
    return box;
}

static void setup_main_window() {
    GtkUIManager *uimanager = gtk_ui_manager_new ();
	
	g_signal_connect (uimanager, "connect_proxy",
					  G_CALLBACK (connect_proxy_cb), NULL);
	g_signal_connect (uimanager, "disconnect_proxy",
					  G_CALLBACK (disconnect_proxy_cb), NULL);
	
	
	
    add_action_groups (uimanager);
    parse_ui (uimanager);
	
    GtkAccelGroup *accels = gtk_ui_manager_get_accel_group (uimanager);
    assert (accels);
	
	GtkWidget *hbox = gtk_hbox_new(FALSE, 0) ;

    GtkWidget *vbox1 = setup_layout (uimanager);
    //	gtk_window_set_default_size((GtkWindow*)vbox1, 600, -1);
	gtk_box_pack_start((GtkBox*)hbox, vbox1, TRUE, TRUE, 0);
		
	GtkWidget *vbox2 = gtk_vbox_new(FALSE, 0);
	sub_label = gtk_label_new(NULL);
	gtk_box_pack_start((GtkBox*)vbox2, sub_label, FALSE, FALSE, 0);

    selected_label = gtk_label_new("\nChoose Particle Size:");
    gtk_box_pack_start((GtkBox*)vbox2, selected_label, FALSE, FALSE, 0);

	selected_size = gtk_hscale_new_with_range(30, 100, 2);
	gtk_box_pack_start((GtkBox*)vbox2, selected_size, FALSE, FALSE, 0);

	selected_view = (GtkWidget*)GTK_IMAGE_VIEW(gtk_image_view_new());
	gtk_image_view_set_fitting((GtkImageView*)selected_view, FALSE);
	gtk_image_view_set_zoom (GTK_IMAGE_VIEW (selected_view), SELECTED_ZOOM);

    sub_scroll_win = gtk_image_scroll_win_new (GTK_IMAGE_VIEW (selected_view));
    //sub_scroll_win = gtk_scrolled_window_new(NULL, NULL);
    //    gtk_scrolled_window_add_with_viewport((GtkScrolledWindow*)sub_scroll_win, selected_view);
    //    gtk_window_set_default_size((GtkWindow*)sub_scroll_win, SELECTED_ZOOM*SELECTED_PER_ROW*MAX_W, 600);
    //g_signal_connect (G_OBJECT (sub_scroll_win), "button_press_event",
    g_signal_connect(G_OBJECT(selected_view), "button_press_event",
					  G_CALLBACK (button_press_callback), NULL);
	gtk_box_pack_start((GtkBox*)vbox2, sub_scroll_win, TRUE, TRUE, 0);
	
	gtk_widget_set_usize(vbox2, SELECTED_ZOOM*SELECTED_PER_ROW*MAX_W+35, -1);

	gtk_box_pack_start((GtkBox*)hbox, vbox2, FALSE, FALSE, 0);

    good_particle_button = gtk_button_new_with_label("Regular particle <-> Good particle");
    g_signal_connect(G_OBJECT(good_particle_button), "clicked",
                     G_CALLBACK(good_particle_callback), NULL);
    gtk_box_pack_start((GtkBox*)vbox2, good_particle_button, FALSE, FALSE, 0);

    delete_particle_button = gtk_button_new_with_label("Delete selected particle");
    g_signal_connect(G_OBJECT(delete_particle_button), "clicked",
                     G_CALLBACK(delete_particle_callback), NULL);
    gtk_box_pack_start((GtkBox*)vbox2, delete_particle_button, FALSE, FALSE, 0);
	
	
	
	
	
    main_window = main_window_new (hbox, 900, 600);
    gtk_window_add_accel_group (main_window, accels);
	
	gtk_widget_grab_focus (GTK_WIDGET (view));
	
	// Setup context ID:s
	help_msg_cid = gtk_statusbar_get_context_id (GTK_STATUSBAR (statusbar),
												 "help_msg");
	image_info_cid = gtk_statusbar_get_context_id (GTK_STATUSBAR (statusbar),
												   "image_info");
	g_object_unref (uimanager);
}

int main (int argc, char *argv[]) {
    //    convert_resize_image("/work/lim/RNA/images/RNA_001.tif", "/home/eecs/lim/test.pcap", 8);
    //    exit(1);

	char **filenames = NULL;
	GOptionEntry options[] = {
		{
			G_OPTION_REMAINING, '\0', 0, G_OPTION_ARG_FILENAME_ARRAY,
			&filenames, NULL, "[FILE...]"
		},
		{NULL}
	};
	GOptionContext *ctx = g_option_context_new ("Sample image viewer");
	g_option_context_add_main_entries (ctx, options, "example1");
	g_option_context_parse (ctx, &argc, &argv, NULL);
	g_option_context_free (ctx);
	
	gtk_init (&argc, &argv);

	view = (GtkWidget*)GTK_IMAGE_VIEW (gtk_image_view_new ());

    dragger = gtk_image_tool_dragger_new (GTK_IMAGE_VIEW (view));
	painter = gtk_image_tool_painter_new(GTK_IMAGE_VIEW(view));
    selector = gtk_image_tool_selector_new (GTK_IMAGE_VIEW (view));

    //    gtk_signal_connect( GTK_OBJECT(button), "button_press_event",
    //                        GTK_SIGNAL_FUNC(button_press_callback), 
    //                        NULL);

    CUR_APP_PATH = argv[0];

	setup_main_window ();
	
	if (filenames)
		load_filename (filenames[0]);

	gtk_widget_show_all (GTK_WIDGET (main_window));
	
    gtk_main ();

    return 0;
}

