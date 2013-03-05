/* -*- Mode: C; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; coding: utf-8 -*- */

#include <glib.h>
#include <gtk/gtk.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <dirent.h>

//#include <painter.h>
//#include "selector.h"
#include "TEXTONSVM.h"


// Defines for backwards compatibility with GTK+ 2.6
#if !GTK_CHECK_VERSION(2,8,0)
#define GTK_STOCK_FULLSCREEN ""
#endif


#define number_template ("Detected Particles (Shown / Total) = (%d/%d)")

// TODO:
// (1) Load/Save function (data + dontcare)
// (3) Dynamic range
// (2) When threshold is applied, should I reduce particles in the list?
// (4) Overall threshold graph
// (5) Perf range between 0-100?

extern GdkPixbuf* cur_buf;

//////////////////////////////////////////////////////////////////////
///// Global data ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static GtkFileChooserDialog *open_dialog = NULL;
GtkFileFilter *imgFilter = NULL;
static GtkWindow *main_window;
static GtkActionGroup *default_group = NULL;
static GtkActionGroup *image_group = NULL;
static GtkActionGroup *transform_group = NULL;
static gboolean is_fullscreen = FALSE;
static GtkWidget *statusbar = NULL;
GtkImageView *view = NULL;
GtkWidget *sub_label = NULL, *selected_size = NULL;
GdkPixbuf *original_img = NULL, *detected_img = NULL, *detected_th_img = NULL, *current_img = NULL;
GdkPixmap *pixmap = NULL;
GtkWidget *plot_image_view = NULL, *plot_overall_view = NULL;
GtkWidget *det_scroll_win = NULL, *det_center_box = NULL;
GtkWidget *number_label = NULL;

// Label that displays the active selection.
static GtkWidget *sel_info_label = NULL;

// Context ID:s for the Statusbar
int help_msg_cid = -1;
int image_info_cid = -1;
int currentTool = -1;

int imgInd = 0;
vector<string> imgList;

char* filepath = NULL;
string CUR_APP_PATH;

string current_img_name;

extern short color_map[64*3];

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

    imgFilter = gtk_file_filter_new();
    gtk_file_filter_set_name(imgFilter, "NORM");
    gtk_file_filter_add_pattern(imgFilter, "*.NORM");
    gtk_file_filter_add_pattern(imgFilter, "*.norm");

    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(open_dialog), imgFilter);
}

/*
//////////////////////////////////////////////////////////////////////
///// Loader dialog //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void init_load_dialog ()
{
	load_dialog = (GtkFileChooserDialog *)
	gtk_file_chooser_dialog_new ("Load Data",
								 main_window,
								 GTK_FILE_CHOOSER_ACTION_OPEN,
								 GTK_STOCK_CANCEL,
								 GTK_RESPONSE_CANCEL,
								 GTK_STOCK_OPEN,
								 GTK_RESPONSE_ACCEPT,
								 NULL);
}
*/

 /*
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
}
 */

  /*
///////
// mouse callback
///////
static gint gt_button_press_callback(GtkWidget *widget,
                                  GdkEventButton *event,
                                  gpointer dat)
{
    GdkRectangle image_rect;

    GdkRectangle viewport, draw_rect;
    //if (!gtk_image_view_get_viewport((GtkImageView*)gt_center_view, &viewport))
    //        return 0;
    //if (!gtk_image_view_get_draw_rect((GtkImageView*)gt_center_view, &draw_rect))
    //        return 0;

    image_rect.x = (int)((viewport.x - draw_rect.x + event->x)/SELECTED_ZOOM);
    image_rect.y = (int)((viewport.y - draw_rect.y + event->y)/SELECTED_ZOOM);

    g_print("%d %d\n", image_rect.x, image_rect.y);

    //MarkGTPoint(image_rect.x, image_rect.y);

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
  */


void update_number_label() {
    char tmp_buffer[200];    
    sprintf(tmp_buffer, number_template, TotalShown(), TotalDetected());
    gtk_label_set_text((GtkLabel*)number_label, tmp_buffer);
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


static void load_filename (const char *path) {
    printf("%s\n", path);
    
    if (filepath) {
        free(filepath);
    }
    filepath = (char*)malloc(sizeof(char)*(strlen(path)+1));
    strcpy(filepath, path);

	original_img = LoadImage(string(path));
    detected_th_img = LoadDetection(path);
    detected_img = LoadConfidence(path);
	LoadCenters(path);

    update_number_label();
	
    compute_perf(path);
    gtk_image_set_from_pixbuf((GtkImage*)plot_image_view, PlotImagePR());
    gtk_image_set_from_pixbuf((GtkImage*)plot_overall_view, PlotOverallPR());

	currentTool = 10;
	current_img = gdk_pixbuf_copy(original_img);
    PlotGTCenters(current_img, false);
    PlotDetCenters(current_img);

    /*
    if ((!original_img) || (!detected_img) || (!detected_th_img)) {
        printf ("No anim!\n");
        return;
    }
    */
    
    if (cur_buf) {
        g_object_unref(cur_buf);
    }
    cur_buf = gdk_pixbuf_copy(original_img);
    ShowDetectedCenter();
    //ShowGTCenter();

	gtk_image_view_set_pixbuf(view, current_img, TRUE);
	
    char *basename = g_path_get_basename (path);
    gtk_window_set_title (main_window, basename);
    push_image_info (basename, original_img);
    g_free (basename);

    gtk_action_group_set_sensitive (image_group, TRUE);

    /* Only active the transform_group if the loaded object is a single
       image -- transformations cannot be applied to animations. */
    gboolean is_image = TRUE;
    gtk_action_group_set_sensitive (transform_group, is_image);	

    current_img_name = path;

    //    calculate_PR_pcap(1733, 1250, 4);

    int tmp_th = 0;
    string filename = path;
    FILE *fin = fopen((filename.substr(0,filename.find_last_of(".")) + ".th").c_str(), "r");
    if (fin) {
        fscanf(fin, "%d\n", &tmp_th);
        fclose(fin);
    }
    gtk_range_set_value((GtkRange*)selected_size, tmp_th);
    UpdateConfidence();
    //UpdateParticles();
}

//////////////////////////////////////////////////////////////////////
///// Callbacks //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static void change_image_tool_cb (GtkAction      *action,
                      GtkRadioAction *current) {	
    int value = gtk_radio_action_get_current_value (current);
	currentTool = value;
	
	if (value == 10) {
		current_img = gdk_pixbuf_copy(original_img);
		PlotGTCenters(current_img, false);
		PlotDetCenters(current_img);
        if (cur_buf) {
            g_object_unref(cur_buf);
        }
        cur_buf = gdk_pixbuf_copy(original_img);
		gtk_image_view_set_pixbuf(view, current_img, FALSE);
	}
	else if (value == 20) {
		gtk_image_view_set_pixbuf(view, detected_img, FALSE);
	}
	else if (value == 30) {
		current_img = gdk_pixbuf_copy(detected_th_img);
		PlotGTCenters(current_img, true);
		PlotDetCenters(current_img);
        if (cur_buf) {
            g_object_unref(cur_buf);
        }
        cur_buf = gdk_pixbuf_copy(detected_th_img);
		gtk_image_view_set_pixbuf(view, current_img, FALSE);
	}
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
	if (!open_dialog) {
		init_open_dialog ();
    }

	if (gtk_dialog_run (GTK_DIALOG (open_dialog)) == GTK_RESPONSE_ACCEPT)
	{
		char *fname;
		fname = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (open_dialog));
        

        string path = fname;
        string dir = path.substr(0, path.find_last_of("/")+1);

        imgList.clear();
        DIR *dirp;
        struct dirent *entry;
        
        if ((dirp = opendir(dir.c_str()))) {
            while ((entry = readdir(dirp))) {
                if (entry->d_type == 8) {
                    string filename = entry->d_name;

                    if (filename.substr(filename.find_last_of(".")+1) == "norm") {
                        imgList.push_back(dir + filename);
                    }
                    if ((imgList.size() > 0) && (imgList.back().compare(fname) == 0)) {
                        imgInd = imgList.size()-1;
                    }
                }
            }
        }
        closedir(dirp);

		load_filename (fname);
		g_free (fname);
	}
    gtk_widget_destroy(GTK_WIDGET (open_dialog));
    open_dialog = NULL;
}

/*
static void save_data_cb(GtkAction *action) {
	if (!save_dialog)
		init_save_dialog ();
	if (gtk_dialog_run (GTK_DIALOG (save_dialog)) == GTK_RESPONSE_ACCEPT)
	{
		char *fname;
		fname = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (save_dialog));
		g_print("%s\n", fname);
		g_free (fname);
	}
	gtk_widget_hide (GTK_WIDGET (save_dialog));	
}
*/

static void fullscreen_cb ()
{
	// I do not have the patience to implement all things you do to
	// fullscreen for real. This is a faked approximation.
	is_fullscreen = !is_fullscreen;
	if (is_fullscreen)
		gtk_window_fullscreen (main_window);
	else
		gtk_window_unfullscreen (main_window);
    
    gtk_image_view_set_show_cursor (GTK_IMAGE_VIEW (view), !is_fullscreen);
	gtk_image_view_set_show_frame (GTK_IMAGE_VIEW (view), !is_fullscreen);
	gtk_image_view_set_black_bg (GTK_IMAGE_VIEW (view), is_fullscreen);
}

static void
change_zoom_quality_cb (GtkAction      *action,
						GtkRadioAction *current)
{
	if (gtk_radio_action_get_current_value (current))
		gtk_image_view_set_interpolation (GTK_IMAGE_VIEW (view),
                                          GDK_INTERP_BILINEAR);
	else
		gtk_image_view_set_interpolation (GTK_IMAGE_VIEW (view),
                                          GDK_INTERP_NEAREST);
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
    /*
	{
		"Load",
		GTK_STOCK_OPEN,
		"_Load data",
		NULL,
		"Load data",
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
    */
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
		"Fullscreen",
		GTK_STOCK_FULLSCREEN,
		"_Fullscreen Mode",
		"F11",
		"View image in fullscreen",
		G_CALLBACK (fullscreen_cb)
	}
};

static const GtkRadioActionEntry image_tools[] = {
    {
        "DraggerTool",
        GTK_STOCK_REFRESH,
        "_Original",
        NULL,
        "Use the hand tool",
        10
    },
    {
        "SelectorTool",
        GTK_STOCK_MEDIA_PAUSE,
        "_Detected with Th",
        NULL,
        "Use the rectangular selection tool",
        20
    }, 
    {
        "PainterTool",
        GTK_STOCK_MEDIA_PLAY,
        "_Confidence Map",
        NULL,
        "Use the painter tool",
        30
    }
};

gchar *ui_info =
    "<ui>"
    "  <menubar name = 'MenuBar'>"
    "    <menu action = 'FileMenu'>"
    "      <menuitem action = 'Open'/>"
    "      <separator/>" 
    "      <menuitem action = 'Quit'/>"
    "    </menu>"
    "    <menu action = 'EditMenu'>"
    //    "      <menuitem action = 'Transform'/>"
    //    "      <separator/>" 
    "      <menuitem action = 'DraggerTool'/>"
    "      <menuitem action = 'SelectorTool'/>"
    "      <menuitem action = 'PainterTool'/>"
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
    "    <toolitem action='Quit'/>"
    "    <toolitem action='Open'/>"
    "    <separator/>"
    "    <toolitem action='DraggerTool'/>"
    "    <toolitem action='SelectorTool'/>"
    "    <toolitem action='PainterTool'/>"
    "    <separator/>"
    "    <toolitem action='ZoomIn'/>"
    "    <toolitem action='ZoomOut'/>"
	"    <toolitem action='ZoomNormal'/>"
	"    <toolitem action='ZoomFit'/>"
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

void UpdateConfidence() {
    ApplyThreshold();

    update_number_label();
    
	//GdkPixbuf *anim = ApplyThreshold();	
    //    detected_img = anim;

    gtk_image_set_from_pixbuf((GtkImage*)plot_image_view, PlotImagePR());
    gtk_image_set_from_pixbuf((GtkImage*)plot_overall_view, PlotOverallPR());
    //    gtk_image_view_set_pixbuf((GtkImageView*)plot_image_view, PlotPR(filepath), TRUE);
	
	if (currentTool== 10) {
		current_img = gdk_pixbuf_copy(original_img);
		PlotGTCenters(current_img, false);  // normal green
		PlotDetCenters(current_img);
		gtk_image_view_set_pixbuf(view, current_img, FALSE);
	}
    //	else if (currentTool == 20) {
    //		gtk_image_view_set_pixbuf(view, detected_img, FALSE);
    //	}
	else if (currentTool == 30) {
		current_img = gdk_pixbuf_copy(detected_th_img);
		PlotGTCenters(current_img, true); // darkgreen
		PlotDetCenters(current_img);
		gtk_image_view_set_pixbuf(view, current_img, FALSE);
	}
}

void LoadAnotherImage(int mv) {
    int imgN = imgList.size();
    if (imgN == 1) {
        // Stay
    } else {
        // Move
        imgInd = (imgInd + mv) % imgN;
        if (imgInd < 0) {
            imgInd = imgInd + imgN;
        }
		load_filename(imgList[imgInd].c_str());
    }
}

void LoadPrevImage() {
    LoadAnotherImage(-1);
}
void LoadNextImage() {
    LoadAnotherImage(1);
}
void UpdateParticles() {
    const int threshold = gtk_range_get_value((GtkRange*)selected_size) * 255.0 / 100;
    printf("Threshold = %d\n", threshold);

    update_box_file(current_img_name, threshold);
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
	//gtk_window_set_default_size((GtkWindow*)vbox1, 600, -1);
	gtk_box_pack_start((GtkBox*)hbox, vbox1, TRUE, TRUE, 0);
		
	GtkWidget *vbox2 = gtk_vbox_new(FALSE, 0);

    GtkWidget *vbox_threshold = gtk_vbox_new(FALSE, 0);
	sub_label = gtk_label_new(NULL);
	gtk_box_pack_start((GtkBox*)vbox_threshold, sub_label, FALSE, FALSE, 0);
	gtk_label_set_text((GtkLabel*)sub_label, "Choose threshold");
	
	selected_size = gtk_hscale_new_with_range(0, 100, 4);
	gtk_box_pack_start((GtkBox*)vbox_threshold, selected_size, FALSE, FALSE, 0);
	
	GtkWidget *threshold_button = gtk_button_new_with_label("Apply Threshold");
	gtk_box_pack_start((GtkBox*)vbox_threshold, threshold_button, FALSE, FALSE, 0);
	gtk_signal_connect((GtkObject*)threshold_button, "clicked",
					   GTK_SIGNAL_FUNC(UpdateConfidence), NULL);
    
    GtkWidget *vbox_nextimage = gtk_vbox_new(FALSE, 0);
    GtkWidget *prev_image_button = gtk_button_new_with_label("Prev Image");
    GtkWidget *next_image_button = gtk_button_new_with_label("Next Image");
    GtkWidget *update_particle_button = gtk_button_new_with_label("Update particles");
	gtk_box_pack_start((GtkBox*)vbox_nextimage, prev_image_button, TRUE, TRUE, 0);
	gtk_box_pack_start((GtkBox*)vbox_nextimage, next_image_button, TRUE, TRUE, 0);
	gtk_box_pack_start((GtkBox*)vbox_nextimage, update_particle_button, TRUE, TRUE, 0);
	gtk_signal_connect((GtkObject*)prev_image_button, "clicked",
					   GTK_SIGNAL_FUNC(LoadPrevImage), NULL);
	gtk_signal_connect((GtkObject*)next_image_button, "clicked",
					   GTK_SIGNAL_FUNC(LoadNextImage), NULL);
	gtk_signal_connect((GtkObject*)update_particle_button, "clicked",
					   GTK_SIGNAL_FUNC(UpdateParticles), NULL);

    GtkWidget *hbox_top = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start((GtkBox*)hbox_top, vbox_threshold, TRUE, TRUE, 1);
    gtk_box_pack_start((GtkBox*)hbox_top, vbox_nextimage, FALSE, FALSE, 1);

    gtk_box_pack_start((GtkBox*)vbox2, hbox_top, FALSE, FALSE, 0);
    gtk_widget_set_usize(vbox2, 280, -1);

    char tmp_buffer[200];
    sprintf(tmp_buffer, number_template, 0, 0);
    number_label = gtk_label_new(tmp_buffer);
    gtk_box_pack_start((GtkBox*)vbox2, number_label, FALSE, FALSE, 0);
    
    /*	
	GdkColormap *cmap = gdk_colormap_get_system(); 
	pixmap = gdk_pixmap_new(NULL, 280, 280, gdk_colormap_get_visual(cmap)->depth);
	cairo_t *cr = gdk_cairo_create(GDK_DRAWABLE(pixmap)); 
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_rectangle(cr, 0, 0, 280, 280);
    cairo_fill(cr); 
    */

    // plot imageview
    GtkWidget *vbox_tmp1 = gtk_vbox_new(FALSE, 0);
    GtkWidget *hbox_tmp1 = gtk_hbox_new(FALSE, 0);
    GtkWidget *hbox_tmp2 = gtk_hbox_new(FALSE, 0);
    GtkWidget *prec_label = gtk_label_new("precision");
    GtkWidget *recall_label = gtk_label_new("recall (image)");
    gtk_label_set_angle((GtkLabel*)prec_label, 90);
    plot_image_view = gtk_image_new();
    gtk_image_set_from_pixbuf((GtkImage*)plot_image_view, PlotImagePR());
    gtk_box_pack_start((GtkBox*)hbox_tmp1, prec_label, FALSE, FALSE, 0);
    gtk_box_pack_start((GtkBox*)hbox_tmp1, plot_image_view, FALSE, FALSE, 0);
    gtk_box_pack_start((GtkBox*)hbox_tmp2, hbox_tmp1, TRUE, FALSE, 0);
    gtk_box_pack_start((GtkBox*)vbox_tmp1, hbox_tmp2, FALSE, FALSE, 0);
    gtk_box_pack_start((GtkBox*)vbox_tmp1, recall_label, FALSE, FALSE, 0);

    GtkWidget *vbox_tmp2 = gtk_vbox_new(FALSE, 0);
    GtkWidget *hbox_tmp3 = gtk_hbox_new(FALSE, 0);
    GtkWidget *hbox_tmp4 = gtk_hbox_new(FALSE, 0);
    GtkWidget *prec_label2 = gtk_label_new("precision");
    GtkWidget *recall_label2 = gtk_label_new("recall (overall)");
    gtk_label_set_angle((GtkLabel*)prec_label2, 90);
    plot_overall_view = gtk_image_new();
    gtk_image_set_from_pixbuf((GtkImage*)plot_overall_view, PlotOverallPR());
    gtk_box_pack_start((GtkBox*)hbox_tmp3, prec_label2, FALSE, FALSE, 0);
    gtk_box_pack_start((GtkBox*)hbox_tmp3, plot_overall_view, FALSE, FALSE, 0);
    gtk_box_pack_start((GtkBox*)hbox_tmp4, hbox_tmp3, TRUE, FALSE, 0);
    gtk_box_pack_start((GtkBox*)vbox_tmp2, hbox_tmp4, FALSE, FALSE, 0);
    gtk_box_pack_start((GtkBox*)vbox_tmp2, recall_label2, FALSE, FALSE, 0);

    GtkWidget *hbox_tmp5 = gtk_hbox_new(FALSE, 0);    
    gtk_box_pack_start((GtkBox*)hbox_tmp5, vbox_tmp1, FALSE, FALSE, 0);
    gtk_box_pack_start((GtkBox*)hbox_tmp5, vbox_tmp2, FALSE, FALSE, 0);
    gtk_box_pack_start((GtkBox*)vbox2, hbox_tmp5, FALSE, FALSE, 0);
    
    

    // det_cen win
	sub_label = gtk_label_new(NULL);
	gtk_box_pack_start((GtkBox*)vbox2, sub_label, FALSE, FALSE, 0);
	gtk_label_set_text((GtkLabel*)sub_label, "Detected Particles");

    /*
    det_center_view = (GtkWidget*)GTK_IMAGE_VIEW(gtk_image_view_new());
    g_signal_connect(G_OBJECT(det_center_view), "button_press_event",
					  G_CALLBACK (det_button_press_callback), NULL);
	gtk_image_view_set_fitting((GtkImageView*)det_center_view, FALSE);
	gtk_image_view_set_zoom (GTK_IMAGE_VIEW (det_center_view), SELECTED_ZOOM);
    */

    //sub_scroll_win1 = gtk_image_scroll_win_new (GTK_IMAGE_VIEW (det_center_view));
    det_center_box = gtk_vbox_new(FALSE, 0);
    det_scroll_win = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (det_scroll_win),
                                    GTK_POLICY_AUTOMATIC, GTK_POLICY_ALWAYS);
    gtk_scrolled_window_add_with_viewport((GtkScrolledWindow*)det_scroll_win, det_center_box);
    //gtk_window_set_default_size((GtkWindow*)sub_scroll_win1, SELECTED_ZOOM*SELECTED_PER_ROW*MAX_W, 400);
	gtk_box_pack_start((GtkBox*)vbox2, det_scroll_win, TRUE, TRUE, 0);


    /*
    // gt_cen_win
	sub_label = gtk_label_new(NULL);
	gtk_box_pack_start((GtkBox*)vbox2, sub_label, FALSE, FALSE, 0);
	gtk_label_set_text((GtkLabel*)sub_label, "Ground Truth");

    gt_center_view = (GtkWidget*)GTK_IMAGE_VIEW(gtk_image_view_new());
    g_signal_connect(G_OBJECT(gt_center_view), "button_press_event",
					  G_CALLBACK (gt_button_press_callback), NULL);
	gtk_image_view_set_fitting((GtkImageView*)gt_center_view, FALSE);
	gtk_image_view_set_zoom (GTK_IMAGE_VIEW (gt_center_view), SELECTED_ZOOM);

    sub_scroll_win2 = gtk_image_scroll_win_new (GTK_IMAGE_VIEW (gt_center_view));
    //gtk_window_set_default_size((GtkWindow*)sub_scroll_win2, SELECTED_ZOOM*SELECTED_PER_ROW*MAX_W, 400);
    //    g_signal_connect(G_OBJECT(gt_center_view), "button_press_event",
    //					  G_CALLBACK (button_press_callback), NULL);
	gtk_box_pack_start((GtkBox*)vbox2, sub_scroll_win2, TRUE, TRUE, 0);
    */


	
	gtk_widget_set_usize(vbox2, SELECTED_ZOOM*SELECTED_PER_ROW*MAX_W+35, -1);

	gtk_box_pack_start((GtkBox*)hbox, vbox2, FALSE, FALSE, 0);
	
	
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

	view = GTK_IMAGE_VIEW (gtk_image_view_new ());

    dragger = gtk_image_tool_dragger_new (GTK_IMAGE_VIEW (view));

	setup_main_window ();
	
	gtk_widget_show_all (GTK_WIDGET (main_window));
	
    CUR_APP_PATH = argv[0];
    string colormap_path = CUR_APP_PATH.substr(0,CUR_APP_PATH.find_last_of("/"))+"/color.map";

    // Loading colormap
	FILE *fin = fopen(colormap_path.c_str(), "rt");
	fread(color_map, sizeof(short), 64*3, fin);
    /* 
	int i;
	for (i = 0; i < 64*3; i++)
		g_print("%d ", color_map[i]);
    */
	fclose(fin);

    gtk_main ();
	
    return 0;
}
