/* -*- Mode: C; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; coding: utf-8 -*- 
 *
 * Copyright © 2007-2008 Björn Lindqvist <bjourne@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2, or
 * (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */
#include <math.h>
#include "selector.h"
#include "TEXTONSVM.h"

#include <gdk/gdkkeysyms.h>

//#define DEBUG_MODE

extern GtkWidget *view;
extern GtkWidget *selected_size;
extern GtkWidget *sub_scroll_win;
extern bool IMG_FLIPPED;

/*************************************************************/
/***** Static stuff ******************************************/
/*************************************************************/
static gboolean
gtk_image_view_widget_to_image_rect (GtkImageView *view,
                                     GdkRectangle *rect_in,
                                     GdkRectangle *rect_out) {
    GdkRectangle viewport, draw_rect;
    if (!gtk_image_view_get_viewport (view, &viewport))
        return FALSE;
    if (!gtk_image_view_get_draw_rect (view, &draw_rect))
        return FALSE;

    rect_out->x = viewport.x - draw_rect.x + rect_in->x;
    rect_out->y = viewport.y - draw_rect.y + rect_in->y;
    gdouble zoom = gtk_image_view_get_zoom (view);
    //    rect_out->x = (int) ((gdouble) rect_out->x / zoom - rect_in->width/2);
    //    rect_out->y = (int) ((gdouble) rect_out->y / zoom - rect_in->height/2);
    rect_out->x = (int) ((gdouble) rect_out->x / zoom);
    rect_out->y = (int) ((gdouble) rect_out->y / zoom);

	rect_out->width = rect_in->width;
	rect_out->height = rect_in->height;

    // Clip it to the pixbufs area.
    GdkPixbuf *pixbuf = gtk_image_view_get_pixbuf (view);
    GdkRectangle pb_rect = {0, 0,
                            gdk_pixbuf_get_width (pixbuf),
                            gdk_pixbuf_get_height (pixbuf)};
    gdk_rectangle_intersect (&pb_rect, rect_out, rect_out);
    if (!rect_out->width || !rect_out->height)
        return FALSE;

    return TRUE;
}

static void gtk_image_tool_selector_paint (GtkImageToolSelector *selector,
                              GdkRectangle        *rect) {
    GtkImageView *view = selector->view;
    GdkPixbuf *pixbuf = gtk_image_view_get_pixbuf (view);
    guchar *pixels = gdk_pixbuf_get_pixels (pixbuf);
    int stride = gdk_pixbuf_get_rowstride (pixbuf);
    int n_chans = gdk_pixbuf_get_n_channels (pixbuf);
	int y, x, ofs;

    int size = gtk_range_get_value((GtkRange*)selected_size);
	
    short R = 0xFF, G = 0x00, B = 0x00;
    for (y = rect->y - size/2; y < rect->y + size/2; y++) {
		for (x = rect->x - size/2; x < rect->x - size/2 + 5; x++) {
			ofs = y * stride + x * n_chans;
			pixels[ofs] = R;
			pixels[ofs+1] = G;
			pixels[ofs+2] = B;
		}

		for (x = rect->x + size/2 - 5; x < rect->x + size/2; x++) {
			ofs = y * stride + x * n_chans;
			pixels[ofs] = R;
			pixels[ofs+1] = G;
			pixels[ofs+2] = B;
		}
	}
	
	for (x = rect->x - size/2; x < rect->x + size/2; x++) {
		for (y = rect->y - size/2; y < rect->y - size/2 + 5; y++) {
			ofs = y * stride + x * n_chans;
			pixels[ofs] = R;
			pixels[ofs+1] = G;
			pixels[ofs+2] = B;
		}
		
		for (y = rect->y + size/2 - 5; y < rect->y + size/2; y++) {
			ofs = y * stride + x * n_chans;
			pixels[ofs] = R;
			pixels[ofs+1] = G;
			pixels[ofs+2] = B;
		}		
	}
}

static void gtk_image_tool_selector_paint_at (GtkImageToolSelector *selector, int  wx, int wy) {
    GtkImageView *view = selector->view;

    int size = gtk_range_get_value((GtkRange*)selected_size);
    //    GdkRectangle wid_rect = {wx+size/2, wy+size/2, size, size};
    GdkRectangle wid_rect = {wx, wy, size, size};
    GdkRectangle image_rect;
    if (!gtk_image_view_widget_to_image_rect (view, &wid_rect, &image_rect))
        return;

	AddSelectedPoint(image_rect.x, image_rect.y, image_rect.width, image_rect.height, false, IMG_FLIPPED);
	ShowSelected();

    //GtkAdjustment *adj = gtk_scrolled_window_get_vadjustment((GtkScrolledWindow*)sub_scroll_win);
    //GtkAdjustment *adj = ((GtkImageScrollWin*)sub_scroll_win)->vscroll;
    //g_print("%d        %d\n", adj->upper, adj->page_size);
    //    gtk_adjustment_set_value(adj, adj->upper - adj->page_size);
    //    gtk_scrolled_window_set_placement((GtkScrolledWindow*)sub_scroll_win, GTK_CORNER_TOP_LEFT);
	
    gtk_image_tool_selector_paint (selector, &image_rect);
    gtk_image_view_damage_pixels (view, &image_rect);
}

/*************************************************************/
/***** Implementation of the GtkIImageTool interface *********/
/*************************************************************/
static GdkCursor*
cursor_at_point (GtkIImageTool *tool,
                 int            x,
                 int            y)
{
    GtkImageToolSelector *selector = GTK_IMAGE_TOOL_SELECTOR (tool);
    return selector->crosshair;
}

static gboolean
button_press (GtkIImageTool  *tool,
              GdkEventButton *ev)
{
    GtkImageToolSelector *selector = GTK_IMAGE_TOOL_SELECTOR (tool);
    if (ev->button != 1)
        return FALSE;

	gtk_widget_grab_focus (GTK_WIDGET (view));
    gtk_image_tool_selector_paint_at (selector, ev->x, ev->y);

    return mouse_handler_button_press (selector->mouse_handler, ev);
}

static gboolean
button_release (GtkIImageTool  *tool,
                GdkEventButton *ev)
{
    GtkImageToolSelector *selector = GTK_IMAGE_TOOL_SELECTOR (tool);
    return mouse_handler_button_release (selector->mouse_handler, ev);
}

static gboolean
motion_notify (GtkIImageTool  *tool,
               GdkEventMotion *ev)
{
    return FALSE;
}

static void
pixbuf_changed (GtkIImageTool *tool,
                gboolean       reset_fit,
                GdkRectangle  *rect)
{
    GtkImageToolSelector *selector = GTK_IMAGE_TOOL_SELECTOR (tool);
    gdk_pixbuf_draw_cache_invalidate (selector->cache);
}

static void paint_image (GtkIImageTool     *tool,
             GdkPixbufDrawOpts *opts,
             GdkDrawable       *drawable) {
    GtkImageToolSelector *selector = GTK_IMAGE_TOOL_SELECTOR (tool);
    gdk_pixbuf_draw_cache_draw (selector->cache, opts, drawable);
}

/*************************************************************/
/***** Stuff that deals with the type ************************/
/*************************************************************/
static void
gtk_iimage_tool_interface_init (gpointer g_iface,
                                gpointer iface_data)
{
    GtkIImageToolClass *klass = (GtkIImageToolClass *) g_iface;
    klass->cursor_at_point = cursor_at_point;
    klass->button_press = button_press;
    klass->button_release = button_release;
    klass->motion_notify = motion_notify;
    klass->pixbuf_changed = pixbuf_changed;
    klass->paint_image = paint_image;
}

G_DEFINE_TYPE_EXTENDED (GtkImageToolSelector,
                        gtk_image_tool_selector,
                        G_TYPE_OBJECT,
                        0,
                        G_IMPLEMENT_INTERFACE (GTK_TYPE_IIMAGE_TOOL,
                                               gtk_iimage_tool_interface_init));

static void
gtk_image_tool_selector_finalize (GObject *object)
{
    GtkImageToolSelector *selector = GTK_IMAGE_TOOL_SELECTOR (object);
    gdk_pixbuf_draw_cache_free (selector->cache);
    gdk_cursor_unref (selector->crosshair);

    /* Chain up */
    G_OBJECT_CLASS (gtk_image_tool_selector_parent_class)->finalize (object);
}

static void
gtk_image_tool_selector_class_init (GtkImageToolSelectorClass *klass)
{
    GObjectClass *object_class = (GObjectClass *) klass;
    object_class->finalize = gtk_image_tool_selector_finalize;
}

static void
gtk_image_tool_selector_init (GtkImageToolSelector *tool)
{
    tool->crosshair = gdk_cursor_new (GDK_CROSSHAIR);
    tool->cache = gdk_pixbuf_draw_cache_new ();
    tool->mouse_handler = mouse_handler_new (tool->crosshair);
}

static void move_box(GtkWidget* widget, GdkEventKey *event) {
	gint k = event->keyval;

	if ((k == 'a') || (k == 'A')) {
#ifdef DEBUG_MODE
		g_print("left");
#endif
		MoveSelected(-5, 0);
        ShowSelected();
	} else if ((k == 'd') || (k == 'D')) {
#ifdef DEBUG_MODE
		g_print("right");
#endif
		MoveSelected(5, 0);		
        ShowSelected();
	} else if ((k == 'w') || (k == 'W')) {
#ifdef DEBUG_MODE
		g_print("up");
#endif
		MoveSelected(0, -5);		
        ShowSelected();
	} else if ((k == 'x') || (k == 'X')) {
#ifdef DEBUG_MODE
		g_print("down");
#endif
		MoveSelected(0, 5);		
        ShowSelected();
	}
	
	
}

GtkIImageTool*
gtk_image_tool_selector_new (GtkImageView *view)
{
	g_signal_connect(view, "key_press_event", G_CALLBACK(move_box), NULL);
	
    g_return_val_if_fail (view, NULL);
    GtkImageToolSelector *selector = (GtkImageToolSelector *)
        g_object_new (GTK_TYPE_IMAGE_TOOL_SELECTOR, NULL);
    selector->view = view;
    return GTK_IIMAGE_TOOL (selector);
}
