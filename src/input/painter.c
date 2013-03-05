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
extern "C" {
#include <gtkimageview/gtkimagetoolpainter.h>
}

#include "TEXTONSVM.h"

//extern int* dontcare;

/*************************************************************/
/***** Static stuff ******************************************/
/*************************************************************/
static gboolean
gtk_image_view_widget_to_image_rect (GtkImageView *view,
                                     GdkRectangle *rect_in,
                                     GdkRectangle *rect_out)
{
    GdkRectangle viewport, draw_rect;
    if (!gtk_image_view_get_viewport (view, &viewport))
        return FALSE;
    if (!gtk_image_view_get_draw_rect (view, &draw_rect))
        return FALSE;

    rect_out->x = viewport.x - draw_rect.x + rect_in->x;
    rect_out->y = viewport.y - draw_rect.y + rect_in->y;
    gdouble zoom = gtk_image_view_get_zoom (view);
    rect_out->x = (int) ((gdouble) rect_out->x / zoom);
    rect_out->y = (int) ((gdouble) rect_out->y / zoom);
	
	rect_out->x -= rect_in->width/2;
	rect_out->y -= rect_in->height/2;

	rect_out->width = rect_in->width;
	rect_out->height = rect_in->height;
//    rect_out->width = (int) MAX((gdouble) rect_in->width / zoom, 1);
//    rect_out->height = (int) MAX((gdouble) rect_in->height / zoom, 1);

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

static void
gtk_image_tool_painter_paint (GtkImageToolPainter *painter,
                              GdkRectangle        *rect)
{
    GtkImageView *view = painter->view;
    GdkPixbuf *pixbuf = gtk_image_view_get_pixbuf (view);
    guchar *pixels = gdk_pixbuf_get_pixels (pixbuf);
    int stride = gdk_pixbuf_get_rowstride (pixbuf);
    int n_chans = gdk_pixbuf_get_n_channels (pixbuf);
	int y, x; //, n;


    for (y = rect->y; y < rect->y + rect->height; y++)
        for (x = rect->x; x < rect->x + rect->width; x++)
        {
            int ofs = y * stride + x * n_chans;
			pixels[ofs] = 0x00;
			pixels[ofs+1] = 0xFF;
			pixels[ofs+2] = 0x00;
			
			SetDontCare(x, y, 1);
        }
}

static void
gtk_image_tool_painter_paint_at (GtkImageToolPainter *painter,
                                 int                  wx,
                                 int                  wy)
{
    GtkImageView *view = painter->view;

    GdkPixbuf *pixbuf = gtk_image_view_get_pixbuf(view);

    int height = gdk_pixbuf_get_height(pixbuf);
    int width = gdk_pixbuf_get_width(pixbuf);
    int wh;
    if (height > width)
        wh = width / 15;
    else
        wh = height / 15;

    GdkRectangle wid_rect = {wx, wy, wh, wh};
    GdkRectangle image_rect;
    if (!gtk_image_view_widget_to_image_rect (view, &wid_rect, &image_rect))
        return;

    gtk_image_tool_painter_paint (painter, &image_rect);
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
    GtkImageToolPainter *painter = GTK_IMAGE_TOOL_PAINTER (tool);
    return painter->crosshair;
}

static gboolean
button_press (GtkIImageTool  *tool,
              GdkEventButton *ev)
{
    GtkImageToolPainter *painter = GTK_IMAGE_TOOL_PAINTER (tool);
    if (ev->button != 1)
        return FALSE;

    gtk_image_tool_painter_paint_at (painter, ev->x, ev->y);

    return mouse_handler_button_press (painter->mouse_handler, ev);
}

static gboolean
button_release (GtkIImageTool  *tool,
                GdkEventButton *ev)
{
    GtkImageToolPainter *painter = GTK_IMAGE_TOOL_PAINTER (tool);
    return mouse_handler_button_release (painter->mouse_handler, ev);
}

static gboolean
motion_notify (GtkIImageTool  *tool,
               GdkEventMotion *ev)
{
    GtkImageToolPainter *painter = GTK_IMAGE_TOOL_PAINTER (tool);
    mouse_handler_motion_notify (painter->mouse_handler, ev);
    if (!painter->mouse_handler->dragging)
        return FALSE;

    gtk_image_tool_painter_paint_at (painter, ev->x, ev->y);

    return FALSE;
}

static void
pixbuf_changed (GtkIImageTool *tool,
                gboolean       reset_fit,
                GdkRectangle  *rect)
{
    GtkImageToolPainter *painter = GTK_IMAGE_TOOL_PAINTER (tool);
    gdk_pixbuf_draw_cache_invalidate (painter->cache);
}

static void
paint_image (GtkIImageTool     *tool,
             GdkPixbufDrawOpts *opts,
             GdkDrawable       *drawable)
{
    GtkImageToolPainter *painter = GTK_IMAGE_TOOL_PAINTER (tool);
    gdk_pixbuf_draw_cache_draw (painter->cache, opts, drawable);
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

G_DEFINE_TYPE_EXTENDED (GtkImageToolPainter,
                        gtk_image_tool_painter,
                        G_TYPE_OBJECT,
                        0,
                        G_IMPLEMENT_INTERFACE (GTK_TYPE_IIMAGE_TOOL,
                                               gtk_iimage_tool_interface_init));

static void
gtk_image_tool_painter_finalize (GObject *object)
{
    GtkImageToolPainter *painter = GTK_IMAGE_TOOL_PAINTER (object);
    gdk_pixbuf_draw_cache_free (painter->cache);
    gdk_cursor_unref (painter->crosshair);

    /* Chain up */
    G_OBJECT_CLASS (gtk_image_tool_painter_parent_class)->finalize (object);
}

static void
gtk_image_tool_painter_class_init (GtkImageToolPainterClass *klass)
{
    GObjectClass *object_class = (GObjectClass *) klass;
    object_class->finalize = gtk_image_tool_painter_finalize;
}

static void
gtk_image_tool_painter_init (GtkImageToolPainter *tool)
{
    tool->crosshair = gdk_cursor_new (GDK_CROSSHAIR);
    tool->cache = gdk_pixbuf_draw_cache_new ();
    tool->mouse_handler = mouse_handler_new (tool->crosshair);
}

GtkIImageTool*
gtk_image_tool_painter_new (GtkImageView *view)
{
    g_return_val_if_fail (view, NULL);
    GtkImageToolPainter *painter = (GtkImageToolPainter *)
        g_object_new (GTK_TYPE_IMAGE_TOOL_PAINTER, NULL);
    painter->view = view;
    return GTK_IIMAGE_TOOL (painter);
}
