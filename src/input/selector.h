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
#ifndef __GTKIMAGETOOLSELECTOR_H__
#define __GTKIMAGETOOLSELECTOR_H__

extern "C" {
#include <gtkimageview/gtkiimagetool.h>
#include <gtkimageview/gtkimageview.h>
#include <gtkimageview/mouse_handler.h>
}

G_BEGIN_DECLS

#define GTK_TYPE_IMAGE_TOOL_SELECTOR            (gtk_image_tool_selector_get_type ())
#define GTK_IMAGE_TOOL_SELECTOR(obj)            (G_TYPE_CHECK_INSTANCE_CAST ((obj), GTK_TYPE_IMAGE_TOOL_SELECTOR, GtkImageToolSelector))
#define GTK_IMAGE_TOOL_SELECTOR_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST ((klass), GTK_TYPE_IMAGE_TOOL_SELECTOR, GtkImageToolSelectorClass))
#define GTK_IS_IMAGE_TOOL_SELECTOR(obj)         (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GTK_TYPE_IMAGE_TOOL_SELECTOR))
#define GTK_IS_IMAGE_TOOL_SELECTOR_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), GTK_TYPE_IMAGE_TOOL_SELECTOR))
#define GTK_IMAGE_TOOL_SELECTOR_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj), GTK_TYPE_IMAGE_TOOL_SELECTOR, GtkImageToolSelectorClass))

typedef struct _GtkImageToolSelector GtkImageToolSelector;
typedef struct _GtkImageToolSelectorClass GtkImageToolSelectorClass;

struct _GtkImageToolSelector
{
    GObject             parent;
    GtkImageView       *view;

    /* Cursor to use */
    GdkCursor          *crosshair;

    /* Cache for the image */
    GdkPixbufDrawCache *cache;

    MouseHandler       *mouse_handler;
};

struct _GtkImageToolSelectorClass
{
    GObjectClass parent;
};

GType         gtk_image_tool_selector_get_type     (void);

/* Constructors */
GtkIImageTool *gtk_image_tool_selector_new         (GtkImageView *view);


G_END_DECLS

#endif

