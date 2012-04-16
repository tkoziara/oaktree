/*
 * viewer.h
 * --------
 */

#include "oaktree.h"

#ifndef __viewer__
#define __viewer__

typedef int  (*View_Menu) (char***, int**); /* menu set up (CALLED FIRST) */
typedef void (*View_Init) (void); /* initialisation (CALLED SECOND) */
typedef int  (*View_Idle) (void); /* idle callback => return != 0 to cose view update */
typedef void (*View_Quit) (void); /* on viewer exit */
typedef void (*View_Render) (void); /* window rendering */
typedef void (*View_Key) (int key, int x, int y); /* keyboard */
typedef void (*View_Keyspec) (int key, int x, int y); /* special keyboard */
typedef void (*View_Mouse) (int button, int state, int x, int y); /* mouse */
typedef void (*View_Motion) (int x, int y); /* mouse motion */
typedef void (*View_Passive_Motion) (int x, int y); /* passive mouse motion */

/* create main window */
void GLV (
  int *argc,
  char **argv,
  char *title,
  int width,
  int height,
  REAL *extents,
  View_Menu menu,
  View_Init init,
  View_Idle idle,
  View_Quit quit,
  View_Render render,
  View_Key key,
  View_Keyspec keyspec,
  View_Mouse mouse,
  View_Motion motion,
  View_Passive_Motion passive);


/* redraw all */
void GLV_Redraw_All (void);

/* reset scene extents */
void GLV_Reset_Extents (REAL *extents);

/* update scene extents */
void GLV_Update_Extents (REAL *extents);

/* get minimal view extent */
REAL GLV_Minimal_Extent ();

/* get main window sizes */
void GLV_Sizes (int *w, int *h);

/* open a 2D window with 1-to-1
 * coordinate-to-pixel mapping;
 * window's content is not exported
 * in BMP and AVI format */
int GLV_Open_Window (
  int x,
  int y,
  int w,
  int h,
  View_Render render);

/* close window */
void GLV_Close_Window (int window);

/* open a viewport whose content
 * is exported in BMP and AVI format */
int GLV_Open_Viewport (
  int x, /* negtive => maintain relative position when resizing */
  int y, /* -||- */
  int w, /* negative => maintain aspect ratio when resizing */
  int h, /* -||- */
  int is3D,
  View_Render render);

/* move viewport */
void GLV_Move_Viewport (int viewport, int x, int y, int w, int h);

/* resize viewport */
void GLV_Resize_Viewport (int viewport, int w, int h);

/* close viewport */
void GLV_Close_Viewport (int viewport);

/* show tiled text intput window 
 * and return read text by callback */
void GLV_Read_Text (char *title, void (*done) (char *text));

/* check whether the last text reading is still active */
int GLV_Reading_Text ();

/* output text at specified coordinates */
enum {GLV_FONT_8_BY_13 = 8, GLV_FONT_9_BY_15 = 9,
      GLV_FONT_10 = 10, GLV_FONT_12 = 12, GLV_FONT_18 = 18};
void GLV_Print (REAL x, REAL y, REAL z, int font, char *fmt, ...);

/* get width of printed text in pixels */
int  GLV_Print_Width (int font, char *fmt, ...);

/* output screen shot bitmap */
void GLV_Screen_Bitmap (char *path);

/* take over mouse */
void GLV_Hold_Mouse ();

/* release mouse takeover */
void GLV_Release_Mouse ();

/* viewer specific projection matrix (useful for picking) */
void GLV_SetProjectionMatrix (int w, int h);

/* enable drawing rectangle in screen coordinates */
void GLV_Rectangle_On (int x1, int y1, int x2, int y2);

/* disable drawing rectangle */
void GLV_Rectangle_Off ();

/* stop filming */
void GLV_AVI_Stop ();

/* set window title */
void GLV_Window_Title (char *fmt, ...);

/* set trackball center */
void GLV_Trackball_Center (REAL *point);

#endif
