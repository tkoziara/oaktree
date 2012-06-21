/*
 * viewer.h
 * --------
 */

#ifndef __viewer__
#define __viewer__

typedef int  (*VIEWER_MENU) (char***, int**); /* menu set up (called first) */
typedef void (*VIEWER_INIT) (void); /* initialisation (called second) */
typedef int  (*VIEWER_IDLE) (void); /* idle callback => return != 0 to cose view update */
typedef void (*VIEWER_QUIT) (void); /* on viewer exit */
typedef void (*VIEWER_RENDER) (void); /* window rendering */
typedef void (*VIEWER_KEY) (int key, int x, int y); /* keyboard */
typedef void (*VIEWER_KEYSPEC) (int key, int x, int y); /* special keyboard */
typedef void (*VIEWER_MOUSE) (int button, int state, int x, int y); /* mouse */
typedef void (*VIEWER_MOTION) (int x, int y); /* mouse motion */
typedef void (*VIEWER_PASSIVE_MOTION) (int x, int y); /* passive mouse motion */

/* create main window */
void viewer (
  int *argc,
  char **argv,
  char *title,
  int width,
  int height,
  REAL *extents,
  VIEWER_MENU menu,
  VIEWER_INIT init,
  VIEWER_IDLE idle,
  VIEWER_QUIT quit,
  VIEWER_RENDER render,
  VIEWER_KEY key,
  VIEWER_KEYSPEC keyspec,
  VIEWER_MOUSE mouse,
  VIEWER_MOTION motion,
  VIEWER_PASSIVE_MOTION passive);

/* redraw all */
void viewer_redraw_all (void);

/* reset scene extents */
void viewer_reset_extents (REAL *extents);

/* update scene extents */
void viewer_update_extents (REAL *extents);

/* get minimal view extent */
REAL viewer_minimal_extent ();

/* get main window sizes */
void viewer_sizes (int *w, int *h);

/* open a 2d window with 1-to-1
 * coordinate-to-pixel mapping;
 * window's content is not exported
 * in bmp and avi format */
int viewer_open_window (
  int x,
  int y,
  int w,
  int h,
  VIEWER_RENDER render);

/* close window */
void viewer_close_window (int window);

/* open a viewport whose content
 * is exported in bmp and avi format */
int viewer_open_viewport (
  int x, /* negtive => maintain relative position when resizing */
  int y, /* -||- */
  int w, /* negative => maintain aspect ratio when resizing */
  int h, /* -||- */
  int is3d,
  VIEWER_RENDER render);

/* move viewport */
void viewer_move_viewport (int viewport, int x, int y, int w, int h);

/* resize viewport */
void viewer_resize_viewport (int viewport, int w, int h);

/* close viewport */
void viewer_close_viewport (int viewport);

/* show tiled text intput window 
 * and return read text by callback */
void viewer_read_text (char *title, void (*done) (char *text));

/* check whether the last text reading is still active */
int viewer_reading_text ();

/* output text at specified coordinates */
enum {FONT_8_BY_13 = 8, FONT_9_BY_15 = 9, FONT_10 = 10, FONT_12 = 12, FONT_18 = 18};
void viewer_print (REAL x, REAL y, REAL z, int font, char *fmt, ...);

/* get width of printed text in pixels */
int viewer_print_width (int font, char *fmt, ...);

/* output screen shot bitmap */
void viewer_screen_bitmap (char *path);

/* take over mouse */
void viewer_hold_mouse ();

/* release mouse takeover */
void viewer_release_mouse ();

/* viewer specific projection matrix (useful for picking) */
void viewer_set_projection_matrix (int w, int h);

/* enable drawing rectangle in screen coordinates */
void viewer_rectangle_on (int x1, int y1, int x2, int y2);

/* disable drawing rectangle */
void viewer_rectangle_off ();

/* stop filming */
void viewer_avi_stop ();

/* set window title */
void viewer_window_title (char *fmt, ...);

/* set trackball center */
void viewer_trackball_center (REAL *point);

#endif
