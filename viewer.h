/*
 * viewer.h
 * --------
 */

#ifndef __viewer__
#define __viewer__

typedef int  (*View_Menu) (char***, int**); /* menu set up (called first) */
typedef void (*View_Init) (void); /* initialisation (called second) */
typedef int  (*View_Idle) (void); /* idle callback => return != 0 to cose view update */
typedef void (*View_Quit) (void); /* on viewer exit */
typedef void (*View_Render) (void); /* window rendering */
typedef void (*View_Key) (int key, int x, int y); /* keyboard */
typedef void (*View_Keyspec) (int key, int x, int y); /* special keyboard */
typedef void (*View_Mouse) (int button, int state, int x, int y); /* mouse */
typedef void (*View_Motion) (int x, int y); /* mouse motion */
typedef void (*View_Passive_Motion) (int x, int y); /* passive mouse motion */

/* create main window */
void viewer (
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
  View_Render render);

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
  View_Render render);

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
int  viewer_print_width (int font, char *fmt, ...);

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
