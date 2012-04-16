/*
 * bmp.h
 * -----
 */

#ifndef __bmp__
#define __bmp__

/* allocate RGB buffer */
void* RGB_Alloc (int width, int height);

/* free the buffer */
void RGB_Free (void *rgb);

/* output an RGB bitmap file (24 bits per pixel) */
void BMP_Output (int width, int height, void *rgb, const char *path);

/* return expected duration */
double AVI_Duration (int frames, int fps);

/* open an RGB avi file (24 bits per pixel, 'fps' frames per second) */
void* AVI_Open (int width, int height, int fps, const char *path);

/* output one frame  */
void AVI_Frame (void *avi, void *rgb);

/* close the avi file */
void AVI_Close (void *avi);

#endif
