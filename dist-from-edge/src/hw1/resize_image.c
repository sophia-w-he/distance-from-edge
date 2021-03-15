// Sophia He
// CSE 455 Computer Vision
// Homework 1
// April 22, 2020

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"


/***********************************************************************
  We've been talking a lot about resizing and interpolation in class,
  now's your time to do it!
  In order to make grading easier, please only edit the files we mention to submit.
  You will submit the resize_image.c file on Canvas.
************************************************************************/


/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

/***********************************************************************
  This function performs nearest-neighbor interpolation on image "im"
  given a floating column value "x", row value "y" and integer channel "c",
  it interpolates and returns the interpolated value.
  Remember to use the closest "int", not just type-cast because in C that
  will truncate towards zero.
************************************************************************/
float nn_interpolate(image im, float x, float y, int c)
{
    // round x, y to get nearest neighbor pixel
    float f = get_pixel(im, round(x), round(y), c);
    // return interpolated value
    return f;
}

/***********************************************************************
  This function uses nearest-neighbor interpolation on image "im" to
  construct a new image of size "w x h". Hint:
  - Create a new image that is "w x h" and the same number of channels as "im"
  - Loop over the pixels and map back to the old coordinates.
  - Use nearest-neighbor interpolate to fill in the image.
************************************************************************/
image nn_resize(image im, int w, int h)
{
    // Create a new image that is "wxh", same number of channels as "im"
    image resized_im = make_image(w, h, im.c);
    
    // map coordinates using formula aX + b = Y
    // X is new resized image, Y is the old image im
    float y1, y2, y3, x1, x2, x3;
    x1 = -0.5; // new upper left corner
    y1 = -0.5; // old upper left corner
    
    x2 = (float)resized_im.w - 0.5; // new width
    y2 = (float)im.w - 0.5; // old width
    
    x3 = (float)resized_im.h - 0.5; // new height
    y3 = (float)im.h - 0.5; // old height
    
    float a1 = (y1 - y2)/(x1 - x2); //mapping variables for width
    float b1 = y1 - (a1*x1);
    
    float a2 = (y1 - y3)/(x1 - x3); //mapping variables for height
    float b2 = y1 - (a2*x1);
    
    // Loop over the pixels for new and map back to the old coordinates.
    float i, j, k;
    for (k = 0; k < resized_im.c; ++k){
        for(j = 0; j < resized_im.h; ++j){
            for(i = 0; i < resized_im.w; ++i){
                float x = (a1*i + b1); // get old col
                float y = (a2*j + b2); // get old row
                // use nearest neighbor interpolation to get new pixel value
                float nn = nn_interpolate(im, x, y, k);
                // set new image pixel value
                set_pixel(resized_im, i, j, k, nn);
            }
        }
    }
    return resized_im; // return final resized image
}

/***********************************************************************
  This function performs bilinear interpolation on image "im" given
  a floating column value "x", row value "y" and integer channel "c".
  It interpolates and returns the interpolated value.
************************************************************************/
float bilinear_interpolate(image im, float x, float y, int c)
{
    // get the closest upper left, right and lower left, right pixels
    // to value x,y
    float V1 = get_pixel(im, floor(x), floor(y), c);
    float V2 = get_pixel(im, ceil(x), floor(y), c);
    float V3 = get_pixel(im, floor(x), ceil(y), c);
    float V4 = get_pixel(im, ceil(x), ceil(y), c);
    
    // calculate the area of each quadrant relative to x, y
    float A4 = (x - (float)floor(x))*(y - (float)floor(y));
    float A3 = ((float)ceil(x) - x)*(y - (float)floor(y));
    float A2 = (x - (float)floor(x))*((float)ceil(y) - y);
    float A1 = ((float)ceil(x) - x)*((float)ceil(y) - y);
    
    // use bilinear interpolation formula
    float Q = V1*A1 + V2*A2 + V3*A3 + V4*A4;
    
    return Q; // return final value
}

/***********************************************************************
  This function uses bilinear interpolation on image "im" to construct
  a new image of size "w x h". Hint:
  - Create a new image that is "w x h" and the same number of channels as "im".
  - Loop over the pixels and map back to the old coordinates.
  - Use bilinear interpolate to fill in the image.
************************************************************************/
image bilinear_resize(image im, int w, int h)
{
    // Create a new image that is "wxh", same number of channels as "im"
    image resized_im = make_image(w, h, im.c);
    // map coordinates using formula aX + b = Y
    // X is new resized image, Y is the old image im
    float y1, y2, y3, x1, x2, x3;
    x1 = -0.5; // new upper left corner
    y1 = -0.5; // old upper left corner
    
    x2 = (float)resized_im.w - 0.5; // new width
    y2 = (float)im.w - 0.5; // old width
    
    x3 = (float)resized_im.h - 0.5; // new height
    y3 = (float)im.h - 0.5; // old height
    
    float a1 = (y1 - y2)/(x1 - x2); //mapping variables for width
    float b1 = y1 - (a1*x1);
    
    float a2 = (y1 - y3)/(x1 - x3); //mapping variables for height
    float b2 = y1 - (a2*x1);
    
    // Loop over the pixels for new and map back to the old coordinates.
    float i, j, k;
    for (k = 0; k < resized_im.c; ++k){
        for(j = 0; j < resized_im.h; ++j){
            for(i = 0; i < resized_im.w; ++i){
                float x = (a1*i + b1); // get old col
                float y = (a2*j + b2); // get old row
                // use bilinear interpolation to get new pixel val
                float bl = bilinear_interpolate(im, x, y, k);
                // set new image pixel val
                set_pixel(resized_im, i, j, k, bl);
            }
        }
    }
    
    return resized_im; // return final resized image
}

