// Sophia He
// CSE 455 Spring 2020
// Homework 0
// April 14, 2020

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "image.h"

// Returns the pixel value at column x, row y, and channel c for image im
float get_pixel(image im, int x, int y, int c)
{
    // initialize variables for padding
    int X, Y, Z;
    X = x;
    Y = y;
    Z = c;
    
    // clamp padding: if col, row, or channel value < 0, use 0
    // if greater than image dimension use (dimension - 1)
    if (x >= im.w) X = im.w - 1;
    if (x < 0) X = 0;
    
    if (y >= im.h) Y = im.h - 1;
    if (y < 0) Y = 0;
    
    if (c >= im.c) Z = im.c - 1;
    if (c < 0) Z = 0;
    
    // find the index using new X, Y, Z values
    int index = X + im.w*Y + Z*im.w*im.h ;

    // return the value in the data array
    return im.data[index];
}

// sets the pixel at coordinates x, y, z to the value v in image im
void set_pixel(image im, int x, int y, int c, float v)
{
    // find the index of the pixel
    int index = (x + y*im.w + c*im.w*im.h);
    // if within dimensions
    if (x < im.w && y < im.h && c < im.c) {
        if (x >= 0 && y >= 0 && c >= 0) {
            // set value in data array to v
            im.data[index] = v;
        }
    }
}

// returns a new copy of image im
image copy_image(image im)
{
    // initialize image
    image copy = make_image(im.w, im.h, im.c);
    // loop through pixels for image im
    // set same pixels for copy image
    int i, j, k;
    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                set_pixel(copy, i, j, k, get_pixel(im, i, j, k));
            }
        }
    }
    // return the final copy
    return copy;
}

// converts colored image im with three channels to grayscale
// returns new image with singular channel
image rgb_to_grayscale(image im)
{
    // assert 3 channels (r, g, b)
    assert(im.c == 3);
    // initialize image with single channel
    image gray = make_image(im.w, im.h, 1);
    
    // loop and get r, g, b pixels from im
    float r, g, b;
    int i, j;
    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            r = get_pixel(im, i, j, 0);
            g = get_pixel(im, i, j, 1);
            b = get_pixel(im, i, j, 2);
            // weighted sum
            float Y = 0.299*r + 0.587*g + 0.114*b;
            // set gray image pixel to weighted sum value
            set_pixel(gray, i, j, 0, Y);
        }
    }
    // return final gray image
    return gray;
}

// adds constant factor v to channel c in image im
void shift_image(image im, int c, float v)
{
    // loop through all row, col pixels for channel c
    int i, j;
    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            // set each pixel to its prior val + v
            set_pixel(im, i, j, c, (get_pixel(im, i, j, c) + v));
        }
    }
}

// makes sure the pixel values in the image stay between 0 and 1.
// Implements clamping on the image so that any value below 0 gets set to 0
// and any value above 1 gets set to 1.
void clamp_image(image im)
{
    // loop through every pixel value in image im
    int i, j, k;
    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                // if pixel val is greater than 1, set it to 1
                if (get_pixel(im, i, j, k) > 1) {
                    set_pixel(im, i, j, k, 1);
                }
                // if pixel val is less than 0, set it to 0
                if (get_pixel(im, i, j, k) < 0) {
                    set_pixel(im, i, j, k, 0);
                }
            }
        }
    }
}

// returns the max of a, b, and c
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

// returns the min of a, b, and c
float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

// converts an image in red, green, blue colorspace to hue, saturation, value
void rgb_to_hsv(image im)
{
    // assert 3 channels (r, g, b)
    assert(im.c == 3);
    // loop through all row, col pixels in image
    float r, g, b;
    int i, j;
    for(j = 0; j < im.h; j++){
        for(i = 0; i < im.w; i++){
            // get the red, green, blue from the three channels
            r = get_pixel(im, i, j, 0);
            g = get_pixel(im, i, j, 1);
            b = get_pixel(im, i, j, 2);
            // find the max of the three pixels
            // set value V equal to max
            float V = three_way_max(r, g, b);
            // find the min of the three pixels
            float m = three_way_min(r, g, b);
            // dif between the min and max
            float C = V - m;
            // set saturation value
            float S;
            if (V == 0) { // avoid divide by 0
                S = 0;
            } else {
                S = C / V;
            }
            
            // calculate hue
            float H = 0.0;
            // intermediate formulas
            float H1 = 0.0; //H'
            if (C == 0) {
                H = 0.0;
            } else if (V == r) {
                H1 = (g - b)/C;
            } else if (V == g) {
                H1 = (b - r)/C + 2.0;
            } else if (V == b) {
                H1 = (r - g)/C + 4.0;
            }
            // set hue
            if (H1 < 0) {
                H = (H1/6.0) + 1.0;
            } else {
                H = H1/6.0;
            }
            // H = [0,1) - circle around if not in range
            if (H >= 1) {
                H = H - 1.0;
            }
            if (H < 0) {
                H = H + 1.0;
            }
            
            // set hue, saturation, value
            set_pixel(im, i, j, 0, H);
            set_pixel(im, i, j, 1, S);
            set_pixel(im, i, j, 2, V);
        }
    }
}

// converts an image in r hue, saturation, value colorspace to red, green, blue
void hsv_to_rgb(image im)
{
    // assert 3 channels (r, g, b)
    assert(im.c == 3);
    // loop through all row, col pixels in image
    float h, s, v;
    float p, q, t;
    int i, j;
    for(j = 0; j < im.h; j++){
        for(i = 0; i < im.w; i++){
            // get hue, saturation, value
            h = get_pixel(im, i, j, 0);
            s = get_pixel(im, i, j, 1);
            v = get_pixel(im, i, j, 2);
            
            // intermediate conversion formulas
            float H = h*6.0;
            int h_i = floor(H);
            float f = H - ((float)(h_i));
            p = v*(1.0 - s);
            q = v*(1.0 - f*s);
            t = v*(1.0 - ((1.0 - f)*s));
                        
            // initialize red, green, blue
            float r = 0.0;
            float g = 0.0;
            float b = 0.0;
            
            // 2D conversion array
            // the red, green, blue variables for
            // each potential value of h_i
            float h_i_formulas[6][3] = {
                {v, t, p},
                {q, v, p},
                {p, v, t},
                {p, q, v},
                {t, p, v},
                {v, p, q}
            };
            
            // match red, green, blue to appropriate variable
            // using the 2D conversion array
            for (int index = 0; index < 6; index++) {
                if (h_i == index) {
                    r = h_i_formulas[index][0];
                    g = h_i_formulas[index][1];
                    b = h_i_formulas[index][2];
                }
            }
            
            // set red, green, blue
            set_pixel(im, i, j, 0, r);
            set_pixel(im, i, j, 1, g);
            set_pixel(im, i, j, 2, b);
            
        }
    }
}

// scales channel c in image im by constant factor v
void scale_image(image im, int c, float v) {
    // loop through all row, col pixels for channel c
    int i, j;
    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            // set each pixel to its prior val * v
            set_pixel(im, i, j, c, (get_pixel(im, i, j, c) * v));
        }
    }
}
