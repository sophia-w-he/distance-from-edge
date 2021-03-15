// Sophia He
// CSE 455 HW 4
// May 21, 2020
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"

void draw_line(image im, float x, float y, float dx, float dy)
{
    /***********************************************************************
      This function draws a line on an image with color corresponding to the direction of line.
      image im: image to draw line on
      float x, y: starting point of line
      float dx, dy: vector corresponding to line angle and magnitude
    ************************************************************************/
    assert(im.c == 3);
    float angle = 6*(atan2(dy, dx) / TWOPI + .5);
    int index = floor(angle);
    float f = angle - index;
    float r, g, b;
    if(index == 0){
        r = 1; g = f; b = 0;
    } else if(index == 1){
        r = 1-f; g = 1; b = 0;
    } else if(index == 2){
        r = 0; g = 1; b = f;
    } else if(index == 3){
        r = 0; g = 1-f; b = 1;
    } else if(index == 4){
        r = f; g = 0; b = 1;
    } else {
        r = 1; g = 0; b = 1-f;
    }
    float i;
    float d = sqrt(dx*dx + dy*dy);
    for(i = 0; i < d; i += 1){
        int xi = x + dx*i/d;
        int yi = y + dy*i/d;
        set_pixel(im, xi, yi, 0, r);
        set_pixel(im, xi, yi, 1, g);
        set_pixel(im, xi, yi, 2, b);
    }
}

image make_integral_image(image im)
{
    image integ = make_image(im.w, im.h, im.c);
    /***********************************************************************
      This function makes an integral image or summed area table from an image.
      image im: image to process
      returns: image I such that I[x,y] = sum{i<=x, j<=y}(im[i,j])
    ************************************************************************/
    /*first get im(x,y), then I(x-1,y), then I(x,y-1) and then I(x-1,y-1) in that order,
     where x is the inner loop over width and y is the outer loop over height, using zero padding*/
    float s, i, s1, s2, s3;
    // loop over image
    for (int z = 0; z < im.c; z++) {
        for (int y = 0; y < im.h; y++) {
            for (int x = 0; x < im.w; x++) {
                // get i(x,y)
                i = get_pixel(im, x, y, z);
                // get I(x-1, y)
                s1 = (x != 0)? get_pixel(integ, x - 1, y, z) : 0;
                //get I(x, y-1)
                s2 = (y != 0)? get_pixel(integ, x, y - 1, z) : 0;
                // get I(x-1, y-1)
                s3 = (x == 0 || y == 0)? 0 : get_pixel(integ, x - 1, y - 1, z);
                
                // set s(x,y) = i(x,y) + I(x-1, y) + I(x, y-1) - I(x-1, y-1)
                s = i + s1 + s2 - s3;
                
                set_pixel(integ, x, y, z, s);
            }
        }
    }
    return integ;
}

image box_filter_image(image im, int s)
{
    /***********************************************************************
      This function applies a box filter to an image using an integral image for speed.
      image im: image to smooth
      int s: window size for box filter
      returns: smoothed image
    ************************************************************************/
    int i,j,k;
    image integ = make_integral_image(im);
    image S = make_image(im.w, im.h, im.c);
    // w = 1/2 of window size
    int w = s/2;
    //loop over integral image
    for (k = 0; k < integ.c; k++) {
        for (j = 0; j < integ.h; j++) {
            for (i = 0; i < integ.w; i++) {
                
                float sum = 0.0; //initialize sum
                // get lower right, upper right, lower left, and upper left pixels
                float lower_right = get_pixel(integ, i + w, j + w, k);
                float upper_right = get_pixel(integ, i + w, j - w - 1, k);
                float lower_left = get_pixel(integ, i - w - 1, j + w, k);
                float upper_left = get_pixel(integ, i - w - 1, j - w - 1, k);
                // sum = (lower right - upper right) - (lower left - upper left)
                sum = (lower_right - upper_right) - (lower_left - upper_left);
                // divide by window to get average
                sum /= s*s;
                set_pixel(S, i, j, k, sum); //set new value in S
                
            }
        }
    }
    return S;
}

image time_structure_matrix(image im, image prev, int s)
{
    /***********************************************************************
      This function calculates the time-structure matrix of an image pair.
      image im: the input image.
      image prev: the previous image in sequence.
      int s: window size for smoothing.
      returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
               3rd channel is IxIy, 4th channel is IxIt, 5th channel is IyIt.
    ************************************************************************/
    int converted = 0;
    if(im.c == 3){
        converted = 1;
        im = rgb_to_grayscale(im);
        prev = rgb_to_grayscale(prev);
    }

    image S = make_image(im.w, im.h, 5);
    // mage gx and gy filters
    image Ix_filt = make_gx_filter();
    image Iy_filt = make_gy_filter();
    // apply convolution to image with gy and gx to get Ix image and Iy image
    image Ix_im = convolve_image(im, Ix_filt, 0); //Ix
    image Iy_im = convolve_image(im, Iy_filt, 0); //Iy
    // subtract prev image from current im to get the It image
    image It_im = sub_image(im, prev);
    
    float Ix_squared, Iy_squared, Ix, Iy, It, IxIy, IxIt, IyIt;
    for(int y = 0; y < S.h; ++y) { //loop image
        for(int x = 0; x < S.w; ++x) {
            // Ix and Ix^2
            Ix = get_pixel(Ix_im, x, y, 0);
            Ix_squared = Ix*Ix;
            // Iy and Iy^2
            Iy = get_pixel(Iy_im, x, y, 0);
            Iy_squared = Iy*Iy;
            // IxIy
            IxIy = Ix*Iy;
            // It
            It = get_pixel(It_im, x, y, 0);
            // multiply to get IxIt and IyIt
            IxIt = Ix*It;
            IyIt = Iy*It;
            
            // set structure matrix with time components
            set_pixel(S, x, y, 0, Ix_squared); //Ix^2
            set_pixel(S, x, y, 1, Iy_squared); //Iy^2
            set_pixel(S, x, y, 2, IxIy); //IxIy
            set_pixel(S, x, y, 3, IxIt);
            set_pixel(S, x, y, 4, IyIt);
        }
    }
    
    // smooth the image using box filter image
    image s_smooth = box_filter_image(S, s);
    
    // Free the temporary variables
    free_image(Ix_filt);
    free_image(Iy_filt);
    free_image(Ix_im);
    free_image(Iy_im);
    free_image(S);
    
    return s_smooth; // return final image
}


image velocity_image(image S, int stride)
{
    /***********************************************************************
      This function calculates the velocity given a structure image.
      image S: time-structure image
      int stride: only calculate subset of pixels for speed
      returns: velocity of structure image.
    ************************************************************************/
    image v = make_image(S.w/stride, S.h/stride, 3);
    int i, j;
    matrix M = make_matrix(2,2);
    for(j = (stride-1)/2; j < S.h; j += stride){
        for(i = (stride-1)/2; i < S.w; i += stride){
            float Ixx = S.data[i + S.w*j + 0*S.w*S.h];
            float Iyy = S.data[i + S.w*j + 1*S.w*S.h];
            float Ixy = S.data[i + S.w*j + 2*S.w*S.h];
            float Ixt = S.data[i + S.w*j + 3*S.w*S.h];
            float Iyt = S.data[i + S.w*j + 4*S.w*S.h];

            // calculate vx and vy using the flow equation
            // invert and negate matrix M: -M^-1
            // multiply appropriate matrix components by Ixt and Iyt
            float vx = (-Iyy*Ixt + Ixy*Iyt)/(Ixx*Iyy - Ixy*Ixy);
            float vy = (Ixy*Ixt - Ixx*Iyt)/(Ixx*Iyy - Ixy*Ixy);

            set_pixel(v, i/stride, j/stride, 0, vx);
            set_pixel(v, i/stride, j/stride, 1, vy);
        }
    }
    free_matrix(M);
    return v;
}


void draw_flow(image im, image v, float scale)
{
    /***********************************************************************
      This function draws lines on an image given the velocity.
      image im: image to draw on
      image v: velocity of each pixel
      float scale: scalar to multiply velocity by for drawing
    ************************************************************************/
    int stride = im.w / v.w;
    int i,j;
    for (j = (stride-1)/2; j < im.h; j += stride) {
        for (i = (stride-1)/2; i < im.w; i += stride) {
            float dx = scale*get_pixel(v, i/stride, j/stride, 0);
            float dy = scale*get_pixel(v, i/stride, j/stride, 1);
            if(fabs(dx) > im.w) dx = 0;
            if(fabs(dy) > im.h) dy = 0;
            draw_line(im, i, j, dx, dy);
        }
    }
}


// Constrain the absolute value of each image pixel
// image im: image to constrain
// float v: each pixel will be in range [-v, v]
void constrain_image(image im, float v)
{
    int i;
    for(i = 0; i < im.w*im.h*im.c; ++i){
        if (im.data[i] < -v) im.data[i] = -v;
        if (im.data[i] >  v) im.data[i] =  v;
    }
}

image smooth_image(image im, float sigma)
{
    if(1){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // TODO: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.
        return copy_image(im);
    }
}
// Calculate the optical flow between two images
// image im: current image
// image prev: previous image
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// returns: velocity matrix
image optical_flow_images(image im, image prev, int smooth, int stride)
{
    image S = time_structure_matrix(im, prev, smooth);
    image v = velocity_image(S, stride);
    constrain_image(v, 6);
    image vs = smooth_image(v, 2);
    free_image(v);
    free_image(S);
    return vs;
}

// Run optical flow demo on webcam
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// int div: downsampling factor for images from webcam
void optical_flow_webcam(int smooth, int stride, int div)
{
#ifdef OPENCV
    void * cap;
    cap = open_video_stream(0, 0, 1280, 720, 30);
    image prev = get_image_from_stream(cap);
    image prev_c = nn_resize(prev, prev.w/div, prev.h/div);
    image im = get_image_from_stream(cap);
    image im_c = nn_resize(im, im.w/div, im.h/div);
    while(im.data){
        image copy = copy_image(im);
        image v = optical_flow_images(im_c, prev_c, smooth, stride);
        draw_flow(copy, v, smooth*div);
        int key = show_image(copy, "flow", 5);
        free_image(v);
        free_image(copy);
        free_image(prev);
        free_image(prev_c);
        prev = im;
        prev_c = im_c;
        if(key != -1) {
            key = key % 256;
            printf("%d\n", key);
            if (key == 27) break;
        }
        im = get_image_from_stream(cap);
        im_c = nn_resize(im, im.w/div, im.h/div);
    }
#else
    fprintf(stderr, "Must compile with OpenCV\n");
#endif
}
