// Sophia He
// CSE 455 HW 2
// April 28, 2020

// updated March 14, 2020
// to include get_smallest_dist_from_edge

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

/***************************** Box filter *******************************
  We want to create a box filter. We will only use square box filters.
  One way to do this is make an image,fill it in with all 1s, and then
  normalize it.That's what we'll do because the normalization function may
  be useful in the future!
************************************************************************/

/***********************************************************************
  This function divides each value in an image "im" by the sum of all the
  values in the image and modifies the image in place.
************************************************************************/
void l1_normalize(image im)
{
    // initialize sum as 0
    float sum = 0.0;
    // loop through data and add to sum
    int i, j, k;
    for (k = 0; k < im.c; ++k) {
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                float f = get_pixel(im, i, j, k);
                sum += f;
            }
        }
    }
    
    int x, y, z;
    // loop through data and divide each by sum
    for (z = 0; z < im.c; ++z) {
        for(y = 0; y < im.h; ++y){
            for(x = 0; x < im.w; ++x){
                float g = get_pixel(im, x, y, z);
                float div = g/sum;
                set_pixel(im, x, y, z, div);
            }
        }
    }
}

/***********************************************************************
  This function makes a square filter of size "w x w". Hint:
  Change the "make_image" arguments so that it is a square image of
  width = height = w and number of channels = 1, with all entries equal
  to 1. Then use "l1_normalize" to normalize your filter.
************************************************************************/
image make_box_filter(int w)
{
    // make image with width = height = w
    image im = make_image(w, w, 1);
    // loop through data and set to 1
    for(int j = 0; j < w; ++j){
        for(int i = 0; i < w; ++i){
            set_pixel(im, i, j, 0, 1.0);
        }
    }
    
    // normalize the image
    l1_normalize(im);

    return im; // return final box filter
}

/***********************************************************************
im: an image with shape "h x w x c"
filter: a convolution filter with shape "k1 x k2 x k3".
Preserve: an integer, which is usually either 0 or 1.

- If `filter` and `im` have the same number of channels then it's just a normal
convolution. We sum over spatial and channel dimensions and produce a 1 channel image.
UNLESS:
    If `preserve` is set to 1 we should produce an image with the same number of
    channels as the input. This is useful if, for example, we want to run a box
    filter over an RGB image and get out an RGB image. This means each channel in
    the image will be filtered by the corresponding channel in the filter.
If the `filter` only has one channel but `im` has multiple channels we want to
apply the filter to each of those channels. Then we either sum between channels
or not depending on if `preserve` is set.

Also, `filter` better have either the same number of channels as `im` or have one channel.
I check this with an `assert`.
 
Assumptions: l and w of filter are odd
************************************************************************/
image convolve_image(image im, image filter, int preserve)
{
    assert(im.c == filter.c || filter.c == 1);
    
    // initialize filtered image
    image filt_im = make_image(im.w, im.h, im.c);
    
    // initialize start index to center filter around pixel
    int filt_start_h = 0 - (filter.h/2);
    int filt_start_w = 0 - (filter.w/2);
    
    // perform cross correlation over each single channel
    float sum; // initalize sum
    float f = 0.0; // initialize filter pixel
    for (int k = 0; k < im.c; ++k) { // loop image
        for(int j = 0; j < im.h; ++j){
            for(int i = 0; i < im.w; ++i){
                sum = 0.0; // reset sum to 0 each loop
                
                for (int m = 0; m < filter.h; m++) { //loop filter
                    for (int n = 0; n < filter.w; n++){
                        // index to map filter to image
                        int im_index_w = (i + filt_start_w + n);
                        int im_index_h = (j + filt_start_h + m);
                        
                        // get filter pixel
                        if (im.c == filter.c) {
                            // normal convolution
                            f = (get_pixel(filter, n, m, k));
                        } else if (filter.c == 1 && im.c > 1) {
                            // filter each channel
                            f = (get_pixel(filter, n, m, 0));
                        }
                        
                        // multiply im pixel x filter pixel
                        f *= (get_pixel(im, im_index_w, im_index_h, k));
                        sum += f; // add to sum
                    }
                }
                
                // set final sum
                set_pixel(filt_im, i, j, k, sum);
            }
        }
    }
    
    if (preserve == 0 && filt_im.c > 1) {
        // sum to produce 1 channel image
        image sum_im = make_image(im.w, im.h, 1);
        float sum;
        //loop and sum each pixel for each respective channel
        for(int y = 0; y < im.h; ++y){
            for(int x = 0; x < im.w; ++x){
                sum = 0.0;
                for (int z = 0; z < im.c; ++z) {
                    sum += get_pixel(filt_im, x, y, z);
                }
                set_pixel(sum_im, x, y, 0, sum);
            }
        }
        free_image(filt_im);
        return sum_im;
    }
    // otherwise return filt_im
    return filt_im;
}

/***********************************************************************
 //highpass:
 // 0, -1, 0
 // -1, 4, -1
 // 0, -1, 0
************************************************************************/
image make_highpass_filter()
{
    // create 3x3x1
    image high_pass = make_image(3,3,1);
    // set -1's
    set_pixel(high_pass, 0, 1, 0, -1.0);
    set_pixel(high_pass, 1, 0, 0, -1.0);
    set_pixel(high_pass, 2, 1, 0, -1.0);
    set_pixel(high_pass, 1, 2, 0, -1.0);
    // set center to 4
    set_pixel(high_pass, 1, 1, 0, 4.0);
    return high_pass;
}

/***********************************************************************
 // sharpen:
 // 0, -1, 0
 // -1, 5, -1
 // 0, -1, 0
************************************************************************/
image make_sharpen_filter()
{
    // create 3x3x1
    image sharpen = make_image(3,3,1);
    // set -1's
    set_pixel(sharpen, 0, 1, 0, -1.0);
    set_pixel(sharpen, 1, 0, 0, -1.0);
    set_pixel(sharpen, 2, 1, 0, -1.0);
    set_pixel(sharpen, 1, 2, 0, -1.0);
    // set center to 5
    set_pixel(sharpen, 1, 1, 0, 5.0);
    return sharpen;
}

/***********************************************************************
 // emboss:
 // -2, -1, 0
 // -1, 1, 1
 // 0, 1, 2
************************************************************************/
image make_emboss_filter()
{
    // create 3x3x1
    image emboss = make_image(3,3,1);
    // first row
    set_pixel(emboss, 0, 0, 0, -2.0);
    set_pixel(emboss, 1, 0, 0, -1.0);
    // second row
    set_pixel(emboss, 0, 1, 0, -1.0);
    set_pixel(emboss, 1, 2, 0, 1.0);
    set_pixel(emboss, 2, 1, 0, 1.0);
    // third row
    set_pixel(emboss, 1, 1, 0, 1.0);
    set_pixel(emboss, 2, 2, 0, 2.0);
    return emboss;
}

// Question 2.2.1: Which of these filters should we use `preserve = 1` when we run our convolution and which ones should we not? Why?
// Answer: The highpass kernel should have preserve = 0 because it should be applied to a graytone image
//        The sharpen filter preserves rgb color so the convolution should use preserve = 1 so that it
//        is applied to all 3 bands;
//        The emboss filter  also preserves rgb color so it should also use preserve = 1;

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: The three above filters all need post-processing. They need to be clamped to ensure that
//         The values remain between 0 and 1. Because all three amplify details in the image, this
//         range restriction makes it so noise is not amplified


/***********************************************************************
sigma: a float number for the Gaussian.

Create a Gaussian filter with the given sigma. Note that the kernel size
is the next highest odd integer from 6x sigma.
 
The Gaussian filter is a square filter with a single channel.
The length and width are the kernel size mentioned above.
 
 x and y are horiz and vert distances from the center of the filter

Return the Gaussian filter.
************************************************************************/
image make_gaussian_filter(float sigma)
{
    // kernel size is next highest odd integer from 6*sigma
    float scale = 6.0*sigma;
    int six_sigma = ceil(scale);
    if (six_sigma % 2 == 0) (six_sigma += 1);
    
    // center of filter
    int center = six_sigma/2;
    
    // make square filter
    image gaussian = make_image(six_sigma, six_sigma, 1);
    float G; // G(x, y)
    for(int j = 0; j < gaussian.h; ++j) { //loop filter
        for(int i = 0; i < gaussian.w; ++i) {
            float x = i - center; // get x and y
            float y = j - center;
            // set squared variables
            float sig_squared = sigma*sigma;
            float x_squared = x*x;
            float y_squared = y*y;
            // G(x, y) = 1/(2*pi*sigma^2)*e^(x^2+y^2/(-2*sigma^2))
            float in_exp = (x_squared + y_squared)/(-2.0*sig_squared);
            G = exp(in_exp)/(2.0*M_PI*sig_squared);
            // set pixel to G
            set_pixel(gaussian, i, j, 0, G);
        }
    }
    l1_normalize(gaussian); // normalize so that sum to 1

    return gaussian;
}

/***********************************************************************
The input images a and image b have the same height, width, and channels.
Sum the given two images and return the result.
The result image should also have the same height, width, and channels as the inputs.
************************************************************************/
image add_image(image a, image b)
{
    // a and b have same height, width, channels
    assert (a.w == b.w && a.h == b.h && a.c == b.c);
    // initialize image
    image sum_im = make_image(a.w, a.h, a.c);
    // initialize sum
    float sum = 0.0;
    for (int z = 0; z < a.c; ++z) { //loop through images
        for(int y = 0; y < a.h; ++y) {
            for(int x = 0; x < a.w; ++x) {
                // sum pixels a and b an d set sum_im pixel
                sum = get_pixel(a, x, y, z) + get_pixel(b, x, y, z);
                set_pixel(sum_im, x, y, z, sum);
            }
        }
    }
    return sum_im; // return final summed image
}

/***********************************************************************
The input image a and image b have the same height, width, and channels.
Subtract the given two images (a - b) and return the result.
The result image should also have the same height, width, and channels as the inputs.
************************************************************************/
image sub_image(image a, image b)
{
    // a and b have same height, width, channels
    assert (a.w == b.w && a.h == b.h && a.c == b.c);
    image sub_im = make_image(a.w, a.h, a.c); //initialize subtraction
    float sub = 0.0;
    for (int z = 0; z < a.c; ++z) { //loop images
        for(int y = 0; y < a.h; ++y) {
            for(int x = 0; x < a.w; ++x) {
                // subtract b from a, set sub_im pixel
                sub = get_pixel(a, x, y, z) - get_pixel(b, x, y, z);
                set_pixel(sub_im, x, y, z, sub);
            }
        }
    }
    return sub_im; // return final subtracted image
}

/***********************************************************************
Create a 3x3 Sobel Gx filter and return it
************************************************************************/
image make_gx_filter()
{
    // initialize 3x3x1
    image gx = make_image(3,3,1);
    // first row
    set_pixel(gx, 0, 0, 0, -1.0);
    set_pixel(gx, 1, 0, 0, 0.0);
    set_pixel(gx, 2, 0, 0, 1.0);
    // second row
    set_pixel(gx, 0, 1, 0, -2.0);
    set_pixel(gx, 1, 1, 0, 0.0);
    set_pixel(gx, 2, 1, 0, 2.0);
    // third row
    set_pixel(gx, 0, 2, 0, -1.0);
    set_pixel(gx, 1, 2, 0, 0.0);
    set_pixel(gx, 2, 2, 0, 1.0);
    
    return gx;
}

/***********************************************************************
Create a 3x3 Sobel Gy filter and return it
************************************************************************/
image make_gy_filter()
{
    // initialize 3x3x1
    image gy = make_image(3,3,1);
    // first row
    set_pixel(gy, 0, 0, 0, -1.0);
    set_pixel(gy, 1, 0, 0, -2.0);
    set_pixel(gy, 2, 0, 0, -1.0);
    // second row
    set_pixel(gy, 0, 1, 0, 0.0);
    set_pixel(gy, 1, 1, 0, 0.0);
    set_pixel(gy, 2, 1, 0, 0.0);
    // third row
    set_pixel(gy, 0, 2, 0, 1.0);
    set_pixel(gy, 1, 2, 0, 2.0);
    set_pixel(gy, 2, 2, 0, 1.0);
    
    return gy;
}

/***********************************************************************
im is the input image with "h x w x 1".
Apply Sobel filter to the given image, get the magnitude and gradient,
and return the result.

Hint: the "calloc" function can allocate the memory for your output. You can
assess the first image (magnitute) by calling rst[0] and the second image
by calling rst[1]
************************************************************************/
image *sobel_image(image im)
{
    // init array with magnitude and orientation(aka theta)
    image *rst = calloc(2, sizeof(image));
    // mage gx and gy filters
    image gx_filt = make_gx_filter();
    image gy_filt = make_gy_filter();
    // apply convolution to image with gy and gx
    image gx_im = convolve_image(im, gx_filt, 0);
    image gy_im = convolve_image(im, gy_filt, 0);
    // initialize magnitude, orientation for rst[0], rst[1]
    image magnitude = make_image(im.w, im.h, 1);
    image orientation = make_image(im.w, im.h, 1);
    
    // initialize variables to be used in mag and theta formulas
    float gx_squared, gy_squared, gx, gy;
    float mag, theta;
    for(int y = 0; y < im.h; ++y) { //loop image
        for(int x = 0; x < im.w; ++x) {
            // gx and gx^2
            gx = get_pixel(gx_im, x, y, 0);
            gx_squared = gx*gx;
            //gy and gy^2
            gy = get_pixel(gy_im, x, y, 0);
            gy_squared = gy*gy;
            // magnitude = sqare root of gx^2 + gy^2
            mag = sqrt(gx_squared + gy_squared);
            // theta = arctan of gy/gx
            theta = atan2(gy, gx);
            // set magnitude and orientation
            set_pixel(magnitude, x, y, 0, mag);
            set_pixel(orientation, x, y, 0, theta);
        }
    }
    
    normalize_image(magnitude);
    normalize_image(orientation);
    
    rst[0] = magnitude; // assign final images to rst
    rst[1] = orientation;
    
    free_image(gx_im); // free images
    free_image(gy_im);

    return rst; // return final rst
}

/***********************************************************************
Calculate minimum and maximum pixel values. Normalize the image by
subtracting the minimum and dividing by the max-min difference.
This is a helper function to visualize the sobel magnitude image better.
No TODO here :)
***********************************************************************/
void normalize_image(image im)
{
    int i;
    float min = im.data[0];
    float max = im.data[0];
    for(i = 0; i < im.w*im.h*im.c; ++i){
        if(im.data[i] > max) max = im.data[i];
        if(im.data[i] < min) min = im.data[i];
    }
    for(i = 0; i < im.w*im.h*im.c; ++i){
        im.data[i] = (im.data[i] - min)/(max-min);
    }
}

/*
 Applies the sobel edge detection filter to image im
 Finds the edge point closest to point x, y on image im
 */
image get_smallest_dist_from_edge(image im, int x, int y)
{
    image *rst = sobel_image(im);
    im = rst[0];
    float min_distance = 99999999;
    int best_x = 0;
    int best_y = 0;
    
    // loop through images
    for(int j = 0; j < im.h; ++j) {
        for(int i = 0; i < im.w; ++i) {
            // get pixel
            float pr = get_pixel(rst[0], i, j, 0);
            // if pixel is greater than threshold
            if (pr > .10) {
                // calculate distance between pixel coordinates and x,y
                float dist = sqrt(pow(i - x, 2) + pow(j - y, 2));
                // if distance is less than current minimum
                if (dist < min_distance) {
                    // and update min_distance and best index as necessary
                    min_distance = dist;
                    best_x = i;
                    best_y = j;
                }
                
            } else {
            }
        }
    }
    
    printf("\nmin dist \n");
    printf("%f", min_distance);
    printf("\nbest edge x \n");
    printf("%d", best_x);
    printf("\nbest edge y \n");
    printf("%d", best_y);
    printf("\n");
    
    for(int m = -5; m <= 5; ++m){
        // Set R,G,B values of the horizontal neighbors
        set_pixel(im, x+m, y, 0, 1);
        set_pixel(im, x+m, y, 1, 0);
        set_pixel(im, x+m, y, 2, 1);
        // Set R,G,B values of the vertical neighbors
        set_pixel(im, x, y+m, 0, 1);
        set_pixel(im, x, y+m, 1, 0);
        set_pixel(im, x, y+m, 2, 1);
    }
    
    for(int m = -9; m <= 9; ++m){
        // Set R,G,B values of the horizontal neighbors
        set_pixel(im, best_x+m, best_y, 0, 1);
        set_pixel(im, best_x+m, best_y, 1, 0);
        set_pixel(im, best_x+m, best_y, 2, 1);
        // Set R,G,B values of the vertical neighbors
        set_pixel(im, best_x, best_y+m, 0, 1);
        set_pixel(im, best_x, best_y+m, 1, 0);
        set_pixel(im, best_x, best_y+m, 2, 1);
    }
    
    return im;
    
}

// EXTRA CREDITS BELOW
int compare_float(const void * a, const void * b)
{
    // This function is provided for your convenience
    float fa = *(const float*) a;
    float fb = *(const float*) b;
    return (fa > fb) - (fa < fb);
}

/***********************************************************************
im is the input image.
kernel_size is a positive odd number.

We assume a median filter is a square, with the same height and width.
The kernel size is always a positive odd number. We use "clamp" padding
for borders and corners. The output image should have the same width,
height, and channels as the input image. You should apply median filter
to each channel of the input image `im`.

Hint: use the qsort() function to sort an array. Make use of compare_float() as needed.
************************************************************************/
image apply_median_filter(image im, int kernel_size)
{
    image out = make_image(im.w, im.h, im.c);

    // for centering filter, assume square filter
    int filt_start_w = 0 - (kernel_size/2);
    // loop image
    for (int k = 0; k < im.c; ++k) {
        for(int j = 0; j < im.h; ++j){
            for(int i = 0; i < im.w; ++i){
                // initialize array for sorting
                float med_array[kernel_size*kernel_size];
                int index = 0; // to index array
                // loop filter centered around pixel i,j
                for (int m = 0; m < kernel_size; m++) {
                    for (int n = 0; n < kernel_size; n++){
                        // index for image w and h to map to filter
                        int im_index_w = (i + filt_start_w + n);
                        int im_index_h = (j + filt_start_w + m);
                        // get each pixel from image and append to array
                        float f = get_pixel(im, im_index_w, im_index_h, k);
                        med_array[index] = f;
                        index++; // increase array index
                    }
                }
                // sort array
                qsort(med_array, kernel_size*kernel_size, sizeof(float), compare_float);
                // get median
                float median = med_array[kernel_size*kernel_size/2];
                // set pixel for final image to median
                set_pixel(out, i, j, k, median);
            }
        }
    }

    return out;
}
