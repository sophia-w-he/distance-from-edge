#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "image.h"
#include "matrix.h"


image make_1d_gaussian(float sigma)
/**************************************************
Create a separable (1-D) Gaussian filter of size given by input "sigma"
TODO:
Fill this for extra credit
**************************************************/
{
    return make_image(1,1,1);
}


image structure_matrix(image im, float sigma)
/**************************************************
Create the structure matrix for the image "im".
 **************************************************/
{
    // mage gx and gy filters
    image Ix_filt = make_gx_filter();
    image Iy_filt = make_gy_filter();
    // apply convolution to image with gy and gx
    image Ix_im = convolve_image(im, Ix_filt, 0); //Ix
    image Iy_im = convolve_image(im, Iy_filt, 0); //Iy
    
    // Create an empty image, "S" of size im.w x im.h x 3
    image S = make_image(im.w, im.h, 3);
    
    float Ix_squared, Iy_squared, Ix, Iy, IxIy;
    for(int y = 0; y < S.h; ++y) { //loop image
        for(int x = 0; x < S.w; ++x) {
            // gx and gx^2
            Ix = get_pixel(Ix_im, x, y, 0);
            Ix_squared = Ix*Ix;
            //gy and gy^2
            Iy = get_pixel(Iy_im, x, y, 0);
            Iy_squared = Iy*Iy;
            
            IxIy = Ix*Iy;
            
            // set magnitude and orientation
            set_pixel(S, x, y, 0, Ix_squared); //Ix^2
            set_pixel(S, x, y, 1, Iy_squared); //Iy^2
            set_pixel(S, x, y, 2, IxIy); //IxIy
        }
    }
    
    // make gaussian filter
    image gauss_filt = make_gaussian_filter(sigma);
    // convolve image
    image s_smooth = convolve_image(S, gauss_filt, 1);
    
    // Free the temporary variables
    free_image(Ix_filt);
    free_image(Iy_filt);
    free_image(Ix_im);
    free_image(Iy_im);
    free_image(gauss_filt);
    free_image(S);
    
    //Return the smoothed image
    return s_smooth;
}


image cornerness_response(image S)
/**************************************************
Estimate the cornerness response "R" of each pixel given a structure matrix "S".
Return response map of the same size as structure matrix
**************************************************/
{
    
    image R = make_image(S.w, S.h, 1);
    float determinant, trace, R_val;
    float Ix_squared, Iy_squared, IxIy;
    // Loop over each pixel position in "S"
    for(int y = 0; y < S.h; ++y) { //loop image
        for(int x = 0; x < S.w; ++x) {
            // gx and gx^2
            //Ix = get_pixel(S, x, y, 0);
            Ix_squared = get_pixel(S, x, y, 0);
            //gy and gy^2
            //Iy = get_pixel(Iy_im, x, y, 0);
            Iy_squared = get_pixel(S, x, y, 1);
            
            IxIy = get_pixel(S, x, y, 2);
            
            // compute determinant = Ix^2 * Iy^2 - IxIy * IxIy
            determinant = (Ix_squared*Iy_squared) - (IxIy*IxIy);
            // compute trace = Ix^2 + Iy^2
            trace = Ix_squared + Iy_squared;
            // compute R = determinant - alpha * trace^2, alpha = .06
            R_val = determinant - (0.06*trace*trace);
            
            // set magnitude and orientation
            set_pixel(R, x, y, 0, R_val);
        }
    }
    
    return R;
}


image nms_image(image R, int w)
/**************************************************
Perform non-max supression on the response map "R".
**************************************************/
{
    
    image r = copy_image(R);
    
    // for every pixel in "R":
    float f, r_pix; // initialize filter pixel
    for (int k = 0; k < r.c; ++k) { // loop image
        for(int j = 0; j < r.h; ++j){
            for(int i = 0; i < r.w; ++i){
                
                r_pix = get_pixel(R, i, j, k);
                
                // for every neighboring pixel within "w":
                for (int m = (j - w); m <= (j + w); m++) { //loop filter
                    for (int n = (i - w); n <= (i + w); n++){
                        
                        // get filter pixel
                        f = get_pixel(R, n, m, k);
                        // if neighboring pixel value > current pixel value:
                        if (f > r_pix) {
                            // set current pixel value to be very low (= -999999)
                            set_pixel(r, i, j, k, -999999);
                        }
                        
                    }
                }
                
            }
        }
    }
    
    return r;
}


descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
/**************************************************
Perform Harris corner detection by computing a feature descriptor for each corner in image "im".
float sigma: std. dev for gaussian filter used in structure matrix
float thresh: threshold for detecting corners from response map
int nms: size of the neighborhood to look for local-maxes in response map
int *n: pointer to number of corners detected
Return array d of descriptors of corners in the input image.
**************************************************/
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);
    // Estimate cornerness
    image R = cornerness_response(S);
    // Run NMS on the responses
    image Rnms = nms_image(R, nms);

    // Initialize the count to 0
    int count = 0; // the number of corners
    float val;
    // For each pixel in Rnms:
    for (int k = 0; k < Rnms.c; ++k) { // loop image
        for(int j = 0; j < Rnms.h; ++j){
            for(int i = 0; i < Rnms.w; ++i){
                // if pixel value > "thresh":
                val = get_pixel(Rnms, i, j, k);
                if (val > thresh) count++; // increase value of count
            }
        }
    }

    // array "d" of descriptors of corners
    descriptor *d = calloc(count, sizeof(descriptor));
    int index = 0;

    // For each pixel in Rnms:
    for (int k = 0; k < Rnms.c; ++k) { // loop image
        for(int j = 0; j < Rnms.h; ++j){
            for(int i = 0; i < Rnms.w; ++i){
                val = get_pixel(Rnms, i, j, k);
                // if pixel value > "thresh":
                if (val > thresh) {
                    //  get descriptor for the current pixel
                    descriptor des = make_descriptor(im, (i + Rnms.w*j + Rnms.w*Rnms.h*k));
                    // update the array "d" with this descriptor
                    d[index] = des;
                    index++;
                }
            }
        }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);
    // set *n equal to number of corners in image "im"
    *n = count;
    return d;
}


match *match_descriptors(descriptor *a, int an, descriptor *b, int bn, int *mn)
/**************************************************
Finds best matches between descriptors of two images.
descriptor *a, *b: array of descriptors in two images
int an, bn: number of descriptors in arrays "a" and "b"
int *mn: pointer to number of matches found
Return an array of matches. Details inline.
The "match" structure is defined in "image.h"
**************************************************/
{
    int i,j;

    // Initialize array of matches. We will have at most "an" matches
    *mn = an;
    match *m = calloc(an, sizeof(match));

    // For each descriptor in a..
    for(j = 0; j < an; ++j){
        // Initialize best index having minimum distance
        float min_distance = 99999999;
        int best_index = 0;
        // For each descriptor in b..
        for(i = 0; i < bn; ++i){
            // Compute L1 distance between the two descriptors
            float dist = 0.0;
            for (int x = 0; x < a[j].n; x++) {
                //L1 distance = total sum of abs(a_i - b_i)
                dist += fabs(a[j].data[x] - b[i].data[x]);
            }
            if (dist < min_distance) {
                // and update min_distance and best_index as necessary
                min_distance = dist;
                best_index = i;
            }
        }
        // Update the match
        m[j].ai = j;
        m[j].bi = best_index;
        m[j].p = a[j].p;
        m[j].q = b[best_index].p;
        m[j].distance = min_distance;
    }

    int count = 0;

    // Initialize an array "seen" to all zeros
    // It keeps track of whether a descriptor in "*b" is already matched or not
    int *seen = calloc(bn, sizeof(int));

    // Sort matches based on distance
    qsort(m, *mn, sizeof(match), match_compare);

    // Find one-to-one matches.
    // Each point is only be a part of one match. Some points will not be in a match.
    // Loop over all the *mn matches
    for(j = 0; j < *mn; ++j){
       // if the b-index of the current match is not seen:
        if (seen[m[j].bi] == 0) {
           // set it to seen (make the corresponding value in "seen" to 1)
            seen[m[j].bi] = 1;
           // assign the current match to m[count] and then update count
            m[count] = m[j];
            count++;
        }
    }
    
    // Update the number of final matches
    *mn = count;
    free(seen);
    return m; // Return an array of matches.
}


point project_point(matrix H, point p)
/**************************************************
Apply a projective transformation using homography "H" to the point "p".
The matrix functions are defined in src/matrix.c
Return the projected point "q".
**************************************************/
{
    // Create c matrix using point "p"
    matrix c = make_matrix(3, 1);
    // Fill in "c" matrix with x-coordinate of p, y-coordinate of p, and 1
    c.data[0][0] = p.x; // x-coordinate of p
    c.data[1][0] = p.y; // y-coordinate of p
    c.data[2][0] = 1;   // 1
    
    // Multiply "H" with the created c matrix
    matrix m = matrix_mult_matrix(H, c);

    point q;
    // Assign x,y coordinates of "q" using homogenous coordinates from "m"
    q.x = m.data[0][0]/m.data[2][0];
    q.y = m.data[1][0]/m.data[2][0];

    free_matrix(c);
    free_matrix(m);
    return q; // Return the projected point "q".
}


int model_inliers(matrix H, match *m, int n, float thresh)
/**************************************************
Count number of inliers in an array "m" of "n" matches.
Also bring inliers to the front of array "m".
matrix H: homography between coordinate systems
float thresh: threshold for a match to be an inlier
Note: In this way, you are sorting the matches so that the inliers are the first #'count' elements.
**************************************************/
{
    int count = 0;
    //initialize count to 0
    //Loop over all matches starting from the end (since matches are already sorted by distance)
    for(int i = n - 1; i >= count; i--) {
        //Project the point "p" in the current match using "H" (use project_point())
        point p = project_point(H, m[i].p);
        //compute L2 distance between the projected point and the point "q" in the current match
        float dist = sqrt(pow(m[i].q.x - p.x, 2) + pow(m[i].q.y - p.y, 2));
        //if L2 distance < thresh:
        if (dist < thresh) {
            //swap m[count] with the current match
            match temp = m[count];
            m[count] = m[i];
            m[i] = temp;
            //update count
            count++;
            i++;
        }
    }
    return count; // Return count (i.e. number of matches that are inliers)
}


matrix RANSAC(match *m, int n, float thresh, int iter, int cutoff)
/**************************************************
Performs RANdom SAmple Consensus to calculate homography for noisy matches.
Returns matrix representing most common homography between matches.
Inputs:
array "m" of "n" matches
"thresh": threshold for inlier modeling
"iter": number of iterations
"cutoff": inlier cutoff to exit early
return the best homography
**************************************************/
{
    // Initializations
    int max_inliers = 0;
    matrix best_H = make_translation_homography(256, 0);
    // for "iter" iterations
    for (int i = 0; i < iter; i++) {
        //compute a homography with 4 matches using compute_homography() defined in panorama_helpers.c
        //(Computes homography between two images using "n" random matches from an array "matches" of "mn" matches.)
        matrix H = compute_homography(m, n, 4);
        //if no homography, continue, else
        if (H.data) {
        //compute #inliers in the matches using model_inliers() and the computed homography
            int in = model_inliers(H, m, n, thresh);
        //if #inliers > max_inliers (i.e. new homography is better than old):
            if (in > max_inliers) {
            //compute updated homography using all inliers
                matrix H_update = compute_homography(m, n, in);
                best_H = H_update;
            //if no homography, continue, else
                if (H_update.data) {
            //compute updated #inliers using updated homography and assign it to max_inliers
                    int in_update = model_inliers(H_update, m, n, thresh);
                    max_inliers = in_update;
            //if updated #inliers > cutoff:
                    if (max_inliers > cutoff) {
                //return the best homography
                        return best_H;
                    }
                }
            }
        }
    }
    
    return best_H; // return best homography
}


image combine_images(image a, image b, matrix H)
/**************************************************
Stitches two images "a" and "b" together using a homography "H" and returns the stitched image.
Paste image a and projected image b into a combined_image and return it.
**************************************************/
{
    matrix Hinv = matrix_invert(H);

    // Project the corners of image b into image a coordinates.
    point c1 = project_point(Hinv, make_point(0,0));
    point c2 = project_point(Hinv, make_point(b.w-1, 0));
    point c3 = project_point(Hinv, make_point(0, b.h-1));
    point c4 = project_point(Hinv, make_point(b.w-1, b.h-1));

    // Find top left and bottom right corners of image b warped into image a.
    point topleft, botright;
    botright.x = MAX(c1.x, MAX(c2.x, MAX(c3.x, c4.x)));
    botright.y = MAX(c1.y, MAX(c2.y, MAX(c3.y, c4.y)));
    topleft.x = MIN(c1.x, MIN(c2.x, MIN(c3.x, c4.x)));
    topleft.y = MIN(c1.y, MIN(c2.y, MIN(c3.y, c4.y)));

    // Find how big our new image should be and the offsets from image a.
    int dx = MIN(0, topleft.x);
    int dy = MIN(0, topleft.y);
    int w = MAX(a.w, botright.x) - dx;
    int h = MAX(a.h, botright.y) - dy;

    // Can disable this if you are making very big panoramas.
    // Usually this means there was an error in calculating H.
    if(w > 7000 || h > 7000){
        fprintf(stderr, "output too big, stopping\n");
        return copy_image(a);
    }

    int i,j,k;
    image combined_img = make_image(w, h, a.c);
    
    // Loop over all pixels in image a
    for (k = 0; k < a.c; ++k) {
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                // get current pixel value
                float f = get_pixel(a, i, j, k);
                // assign it to corresponding pixel offset by "dx" and "dy" in "combined_img"
                set_pixel(combined_img, i - dx, j - dy, k, f);
                
            }
        }
    }

    // Loop over all pixels lying between "topleft" and "botright"
    for (int loop3 = 0; loop3 < b.c; loop3++) {
        for (int loop2 = topleft.y; loop2 <= botright.y; loop2++)  {
            for (int loop1 = topleft.x; loop1 <= botright.x; loop1++) {
            //create a point with the current x,y coordinates using make_point() defined in panorama_helpers.c
                point p = make_point(loop1, loop2);
                //project the point to image b coordinate system using homography "H"
                point proj_p = project_point(H, p);
                // if the projected point's coordinates lie within the bounds for image b:
                if (0 <= proj_p.x && proj_p.x < b.w) {
                    if (0 <= proj_p.y && proj_p.y < b.h) {
                        // estimate the value at the projected point location using bilinear_interpolate()
                        float f = bilinear_interpolate(b, proj_p.x, proj_p.y, loop3);
                        // assign it to the corresponding pixel offset by "dx" and "dy" in "combined_img"
                        set_pixel(combined_img, loop1 - dx, loop2 - dy, loop3, f);
                    }
                }
            }
        }
    }
    

    return combined_img;
}


image cylindrical_project(image im, float f)
/**************************************************
Project an image "im" onto a cylinder using focal length "f"
The formulas are given in the lecture slides
Attempt for extra credit. Return the projected image.
**************************************************/
{
    image c = copy_image(im);
    int x,y,z;
    int xc = im.w/2;
    int yc = im.h/2;
    // Loop over all pixels in image im
    for (z = 0; z < im.c; ++z) {
        for(y = 0; y < im.h; ++y){
            for(x = 0; x < im.w; ++x){
                // theta and h
                int theta = (x - xc)/f;
                int h = (y - yc)/f;
                // X', Y', Z'
                int X_p = sin(theta);
                int Y_p = h;
                int Z_p = cos(theta);
                // x', y'
                int x_p = f*(X_p/Z_p) + xc;
                int y_p = f*(Y_p/Z_p) + yc;
                // get current pixel value
                float pix = get_pixel(im, x, y, z);
                
                if (0 <= x_p && x_p < c.w) {
                    if (0 <= y_p && y_p < c.h) {
                        set_pixel(c, x_p, y_p, z, pix);
                    }
                }
                
            }
        }
    }
    return c;
}
