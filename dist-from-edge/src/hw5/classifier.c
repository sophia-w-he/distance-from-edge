// Sophia He
// CSE 455 Spring 2020
// HW 5
// May 29, 2020
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"


void activate_matrix(matrix m, ACTIVATION a)
{    
    /***********************************************************************
    Run an activation function on each element in a matrix,
    modifies the matrix in place
    matrix m: Input to activation function
    ACTIVATION a: function to run

    For each element "x" in the matrix:
    - If the activation is logistic, then the result is 1 / (1 + pow(M_E, -x))   
    - If the activation is relu, then the result is max(0, x).
    - If the activation is LRELU (i.e. Leaky Relu), then the result is max(0.1 * x, x)
    - if the activation is Softmax, then the result is pow(M_E, x) / \Sum_{i=0}^n(pow(M_E, x_i)), 
    where x_i is the i-the element in the row.
    ************************************************************************/
    int i, j;
    for(i = 0; i < m.rows; ++i){
        double sum = 0;
        for(j = 0; j < m.cols; ++j){
            double x = m.data[i][j];
            if(a == LOGISTIC){
                // result is 1 / (1 + pow(M_E, -x))
                m.data[i][j] = 1 / (1 + pow(M_E, -x));
            } else if (a == RELU){
                // result is max(0, x)
                if (0 > x) m.data[i][j] = 0 ;
            } else if (a == LRELU){
                // result is max (0.1*x, x)
                if (0.1 * x > x) m.data[i][j] = 0.1 * x ;
            } else if (a == SOFTMAX){
                // get pow (M_E, x) to divide by sum later
                m.data[i][j] = pow(M_E, x);
            }
            sum += m.data[i][j];
        }
        if (a == SOFTMAX) {
            // result is pow(M_E, x) normalized by sum
            for (j = 0; j < m.cols; ++j) {
                m.data[i][j] = m.data[i][j]/sum;
            }
        }
    }
}


void gradient_matrix(matrix m, ACTIVATION a, matrix d)
{
    /***********************************************************************
    Calculates the gradient of an activation function and multiplies it into
    the delta for a layer
    matrix m: an activated layer output
    ACTIVATION a: activation function for a layer
    matrix d: delta before activation gradient

    Modify matrix d in place. The "delta" ("gradient") is stored in the matrix d.

    Note that the shape for m and d are the same. Each element in m should has
    the corresponding gradient stored in d.

    For each element x in the matrix m:
    - If the activation is softmax or linear, then the gradient in d should not change.
    - If the activation is LOGISTIC, d.data[i][j] should times x * (1.0 - x);
    - If the activation is RELU:
        if x <= 0, then the gradient in d is zero. 
        otherwise, the gradient should stay the same. 
    - If the activation is LRELU:
        if x <= 0, then the gradient in d should divide 10.
        otherwise, the gradient should stay the same. 
    ************************************************************************/
    assert(m.rows == d.rows);
    assert(m.cols == d.cols);
 
    int i, j;
    for(i = 0; i < m.rows; ++i){
        for(j = 0; j < m.cols; ++j){
            double x = m.data[i][j];
            // multiply the correct element of d by the gradient
            //if (a == SOFTMAX || a == LINEAR)
                // data does not change
            if (a == LOGISTIC){
                d.data[i][j] *= x * (1.0 - x);
            } else if (a == RELU && x <= 0) {
                //if x <= 0, then the gradient in d is zero.
                d.data[i][j] = 0;
            } else if (a == LRELU && x <= 0) {
                // if x <= 0, then the gradient in d should divide 10.
                d.data[i][j] /= 10;
            }
        }
    }
}


matrix forward_layer(layer *l, matrix in)
{
    /***********************************************************************
    Forward propagate information through a layer
    layer *l: pointer to the layer, where l->w is the weight,
              and l->activation is the activation.
    matrix in: input to layer
    returns: matrix that is output of the layer

    calculate the output matrix, which is "activation(in \times l->w)". 
    You can simply multiply the matrix in with l->w, and then apply the 
    activate_matrix function.    
    ************************************************************************/

    l->in = in;  // Save the input for backpropagation


    // multiply input by weights and apply activation function.
    matrix out = make_matrix(in.rows, l->w.cols);
    
    out = matrix_mult_matrix(in, l->w);
    
    // activate matrix
    activate_matrix(out, l->activation);

    free_matrix(l->out);// free the old output
    l->out = out;       // Save the current output for gradient calculation
    return out;
}

matrix backward_layer(layer *l, matrix delta)
{
    /***********************************************************************
    Backward propagate derivatives through a layer
    layer *l: pointer to the layer
    matrix delta: partial derivative of loss w.r.t. output of layer
    returns: matrix, partial derivative of loss w.r.t. input to layer

    1.4.1 Apply the gradient_matrix function to add the activation into consideration
    for matrix delta.

    1.4.2 Calculate the l->dw, the delta (gradient) of weight.
    l->dw = x^T \times delta, where x is the input data (i.e. l->in), 
    and delta is calculated from 1.4.1\
    
    1.4.3 Calculate dx and return it.
    dx = delta \times w^T, where w is l->w. 
    ************************************************************************/

    // 1.4.1
    // delta is dL/dy
    // modify it in place to be dL/d(xw)
    gradient_matrix(l->out, l->activation, delta);


    // 1.4.2
    // then calculate dL/dw and save it in l->dw
    free_matrix(l->dw);
    matrix xt = transpose_matrix(l->in); // transpose
    // multiply by transposed matrix
    matrix dw = matrix_mult_matrix(xt, delta);
    free_matrix(xt);
    l->dw = dw;

    
    // 1.4.3
    // finally, calculate dL/dx and return it.
    matrix w_t = transpose_matrix(l->w); // transpose
    matrix dx = make_matrix(l->in.rows, l->in.cols);
    // multiply by transposed matrix
    dx = matrix_mult_matrix(delta, w_t);
    free_matrix(w_t);
    return dx;
}


void update_layer(layer *l, double rate, double momentum, double decay)
{
    /***********************************************************************
    Update the weights at layer l
    layer *l: pointer to the layer
    double rate: learning rate
    double momentum: amount of momentum to use
    double decay: value for weight decay

    1. Calculate Calculate Δw_t and update l->v
        Δw_t = dL/dw_t - λw_t + mΔw_{t-1}, where 
    m is the momentum, λ is the decay, 
    Δw_{t-1} is the gradient update in previous batch (i.e. l->v),
    w_t is the current weight (i.e. l->w), 
    dL/dw_t is the current gradient (i.e. l->dw).

    Hint: use axpy_matrix for convenience. 
    
    2. Update l->w, the weight for this layer
    New weight = learning rate * Δw_t + previous weight
    ************************************************************************/
    // begin calculation Δw_t = dL/dw_t - λw_t + mΔw_{t-1}
    matrix lambda_wt = copy_matrix(l->w);
    int i, j;
    // scale copy of l->w with decay: λw_t
    for(i = 0; i < lambda_wt.rows; ++i){
        for(j =0 ; j < lambda_wt.cols; ++j){
            lambda_wt.data[i][j] *= decay;
        }
    }
    matrix p = make_matrix(l->dw.rows, l->dw.cols);
    // subtract from l->dw: dL/dw_t - λw_t
    for(i = 0; i < p.rows; ++i){
        for(j = 0; j < p.cols; ++j){
            p.data[i][j] = l->dw.data[i][j] - lambda_wt.data[i][j];
        }
    }
    // scale and sum to achieve:
    // Δw_t = dL/dw_t - λw_t + mΔw_{t-1}
    // save it to l->v
    l->v = axpy_matrix(momentum, l->v, p);
    free_matrix(lambda_wt);
    free_matrix(p);


    // Update l->w
    // New weight = learning rate * Δw_t + previous weight
    matrix new_w = axpy_matrix(rate, l->v, l->w);
    l->w = new_w;

}

// Make a new layer for our model
// int input: number of inputs to the layer
// int output: number of outputs from the layer
// ACTIVATION activation: the activation function to use
layer make_layer(int input, int output, ACTIVATION activation)
{
    layer l;
    l.in  = make_matrix(1,1);
    l.out = make_matrix(1,1);
    l.w   = random_matrix(input, output, sqrt(2./input));
    l.v   = make_matrix(input, output);
    l.dw  = make_matrix(input, output);
    l.activation = activation;
    return l;
}

// Run a model on input X
// model m: model to run
// matrix X: input to model
// returns: result matrix
matrix forward_model(model m, matrix X)
{
    int i;
    for(i = 0; i < m.n; ++i){
        X = forward_layer(m.layers + i, X);
    }
    return X;
}

// Run a model backward given gradient dL
// model m: model to run
// matrix dL: partial derivative of loss w.r.t. model output dL/dy
void backward_model(model m, matrix dL)
{
    matrix d = copy_matrix(dL);
    int i;
    for(i = m.n-1; i >= 0; --i){
        matrix prev = backward_layer(m.layers + i, d);
        free_matrix(d);
        d = prev;
    }
    free_matrix(d);
}

// Update the model weights
// model m: model to update
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_model(model m, double rate, double momentum, double decay)
{
    int i;
    for(i = 0; i < m.n; ++i){
        update_layer(m.layers + i, rate, momentum, decay);
    }
}

// Find the index of the maximum element in an array
// double *a: array
// int n: size of a, |a|
// returns: index of maximum element
int max_index(double *a, int n)
{
    if(n <= 0) return -1;
    int i;
    int max_i = 0;
    double max = a[0];
    for (i = 1; i < n; ++i) {
        if (a[i] > max){
            max = a[i];
            max_i = i;
        }
    }
    return max_i;
}

// Calculate the accuracy of a model on some data d
// model m: model to run
// data d: data to run on
// returns: accuracy, number correct / total
double accuracy_model(model m, data d)
{
    matrix p = forward_model(m, d.X);
    int i;
    int correct = 0;
    for(i = 0; i < d.y.rows; ++i){
        if(max_index(d.y.data[i], d.y.cols) == max_index(p.data[i], p.cols)) ++correct;
    }
    return (double)correct / d.y.rows;
}

// Calculate the cross-entropy loss for a set of predictions
// matrix y: the correct values
// matrix p: the predictions
// returns: average cross-entropy loss over data points, 1/n Σ(-ylog(p))
double cross_entropy_loss(matrix y, matrix p)
{
    int i, j;
    double sum = 0;
    for(i = 0; i < y.rows; ++i){
        for(j = 0; j < y.cols; ++j){
            sum += -y.data[i][j]*log(p.data[i][j]);
        }
    }
    return sum/y.rows;
}


// Train a model on a dataset using SGD
// model m: model to train
// data d: dataset to train on
// int batch: batch size for SGD
// int iters: number of iterations of SGD to run (i.e. how many batches)
// double rate: learning rate
// double momentum: momentum
// double decay: weight decay
void train_model(model m, data d, int batch, int iters, double rate, double momentum, double decay)
{
    int e;
    for(e = 0; e < iters; ++e){
        data b = random_batch(d, batch);
        matrix p = forward_model(m, b.X);
        fprintf(stderr, "%06d: Loss: %f\n", e, cross_entropy_loss(b.y, p));
        matrix dL = axpy_matrix(-1, p, b.y); // partial derivative of loss dL/dy
        backward_model(m, dL);
        update_model(m, rate/batch, momentum, decay);
        free_matrix(dL);
        free_data(b);
    }
}


// Questions 
// (IMPORTANT. PLEASE ANSWER ALL QUESTIONS in the README.md FILE and SAVE TO a PDF FILE FOR GRADING)
// 
