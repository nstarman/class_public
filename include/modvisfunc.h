/**
 * definitions for modifying the visibility function
 */

#ifndef __MODVISFUNC__
#define __MODVISFUNC__

#include "common.h"

// #define _SPLINE_NATURAL_ 0 /**< natural spline: ddy0=ddyn=0 */
// #define _SPLINE_EST_DERIV_ 1 /**< spline with estimation of first derivative on both edges */

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

// integrate g w/out using its derivatives
// return integral value
double integrate_g_simple(double * x_array,
				          double * y_array,
                          int n_lines,
				          int n_columns,
				          int index_y);

// integrate g w/out using its derivatives
// store integral array
int array_integrate_g_simple(double * x_array,
					         double * y_array,
                             int n_lines,
					         int n_columns,
					         int index_y,
					         int index_inty,
					         ErrorMsg errmsg);

// integrate g using its derivatives
// store integral array
int array_integrate_g(double * x_array,
			          double * y_array,
                      int n_lines,
			          int n_columns,
			          int index_y,
			          int index_inty,
			          ErrorMsg errmsg);

// sample from the original visibility function
// for inverse CDF method
double * visfunc_sample(double * x_array,
                      double * y_array,
                      int index_cdf,
                      int n_lines,
                      int n_columns,
                      int * tau_vis_size,
                      double x_low,
                      double x_up
                      );

// double * visfunc_subsample(struct thermo * pth,
//                            double * x_array,
//                            int n_lines,
//                            int x_low,
//                            int x_up);

// average of visibility function
// TODO find a boilerplate function that does this
double visfunc_avg(double * x_array,
                   int n_lines);

// standard deviation of visibility function
double visfunc_sigma(double * x_array,
                    double xavg,
                    int n_lines);

double visfunc_third_moment(double * x_array,
                            double xavg,
                            int n_lines);

double visfunc_skew_param(int n_lines,
                          double * x_array,
                          double avg,
                          double sigma,
                          double * array);

double visfunc_sigmaz(double muz);

double visfunc_mode(double skew_param);


/** The Parameterizations*/
// TODO, if possible to use the thermo struct


#ifdef __cplusplus
}
#endif

#endif
