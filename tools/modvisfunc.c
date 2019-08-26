/**
 * module with tools for manipulating the visibility function
 * Nathaniel Starkman, Aug 25, 2019
 */

#include "modvisfunc.h"


double integrate_g_simple(double * x_array,
                          double * y_array,
                          int n_lines,
                          int n_columns,
                          int index_y) {
  /** integrate the visibility function w/out derivatives
  returns just the final integral value
  
  INPUT:
  tau_table,
  pth->thermodynamics_table,
  pth->tt_size,
  pth->th_size,
  pth->index_th_g,
  pth->error_message

  local variables
  double intg // the resulting integral value
  int i;      // index
  double h;   // step size
  */

  /** - define local variables */
  double intg;
  int i;
  double h;

  /** - integrating */
  intg  = 0.;

  // for (i=0; i<=n_lines-1; i++) {
  for (i=n_lines-1; i>0; i--) {

    h = (x_array[i-1] - x_array[i]);

    intg += (y_array[(i-1)*n_columns+index_y]+y_array[i*n_columns+index_y])*h/2.;

  }

  return intg;
}


int array_integrate_g_simple(double * x_array,
                             double * y_array,
                             int n_lines,
                             int n_columns,
                             int index_y,
                             int index_inty,
                             ErrorMsg errmsg) {
  // integrate the visibility function (g)
  // without using derivatives of g

  int i;
  double h;

  *(y_array+0*n_columns+index_inty)  = 0.;

  for (i=0; i<n_lines-1; i++) {

    h = (x_array[i] - x_array[i+1]);

    *(y_array+(i+1)*n_columns+index_inty) = *(y_array+i*n_columns+index_inty) +
      (y_array[i*n_columns+index_y]+y_array[(i+1)*n_columns+index_y])*h/2.;

  }

  return _SUCCESS_;
}


int array_integrate_g(double * x_array,
                      double * y_array,
                      int n_lines,
                      int n_columns,
                      int index_y,
                      int index_inty,
                      ErrorMsg errmsg) {

  int i;
  double h;

  *(y_array+(n_lines-1)*n_columns+index_inty)  = 0.;

  for (i=n_lines-1; i>0; i--) {

    h = (x_array[i-1] - x_array[i]);

    *(y_array+(i-1)*n_columns+index_inty) = *(y_array+i*n_columns+index_inty) +
      (y_array[(i-1)*n_columns+index_y]+y_array[i*n_columns+index_y])*h/2.;

  }

  return _SUCCESS_;
}


double * visfunc_sample(double * x_array,
                        double * y_array,
                        int index_cdf,
                        int n_lines,
                        int n_columns,
                        int* tau_vis_size,
                        double x_low,
                        double x_up
                        ) {
  // sample a PDF by inverse CDF sampling
  // between two values

  /** - define local variables */
  int i;        // index through the cdf, and pdf
  int j;        // index through the uniform sample
  double *xx;   // x sample
  double *x;    // x subsample
  double u;     // uniform sample
  double N;     // double form of n_lines
  // double u0; // uniform sample

  /** SAMPLE X */
  xx=malloc((n_lines) * sizeof(double));
  i=0;
  u=0.;
  N=(double)n_lines;
  // u0 = 0.; // can start at any tiny number

  for (j=0; j<n_lines; j++) {
      u = j/N;  // 0<uj<1
      while (u > *(y_array+i*n_columns+index_cdf)) {
          i++;
      }
      *(xx+j)=*(x_array+i);
  }

  /** SUBSAMPLE X */
  // subsampling because it always includes a few samples from tau=14000
  // which throws everything off. Restricting the domain to tau<500
  /** getting number of elements */
  i=0; // reusing as main index
  j=0; // reusing as size index
  for (i=0; i<n_lines; i++) {
      if ((xx[i] > x_low) && (xx[i] < x_up)) {
          j++;
      }
  }

  x=malloc((j) * sizeof(double));
  *tau_vis_size=j;  // assigning size

  /** now that have size, getting subsample */
  i=0; // reusing as main index
  j=0; // reusing as size index
  for (i=0; i<n_lines; i++) {
      if ((xx[i] > x_low) && (xx[i] < x_up)) {
          x[j]=xx[i];
          j++;
      }
  }

  return x;
}


double visfunc_avg(double * x_array,
                   int n_lines
                   ) {
  // 

  /** - define local variables */
  int i;       // index through the cdf, and pdf
  double avg;  // avg of x sample

  // find x_avg
  avg=0.0;
  for (i=0; i<n_lines; i++) {
      avg+=x_array[i];
  }
  avg/=(double)n_lines;
  // printf("visfunc_avg: %f\n", avg);
  return avg;
}


double visfunc_sigma(double * x_array,
                     double avg,
                     int n_lines) {
  /** - define local variables */
  int i;        // index through the cdf, and pdf
  double E2;

  E2=0;
  for (i=0; i<n_lines; i++) {
      E2+=pow(x_array[i]-avg, 2.);
  }
  E2=E2/n_lines;

  return sqrt(E2);
}


double visfunc_third_moment(double * x_array,
                            double xavg,
                            int n_lines) {
  // central third moment
  /** - define local variables */
  int i;        // index through the cdf, and pdf
  double E3;

  E3=0;
  for (i=0; i<n_lines; i++) {
      E3+=pow(x_array[i]-xavg, 3.);
  }
  E3/=n_lines;
  return E3;
}


double visfunc_skew_param(int n_lines,
                          double * x_array,
                          double avg,
                          double sigma,
                          double * array) {
  /** - define local variables */
  int i;        // index through the cdf, and pdf
  // double *x;    // x sample
  double E3;    // 3rd moment
  double skew;   // sample skewness
  double delta; //
  double param;  // skewness parameterization
  // have cdf
  E3=visfunc_third_moment(x_array, avg, n_lines);
  // printf("E3: %f\n", E3);

  // TODO compare definition against Ashour and Abdel-Hamid (2010)
  // skew=sqrt(n_lines*(n_lines-1))*E3/pow(sigma,3.)/(n_lines-2);
  skew=E3/pow(sigma,3.);
  // printf("skew: %f\n", skew);

  // // delta from wikipedia
  // delta=sqrt((_PI_/2.)*(pow(skew,2./3.) / (pow(skew,2./3.) + pow((4.-_PI_)/2., 2./3.))));
  // // printf("delta: %f\n", delta);
  //
  // // skewness paramter
  // param=delta/sqrt(1-delta*delta);
  // // printf("param: %f\n", param);

  return skew;
}

double visfunc_sigmaz(double muz){
    return sqrt(1 - muz*muz);
}

double visfunc_mode(double skew_param){
    // wikipedia numerical approximation

    /** - define local variables */
    double delta;
    double muz;
    double sigz;
    double skewness;
    double mo;

    delta=skew_param/sqrt(1 + skew_param*skew_param);
    muz=sqrt(2./_PI_)*delta;
    sigz=sqrt(1 - muz*muz);

    skewness=((4-_PI_)/2) * pow(delta*sqrt(2/_PI_), 3.) / pow(1-2*delta*delta/_PI_,3./2.);

    // mo(s) \approx \mu_z - gam_1*sigz/2 - sgn(s)/2 Exp(-2pi/|s|])
    mo=muz - (skewness*sigz/2) - exp(-2*_PI_/sqrt(skew_param*skew_param))/2;

    return mo;
}
