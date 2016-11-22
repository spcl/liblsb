/* Source: GNU Scientific Library: https://www.gnu.org/software/gsl/ */



#include "liblsb_internal.hpp"

#include <stdlib.h>

double lsb_median_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n){
  double median ;
  const size_t lhs = (n - 1) / 2 ;
  const size_t rhs = n / 2 ;
  
  if (n == 0)
    return 0.0 ;

  if (lhs == rhs)
    {
      median = sorted_data[lhs * stride] ;
    }
  else 
    {
      median = (sorted_data[lhs * stride] + sorted_data[rhs * stride])/2.0 ;
    }

  return median ;
}


double lsb_quantile_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n, const double f){
  const double index = f * (n - 1) ;
  const size_t lhs = (int)index ;
  const double delta = index - lhs ;
  double result;

  if (n == 0)
    return 0.0 ;

  if (lhs == n - 1)
    {
      result = sorted_data[lhs * stride] ;
    }
  else 
    {
      result = (1 - delta) * sorted_data[lhs * stride] + delta * sorted_data[(lhs + 1) * stride] ;
    }

  return result ;
}


int lsb_fit_linear (const double *x, const size_t xstride, const double *y, const size_t ystride, const size_t n, double *c0, double *c1, double *cov_00, double *cov_01, double *cov_11, double *sumsq){
  double m_x = 0, m_y = 0, m_dx2 = 0, m_dxdy = 0;

  size_t i;

  for (i = 0; i < n; i++)
    {
      m_x += (x[i * xstride] - m_x) / (i + 1.0);
      m_y += (y[i * ystride] - m_y) / (i + 1.0);
    }

  for (i = 0; i < n; i++)
    {
      const double dx = x[i * xstride] - m_x;
      const double dy = y[i * ystride] - m_y;

      m_dx2 += (dx * dx - m_dx2) / (i + 1.0);
      m_dxdy += (dx * dy - m_dxdy) / (i + 1.0);
    }

  /* In terms of y = a + b x */

  {
    double s2 = 0, d2 = 0;
    double b = m_dxdy / m_dx2;
    double a = m_y - m_x * b;

    *c0 = a;
    *c1 = b;

    /* Compute chi^2 = \sum (y_i - (a + b * x_i))^2 */

    for (i = 0; i < n; i++)
      {
        const double dx = x[i * xstride] - m_x;
        const double dy = y[i * ystride] - m_y;
        const double d = dy - b * dx;
        d2 += d * d;
      }

    s2 = d2 / (n - 2.0);        /* chisq per degree of freedom */

    *cov_00 = s2 * (1.0 / n) * (1 + m_x * m_x / m_dx2);
    *cov_11 = s2 * 1.0 / (n * m_dx2);

    *cov_01 = s2 * (-m_x) / (n * m_dx2);

    *sumsq = d2;
  }

  return 0;
}




double lsb_mean (const double data[], const size_t stride, const size_t size){
  /* Compute the arithmetic mean of a dataset using the recurrence relation 
     mean_(n) = mean(n-1) + (data[n] - mean(n-1))/(n+1)   */

  long double mean = 0;
  size_t i;

  for (i = 0; i < size; i++)
    {
      mean += (data[i * stride] - mean) / (i + 1);
    }

  return mean;
}

double min(double a, double b) {
    if (a < b) return a;
    else return b;
}

double max(double a, double b) {
    if (a > b) return a;
    else return b;
}

