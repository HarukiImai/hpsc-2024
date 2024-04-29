#include <cstdio>
#include <cstdlib>
#include <cmath>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N], is[N], buf[i];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
    is[i] = i;
    buf[i] = 0;
  }
  __m256 xvec = _mm256_load_ps(x);
  __m256 yvec = _mm256_load_ps(y);
  __m256 mvec = _mm256_load_ps(m);
  __m256 isvec = _mm256_load_ps(is);
  __m256 zeros = _mm256_setzero_ps();

  for (int i = 0; i < N; i++) {
   __m256 xi = _mm256_set1_ps(x[i]);
   __m256 yi = _mm256_set1_ps(y[i]);
   __m256 ivec = _mm256_set1_ps(i);
   __m256 mask = _mm256_cmp_ps(ivec, isvec, _CMP_NEQ_OQ);

   __m256 rx = _mm256_sub_ps(xi, xvec);
   __m256 ry = _mm256_sub_ps(yi, yvec);
   __m256 rsq = _mm256_add_ps(_mm_mul_ps(rx,rx), _mm256_mul_ps(ry,ry));
   __m256 rr = _mm256_blendv_ps(zeros, _mm256_rsqrt_ps(rsq), mask);
   __m256 fxi = _mm256_mul_ps(_mm256_mul_ps(rx,mvec), _mm256_mul_ps(_mm256_mul_ps(rr, rr), rr);
   __m256 fyi = _mm256_mul_ps(_mm256_mul_ps(ry,mvec), _mm256_mul_ps(_mm256_mul_ps(rr, rr), rr);

   fxi = _mm256_hadd_ps( fxi , _mm256_permute2f128_ps(fxi, fxi, 1));
   fxi = _mm256_hadd_ps(fxi, fxi);
   fxi = _mm256_hadd_ps(fxi, fxi);
   _mm256_store_ps(buf,fxi);
   fx[i] = -1 * buf[0];

   fyi = _mm256_hadd_ps( fyi , _mm256_permute2f128_ps(fyi, fyi, 1));
   fyi = _mm256_hadd_ps(fyi, fyi);
   fyi = _mm256_hadd_ps(fyi, fyi);
   _mm256_store_ps(buf,fyi);
   fy[i] = -1 * buf[0];

   printf("%d %g %g\n", i , fx[i], fy[i]);
   }
//  for(int i=0; i<N; i++) {
//    for(int j=0; j<N; j++) {
//      ifa(i != j) {
//        float rx = x[i] - x[j];
//        float ry = y[i] - y[j];
//        float r = std::sqrt(rx * rx + ry * ry);
//        fx[i] -= rx * m[j] / (r * r * r);
//        fy[i] -= ry * m[j] / (r * r * r);
//      }
//    }
//    printf("%d %g %g\n",i,fx[i],fy[i]);
//  }
}
