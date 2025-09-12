#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include "../include/macro.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/flavour_matrix.h"
#include "../include/function_pointers.h"
#include "../include/gauge_conf.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/su2_monopoles.h"
#include "../include/sun_monopoles.h"
#include "../include/tens_prod.h"
#include "../include/u1_monopoles.h"

// computation of the plaquette (1/NCOLOR the trace of) in position r and
// positive directions i,j
double plaquettep(Gauge_Conf const *const GC, Geometry const *const geo, long r,
                  int i, int j) {
  GAUGE_GROUP matrix;

#ifdef DEBUG
  if (r >= geo->d_volume) {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume,
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if (j >= STDIM || i >= STDIM) {
    fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j,
            STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
#endif

  //
  //       ^ i
  //       |   (2)
  //       +---<---+
  //       |       |
  //   (3) V       ^ (1)
  //       |       |
  //       +--->---+---> j
  //       r   (4)
  //

  equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
  times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
  times_equal_dag(&matrix, &(GC->lattice[r][i]));
  times_equal(&matrix, &(GC->lattice[r][j]));

  return retr(&matrix);
}

// computation of the plaquette (1/NCOLOR the trace of) in position r and
// positive directions i,j
double complex plaquettep_complex(Gauge_Conf const *const GC,
                                  Geometry const *const geo, long r, int i,
                                  int j) {
  GAUGE_GROUP matrix;

#ifdef DEBUG
  if (r >= geo->d_volume) {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume,
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if (j >= STDIM || i >= STDIM) {
    fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j,
            STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
#endif

  //
  //       ^ i
  //       |   (2)
  //       +---<---+
  //       |       |
  //   (3) V       ^ (1)
  //       |       |
  //       +--->---+---> j
  //       r   (4)
  //

  equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
  times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
  times_equal_dag(&matrix, &(GC->lattice[r][i]));
  times_equal(&matrix, &(GC->lattice[r][j]));

  return retr(&matrix) + I * imtr(&matrix);
}

// computation of the plaquette (matrix) in position r and positive directions
// i,j
void plaquettep_matrix(Gauge_Conf const *const GC, Geometry const *const geo,
                       long r, int i, int j, GAUGE_GROUP *matrix) {
#ifdef DEBUG
  if (r >= geo->d_volume) {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume,
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if (j >= STDIM || i >= STDIM) {
    fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j,
            STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
#endif

  //
  //       ^ j
  //       |   (3)
  //       +---<---+
  //       |       |
  //   (4) V       ^ (2)
  //       |       |
  //       +--->---+---> i
  //       r   (1)
  //

  equal(matrix, &(GC->lattice[r][i]));
  times_equal(matrix, &(GC->lattice[nnp(geo, r, i)][j]));
  times_equal_dag(matrix, &(GC->lattice[nnp(geo, r, j)][i]));
  times_equal_dag(matrix, &(GC->lattice[r][j]));
}

// compute the four-leaf clover in position r, in the plane i,j and save it in M
void clover(Gauge_Conf const *const GC, Geometry const *const geo, long r,
            int i, int j, GAUGE_GROUP *M) {
  GAUGE_GROUP aux;
  long k, p;

#ifdef DEBUG
  if (r >= geo->d_volume) {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume,
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if (i >= STDIM || j >= STDIM) {
    fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j,
            STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
#endif

  zero(M);

  //
  //                   i ^
  //                     |
  //              (7)    |     (2)
  //         +-----<-----++----->-----+
  //         |           ||           |
  //         |           ||           |
  //    (6)  ^       (8) V^ (1)       V (3)
  //         |           ||           |
  //         |   (5)     || r   (4)   |
  //       k +-----<-----++-----<-----+------>   j
  //         +----->-----++----->-----+
  //         |    (12)   ||   (13)    |
  //         |           ||           |
  //    (11) ^       (9) V^ (16)       V (14)
  //         |           ||           |
  //         |           ||           |
  //         +------<----++-----<-----+
  //              (10)   p      (15)
  //
  // avanti-avanti
  equal(&aux, &(GC->lattice[r][i]));                        // 1
  times_equal(&aux, &(GC->lattice[nnp(geo, r, i)][j]));     // 2
  times_equal_dag(&aux, &(GC->lattice[nnp(geo, r, j)][i])); // 3
  times_equal_dag(&aux, &(GC->lattice[r][j]));              // 4
  plus_equal(M, &aux);

  k = nnm(geo, r, j);

  // avanti-indietro
  equal_dag(&aux, &(GC->lattice[k][j]));                // 5
  times_equal(&aux, &(GC->lattice[k][i]));              // 6
  times_equal(&aux, &(GC->lattice[nnp(geo, k, i)][j])); // 7
  times_equal_dag(&aux, &(GC->lattice[r][i]));          // 8
  plus_equal(M, &aux);

  p = nnm(geo, r, i);

  // indietro-indietro
  equal_dag(&aux, &(GC->lattice[p][i]));                    // 9
  times_equal_dag(&aux, &(GC->lattice[nnm(geo, k, i)][j])); // 10
  times_equal(&aux, &(GC->lattice[nnm(geo, k, i)][i]));     // 11
  times_equal(&aux, &(GC->lattice[k][j]));                  // 12
  plus_equal(M, &aux);

  // indietro-avanti
  equal(&aux, &(GC->lattice[r][j]));                        // 13
  times_equal_dag(&aux, &(GC->lattice[nnp(geo, p, j)][i])); // 14
  times_equal_dag(&aux, &(GC->lattice[p][j]));              // 15
  times_equal(&aux, &(GC->lattice[p][i]));                  // 16
  plus_equal(M, &aux);
}

// compute the mean plaquettes (spatial, temporal)
void plaquette(Gauge_Conf const *const GC, Geometry const *const geo,
               double *plaqs, double *plaqt) {
  long r;
  double ps, pt;

  ps = 0.0;
  pt = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt)    \
    reduction(+ : ps)
#endif
  for (r = 0; r < (geo->d_volume); r++) {
    int i, j;
    i = 0;
    for (j = 1; j < STDIM; j++) {
      pt += plaquettep(GC, geo, r, i, j);
    }

    for (i = 1; i < STDIM; i++) {
      for (j = i + 1; j < STDIM; j++) {
        ps += plaquettep(GC, geo, r, i, j);
      }
    }
  }

  if (STDIM > 2) {
    ps *= geo->d_inv_vol;
    ps /= ((double)(STDIM - 1) * (STDIM - 2) / 2);
  } else {
    ps = 0.0;
  }

  pt *= geo->d_inv_vol;
  pt /= ((double)STDIM - 1);

  *plaqs = ps;
  *plaqt = pt;
}

// compute the mean plaquette (spatial, temporal) with obc
void plaquette_obc(Gauge_Conf const *const GC, Geometry const *const geo,
                   double *plaqs, double *plaqt) {
  long r;
  double ps, pt;

  ps = 0.0;
  pt = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt)    \
    reduction(+ : ps)
#endif
  for (r = 0; r < (geo->d_volume); r++) {
    int i, j;
    i = 0;
    for (j = 1; j < STDIM; j++) {
      pt += plaquettep(GC, geo, r, i, j);
    }

    for (i = 1; i < STDIM; i++) {
      for (j = i + 1; j < STDIM; j++) {
        ps += plaquettep(GC, geo, r, i, j);
      }
    }
  }

  if (STDIM > 2) {
    ps /= ((double)geo->d_size[0]);
    ps /= ((double)geo->d_nfaces);
  } else {
    ps = 0.0;
  }

  pt /= ((double)geo->d_nfaces_temp);

  *plaqs = ps;
  *plaqt = pt;
}

// computes the number of temporal faces on a sublattice of dimension (Ns -
// dis)^STDIM * Nt
int compute_nfaces_t_sublat(Geometry const *const geo, int dis) {
  int j, k;
  int aux;
  int nfaces_t = 0;

  for (j = 1; j < STDIM; j++) {
    aux = (geo->d_size[0]) * (geo->d_size[j] - dis - 1);
    for (k = 1; k < STDIM; k++) {
      if (k != j)
        aux *= geo->d_size[k] - dis;
    }
    nfaces_t += aux;
  }

  return nfaces_t;
}

// computes the number of spatial faces on a sublattice of dimension (Ns - dis)
int compute_nfaces_sp_sublat(Geometry const *const geo, int dis) {
  int i, j, k;
  int aux;
  int nfaces_sp = 0;

  for (i = 1; i < STDIM; i++) {
    for (j = i + 1; j < STDIM; j++) {
      aux = (geo->d_size[i] - dis - 1) * (geo->d_size[j] - dis - 1);

      for (k = 1; k < STDIM; k++) {
        if (k != i && k != j) {
          aux *= (geo->d_size[k] - dis);
        }
      }
      nfaces_sp += aux;
    }
  }

  return nfaces_sp;
}

void plaquette_obc_fve(Gauge_Conf const *const GC, Geometry const *const geo,
                       double *plaqs, double *plaqt, int dis) {
  long r;
  int nfaces_t, nfaces_sp;
  int measure_sp, measure_t, measure_t_dir[STDIM - 1];
  int cartcoord[STDIM];

  double ps, pt;

  ps = 0.0;
  pt = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt)    \
    reduction(+ : ps)
#endif
  for (r = 0; r < (geo->d_volume); r++) {
    int k;

    si_to_cart(cartcoord, r, geo);

    measure_sp = 1;
    measure_t = 1;
    for (k = 0; k < STDIM - 1; k++) {
      measure_t_dir[k] = 1;
    }

    for (k = 1; k < STDIM; k++) {
      if (cartcoord[k] < dis || cartcoord[k] > (geo->d_size[k] - dis - 1)) {
        measure_t = 0;
        measure_sp = 0;
      } else if (cartcoord[k] == (geo->d_size[k] - dis - 1)) {
        measure_t_dir[k - 1] = 0;
        measure_sp = 0;
      }
    }

    // // debug
    // for(k=0; k<STDIM; k++){
    //   fprintf(stdout, "%d ", cartcoord[k]);
    // }
    // fprintf(stdout, "%d ", measure_t);
    // for(k=0; k<STDIM-1; k++){
    //   fprintf(stdout, "%d ", measure_t_dir[k]);
    // }
    // fprintf(stdout, "\n");

    int i, j;
    i = 0;

    for (j = 1; j < STDIM; j++) {
      pt += plaquettep(GC, geo, r, i, j) * (double)measure_t *
            (double)measure_t_dir[j - 1];
    }

    for (i = 1; i < STDIM; i++) {
      for (j = i + 1; j < STDIM; j++) {
        ps += plaquettep(GC, geo, r, i, j) * (double)measure_sp;
      }
    }
  }

  // fprintf(stdout, "\n");

  nfaces_t = compute_nfaces_t_sublat(geo, 2 * dis);
  nfaces_sp = compute_nfaces_sp_sublat(geo, 2 * dis);

  if (STDIM > 2) {
    ps /= ((double)geo->d_size[0]);
    ps /= ((double)nfaces_sp);
  } else {
    ps = 0.0;
  }

  pt /= ((double)nfaces_t);

  *plaqs = ps;
  *plaqt = pt;
}

// compute the clover discretization of
// sum_{\mu\nu}  Tr(F_{\mu\nu}F_{\mu\nu})/2
void clover_disc_energy(Gauge_Conf const *const GC, Geometry const *const geo,
                        double *energy) {
  long r;
  double ris;

  ris = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
#endif
  for (r = 0; r < geo->d_volume; r++) {
    int i, j;
    GAUGE_GROUP aux1, aux2;

    for (i = 0; i < STDIM; i++) {
      for (j = i + 1; j < STDIM; j++) {
        clover(GC, geo, r, i, j, &aux1);

        ta(&aux1);
        equal(&aux2, &aux1);
        times_equal(&aux1, &aux2);
        ris += -NCOLOR * retr(&aux1) / 16.0;
      }
    }
  }

  *energy = ris * geo->d_inv_vol;
}

// rectangular Wilson loop of size (wi,wj) in direction (i,j) at point r
double Wilsonp(Gauge_Conf const *const GC, Geometry const *const geo, int i,
               int j, int wi, int wj, long r) {
  int aux;
  GAUGE_GROUP matrix;

#ifdef DEBUG
  if (r >= geo->d_volume) {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume,
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if (j >= STDIM || i >= STDIM) {
    fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j,
            STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
#endif

  //
  //       ^ i
  //       |       r2
  //    r3 +---<---+
  //       |       |
  //       V       ^ wi
  //       |       |
  //       +--->---+---> j
  //       r   wj  r1
  //

  one(&matrix);
  // now we are in r
  for (aux = 0; aux < wj; aux++) {
    times_equal(&matrix, &(GC->lattice[r][j]));
    r = nnp(geo, r, j);
  }
  // now we are in r1
  for (aux = 0; aux < wi; aux++) {
    times_equal(&matrix, &(GC->lattice[r][i]));
    r = nnp(geo, r, i);
  }
  // now we are in r2
  for (aux = 0; aux < wj; aux++) {
    r = nnm(geo, r, j);
    times_equal_dag(&matrix, &(GC->lattice[r][j]));
  }
  // now we are in r3
  for (aux = 0; aux < wi; aux++) {
    r = nnm(geo, r, i);
    times_equal_dag(&matrix, &(GC->lattice[r][i]));
  }
  // now we are in r

  return retr(&matrix);
}

// averaged temporal rectangular Wilson loop of size (wt,ws)
double Wilsont(Gauge_Conf const *const GC, Geometry const *const geo, int wt,
               int ws) {
  long r;
  double ris;

  ris = 0;
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
#endif
  for (r = 0; r < geo->d_volume; r++) {
    int j;
    for (j = 1; j < STDIM; j++) {
      ris += Wilsonp(GC, geo, 0, j, wt, ws, r);
    }
  }

  ris *= geo->d_inv_vol;
  ris /= (STDIM - 1);

  return ris;
}

// multi-step staircase Wilson loop of size sqrt(2)*wjk at point r.
// The i direction has length wi.
double staircase_Wilsonp(Gauge_Conf const *const GC, Geometry const *const geo,
                         int i, int j, int k, int wi, int wjk, long r,
                         int extra) {
  int aux;
  GAUGE_GROUP matrix;

#ifdef DEBUG
  if (r >= geo->d_volume) {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume,
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if (j >= STDIM || i >= STDIM) {
    fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j,
            STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
#endif

  //
  //       ^ j
  //       |     |
  //       |    _|
  //       |  _|
  //       |_|
  //       |
  //       +----------> k
  //       r
  //

  one(&matrix);
  // diagonal steps in jk plane
  for (aux = 0; aux < wjk; aux++) {
    times_equal(&matrix, &(GC->lattice[r][j]));
    r = nnp(geo, r, j);
    times_equal(&matrix, &(GC->lattice[r][k]));
    r = nnp(geo, r, k);
  }
  // extra steps in dir j
  if (extra > 0) {
    for (aux = 0; aux < extra; aux++) {
      times_equal(&matrix, &(GC->lattice[r][j]));
      r = nnp(geo, r, j);
    }
  }
  // linear steps in direction i
  for (aux = 0; aux < wi; aux++) {
    times_equal(&matrix, &(GC->lattice[r][i]));
    r = nnp(geo, r, i);
  }
  // extra steps in dir j backwards
  if (extra > 0) {
    for (aux = 0; aux < extra; aux++) {
      r = nnm(geo, r, j);
      times_equal_dag(&matrix, &(GC->lattice[r][j]));
    }
  }
  // diagonal steps in jk plane backwards
  for (aux = 0; aux < wjk; aux++) {
    r = nnm(geo, r, k);
    times_equal_dag(&matrix, &(GC->lattice[r][k]));
    r = nnm(geo, r, j);
    times_equal_dag(&matrix, &(GC->lattice[r][j]));
  }
  // linear steps backwards in direction i
  for (aux = 0; aux < wi; aux++) {
    r = nnm(geo, r, i);
    times_equal_dag(&matrix, &(GC->lattice[r][i]));
  }

  return retr(&matrix);
}

// single-step staircase Wilson loop of size sqrt(wj**2+wk**2) at point r.
// The i direction has length wi.
double singlestep_Wilsonp(Gauge_Conf const *const GC, Geometry const *const geo,
                          int i, int j, int k, int wi, int wj, int wk, long r) {
  int aux;
  GAUGE_GROUP matrix;

#ifdef DEBUG
  if (r >= geo->d_volume) {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume,
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if (j >= STDIM || i >= STDIM) {
    fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j,
            STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
#endif

  //
  //       ^ j
  //       |
  //       |
  //       |  ___
  //       | |
  //       | |
  //       +----------> k
  //       r
  //

  one(&matrix);
  // linear steps in j direction
  for (aux = 0; aux < wj; aux++) {
    times_equal(&matrix, &(GC->lattice[r][j]));
    r = nnp(geo, r, j);
  }
  // linear steps in k direction
  for (aux = 0; aux < wk; aux++) {
    times_equal(&matrix, &(GC->lattice[r][k]));
    r = nnp(geo, r, k);
  }
  // linear steps in direction i
  for (aux = 0; aux < wi; aux++) {
    times_equal(&matrix, &(GC->lattice[r][i]));
    r = nnp(geo, r, i);
  }
  // linear steps backwards in direction k
  for (aux = 0; aux < wk; aux++) {
    r = nnm(geo, r, k);
    times_equal_dag(&matrix, &(GC->lattice[r][k]));
  }
  // linear steps backwards in direction j
  for (aux = 0; aux < wj; aux++) {
    r = nnm(geo, r, j);
    times_equal_dag(&matrix, &(GC->lattice[r][j]));
  }
  // linear steps backwards in direction i
  for (aux = 0; aux < wi; aux++) {
    r = nnm(geo, r, i);
    times_equal_dag(&matrix, &(GC->lattice[r][i]));
  }

  return retr(&matrix);
}

// averaged temporal multi-step staircase Wilson loop of size (wt,ws*sqrt(2)) in
// the (1,2)=(x,y) plane
double staircase_Wilsont_xy(Gauge_Conf const *const GC,
                            Geometry const *const geo, int wt, int ws,
                            int extra) {
  long r;
  double ris;

  ris = 0;
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
#endif
  for (r = 0; r < geo->d_volume; r++) {
    ris += staircase_Wilsonp(GC, geo, 0, 1, 2, wt, ws, r, extra);
  }

  ris *= geo->d_inv_vol;

  return ris;
}

// temporally averaged temporal rectangular Wilson loop of size (wt,ws),
// starting at spatial point rsp and measuring in dir y
double Wilsont_obc(Gauge_Conf const *const GC, Geometry const *const geo,
                   int wt, int ws, long rsp) {
  long r;
  int t;
  double ris;

  ris = 0;
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
#endif
  for (t = 0; t < geo->d_size[0]; t++) {
    r = sisp_and_t_to_si(geo, rsp, t);
    // ris += Wilsonp(GC, geo, 0, 1, wt, ws, r);
    ris += Wilsonp(GC, geo, 0, 2, wt, ws, r);
  }

  ris /= geo->d_size[0];

  return ris;
}

// averaged temporal rectangular Wilson loop of size (wt,ws)a
double Wilsont_obc_avg(Gauge_Conf const *const GC, Geometry const *const geo,
                       int wt, int ws) {
  int count;
  long r;
  double aux, ris;

  ris = 0;
  count = 0;
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
#endif
  for (r = 0; r < geo->d_volume; r++) {
    int j;
    for (j = 1; j < STDIM; j++) {
      aux = Wilsonp(GC, geo, 0, j, wt, ws, r);
      if (fabs(aux) < 1e-12) {
        ris += Wilsonp(GC, geo, 0, j, wt, ws, r);
        count++;
      }
    }
  }
  ris /= count;

  return ris;
}

// temporally averaged temporal multi-step staircase Wilson loop of size
// (wt,ws*sqrt(2)) in the (1,2)=(x,y) plane, at spatial point rsp
double staircase_Wilsont_xy_obc(Gauge_Conf const *const GC,
                                Geometry const *const geo, int wt, int ws,
                                long rsp, int extra) {
  long r;
  int t;
  double ris;

  ris = 0;
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(t) reduction(+ : ris)
#endif
  for (t = 0; t < geo->d_size[0]; t++) {
    r = sisp_and_t_to_si(geo, rsp, t);
    ris += staircase_Wilsonp(GC, geo, 0, 2, 1, wt, ws, r, extra);
    // ris += staircase_Wilsonp(GC, geo, 0, 1, 2, wt, ws, r);
  }

  ris /= geo->d_size[0];

  return ris;
}

// averaged temporal multi-step staircase Wilson loop of size
// (wt,ws*sqrt(2)) in the (1,2)=(x,y) plane
double staircase_Wilsont_xy_obc_avg(Gauge_Conf const *const GC,
                                    Geometry const *const geo, int wt, int ws,
                                    int extra) {
  int count;
  long r;
  double aux, ris;

  ris = 0;
  count = 0;
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
#endif
  for (r = 0; r < geo->d_volume; r++) {
    aux = staircase_Wilsonp(GC, geo, 0, 2, 1, wt, ws, r, extra);
    if (fabs(aux) < 1e-12) {
      ris += staircase_Wilsonp(GC, geo, 0, 2, 1, wt, ws, r, extra);
      count++;
    }
  }

  ris /= count;

  return ris;
}

// temporally averaged temporal multi-step staircase Wilson loop of size
// (wt,sqrt(wj**2+2k**2)) in the (1,2)=(x,y) plane, at spatial point rsp
double singlestep_Wilsont_xy_obc(Gauge_Conf const *const GC,
                                 Geometry const *const geo, int wt, int wx,
                                 int wy, long rsp) {
  long r;
  int t;
  double ris;

  ris = 0;
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
#endif
  for (t = 0; t < geo->d_size[0]; t++) {
    r = sisp_and_t_to_si(geo, rsp, t);
    ris += singlestep_Wilsonp(GC, geo, 0, 1, 2, wt, wx, wy, r);
    // ris += singlestep_Wilsonp(GC, geo, 0, 2, 1, wt, wx, wy, r);
  }

  ris /= geo->d_size[0];

  return ris;
}

// compute the mean Polyakov loop (the trace of)
void polyakov(Gauge_Conf const *const GC, Geometry const *const geo,
              double *repoly, double *impoly) {
  long rsp;
  double rep, imp;

  rep = 0.0;
  imp = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) \
    reduction(+ : imp)
#endif
  for (rsp = 0; rsp < geo->d_space_vol; rsp++) {
    long r;
    int i;
    GAUGE_GROUP matrix;

    r = sisp_and_t_to_si(geo, rsp, 0);

    one(&matrix);
    for (i = 0; i < geo->d_size[0]; i++) {
      times_equal(&matrix, &(GC->lattice[r][0]));
      r = nnp(geo, r, 0);
    }

    rep += retr(&matrix);
    imp += imtr(&matrix);
  }

  *repoly = rep * geo->d_inv_space_vol;
  *impoly = imp * geo->d_inv_space_vol;
}

// compute the mean Polyakov loop in the adjoint representation (the trace of)
void polyakov_adj(Gauge_Conf const *const GC, Geometry const *const geo,
                  double *repoly, double *impoly) {
  long rsp;
  double rep, imp;
  double complex tr;

  rep = 0.0;
  imp = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) \
    reduction(+ : imp)
#endif
  for (rsp = 0; rsp < geo->d_space_vol; rsp++) {
    long r;
    int i;
    GAUGE_GROUP matrix;

    r = sisp_and_t_to_si(geo, rsp, 0);

    one(&matrix);
    for (i = 0; i < geo->d_size[0]; i++) {
      times_equal(&matrix, &(GC->lattice[r][0]));
      r = nnp(geo, r, 0);
    }
    tr = NCOLOR * retr(&matrix) + NCOLOR * imtr(&matrix) * I;

#if NCOLOR == 1
    (void)tr;
    rep += 0.0;
#else
    rep += (cabs(tr) * cabs(tr) - 1) / (NCOLOR * NCOLOR - 1);
#endif

    imp += 0.0;
  }

  *repoly = rep * geo->d_inv_space_vol;
  *impoly = imp * geo->d_inv_space_vol;
}

// compute the mean Polyakov loop and its powers (trace of) in the presence of
// trace deformation
void polyakov_for_tracedef(Gauge_Conf const *const GC,
                           Geometry const *const geo, double *repoly,
                           double *impoly) {
  long rsp;
  double **rep, **imp;
  int j, err;
  long i;

  for (j = 0; j < (int)floor(NCOLOR / 2); j++) {
    repoly[j] = 0.0;
    impoly[j] = 0.0;
  }

  err = posix_memalign((void **)&rep, (size_t)DOUBLE_ALIGN,
                       (size_t)geo->d_space_vol * sizeof(double *));
  if (err != 0) {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__,
            __LINE__);
    exit(EXIT_FAILURE);
  }
  err = posix_memalign((void **)&imp, (size_t)DOUBLE_ALIGN,
                       (size_t)geo->d_space_vol * sizeof(double *));
  if (err != 0) {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__,
            __LINE__);
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < geo->d_space_vol; i++) {
    err = posix_memalign((void **)&(rep[i]), (size_t)DOUBLE_ALIGN,
                         (size_t)(int)floor(NCOLOR / 2) * sizeof(double));
    if (err != 0) {
      fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__,
              __LINE__);
      exit(EXIT_FAILURE);
    }
    err = posix_memalign((void **)&(imp[i]), (size_t)DOUBLE_ALIGN,
                         (size_t)(int)floor(NCOLOR / 2) * sizeof(double));
    if (err != 0) {
      fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__,
              __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  for (i = 0; i < geo->d_space_vol; i++) {
    for (j = 0; j < (int)floor(NCOLOR / 2); j++) {
      rep[i][j] = 0.0;
      imp[i][j] = 0.0;
    }
  }

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(rsp)
#endif
  for (rsp = 0; rsp < geo->d_space_vol; rsp++) {
    long r;
    int k;
    GAUGE_GROUP matrix, matrix2;

    r = sisp_and_t_to_si(geo, rsp, 0);

    one(&matrix);
    for (k = 0; k < geo->d_size[0]; k++) {
      times_equal(&matrix, &(GC->lattice[r][0]));
      r = nnp(geo, r, 0);
    }

    rep[rsp][0] = retr(&matrix);
    imp[rsp][0] = imtr(&matrix);

    equal(&matrix2, &matrix);

    for (k = 1; k < (int)floor(NCOLOR / 2.0); k++) {
      times_equal(&matrix2, &matrix);
      rep[rsp][k] = retr(&matrix2);
      imp[rsp][k] = imtr(&matrix2);
    }
  }

  for (j = 0; j < (int)floor(NCOLOR / 2); j++) {
    for (i = 0; i < geo->d_space_vol; i++) {
      repoly[j] += rep[i][j];
      impoly[j] += imp[i][j];
    }
  }

  for (j = 0; j < (int)floor(NCOLOR / 2.0); j++) {
    repoly[j] *= geo->d_inv_space_vol;
    impoly[j] *= geo->d_inv_space_vol;
  }

  for (i = 0; i < geo->d_space_vol; i++) {
    free(rep[i]);
    free(imp[i]);
  }
  free(rep);
  free(imp);
}

// computing all possible Polyakov loops for a given lattice
void polyvec(Gauge_Conf const *const GC, Geometry const *const geo,
             double complex *polyvec) {
  int rsp;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(rsp)
#endif
  for (rsp = 0; rsp < geo->d_space_vol; rsp++) {
    long r;
    int i;
    GAUGE_GROUP matrix;

    r = sisp_and_t_to_si(geo, rsp, 0);

    one(&matrix);
    for (i = 0; i < geo->d_size[0]; i++) {
      times_equal(&matrix, &(GC->lattice[r][0]));
      r = nnp(geo, r, 0);
    }

    polyvec[rsp] = retr(&matrix) + I * imtr(&matrix);
  }
}

// Polyakov loops in momentum space for a given value of spatial momentum
// \vec{p}
void polyakov_FT(Geometry const *const geo, double complex const *const polyvec,
                 double complex *polyakov_FT, double *spatial_momentum) {
  int cartcoord[STDIM];
  long rsp;
  double complex polyakov_FT_tmp;

  polyakov_FT_tmp = 0.0 + 0.0 * I;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(rsp, i, r, aux)         \
    reduction(+ : polyakov_FT_tmp)
#endif
  for (rsp = 0; rsp < geo->d_space_vol; rsp++) {
    int i;
    long r;
    double complex aux;

    r = sisp_and_t_to_si(geo, rsp, 0);

    si_to_cart(cartcoord, r, geo);

    aux = 1.0 + I * 0.0;
    for (i = 0; i < STDIM - 1; i++) {
      aux *= cexp(I * spatial_momentum[i] * cartcoord[i + 1]);
    }
    aux *= polyvec[rsp];

    polyakov_FT_tmp += aux;
  }

  *(polyakov_FT) = polyakov_FT_tmp * geo->d_inv_space_vol; // /(STDIM - 1) ????
}

// Compute the Polyakov loop correlator in momentum space for a given momentum
void polyakov_corr_FT(Geometry const *const geo,
                      double complex const *const polyvec, double *G_FT,
                      double *spatial_momentum) {
  double complex A;

  polyakov_FT(geo, polyvec, &A, spatial_momentum);

  *G_FT = creal(A) * creal(A) + cimag(A) * cimag(A);
}

// compute the Polyakov loop correlator up to a fixed distance
void polyakov_corr(Geometry const *const geo, GParam const *const param,
                   double complex const *const polyvec,
                   double complex *polycorr) {
  int k;
  long rsp;
  double rep, imp;

  for (k = 0; k < param->d_poly_corr; k++) {
    imp = 0.0;
    rep = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) \
    reduction(+ : imp)
#endif
    for (rsp = 0; rsp < geo->d_space_vol; rsp++) {
      int j, t_tmp;
      long r, rsp_tmp;
      double complex p1, p2;

      for (j = 1; j < STDIM; j++) {
        int i;

        r = sisp_and_t_to_si(geo, rsp, 0);

        for (i = 0; i < k; i++) {
          r = nnp(geo, r, j);
        }

        si_to_sisp_and_t(&rsp_tmp, &t_tmp, geo, r);

        p1 = polyvec[rsp_tmp];
        p2 = polyvec[rsp];

        rep += creal(conj(p2) * p1);
        imp += cimag(conj(p2) * p1);
      }
    }

    *(polycorr + k) = (rep + I * imp) * geo->d_inv_space_vol / (STDIM - 1);
  }
}

// compute the local topological charge at point r
// see readme for more details
double loc_topcharge(Gauge_Conf const *const GC, Geometry const *const geo,
                     GParam const *const param, long r) {
  if (!(STDIM == 4 && NCOLOR > 1) && !(STDIM == 2 && NCOLOR == 1)) {
    (void)GC;
    (void)geo;
    (void)param;
    (void)r;
    fprintf(stderr,
            "Wrong number of dimensions or number of colors! (%s, %d)\n",
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  double ris = 0.0; // initialized just to avoid compiler warnings

#if (STDIM == 4 && NCOLOR > 1)
  GAUGE_GROUP aux1, aux2, aux3;
  double real1, real2, loc_charge;
  const double chnorm = 1.0 / (128.0 * PI * PI);
  int i, dir[4][3], sign;

  dir[0][0] = 0;
  dir[0][1] = 0;
  dir[0][2] = 0;

  dir[1][0] = 1;
  dir[1][1] = 2;
  dir[1][2] = 3;

  dir[2][0] = 2;
  dir[2][1] = 1;
  dir[2][2] = 1;

  dir[3][0] = 3;
  dir[3][1] = 3;
  dir[3][2] = 2;

  sign = -1;
  loc_charge = 0.0;

  for (i = 0; i < 3; i++) {
    clover(GC, geo, r, dir[0][i], dir[1][i], &aux1);
    clover(GC, geo, r, dir[2][i], dir[3][i], &aux2);

    times_dag2(&aux3, &aux2, &aux1); // aux3=aux2*(aux1^{dag})
    real1 = retr(&aux3) * NCOLOR;

    times(&aux3, &aux2, &aux1); // aux3=aux2*aux1
    real2 = retr(&aux3) * NCOLOR;

    loc_charge += ((double)sign * (real1 - real2));
    sign = -sign;
  }
  ris = (loc_charge * chnorm);
#endif

#if (STDIM == 2 && NCOLOR == 1)
  GAUGE_GROUP u1matrix;
  double angle;

  plaquettep_matrix(GC, geo, param, r, 0, 1, &u1matrix);
  angle = atan2(cimag(u1matrix.comp), creal(u1matrix.comp)) / PI2;

  ris = angle;
#endif

  return ris;
}

// compute the topological charge
// see readme for more details
double topcharge(Gauge_Conf const *const GC, Geometry const *const geo,
                 GParam const *const param) {
  if (!(STDIM == 4 && NCOLOR > 1) && !(STDIM == 2 && NCOLOR == 1)) {
    fprintf(stderr,
            "Wrong number of dimensions or number of colors! (%s, %d)\n",
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  double ris;
  long r;

  ris = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
#endif
  for (r = 0; r < (geo->d_volume); r++) {
    ris += loc_topcharge(GC, geo, param, r);
  }

  return ris;
}

// compute GParam::d_nummeas values of the topological charge after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topcharge_cooling(Gauge_Conf const *const GC, Geometry const *const geo,
                       GParam const *const param, double *charge,
                       double *meanplaq) {
  if (!(STDIM == 4 && NCOLOR > 1) && !(STDIM == 2 && NCOLOR == 1)) {
    fprintf(stderr,
            "Wrong number of dimensions or number of colors! (%s, %d)\n",
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  if (param->d_coolsteps > 0) // if using cooling
  {
    Gauge_Conf helperconf;
    double ris, plaqs, plaqt;
    int iter;

    init_gauge_conf_from_gauge_conf(&helperconf, GC, geo);
    // helperconf is a copy of the configuration

    for (iter = 0; iter < (param->d_coolrepeat); iter++) {
      cooling(&helperconf, geo, param->d_coolsteps);

      ris = topcharge(&helperconf, geo, param);
      charge[iter] = ris;

      plaquette(&helperconf, geo, &plaqs, &plaqt);
#if (STDIM == 4)
      meanplaq[iter] = 0.5 * (plaqs + plaqt);
#else
      meanplaq[iter] = plaqt;
#endif
    }

    free_gauge_conf(&helperconf, geo);
  } else // no cooling
  {
    double ris, plaqs, plaqt;
    int iter;

    ris = topcharge(GC, geo, param);
    plaquette(GC, geo, &plaqs, &plaqt);

    for (iter = 0; iter < (param->d_coolrepeat); iter++) {
      charge[iter] = ris;
#if (STDIM == 4)
      meanplaq[iter] = 0.5 * (plaqs + plaqt);
#else
      meanplaq[iter] = plaqt;
#endif
    }
  }
}

void perform_measures_localobs(Gauge_Conf const *const GC,
                               Geometry const *const geo,
                               GParam const *const param, FILE *datafilep,
                               FILE *monofilep) {
  int i, ws, wt, max_wt, max_ws;
  double plaqs, plaqt, polyre, polyim;

  plaquette(GC, geo, &plaqs, &plaqt);
  polyakov(GC, geo, &polyre, &polyim);

  fprintf(datafilep, "%.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim);

  max_wt = MIN(10, (int)geo->d_size[0] / 4);
  for (i = 1; i < STDIM; i++) {
    max_ws = MIN(10, (int)geo->d_size[i] / 4);
  }

  for (wt = 1; wt <= max_wt; wt++) {
    for (ws = 1; ws <= max_ws; ws++) {
      fprintf(datafilep, "%.12g ", Wilsont(GC, geo, wt, ws));
    }
  }

// topological observables
#if ((STDIM == 4 && NCOLOR > 1) || (STDIM == 2 && NCOLOR == 1))
  int err;
  double *charge, *meanplaq, charge_nocooling;

  charge_nocooling = topcharge(GC, geo, param);
  fprintf(datafilep, " %.12g ", charge_nocooling);

  err = posix_memalign((void **)&charge, (size_t)DOUBLE_ALIGN,
                       (size_t)param->d_coolrepeat * sizeof(double));
  if (err != 0) {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__,
            __LINE__);
    exit(EXIT_FAILURE);
  }
  err = posix_memalign((void **)&meanplaq, (size_t)DOUBLE_ALIGN,
                       (size_t)param->d_coolrepeat * sizeof(double));
  if (err != 0) {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__,
            __LINE__);
    exit(EXIT_FAILURE);
  }

  topcharge_cooling(GC, geo, param, charge, meanplaq);
  for (i = 0; i < param->d_coolrepeat; i++) {
    fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
  }
  fprintf(datafilep, "\n");

  free(charge);
  free(meanplaq);
#else
  fprintf(datafilep, "\n");
#endif
  fflush(datafilep);

  // monopole observables
  if (param->d_mon_meas == 1) {
#if STDIM == 4
    int subg, subgnum;
    Gauge_Conf helperconf;

    init_gauge_conf_from_gauge_conf(&helperconf, GC, geo);
    alloc_diag_proj_stuff(&helperconf, geo);

    // MAG gauge fixing
    max_abelian_gauge_fix(&helperconf, geo);

    // diagonal projection
    diag_projection(&helperconf, geo);

    // loop on all the U(1) subgroups
    if (NCOLOR > 1) {
      subgnum = NCOLOR - 1;
    } else {
      subgnum = 1;
    }
    for (subg = 0; subg < subgnum; subg++) {
      // extract the abelian component subg and save it to GC->u1_subg
      U1_extract(&helperconf, geo, subg);

      // compute monopole observables
      monopoles_obs(&helperconf, geo, param, subg, monofilep);
    }

    free_diag_proj_stuff(&helperconf, geo);
    free_gauge_conf(&helperconf, geo);

    fflush(monofilep);
#else
    (void)monofilep;
#endif
  }
}

void perform_measures_localobs_obc(Gauge_Conf const *const GC,
                                   Geometry const *const geo,
                                   GParam const *const param, FILE *datafilep,
                                   FILE *datafileW, FILE *datafilesW) {
  int dis, wt, ws;
  double plaqs, plaqt;

  for (dis = 0; dis <= param->d_dis_max; dis++) {
    plaquette_obc_fve(GC, geo, &plaqs, &plaqt, dis);
    fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
  }
  fprintf(datafilep, "\n");

  for (ws = 1; ws < (int)geo->d_size[1]; ws++) {
    for (wt = 1; wt < (int)geo->d_size[0]; wt++) {
      fprintf(datafileW, "%.12g ", Wilsont_obc_avg(GC, geo, wt, ws));
    }
  }
  fprintf(datafileW, "\n");

  if ((int)geo->d_size[1] == 3) {
    for (wt = 1; wt < (int)geo->d_size[0]; wt++) {
      fprintf(datafilesW, "%.12g ",
              staircase_Wilsont_xy_obc_avg(GC, geo, wt, 1, 0));
    }
    fprintf(datafilesW, "\n");
  }

  fflush(datafilep);
  fflush(datafileW);
  fflush(datafilesW);
}

// perform local observables in the case of trace deformation, it computes all
// the order parameters
void perform_measures_localobs_with_tracedef(Gauge_Conf const *const GC,
                                             Geometry const *const geo,
                                             GParam const *const param,
                                             FILE *datafilep, FILE *monofilep,
                                             double complex *poly_vec) {
  int i;
  int n_tr_def = (int)floor(NCOLOR / 2);
  double poly_sq[n_tr_def];
  double plaqs, plaqt, polyre[NCOLOR / 2 + 1],
      polyim[NCOLOR / 2 + 1]; // +1 just to avoid warning if NCOLOR=1

  plaquette(GC, geo, &plaqs, &plaqt);
  polyakov_for_tracedef(GC, geo, polyre, polyim);

  fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);

  for (i = 0; i < (int)floor(NCOLOR / 2); i++) {
    fprintf(datafilep, "%.12g %.12g ", polyre[i], polyim[i]);
  }

  polyvec(GC, geo, poly_vec);

  double spatial_momentum[STDIM - 1];
  double G_FT = 0.0;

  for (i = 0; i < STDIM - 1; i++) {
    spatial_momentum[i] = 0.0;
  }

  polyakov_corr_FT(geo, poly_vec, &G_FT, spatial_momentum);

  fprintf(datafilep, "%.12g ", G_FT);

  G_FT = 0.0;

  for (i = 0; i < STDIM - 1; i++) {
    int k;
    double p_min, tmp_G_FT;

    p_min = PI2 / (geo->d_size[i]);
    for (k = 0; k < STDIM - 1; k++) {
      if (k == i)
        spatial_momentum[k] = p_min;
      else
        spatial_momentum[k] = 0.0;
    }

    polyakov_corr_FT(geo, poly_vec, &tmp_G_FT, spatial_momentum);

    G_FT += tmp_G_FT;
  }

  G_FT /= (STDIM - 1);

  fprintf(datafilep, "%.12g ", G_FT);

  GAUGE_GROUP matrix, matrix2;

  for (i = 0; i < n_tr_def; i++) {
    poly_sq[i] = 0.0;
  }

  for (long rsp = 0; rsp < geo->d_space_vol; rsp++) {
    long r;

    r = sisp_and_t_to_si(geo, rsp, 0);

    one(&matrix);
    for (i = 0; i < geo->d_size[0]; i++) {
      times_equal(&matrix, &(GC->lattice[r][0]));
      r = nnp(geo, r, 0);
    }

    one(&matrix2);
    for (i = 0; i < n_tr_def; i++) {
      times_equal(&matrix2, &matrix);
      poly_sq[i] +=
          retr(&matrix2) * retr(&matrix2) + imtr(&matrix2) * imtr(&matrix2);
    }
  }

  for (i = 0; i < n_tr_def; i++) {
    poly_sq[i] *= geo->d_inv_space_vol;
    fprintf(datafilep, "%.12g ", poly_sq[i]);
  }

// topological observables
#if (STDIM == 4 && NCOLOR > 1 || (STDIM == 2 && NCOLOR == 1))
  int err;
  double *charge, *meanplaq, charge_nocooling;

  charge_nocooling = topcharge(GC, geo, param);

  fprintf(datafilep, "%.12g ", charge_nocooling);

  err = posix_memalign((void **)&charge, (size_t)DOUBLE_ALIGN,
                       (size_t)param->d_coolrepeat * sizeof(double));
  if (err != 0) {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__,
            __LINE__);
    exit(EXIT_FAILURE);
  }
  err = posix_memalign((void **)&meanplaq, (size_t)DOUBLE_ALIGN,
                       (size_t)param->d_coolrepeat * sizeof(double));
  if (err != 0) {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__,
            __LINE__);
    exit(EXIT_FAILURE);
  }

  topcharge_cooling(GC, geo, param, charge, meanplaq);
  for (i = 0; i < param->d_coolrepeat; i++) {
    fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
  }
  fprintf(datafilep, "\n");

  free(charge);
  free(meanplaq);
#else
  fprintf(datafilep, "\n");
#endif

  fflush(datafilep);

  // monopole observables
  if (param->d_mon_meas == 1) {
#if (STDIM == 4)
    Gauge_Conf helperconf;
    int subg, subgnum;

    init_gauge_conf_from_gauge_conf(&helperconf, GC, geo);
    alloc_diag_proj_stuff(&helperconf, geo);

    // MAG gauge fixing
    max_abelian_gauge_fix(&helperconf, geo);

    // diagonal projection
    diag_projection(&helperconf, geo);

    // loop on all the U(1) subgroups
    if (NCOLOR > 1) {
      subgnum = NCOLOR - 1;
    } else {
      subgnum = 1;
    }
    for (subg = 0; subg < subgnum; subg++) {
      // extract the abelian component subg and save it to GC->u1_subg
      U1_extract(&helperconf, geo, subg);

      // compute monopole observables
      monopoles_obs(&helperconf, geo, param, subg, monofilep);
    }

    free_diag_proj_stuff(&helperconf, geo);
    free_gauge_conf(&helperconf, geo);

    fflush(monofilep);
#else
    (void)monofilep;
#endif
  }
}

// compute the average value of \sum_{flavours} Re(H_x U_{x,mu} H_{x+mu})
void higgs_interaction(Gauge_Conf const *const GC, Geometry const *const geo,
                       double *he) {
  long r;
  double ris = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
#endif
  for (r = 0; r < (geo->d_volume); r++) {
    int i;
    double aux = 0.0;
    GAUGE_VECS v1;
    GAUGE_GROUP matrix;

    for (i = 0; i < STDIM; i++) {
      equal(&matrix, &(GC->lattice[r][i]));

      matrix_times_vector_all_vecs(&v1, &matrix, &(GC->higgs[nnp(geo, r, i)]));
      aux += re_scal_prod_vecs(&(GC->higgs[r]), &v1);
    }

    ris += aux;
  }

  ris /= (double)STDIM;
  ris *= geo->d_inv_vol;

  *he = ris;
}

// compute flavour related observables
//
// flavour matrices Qh and Dh HAVE TO BE INITIALIZED before calling this
// function
//
// tildeG0=ReTr[(\sum_x Q_x)(\sum_y Q_y)]/volume/NHIGGS
// tildeGminp=ReTr[(\sum_x Q_xe^{ipx})(\sum_y Q_ye^{-ipy)]/volume/NHIGGS
//
// tildeG0 is susceptibility/NHIGGS, tildeGminp is used to compute the 2nd
// momentum correlation function
//
// tildeD0=conj(\sum_x D_x) (\sum_y D_y) / volume
// tildeDminp=(\sum_x D_x e^{ipx}) conj(\sum_y D_y e^{ipy}) /volume
//
// tildeD0 is a U1 susceptibility, tildeDminp is used to compute the 2nd
// momentum correlation function
void compute_flavour_observables(Gauge_Conf const *const GC,
                                 Geometry const *const geo, double *tildeG0,
                                 double *tildeGminp, double *tildeD0,
                                 double *tildeDminp) {
  int coord[STDIM];
  long r;
  const double p = 2.0 * PI / (double)geo->d_size[1];
  double complex D, Dp;
  FMatrix Q, Qp, Qmp, tmp1, tmp2;

  // Q =sum_x Q_x
  // Qp=sum_x e^{ipx}Q_x
  // Qmp=sum_x e^{-ipx}Q_x
  //
  // D, Dp and are the analogous of Q and Qp for D

  D = 0.0 + 0.0 * I;
  Dp = 0.0 + 0.0 * I;

  zero_FMatrix(&Q);
  zero_FMatrix(&Qp);
  zero_FMatrix(&Qmp);
  for (r = 0; r < (geo->d_volume); r++) {
    equal_FMatrix(&tmp1, &(GC->Qh[r]));
    equal_FMatrix(&tmp2, &tmp1);

    plus_equal_FMatrix(&Q, &tmp1);
    D += (GC->Dh[r]);

    si_to_cart(coord, r, geo);

    times_equal_complex_FMatrix(&tmp1, cexp(I * ((double)coord[1]) * p));
    plus_equal_FMatrix(&Qp, &tmp1);
    Dp += ((GC->Dh[r]) * cexp(I * ((double)coord[1]) * p));

    times_equal_complex_FMatrix(&tmp2, cexp(-I * ((double)coord[1]) * p));
    plus_equal_FMatrix(&Qmp, &tmp2);
  }

  equal_FMatrix(&tmp1, &Q);
  times_equal_FMatrix(&tmp1, &Q);

  *tildeG0 = retr_FMatrix(&tmp1) * geo->d_inv_vol;
  *tildeD0 = creal(conj(D) * D) * geo->d_inv_vol;

  equal_FMatrix(&tmp1, &Qp);
  times_equal_FMatrix(&tmp1, &Qmp);
  *tildeGminp = retr_FMatrix(&tmp1) * geo->d_inv_vol;
  *tildeDminp = creal(Dp * conj(Dp)) * geo->d_inv_vol;
}

// compute correlators of flavour observables
//
// flavour matrices Qh and Dh HAVE TO BE INITIALIZED before calling this
// function
//
// corrQQ is the correlato ReTr[Q_x Q_{x+d}]/N_higgs
// corr0string0 is the correlator \sum_f Re[hf^{dag}
// U_{x,1}U_{x+1,1}....Q_{x+d-1,1} hf], where hf is the f-th flavour
// corr0string1 is the correlator Re[h0^{dag} U_{x,1}U_{x+1,1}....U_{x+d-1,1}
// h1], where h1 is the second flavour
void compute_flavour_observables_corr(Gauge_Conf const *const GC,
                                      Geometry const *const geo, double *corrQQ,
                                      double *corr0string0,
                                      double *corr0string1) {
  int dist;
  long r;
  double accumulator1, accumulator2;

  for (dist = 0; dist < geo->d_size[1]; dist++) {
    accumulator1 = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r)                      \
    reduction(+ : accumulator1)
#endif
    for (r = 0; r < (geo->d_volume); r++) {
      int i;
      long r1;
      FMatrix tmp1;

      equal_FMatrix(&tmp1, &(GC->Qh[r]));
      r1 = r;
      for (i = 0; i < dist; i++) {
        r1 = nnp(geo, r1, 1);
      }
      times_equal_FMatrix(&tmp1, &(GC->Qh[r1]));
      accumulator1 += retr_FMatrix(&tmp1);
    }
    accumulator1 *= geo->d_inv_vol;
    corrQQ[dist] = accumulator1;

    accumulator1 = 0.0;
    accumulator2 = 0.0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r)                      \
    reduction(+ : accumulator1) reduction(+ : accumulator2)
#endif
    for (r = 0; r < (geo->d_volume); r++) {
      int i;
      long r1;
      GAUGE_VECS phi1, phi2;
      GAUGE_GROUP U;

      equal_vecs(&phi1, &(GC->higgs[r]));
      r1 = r;
      one(&U);
      for (i = 0; i < dist; i++) {
        times_equal(&U, &(GC->lattice[r1][1]));
        r1 = nnp(geo, r1, 1);
      }
      matrix_times_vector_all_vecs(&phi2, &U, &(GC->higgs[r1]));
      accumulator1 += re_scal_prod_vecs(&phi1, &phi2);
#if NHIGGS > 1
      accumulator2 += re_scal_prod_single_vecs(&phi1, &phi2, 0, 1);
#else
      accumulator2 += 0.0;
#endif
    }
    accumulator1 *= geo->d_inv_vol;
    accumulator2 *= geo->d_inv_vol;

    corr0string0[dist] = accumulator1;
    corr0string1[dist] = accumulator2;
  }
}

void perform_measures_higgs(Gauge_Conf *GC, Geometry const *const geo,
                            FILE *datafilep) {
  double plaqs, plaqt, polyre, polyim, he, tildeG0, tildeGminp, tildeD0,
      tildeDminp;
  long r;

  plaquette(GC, geo, &plaqs, &plaqt);
  polyakov(GC, geo, &polyre, &polyim);
  higgs_interaction(GC, geo, &he);

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r)
#endif
  for (r = 0; r < (geo->d_volume); r++) {
    init_FMatrix_vecs(&(GC->Qh[r]), &(GC->higgs[r]));
    GC->Dh[r] = HiggsU1Obs_vecs(&(GC->higgs[r]));
  }

  compute_flavour_observables(GC, geo, &tildeG0, &tildeGminp, &tildeD0,
                              &tildeDminp);

  fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
  fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
  fprintf(datafilep, "%.12g ", he);
  fprintf(datafilep, "%.12g %.12g ", tildeG0, tildeGminp);
  fprintf(datafilep, "%.12g %.12g ", tildeD0, tildeDminp);

  /*
  // for correlators

  int err, i;
  double *corrQQ, *corr0string0, *corr0string1;
  err=posix_memalign((void**) &(corrQQ), (size_t) DOUBLE_ALIGN, (size_t)
  param->d_size[1] * sizeof(double)); err+=posix_memalign((void**)
  &(corr0string0), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[1] *
  sizeof(double)); err+=posix_memalign((void**) &(corr0string1), (size_t)
  DOUBLE_ALIGN, (size_t) param->d_size[1] * sizeof(double)); if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the correlators! (%s, %d)\n",
  __FILE__, __LINE__); exit(EXIT_FAILURE);
    }

  compute_flavour_observables_corr(GC,
                                   geo,
                                   param,
                                   corrQQ,
                                   corr0string0,
                                   corr0string1);
  for(i=0; i<param->d_size[1]; i++)
     {
     fprintf(datafilep, "%.12g ", corrQQ[i]);
     }
  for(i=0; i<param->d_size[1]; i++)
     {
     fprintf(datafilep, "%.12g ", corr0string0[i]);
     }
   for(i=0; i<param->d_size[1]; i++)
     {
     fprintf(datafilep, "%.12g ", corr0string1[i]);
     }

  free(corrQQ);
  free(corr0string0);
  free(corr0string1);
  */

  fprintf(datafilep, "\n");

  fflush(datafilep);
}

// this is a function to be used just to test some fine points
// most notably TrP^2=TrQ^2+1/NHIGGS
void perform_measures_higgs_for_testing(Gauge_Conf *GC,
                                        Geometry const *const geo,
                                        FILE *datafilep) {
  double plaqs, plaqt, polyre, polyim, he, p2, tildeG0, tildeGminp, tildeD0,
      tildeDminp;
  long r;

  plaquette(GC, geo, &plaqs, &plaqt);
  polyakov(GC, geo, &polyre, &polyim);
  higgs_interaction(GC, geo, &he);

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r)
#endif
  for (r = 0; r < (geo->d_volume); r++) {
    init_FMatrix_vecs(&(GC->Qh[r]), &(GC->higgs[r]));
    GC->Dh[r] = HiggsU1Obs_vecs(&(GC->higgs[r]));
  }

  fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);
  fprintf(datafilep, "%.12g %.12g ", polyre, polyim);
  fprintf(datafilep, "%.12g ", he);

  compute_flavour_observables(GC, geo, &tildeG0, &tildeGminp, &tildeD0,
                              &tildeDminp);

  fprintf(datafilep, "%.12g %.12g ", tildeG0, tildeGminp);
  fprintf(datafilep, "%.12g %.12g ", tildeD0, tildeDminp);

  p2 = 0.0;
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : p2)
#endif
  for (r = 0; r < (geo->d_volume); r++) {
    FMatrix tmp1, tmp2;
    equal_FMatrix(&tmp1, &(GC->Qh[r]));
    equal_FMatrix(&tmp2, &tmp1);
    times_equal_FMatrix(&tmp1, &tmp2);
    p2 += retr_FMatrix(&tmp1) * NHIGGS + 1. / NHIGGS;
  }
  p2 *= geo->d_inv_vol;

  fprintf(datafilep, "%.12g ", p2);

  fprintf(datafilep, "\n");

  fflush(datafilep);
}

// fix maximal abelian gauge
// following the procedure described in
// C. Bonati, M. D'Elia Nuc. Phys. B 877 (2013) 233-259 [ 1308.0302 ]
void max_abelian_gauge_fix(Gauge_Conf *GC, Geometry const *const geo) {
  int i, dir;
  long r;
  double lambda[NCOLOR];
  const double overrelaxparam = 1.85; // 1.0 means no overrelaxation
  const double target = 1.0e-8;
  double nondiag, nondiagaux;

  // inizialize the matrix lambda = diag((N-1)/2, (N-1)/2-1, ..., -(N-1)/2)
  for (i = 0; i < NCOLOR; i++) {
    lambda[i] = ((double)NCOLOR - 1.) / 2. - (double)i;
  }

  nondiag = 1;
  while (nondiag > target) {
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r, dir)
#endif
    for (r = 0; r < geo->d_volume / 2; r++) {
      GAUGE_GROUP G_mag, help,
          X_links[2 * STDIM]; // X_links contains the 2*STDIM links used in
                              // the computation of X(n)

      // initialize X_links[2*STDIM] with the 2*STDIM links surrounding the
      // point r links 0 to (STDIM-1) are forward, while links STDIM to
      // (2*STDIM-1) are backwards.
      for (dir = 0; dir < STDIM; dir++) {
        equal(&(X_links[dir]), &(GC->lattice[r][dir]));
        equal(&(X_links[dir + STDIM]), &(GC->lattice[nnm(geo, r, dir)][dir]));
      }

      comp_MAG_gauge_transformation(X_links, lambda, overrelaxparam, &G_mag);

      // apply the gauge transformation
      for (dir = 0; dir < STDIM; dir++) {
        times(&help, &G_mag, &(GC->lattice[r][dir]));
        equal(&(GC->lattice[r][dir]), &help);

        times_equal_dag(&(GC->lattice[nnm(geo, r, dir)][dir]), &G_mag);
      }
    }

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r, dir)
#endif
    for (r = geo->d_volume / 2; r < geo->d_volume; r++) {
      GAUGE_GROUP G_mag, help,
          X_links[2 * STDIM]; // X_links contains the 2*STDIM links used in
                              // the computation of X(n)

      // initialize X_links[2*STDIM] with the 2*STDIM links surrounding the
      // point r links 0 to (STDIM-1) are forward, while links STDIM to
      // (2*STDIM-1) are backwards.
      for (dir = 0; dir < STDIM; dir++) {
        equal(&(X_links[dir]), &(GC->lattice[r][dir]));
        equal(&(X_links[dir + STDIM]), &(GC->lattice[nnm(geo, r, dir)][dir]));
      }

      comp_MAG_gauge_transformation(X_links, lambda, overrelaxparam, &G_mag);

      // apply the gauge transformation
      for (dir = 0; dir < STDIM; dir++) {
        times(&help, &G_mag, &(GC->lattice[r][dir]));
        equal(&(GC->lattice[r][dir]), &help);

        times_equal_dag(&(GC->lattice[nnm(geo, r, dir)][dir]), &G_mag);
      }
    }

    // nondiagaux is the sum of the squares of the out-diagonal terms
    nondiagaux = 0;

#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r, dir)                 \
    reduction(+ : nondiagaux)
#endif
    for (r = 0; r < geo->d_volume; r++) {
      GAUGE_GROUP X_links[2 * STDIM]; // X_links contains the 2*STDIM links
                                      // used in the computation of X(n)
      double counter;

      for (dir = 0; dir < STDIM; dir++) {
        equal(&(X_links[dir]), &(GC->lattice[r][dir]));
        equal(&(X_links[dir + STDIM]), &(GC->lattice[nnm(geo, r, dir)][dir]));
      }
      comp_outdiagnorm_of_X(X_links, lambda, &counter);
      nondiagaux += counter;
    }

    nondiag = nondiagaux * geo->d_inv_vol / (double)NCOLOR / (double)NCOLOR;

    // printf("%g  %g\n", nondiag, nondiag/target);
    // fflush(stdout);
  }

// unitarize all the links
#ifdef OPENMP_MODE
#pragma omp parallel for num_threads(NTHREADS) private(r, dir)
#endif
  for (r = 0; r < (geo->d_volume); r++) {
    for (dir = 0; dir < STDIM; dir++) {
      unitarize(&(GC->lattice[r][dir]));
    }
  }
}

// extract the diagonal part of the links after gauge fixing.
// the phases are saved in GC->diag_proj but these are NOT the monopole phases
// (see U1_extract)
void diag_projection(Gauge_Conf *GC, Geometry const *const geo) {
  int dir;
  long r;

  for (r = 0; r < geo->d_volume; r++) {
    for (dir = 0; dir < STDIM; dir++) {
      diag_projection_single_site(GC, &(GC->lattice[r][dir]), r, dir);
    }
  }
}

// extract the abelian components of the link
// following the procedure described in
// Bonati, D'Elia https://arxiv.org/abs/1308.0302
// and save them in GC->u1_subg
//
// also intialize GC->uflag to zero
void U1_extract(Gauge_Conf *GC, Geometry const *const geo, int subg) {
  int dir, i;
  long r;

  for (r = 0; r < geo->d_volume; r++) {
    for (dir = 0; dir < STDIM; dir++) {
      GC->u1_subg[r][dir] = 0.0;
      for (i = 0; i <= subg; i++) {
        GC->u1_subg[r][dir] += GC->diag_proj[r][dir][i];
      }

      GC->uflag[r][dir] = 0;
    }
  }
}

// Compute the forward derivative of the abelian part of the plaquette Fjk in
// direction i. the angle is chosen in between -pi and pi.
void Di_Fjk(Gauge_Conf *GC, Geometry const *const geo, long r, int idir,
            int jdir, int kdir, double *DiFjk)

{
  double prpi, pr; // pr -> plaquette at site r, prpi plaquette at site r+idir

  //
  //       ^ k
  //       |
  //       +---<---+
  //       |       |
  //       V       ^         pr
  //       |       |
  //       +--->---+---> j
  //       r
  //

  pr = GC->u1_subg[r][jdir] - GC->u1_subg[r][kdir];
  pr += GC->u1_subg[nnp(geo, r, jdir)][kdir] -
        GC->u1_subg[nnp(geo, r, kdir)][jdir];

  //
  //       ^ k
  //       |   (2)
  //       +---<---+
  //       |       |
  //   (3) V       ^ (1)        prpi
  //       |       |
  //       +--->---+---> j
  //      r+i    (4)
  //

  r = nnp(geo, r, idir);

  prpi = GC->u1_subg[r][jdir] - GC->u1_subg[r][kdir];
  prpi += GC->u1_subg[nnp(geo, r, jdir)][kdir] -
          GC->u1_subg[nnp(geo, r, kdir)][jdir];

  *DiFjk = 2.0 * (atan(tan(prpi / 2.0)) - atan(tan(pr / 2.0)));
}

// compute the DeGrand-DeTar currents
int DeGrand_current(Gauge_Conf *GC, Geometry const *const geo, long r,
                    int dir) {
  if (STDIM != 4) {
    fprintf(stderr, "Wrong number of dimensions! (%s, %d)\n", __FILE__,
            __LINE__);
    exit(EXIT_FAILURE);
  }

  double der1, der2, der3;
  int ris;

  if (dir == 0) {
    Di_Fjk(GC, geo, r, 1, 2, 3, &der1);
    Di_Fjk(GC, geo, r, 3, 1, 2, &der2);
    Di_Fjk(GC, geo, r, 2, 1, 3, &der3);

    ris = (int)round(((der1 + der2 - der3) / PI2));
  } else if (dir == 1) {
    Di_Fjk(GC, geo, r, 3, 2, 0, &der1);
    Di_Fjk(GC, geo, r, 0, 3, 2, &der2);
    Di_Fjk(GC, geo, r, 2, 3, 0, &der3);

    ris = (int)round(((der1 + der2 - der3) / PI2));
  } else if (dir == 2) {
    Di_Fjk(GC, geo, r, 3, 0, 1, &der1);
    Di_Fjk(GC, geo, r, 0, 1, 3, &der2);
    Di_Fjk(GC, geo, r, 1, 0, 3, &der3);

    ris = (int)round(((der1 + der2 - der3) / PI2));
  } else {
    Di_Fjk(GC, geo, r, 0, 2, 1, &der1);
    Di_Fjk(GC, geo, r, 2, 1, 0, &der2);
    Di_Fjk(GC, geo, r, 1, 2, 0, &der3);

    ris = (int)round(((der1 + der2 - der3) / PI2));
  }

  return ris;
}

// search for monopole wrappings passing from r_tback
// this function can be invoked in two different ways
//
// or r=nnp(geo, r_tback, 0) and DeGrand_current(GC, geo, r_tback, 0)!=0
// (forward case) or r=nnm(geo, r_tback, 0) and DeGrand_current(GC, geo, r,
// 0)!=0  (backward case)
//
// GC->uflag[][] is initialized in monopole_obs
//
// nonzero DeGrand_current(GC, geo, nnp(geo, r, dir), dir ) are associated to
// uflag[r][dir]
//
// num_wrap = number of wrappings
void wrap_search(Gauge_Conf *GC, Geometry const *const geo,
                 GParam const *const param, long r, long r_tback,
                 int *num_wrap) {
#if STDIM != 4
  fprintf(stderr, "Wrong number of dimensions! (%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
#endif

  int dir, n_mu;

  if (r == r_tback) {
    return;
  } else {
    // forward case
    for (dir = 0; dir < STDIM; dir++) {
      n_mu = DeGrand_current(GC, geo, nnp(geo, r, dir), dir);

      // if not all the monopole currents have been followed
      if (n_mu > GC->uflag[r][dir]) {
        GC->uflag[r][dir] += 1;

        if ((geo->d_timeslice[r] == geo->d_size[0] - 1) && (dir == 0)) {
          *num_wrap += 1;
        }

        wrap_search(GC, geo, param, nnp(geo, r, dir), r_tback, num_wrap);

        return;
      }
    }

    // backward case
    for (dir = 0; dir < STDIM; dir++) {
      n_mu = DeGrand_current(GC, geo, r, dir);

      if (n_mu < GC->uflag[nnm(geo, r, dir)][dir]) {
        GC->uflag[nnm(geo, r, dir)][dir] -= 1;

        if ((geo->d_timeslice[r] == 0) && (dir == 0)) {
          *num_wrap -= 1;
        }

        wrap_search(GC, geo, param, nnm(geo, r, dir), r_tback, num_wrap);

        return;
      }
    }
  }
}

// GC->uflag[][] has to be initialized to zero before calling this function
// (when GC->uflag is allocated it is also initialized to zero)
void monopoles_obs(Gauge_Conf *GC, Geometry const *const geo,
                   GParam const *const param, int subg, FILE *monofilep) {
  double mean_wrap;
  long r, rsp, r_tback, r_tbackback;
  int n_mu, num_wrap, mono_charge;
  int cartcoord[4];

  mean_wrap = 0.0; // mean value of monopole wraps for unit volume

  for (rsp = 0; rsp < geo->d_space_vol; rsp++) {
    r = sisp_and_t_to_si(geo, rsp, 1);                            // t=1 slice
    r_tback = sisp_and_t_to_si(geo, rsp, 0);                      // t=0 slice
    r_tbackback = sisp_and_t_to_si(geo, rsp, geo->d_size[0] - 1); // t=T-1 slice

    // check the t=1 temporal slice to find monopoles currents
    n_mu = DeGrand_current(GC, geo, r, 0);

    // start following monopole charge in forward direction. Maximum lattice
    // charge is +2 so we try twice
    for (mono_charge = 0; mono_charge < 2; mono_charge++) {
      // nonzero DeGrand_current(GC, geo, nnp(geo, r, dir), dir ) are
      // associated to uflag[r][dir]
      if (n_mu > GC->uflag[r_tback][0]) {
        GC->uflag[r_tback][0] += 1;

        num_wrap = 0;
        wrap_search(GC, geo, param, r, r_tback, &num_wrap);

        mean_wrap += abs(num_wrap);

        lexeo_to_cart(cartcoord, r_tback, geo);
        if (n_mu == 1) {
          fprintf(monofilep, "%ld ", GC->update_index);

          for (int k = 0; k < 4; k++) {
            fprintf(monofilep, "%d ", cartcoord[k]);
          }
          fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
        } else if (GC->uflag[r_tback][0] ==
                   1) // this is to print only once monopole of charge +2
        {
          fprintf(monofilep, "%ld ", GC->update_index);

          for (int k = 0; k < 4; k++) {
            fprintf(monofilep, "%d ", cartcoord[k]);
          }
          fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
        }
      }
    }

    n_mu = DeGrand_current(GC, geo, r_tback, 0);

    // start following monopole charge in backward direction. Maximum lattice
    // charge is +2 so we try twice
    for (mono_charge = 0; mono_charge < 2; mono_charge++) {
      // nonzero DeGrand_current(GC, geo, nnp(geo, r, dir), dir ) are
      // associated to uflag[r][dir]
      if (n_mu < GC->uflag[r_tbackback][0]) {
        GC->uflag[r_tbackback][0] -= 1;

        num_wrap = -1;
        wrap_search(GC, geo, param, r_tbackback, r_tback, &num_wrap);

        lexeo_to_cart(cartcoord, r_tback, geo);
        if (n_mu == -1) {
          fprintf(monofilep, "%ld ", GC->update_index);

          for (int k = 0; k < 4; k++) {
            fprintf(monofilep, "%d ", cartcoord[k]);
          }
          fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
        } else if (GC->uflag[r][0] ==
                   -1) // this is to print only once monopole of charge +2
        {
          fprintf(monofilep, "%ld ", GC->update_index);

          for (int k = 0; k < 4; k++) {
            fprintf(monofilep, "%d ", cartcoord[k]);
          }
          fprintf(monofilep, "%d %d %d\n", subg, n_mu, num_wrap);
        }
      }
    }
  }
}

#endif
