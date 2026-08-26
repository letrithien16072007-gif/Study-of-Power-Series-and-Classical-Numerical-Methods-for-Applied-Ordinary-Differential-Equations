/**
 * ============================================================================
 * ODE SOLVER ACCURACY, TIMING, AND PADÉ APPROXIMATION BENCHMARKS
 * ============================================================================
 * Project: Comparative Study of Power Series and Classical Solvers for Applied
 *          ODEs
 * Department of Mathematics, FPT University, HCMC, Vietnam
 *
 * Cross-platform compatible (Linux / macOS / Windows MinGW / MSVC)
 * ============================================================================
 */

#ifndef _POSIX_C_SOURCE
#define POSIX_C_SOURCE 199309L
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#include <windows.h>
typedef struct {
  LARGE_INTEGER start;
  LARGE_INTEGER end;
} benchmark_timer_t;

static inline void timer_start(benchmark_timer_t *t) {
  QueryPerformanceCounter(&t->start);
}

static inline void timer_end(benchmark_timer_t *t) {
  QueryPerformanceCounter(&t->end);
}

static inline double timer_elapsed_us(benchmark_timer_t *t) {
  LARGE_INTEGER freq;
  QueryPerformanceFrequency(&freq);
  return (double)(t->end.QuadPart - t->start.QuadPart) * 1000000.0 /
         (double)freq.QuadPart;
}

static inline void ensure_data_dir(void) { _mkdir("data"); }
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

typedef struct {
  struct timespec start;
  struct timespec end;
} benchmark_timer_t;

static inline void timer_start(benchmark_timer_t *t) {
  clock_gettime(CLOCK_MONOTONIC, &t->start);
}

static inline void timer_end(benchmark_timer_t *t) {
  clock_gettime(CLOCK_MONOTONIC, &t->end);
}

static inline double timer_elapsed_us(benchmark_timer_t *t) {
  return ((t->end.tv_sec - t->start.tv_sec) * 1000000.0) +
         ((t->end.tv_nsec - t->start.tv_nsec) / 1000.0);
}

static inline void ensure_data_dir(void) { mkdir("data", 0755); }
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SPARSE_PTS 50

/**
 * @brief Benchmark Configuration Parameters (Customisable via CLI args).
 */
typedef struct {
  int runs;
  double x_eval;
  double x0;
  double tol;
  double newton_y0;
  double newton_tm;
  double newton_k;
  double airy_y0;
  double airy_v0;
  double logistic_y0;
  double logistic_k;
  double logistic_r;
  double step_h1;
  double step_h2;
  int series_n;
  int target_part;
} benchmark_config_t;

/**
 * @brief Result structure for Algorithm 1 execution.
 */
typedef struct {
  double alpha;
  double beta;
  double R_conv;
  int M_min;
  double x0_threshold;
  double pade_val;
  double exact_val;
  double error;
  int converged_M;
  int is_taylor;
  char decision_note[256];
} alg1_result_t;

static void set_default_config(benchmark_config_t *cfg) {
  cfg->runs = 10000;
  cfg->x_eval = 5.0;
  cfg->x0 = 0.0;
  cfg->tol = 1e-4;
  cfg->newton_y0 = 100.0;
  cfg->newton_tm = 20.0;
  cfg->newton_k = 0.5;
  cfg->airy_y0 = 1.0;
  cfg->airy_v0 = 0.0;
  cfg->logistic_y0 = 10.0;
  cfg->logistic_k = 100.0;
  cfg->logistic_r = 1.0;
  cfg->step_h1 = 0.1;
  cfg->step_h2 = 0.01;
  cfg->series_n = 20;
  cfg->target_part = 0; /* 0 means run all parts */
}

static void print_usage(const char *prog_name) {
  printf("Usage: %s [options]\n", prog_name);
  printf("Options:\n");
  printf("  --runs <int>           Number of benchmark runs for timing "
         "(default: 10000)\n");
  printf("  --x-eval <double>      Evaluation target point x (default: 5.0)\n");
  printf("  --x0 <double>          Expansion point x0 (default: 0.0)\n");
  printf("  --tol <double>         Algorithm 1 tolerance (default: 1e-4)\n");
  printf(
      "  --newton-y0 <double>   Newton cooling initial y0 (default: 100.0)\n");
  printf("  --newton-tm <double>   Newton cooling ambient temp Tm (default: "
         "20.0)\n");
  printf("  --newton-k <double>    Newton cooling rate k (default: 0.5)\n");
  printf("  --airy-y0 <double>     Airy initial y(0) (default: 1.0)\n");
  printf("  --airy-v0 <double>     Airy initial y'(0) (default: 0.0)\n");
  printf("  --logistic-y0 <double> Logistic initial y0 (default: 10.0)\n");
  printf("  --logistic-k <double>  Logistic capacity K (default: 100.0)\n");
  printf("  --logistic-r <double>  Logistic growth rate r (default: 1.0)\n");
  printf("  --step-h1 <double>     Coarse step size h1 (default: 0.1)\n");
  printf("  --step-h2 <double>     Fine step size h2 (default: 0.01)\n");
  printf("  --series-n <int>       Power series terms N (default: 20)\n");
  printf("  --part <int>           Run specific part [1|2|3|4|0=all] (default: "
         "0)\n");
  printf("  --help, -h             Show this help message\n");
}

static void parse_cli_args(int argc, char **argv, benchmark_config_t *cfg) {
  set_default_config(cfg);
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--runs") == 0 && i + 1 < argc) {
      cfg->runs = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--x-eval") == 0 && i + 1 < argc) {
      cfg->x_eval = atof(argv[++i]);
    } else if (strcmp(argv[i], "--x0") == 0 && i + 1 < argc) {
      cfg->x0 = atof(argv[++i]);
    } else if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc) {
      cfg->tol = atof(argv[++i]);
    } else if (strcmp(argv[i], "--newton-y0") == 0 && i + 1 < argc) {
      cfg->newton_y0 = atof(argv[++i]);
    } else if (strcmp(argv[i], "--newton-tm") == 0 && i + 1 < argc) {
      cfg->newton_tm = atof(argv[++i]);
    } else if (strcmp(argv[i], "--newton-k") == 0 && i + 1 < argc) {
      cfg->newton_k = atof(argv[++i]);
    } else if (strcmp(argv[i], "--airy-y0") == 0 && i + 1 < argc) {
      cfg->airy_y0 = atof(argv[++i]);
    } else if (strcmp(argv[i], "--airy-v0") == 0 && i + 1 < argc) {
      cfg->airy_v0 = atof(argv[++i]);
    } else if (strcmp(argv[i], "--logistic-y0") == 0 && i + 1 < argc) {
      cfg->logistic_y0 = atof(argv[++i]);
    } else if (strcmp(argv[i], "--logistic-k") == 0 && i + 1 < argc) {
      cfg->logistic_k = atof(argv[++i]);
    } else if (strcmp(argv[i], "--logistic-r") == 0 && i + 1 < argc) {
      cfg->logistic_r = atof(argv[++i]);
    } else if (strcmp(argv[i], "--step-h1") == 0 && i + 1 < argc) {
      cfg->step_h1 = atof(argv[++i]);
    } else if (strcmp(argv[i], "--step-h2") == 0 && i + 1 < argc) {
      cfg->step_h2 = atof(argv[++i]);
    } else if (strcmp(argv[i], "--series-n") == 0 && i + 1 < argc) {
      cfg->series_n = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--part") == 0 && i + 1 < argc) {
      cfg->target_part = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
      print_usage(argv[0]);
      exit(0);
    }
  }
}

/**
 * @brief Solves a linear system Ax = b using Gaussian elimination with partial
 * pivoting.
 */
static int solve_linear_system(int n, const double *A, const double *b,
                               double *x) {
  double *M = (double *)malloc(n * (n + 1) * sizeof(double));
  if (!M) {
    return 0;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      M[(i * (n + 1)) + j] = A[(i * n) + j];
    }
    M[(i * (n + 1)) + n] = b[i];
  }

  for (int i = 0; i < n; i++) {
    int pivot_row = i;
    double max_val = fabs(M[(i * (n + 1)) + i]);
    for (int r = i + 1; r < n; r++) {
      double val = fabs(M[(r * (n + 1)) + i]);
      if (val > max_val) {
        max_val = val;
        pivot_row = r;
      }
    }
    if (max_val < 1e-12) {
      free(M);
      return 0;
    }
    if (pivot_row != i) {
      for (int c = i; c <= n; c++) {
        double temp = M[(i * (n + 1)) + c];
        M[(i * (n + 1)) + c] = M[(pivot_row * (n + 1)) + c];
        M[(pivot_row * (n + 1)) + c] = temp;
      }
    }

    for (int r = i + 1; r < n; r++) {
      double factor = M[(r * (n + 1)) + i] / M[(i * (n + 1)) + i];
      for (int c = i; c <= n; c++) {
        M[(r * (n + 1)) + c] -= factor * M[(i * (n + 1)) + c];
      }
    }
  }

  for (int i = n - 1; i >= 0; i--) {
    double sum = M[(i * (n + 1)) + n];
    for (int j = i + 1; j < n; j++) {
      sum -= M[(i * (n + 1)) + j] * x[j];
    }
    x[i] = sum / M[(i * (n + 1)) + i];
  }

  free(M);
  return 1;
}

/**
 * @brief Tikhonov-regularized solver for ill-conditioned Hankel matrices.
 */
static int solve_linear_system_regularized(int n, const double *A,
                                           const double *b, double *x,
                                           double lambda) {
  double *AtA = (double *)calloc(n * n, sizeof(double));
  double *Atb = (double *)calloc(n, sizeof(double));
  if (!AtA || !Atb) {
    free(AtA);
    free(Atb);
    return 0;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0.0;
      for (int k = 0; k < n; k++) {
        sum += A[(k * n) + i] * A[(k * n) + j];
      }
      AtA[(i * n) + j] = sum;
    }
    AtA[(i * n) + i] += lambda;
  }

  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int k = 0; k < n; k++) {
      sum += A[(k * n) + i] * b[k];
    }
    Atb[i] = sum;
  }

  int status = solve_linear_system(n, AtA, Atb, x);
  free(AtA);
  free(Atb);
  return status;
}

/* ============================================================================
 * Case 1: Newton's Law of Cooling (ODE: y' = -k*(y - Tm))
 * ============================================================================
 */

static double newton_exact(double x, double y0, double Tm, double k) {
  return Tm + ((y0 - Tm) * exp(-k * x));
}

static double newton_euler(double y0, double Tm, double k, double h,
                           double x_max) {
  int steps = (int)((x_max / h) + 0.5);
  double y = y0;
  for (int i = 0; i < steps; i++) {
    y += h * (-k * (y - Tm));
  }
  return y;
}

static double newton_rk4(double y0, double Tm, double k, double h,
                         double x_max) {
  int steps = (int)((x_max / h) + 0.5);
  double y = y0;
  for (int i = 0; i < steps; i++) {
    double k1 = -k * (y - Tm);
    double k2 = -k * ((y + (0.5 * h * k1)) - Tm);
    double k3 = -k * ((y + (0.5 * h * k2)) - Tm);
    double k4 = -k * ((y + (h * k3)) - Tm);
    y += (h / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);
  }
  return y;
}

static double newton_power_series(int N, double x, double y0, double Tm,
                                  double k) {
  double *a = (double *)malloc((N + 1) * sizeof(double));
  if (!a) {
    return 0.0;
  }
  a[0] = y0;
  if (N >= 1) {
    a[1] = -k * (a[0] - Tm);
  }
  for (int n = 1; n < N; n++) {
    a[n + 1] = -k * a[n] / (n + 1);
  }

  double sum = a[0];
  double term = 1.0;
  for (int n = 1; n <= N; n++) {
    term *= x;
    sum += a[n] * term;
  }
  free(a);
  return sum;
}

static void newton_power_series_trajectory(int N, int num_pts, double x_max,
                                           double y0, double Tm, double k,
                                           double *out) {
  double *a = (double *)malloc((N + 1) * sizeof(double));
  if (!a) {
    return;
  }
  a[0] = y0;
  if (N >= 1) {
    a[1] = -k * (a[0] - Tm);
  }
  for (int n = 1; n < N; n++) {
    a[n + 1] = -k * a[n] / (n + 1);
  }

  double dx = x_max / (num_pts - 1);
  for (int i = 0; i < num_pts; i++) {
    double xi = i * dx;
    double sum = a[0];
    double term = 1.0;
    for (int n = 1; n <= N; n++) {
      term *= xi;
      sum += a[n] * term;
    }
    out[i] = sum;
  }
  free(a);
}

/* ============================================================================
 * Case 2: Airy Equation (ODE: y'' - x*y = 0, y(0)=y0, y'(0)=v0)
 * ============================================================================
 */

static double airy_euler(double y0, double v0, double h, double x_max) {
  int steps = (int)((x_max / h) + 0.5);
  double y = y0;
  double v = v0;
  double x = 0.0;
  for (int i = 0; i < steps; i++) {
    double dy = v;
    double dv = x * y;
    y += h * dy;
    v += h * dv;
    x += h;
  }
  return y;
}

static double airy_rk4(double y0, double v0, double h, double x_max) {
  int steps = (int)((x_max / h) + 0.5);
  double y = y0;
  double v = v0;
  double x = 0.0;
  for (int i = 0; i < steps; i++) {
    double k1_y = v;
    double k1_v = x * y;

    double k2_y = v + (0.5 * h * k1_v);
    double k2_v = (x + (0.5 * h)) * (y + (0.5 * h * k1_y));

    double k3_y = v + (0.5 * h * k2_v);
    double k3_v = (x + (0.5 * h)) * (y + (0.5 * h * k2_y));

    double k4_y = v + (h * k3_v);
    double k4_v = (x + h) * (y + (h * k3_y));

    y += (h / 6.0) * (k1_y + (2.0 * k2_y) + (2.0 * k3_y) + k4_y);
    v += (h / 6.0) * (k1_v + (2.0 * k2_v) + (2.0 * k3_v) + k4_v);
    x += h;
  }
  return y;
}

static double airy_power_series(int N, double x, double y0, double v0) {
  double *a = (double *)malloc((N + 1) * sizeof(double));
  if (!a) {
    return 0.0;
  }
  a[0] = y0;
  if (N >= 1) {
    a[1] = v0;
  }
  if (N >= 2) {
    a[2] = 0.0;
  }

  for (int n = 1; n <= N - 2; n++) {
    a[n + 2] = a[n - 1] / ((n + 1) * (n + 2));
  }

  double sum = a[0];
  double term = 1.0;
  for (int n = 1; n <= N; n++) {
    term *= x;
    sum += a[n] * term;
  }
  free(a);
  return sum;
}

static void airy_power_series_trajectory(int N, int num_pts, double x_max,
                                         double y0, double v0, double *out) {
  double *a = (double *)malloc((N + 1) * sizeof(double));
  if (!a) {
    return;
  }
  a[0] = y0;
  if (N >= 1) {
    a[1] = v0;
  }
  if (N >= 2) {
    a[2] = 0.0;
  }

  for (int n = 1; n <= N - 2; n++) {
    a[n + 2] = a[n - 1] / ((n + 1) * (n + 2));
  }

  double dx = x_max / (num_pts - 1);
  for (int i = 0; i < num_pts; i++) {
    double xi = i * dx;
    double sum = a[0];
    double term = 1.0;
    for (int n = 1; n <= N; n++) {
      term *= xi;
      sum += a[n] * term;
    }
    out[i] = sum;
  }
  free(a);
}

/* ============================================================================
 * Case 3: Logistic Growth Model & Padé Approximants
 * ============================================================================
 */

static double logistic_exact(double x, double y0, double K, double r) {
  double num = K * y0;
  double den = y0 + ((K - y0) * exp(-r * x));
  return num / den;
}

static double logistic_euler(double y0, double K, double r, double h,
                             double x_max) {
  int steps = (int)((x_max / h) + 0.5);
  double y = y0;
  for (int i = 0; i < steps; i++) {
    y += h * (r * y * (1.0 - (y / K)));
  }
  return y;
}

static double logistic_rk4(double y0, double K, double r, double h,
                           double x_max) {
  int steps = (int)((x_max / h) + 0.5);
  double y = y0;
  for (int i = 0; i < steps; i++) {
    double k1 = r * y * (1.0 - (y / K));
    double y1 = y + (0.5 * h * k1);
    double k2 = r * y1 * (1.0 - (y1 / K));
    double y2 = y + (0.5 * h * k2);
    double k3 = r * y2 * (1.0 - (y2 / K));
    double y3 = y + (h * k3);
    double k4 = r * y3 * (1.0 - (y3 / K));
    y += (h / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);
  }
  return y;
}

static void compute_logistic_coefficients(int N, double y0, double K, double r,
                                          double *a) {
  a[0] = y0;
  for (int n = 0; n < N; n++) {
    double cn = 0.0;
    for (int k = 0; k <= n; k++) {
      cn += a[k] * a[n - k];
    }
    a[n + 1] = (r / (n + 1.0)) * (a[n] - (cn / K));
  }
}

static double logistic_power_series(int N, double x, double y0, double K,
                                    double r) {
  double *a = (double *)malloc((N + 1) * sizeof(double));
  if (!a) {
    return 0.0;
  }
  compute_logistic_coefficients(N, y0, K, r, a);
  double sum = a[0];
  double term = 1.0;
  for (int n = 1; n <= N; n++) {
    term *= x;
    sum += a[n] * term;
  }
  free(a);
  return sum;
}

static void logistic_power_series_trajectory(int N, int num_pts, double x_max,
                                             double y0, double K, double r,
                                             double *out) {
  double *a = (double *)malloc((N + 1) * sizeof(double));
  if (!a) {
    return;
  }
  compute_logistic_coefficients(N, y0, K, r, a);

  double dx = x_max / (num_pts - 1);
  for (int i = 0; i < num_pts; i++) {
    double xi = i * dx;
    double sum = a[0];
    double term = 1.0;
    for (int n = 1; n <= N; n++) {
      term *= xi;
      sum += a[n] * term;
    }
    out[i] = sum;
  }
  free(a);
}

static double evaluate_pade(int L, int M, double x, double x0, double y0,
                            double K, double r) {
  int N = L + M;
  double *a = (double *)malloc((N + 1) * sizeof(double));
  if (!a) {
    return 0.0;
  }

  if (x0 == 0.0) {
    compute_logistic_coefficients(N, y0, K, r, a);
  } else {
    double y_x0 = logistic_rk4(y0, K, r, 0.0001, x0);
    compute_logistic_coefficients(N, y_x0, K, r, a);
  }

  double *q = (double *)malloc((M + 1) * sizeof(double));
  if (!q) {
    free(a);
    return 0.0;
  }
  q[0] = 1.0;

  if (M > 0) {
    double *A_mat = (double *)malloc(M * M * sizeof(double));
    double *b_vec = (double *)malloc(M * sizeof(double));
    double *q_sub = (double *)malloc(M * sizeof(double));

    if (A_mat && b_vec && q_sub) {
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
          int idx = L + 1 + i - (j + 1);
          A_mat[(i * M) + j] = (idx >= 0) ? a[idx] : 0.0;
        }
        b_vec[i] = -a[L + 1 + i];
      }

      int ok = solve_linear_system(M, A_mat, b_vec, q_sub);
      if (!ok) {
        ok = solve_linear_system_regularized(M, A_mat, b_vec, q_sub, 1e-10);
      }

      if (!ok) {
        for (int j = 1; j <= M; j++) {
          q[j] = 0.0;
        }
      } else {
        for (int j = 1; j <= M; j++) {
          q[j] = q_sub[j - 1];
        }
      }
    }

    free(A_mat);
    free(b_vec);
    free(q_sub);
  }

  double *p = (double *)malloc((L + 1) * sizeof(double));
  if (!p) {
    free(a);
    free(q);
    return 0.0;
  }

  for (int k = 0; k <= L; k++) {
    double sum = 0.0;
    for (int j = 0; j <= k && j <= M; j++) {
      sum += a[k - j] * q[j];
    }
    p[k] = sum;
  }

  double dx = x - x0;
  double P_val = 0.0;
  double term = 1.0;
  for (int k = 0; k <= L; k++) {
    P_val += p[k] * term;
    term *= dx;
  }

  double Q_val = 0.0;
  term = 1.0;
  for (int j = 0; j <= M; j++) {
    Q_val += q[j] * term;
    term *= dx;
  }

  free(a);
  free(p);
  free(q);

  if (fabs(Q_val) < 1e-14) {
    return 0.0;
  }
  return P_val / Q_val;
}

/* ============================================================================
 * Algorithm 1 Logic for ALL 3 Models
 * ============================================================================
 */

/**
 * @brief Runs Algorithm 1 decision process for Newton's Law of Cooling.
 */
static alg1_result_t run_algorithm_1_newton(double y0, double Tm, double k,
                                            double x0, double x_eval,
                                            int series_n, int verbose) {
  alg1_result_t res;
  memset(&res, 0, sizeof(res));

  res.exact_val = newton_exact(x_eval, y0, Tm, k);
  res.alpha = 0.0;
  res.beta = 0.0;
  res.R_conv = 1e300; /* Infinite convergence radius in C */
  res.is_taylor = 1;
  res.converged_M = 0;
  res.pade_val = newton_power_series(series_n, x_eval, y0, Tm, k);
  res.error = fabs(res.pade_val - res.exact_val);
  snprintf(res.decision_note, sizeof(res.decision_note),
           "Kept Taylor Series N=%d (R_conv = inf >= |x-x0|=%.2f)", series_n,
           fabs(x_eval - x0));

  if (verbose) {
    printf("\n=== Running Algorithm 1 for Newton's Cooling ===\n");
    printf("Target: x_eval = %.2f, Expansion Point x0 = %.2f\n", x_eval, x0);
    printf("Step 1: Singularity structure: No finite poles in complex plane "
           "(xp = inf)\n");
    printf("Step 2: Convergence radius R_conv = infinity\n");
    printf("Step 3: |x_eval - x0| = %.2f <= R_conv (inf). Taylor series is "
           "valid everywhere!\n",
           fabs(x_eval - x0));
    printf("Decision: Stayed with Power Series N=%d (Value = %.6f, Abs Error = "
           "%e)\n",
           series_n, res.pade_val, res.error);
  }
  return res;
}

/**
 * @brief Runs Algorithm 1 decision process for Airy Equation.
 */
static alg1_result_t run_algorithm_1_airy(double y0, double v0, double x0,
                                          double x_eval, int series_n,
                                          int verbose) {
  alg1_result_t res;
  memset(&res, 0, sizeof(res));

  res.exact_val = 534.854243; /* Reference value at x=5.0 */
  res.alpha = 0.0;
  res.beta = 0.0;
  res.R_conv = 1e300; /* Infinite convergence radius in C */
  res.is_taylor = 1;
  res.converged_M = 0;
  res.pade_val = airy_power_series(series_n, x_eval, y0, v0);
  res.error = fabs(res.pade_val - res.exact_val);
  snprintf(res.decision_note, sizeof(res.decision_note),
           "Kept Taylor Series N=%d (R_conv = inf >= |x-x0|=%.2f)", series_n,
           fabs(x_eval - x0));

  if (verbose) {
    printf("\n=== Running Algorithm 1 for Airy Equation ===\n");
    printf("Target: x_eval = %.2f, Expansion Point x0 = %.2f\n", x_eval, x0);
    printf("Step 1: Singularity structure: Entire function, no finite poles "
           "(xp = inf)\n");
    printf("Step 2: Convergence radius R_conv = infinity\n");
    printf("Step 3: |x_eval - x0| = %.2f <= R_conv (inf). Taylor series is "
           "valid everywhere!\n",
           fabs(x_eval - x0));
    printf("Decision: Stayed with Power Series N=%d (Value = %.6f, Abs Error = "
           "%e)\n",
           series_n, res.pade_val, res.error);
  }
  return res;
}

/**
 * @brief Runs Algorithm 1 decision process for Logistic Growth Model.
 */
static alg1_result_t run_algorithm_1_logistic(double y0, double K, double r,
                                              double x0, double x_eval,
                                              double tol, int verbose) {
  alg1_result_t res;
  memset(&res, 0, sizeof(res));

  res.exact_val = logistic_exact(x_eval, y0, K, r);

  double ratio = (K - y0) / y0;
  res.alpha = (1.0 / r) * log(ratio);
  res.beta = M_PI / r;
  res.R_conv =
      sqrt(((res.alpha - x0) * (res.alpha - x0)) + (res.beta * res.beta));

  if (verbose) {
    printf("\n=== Running Algorithm 1 for Logistic Growth ===\n");
    printf("Target: x_eval = %.2f, Expansion Point x0 = %.2f, Tol = %e\n",
           x_eval, x0, tol);
    printf("Step 1: Singularities at xp = %.4f +/- %.4fi\n", res.alpha,
           res.beta);
    printf("Step 2: Convergence radius R_conv = %.4f\n", res.R_conv);
  }

  double dx = fabs(x_eval - x0);
  if (dx <= res.R_conv) {
    res.is_taylor = 1;
    res.converged_M = 0;
    snprintf(res.decision_note, sizeof(res.decision_note),
             "Kept Taylor Series (|x-x0|=%.2f <= R_conv=%.4f)", dx, res.R_conv);
    if (verbose) {
      printf("Step 3: |x_eval - x0| = %.2f <= R_conv (%.4f). Taylor series is "
             "valid.\n",
             dx, res.R_conv);
    }
    double *a = (double *)malloc(21 * sizeof(double));
    if (a) {
      compute_logistic_coefficients(20, y0, K, r, a);
      double sum = a[0];
      double term = 1.0;
      for (int n = 1; n <= 20; n++) {
        term *= (x_eval - x0);
        sum += a[n] * term;
      }
      free(a);
      res.pade_val = sum;
    }
    res.error = fabs(res.pade_val - res.exact_val);
    return res;
  }

  res.is_taylor = 0;
  if (verbose) {
    printf("Step 3: |x_eval - x0| = %.2f > R_conv (%.4f) -> Taylor series "
           "diverges! Switching to Padé.\n",
           dx, res.R_conv);
  }

  int np = 1;
  res.M_min = 2 * np;
  if (verbose) {
    printf("Step 4: Conjugate pole pairs np = %d => M_min = %d\n", np,
           res.M_min);
  }

  double numerator =
      (x_eval * x_eval) - (res.alpha * res.alpha) - (res.beta * res.beta);
  double denominator = 2.0 * (x_eval - res.alpha);
  res.x0_threshold = numerator / denominator;
  if (verbose) {
    printf("Step 5: Expansion shift condition requires x0 > %.4f\n",
           res.x0_threshold);
  }

  int M = res.M_min;
  double prev_val = evaluate_pade(M, M, x_eval, x0, y0, K, r);
  double curr_val = prev_val;
  if (verbose) {
    printf("  Iter 0: Padé [%d/%d] = %.6f\n", M, M, prev_val);
  }

  FILE *f_alg = fopen("data/algorithm1_results.csv", "w");
  if (f_alg) {
    fprintf(f_alg, "Step,Description,Value\n");
    fprintf(f_alg, "1,Singularity_alpha,%.6f\n", res.alpha);
    fprintf(f_alg, "1,Singularity_beta,%.6f\n", res.beta);
    fprintf(f_alg, "2,R_conv,%.6f\n", res.R_conv);
    fprintf(f_alg, "4,M_min,%d\n", res.M_min);
    fprintf(f_alg, "5,x0_threshold,%.6f\n", res.x0_threshold);
    fprintf(f_alg, "6,Pade_%d_%d,%.6f\n", M, M, prev_val);
  }

  while (M < 10) {
    M++;
    curr_val = evaluate_pade(M, M, x_eval, x0, y0, K, r);
    double diff = fabs(curr_val - prev_val);
    if (verbose) {
      printf("  Iter %d: Padé [%d/%d] = %.6f (Diff = %e)\n", M - res.M_min, M,
             M, curr_val, diff);
    }
    if (f_alg) {
      fprintf(f_alg, "6,Pade_%d_%d,%.6f\n", M, M, curr_val);
    }
    if (diff < tol) {
      res.converged_M = M;
      res.pade_val = curr_val;
      res.error = fabs(res.pade_val - res.exact_val);
      snprintf(res.decision_note, sizeof(res.decision_note),
               "Switched to Padé [%d/%d] (|x-x0|=%.2f > R_conv=%.4f)", M, M, dx,
               res.R_conv);
      if (verbose) {
        printf("Algorithm 1 Converged at Padé [%d/%d] with value = %.6f (Exact "
               "= %.6f)\n",
               M, M, curr_val, res.exact_val);
      }
      if (f_alg) {
        fclose(f_alg);
      }
      return res;
    }
    prev_val = curr_val;
  }

  if (f_alg) {
    fclose(f_alg);
  }
  res.converged_M = M;
  res.pade_val = curr_val;
  res.error = fabs(res.pade_val - res.exact_val);
  snprintf(res.decision_note, sizeof(res.decision_note),
           "Switched to Padé [%d/%d] (|x-x0|=%.2f > R_conv=%.4f)", M, M, dx,
           res.R_conv);
  return res;
}

/* ============================================================================
 * Part Runners
 * ============================================================================
 */

static void run_part_1_accuracy(const benchmark_config_t *cfg) {
  printf("--- Part 1: Accuracy at x = %.2f ---\n", cfg->x_eval);

  double exact_newton =
      newton_exact(cfg->x_eval, cfg->newton_y0, cfg->newton_tm, cfg->newton_k);
  double exact_logistic = logistic_exact(cfg->x_eval, cfg->logistic_y0,
                                         cfg->logistic_k, cfg->logistic_r);
  double ref_airy = 534.854243;

  /* Run Algorithm 1 for all 3 models */
  alg1_result_t alg1_n =
      run_algorithm_1_newton(cfg->newton_y0, cfg->newton_tm, cfg->newton_k,
                             cfg->x0, cfg->x_eval, cfg->series_n, 0);

  alg1_result_t alg1_a = run_algorithm_1_airy(
      cfg->airy_y0, cfg->airy_v0, cfg->x0, cfg->x_eval, cfg->series_n, 0);

  alg1_result_t alg1_l = run_algorithm_1_logistic(
      cfg->logistic_y0, cfg->logistic_k, cfg->logistic_r, cfg->x0, cfg->x_eval,
      cfg->tol, 0);

  /* Newton's Cooling */
  printf("Newton's Cooling (Exact: %.6f)\n", exact_newton);
  double val_n_e1 = newton_euler(cfg->newton_y0, cfg->newton_tm, cfg->newton_k,
                                 cfg->step_h1, cfg->x_eval);
  double val_n_e2 = newton_euler(cfg->newton_y0, cfg->newton_tm, cfg->newton_k,
                                 cfg->step_h2, cfg->x_eval);
  double val_n_r1 = newton_rk4(cfg->newton_y0, cfg->newton_tm, cfg->newton_k,
                               cfg->step_h1, cfg->x_eval);
  double val_n_r2 = newton_rk4(cfg->newton_y0, cfg->newton_tm, cfg->newton_k,
                               cfg->step_h2, cfg->x_eval);
  double val_n_s5 = newton_power_series(5, cfg->x_eval, cfg->newton_y0,
                                        cfg->newton_tm, cfg->newton_k);
  double val_n_s10 = newton_power_series(10, cfg->x_eval, cfg->newton_y0,
                                         cfg->newton_tm, cfg->newton_k);
  double val_n_sn =
      newton_power_series(cfg->series_n, cfg->x_eval, cfg->newton_y0,
                          cfg->newton_tm, cfg->newton_k);

  printf("  Euler h=%.2f : Value = %.6f, Abs Error = %e\n", cfg->step_h1,
         val_n_e1, fabs(val_n_e1 - exact_newton));
  printf("  Euler h=%.2f: Value = %.6f, Abs Error = %e\n", cfg->step_h2,
         val_n_e2, fabs(val_n_e2 - exact_newton));
  printf("  RK4 h=%.2f   : Value = %.6f, Abs Error = %e\n", cfg->step_h1,
         val_n_r1, fabs(val_n_r1 - exact_newton));
  printf("  RK4 h=%.2f  : Value = %.6f, Abs Error = %e\n", cfg->step_h2,
         val_n_r2, fabs(val_n_r2 - exact_newton));
  printf("  Series N=5  : Value = %.6f, Abs Error = %e\n", val_n_s5,
         fabs(val_n_s5 - exact_newton));
  printf("  Series N=10 : Value = %.6f, Abs Error = %e\n", val_n_s10,
         fabs(val_n_s10 - exact_newton));
  printf("  Series N=%d : Value = %.6f, Abs Error = %e\n", cfg->series_n,
         val_n_sn, fabs(val_n_sn - exact_newton));
  printf("  Algorithm 1 : Value = %.6f, Abs Error = %e [%s]\n\n",
         alg1_n.pade_val, alg1_n.error, alg1_n.decision_note);

  /* Airy Equation */
  printf("Airy Equation (Reference N=100: %.6f)\n", ref_airy);
  double val_a_e1 =
      airy_euler(cfg->airy_y0, cfg->airy_v0, cfg->step_h1, cfg->x_eval);
  double val_a_e2 =
      airy_euler(cfg->airy_y0, cfg->airy_v0, cfg->step_h2, cfg->x_eval);
  double val_a_r1 =
      airy_rk4(cfg->airy_y0, cfg->airy_v0, cfg->step_h1, cfg->x_eval);
  double val_a_r2 =
      airy_rk4(cfg->airy_y0, cfg->airy_v0, cfg->step_h2, cfg->x_eval);
  double val_a_s5 =
      airy_power_series(5, cfg->x_eval, cfg->airy_y0, cfg->airy_v0);
  double val_a_s10 =
      airy_power_series(10, cfg->x_eval, cfg->airy_y0, cfg->airy_v0);
  double val_a_sn =
      airy_power_series(cfg->series_n, cfg->x_eval, cfg->airy_y0, cfg->airy_v0);

  printf("  Euler h=%.2f : Value = %.6f, Abs Error = %e\n", cfg->step_h1,
         val_a_e1, fabs(val_a_e1 - ref_airy));
  printf("  Euler h=%.2f: Value = %.6f, Abs Error = %e\n", cfg->step_h2,
         val_a_e2, fabs(val_a_e2 - ref_airy));
  printf("  RK4 h=%.2f   : Value = %.6f, Abs Error = %e\n", cfg->step_h1,
         val_a_r1, fabs(val_a_r1 - ref_airy));
  printf("  RK4 h=%.2f  : Value = %.6f, Abs Error = %e\n", cfg->step_h2,
         val_a_r2, fabs(val_a_r2 - ref_airy));
  printf("  Series N=5  : Value = %.6f, Abs Error = %e\n", val_a_s5,
         fabs(val_a_s5 - ref_airy));
  printf("  Series N=10 : Value = %.6f, Abs Error = %e\n", val_a_s10,
         fabs(val_a_s10 - ref_airy));
  printf("  Series N=%d : Value = %.6f, Abs Error = %e\n", cfg->series_n,
         val_a_sn, fabs(val_a_sn - ref_airy));
  printf("  Algorithm 1 : Value = %.6f, Abs Error = %e [%s]\n\n",
         alg1_a.pade_val, alg1_a.error, alg1_a.decision_note);

  /* Logistic Growth */
  printf("Logistic Growth (Exact: %.6f)\n", exact_logistic);
  double val_l_e1 = logistic_euler(cfg->logistic_y0, cfg->logistic_k,
                                   cfg->logistic_r, cfg->step_h1, cfg->x_eval);
  double val_l_e2 = logistic_euler(cfg->logistic_y0, cfg->logistic_k,
                                   cfg->logistic_r, cfg->step_h2, cfg->x_eval);
  double val_l_r1 = logistic_rk4(cfg->logistic_y0, cfg->logistic_k,
                                 cfg->logistic_r, cfg->step_h1, cfg->x_eval);
  double val_l_r2 = logistic_rk4(cfg->logistic_y0, cfg->logistic_k,
                                 cfg->logistic_r, cfg->step_h2, cfg->x_eval);
  double val_l_s5 = logistic_power_series(5, cfg->x_eval, cfg->logistic_y0,
                                          cfg->logistic_k, cfg->logistic_r);
  double val_l_s10 = logistic_power_series(10, cfg->x_eval, cfg->logistic_y0,
                                           cfg->logistic_k, cfg->logistic_r);
  double val_l_sn =
      logistic_power_series(cfg->series_n, cfg->x_eval, cfg->logistic_y0,
                            cfg->logistic_k, cfg->logistic_r);

  printf("  Euler h=%.2f : Value = %.6f, Abs Error = %e\n", cfg->step_h1,
         val_l_e1, fabs(val_l_e1 - exact_logistic));
  printf("  Euler h=%.2f: Value = %.6f, Abs Error = %e\n", cfg->step_h2,
         val_l_e2, fabs(val_l_e2 - exact_logistic));
  printf("  RK4 h=%.2f   : Value = %.6f, Abs Error = %e\n", cfg->step_h1,
         val_l_r1, fabs(val_l_r1 - exact_logistic));
  printf("  RK4 h=%.2f  : Value = %.6f, Abs Error = %e\n", cfg->step_h2,
         val_l_r2, fabs(val_l_r2 - exact_logistic));
  printf("  Series N=5  : Value = %.6f, Abs Error = %e\n", val_l_s5,
         fabs(val_l_s5 - exact_logistic));
  printf("  Series N=10 : Value = %.6f, Abs Error = %e\n", val_l_s10,
         fabs(val_l_s10 - exact_logistic));
  printf("  Series N=%d : Value = %.6f, Abs Error = %e\n", cfg->series_n,
         val_l_sn, fabs(val_l_sn - exact_logistic));
  printf("  Algorithm 1 : Value = %.6f, Abs Error = %e [%s]\n", alg1_l.pade_val,
         alg1_l.error, alg1_l.decision_note);

  /* Save Accuracy Results to CSV */
  ensure_data_dir();
  FILE *f_acc = fopen("data/accuracy_results.csv", "w");
  if (f_acc) {
    fprintf(f_acc, "Model,Method,Config,Value,AbsError,DecisionNote\n");
    fprintf(f_acc, "Newton,Euler,h=%.2f,%.6f,%e,N/A\n", cfg->step_h1, val_n_e1,
            fabs(val_n_e1 - exact_newton));
    fprintf(f_acc, "Newton,Euler,h=%.2f,%.6f,%e,N/A\n", cfg->step_h2, val_n_e2,
            fabs(val_n_e2 - exact_newton));
    fprintf(f_acc, "Newton,RK4,h=%.2f,%.6f,%e,N/A\n", cfg->step_h1, val_n_r1,
            fabs(val_n_r1 - exact_newton));
    fprintf(f_acc, "Newton,RK4,h=%.2f,%.6f,%e,N/A\n", cfg->step_h2, val_n_r2,
            fabs(val_n_r2 - exact_newton));
    fprintf(f_acc, "Newton,Series,N=%d,%.6f,%e,N/A\n", cfg->series_n, val_n_sn,
            fabs(val_n_sn - exact_newton));
    fprintf(f_acc, "Newton,Algorithm1,Auto,%.6f,%e,\"%s\"\n", alg1_n.pade_val,
            alg1_n.error, alg1_n.decision_note);

    fprintf(f_acc, "Airy,Euler,h=%.2f,%.6f,%e,N/A\n", cfg->step_h2, val_a_e2,
            fabs(val_a_e2 - ref_airy));
    fprintf(f_acc, "Airy,RK4,h=%.2f,%.6f,%e,N/A\n", cfg->step_h2, val_a_r2,
            fabs(val_a_r2 - ref_airy));
    fprintf(f_acc, "Airy,Series,N=%d,%.6f,%e,N/A\n", cfg->series_n, val_a_sn,
            fabs(val_a_sn - ref_airy));
    fprintf(f_acc, "Airy,Algorithm1,Auto,%.6f,%e,\"%s\"\n", alg1_a.pade_val,
            alg1_a.error, alg1_a.decision_note);

    fprintf(f_acc, "Logistic,Euler,h=%.2f,%.6f,%e,N/A\n", cfg->step_h2,
            val_l_e2, fabs(val_l_e2 - exact_logistic));
    fprintf(f_acc, "Logistic,RK4,h=%.2f,%.6f,%e,N/A\n", cfg->step_h2, val_l_r2,
            fabs(val_l_r2 - exact_logistic));
    fprintf(f_acc, "Logistic,Series,N=%d,%.6f,%e,N/A\n", cfg->series_n,
            val_l_sn, fabs(val_l_sn - exact_logistic));
    fprintf(f_acc, "Logistic,Algorithm1,Tol=%e,%.6f,%e,\"%s\"\n", cfg->tol,
            alg1_l.pade_val, alg1_l.error, alg1_l.decision_note);
    fclose(f_acc);
  }
}

static void run_part_2_pade(const benchmark_config_t *cfg) {
  printf("\n--- Part 2: Padé Approximants for Logistic Growth at x=%.2f ---\n",
         cfg->x_eval);
  double exact_logistic = logistic_exact(cfg->x_eval, cfg->logistic_y0,
                                         cfg->logistic_k, cfg->logistic_r);

  int orders[][2] = {{1, 1}, {2, 1}, {2, 2}, {3, 2}, {4, 3},
                     {4, 4}, {5, 5}, {6, 6}, {8, 8}};
  int num_orders = sizeof(orders) / sizeof(orders[0]);

  ensure_data_dir();
  FILE *f_pade = fopen("data/pade_order_results.csv", "w");
  if (f_pade) {
    fprintf(f_pade, "L,M,Value,AbsError\n");
  }

  for (int i = 0; i < num_orders; i++) {
    int L = orders[i][0];
    int M = orders[i][1];
    double val = evaluate_pade(L, M, cfg->x_eval, cfg->x0, cfg->logistic_y0,
                               cfg->logistic_k, cfg->logistic_r);
    double err = fabs(val - exact_logistic);
    printf("  Padé [%d/%d]: Value = %.6f, Abs Error = %e\n", L, M, val, err);
    if (f_pade) {
      fprintf(f_pade, "%d,%d,%.6f,%e\n", L, M, val, err);
    }
  }
  if (f_pade) {
    fclose(f_pade);
  }

  /* Padé Matrix L x M */
  FILE *f_matrix = fopen("data/pade_matrix_results.csv", "w");
  if (f_matrix) {
    fprintf(f_matrix, "L/M,1,2,3,4,5,6,7\n");
    for (int L = 1; L <= 7; L++) {
      fprintf(f_matrix, "%d", L);
      for (int M = 1; M <= 7; M++) {
        double val = evaluate_pade(L, M, cfg->x_eval, cfg->x0, cfg->logistic_y0,
                                   cfg->logistic_k, cfg->logistic_r);
        double err = fabs(val - exact_logistic);
        fprintf(f_matrix, ",%e", err);
      }
      fprintf(f_matrix, "\n");
    }
    fclose(f_matrix);
  }

  /* Expansion Point Shift Verification */
  printf("\nExpansion point shift verification (Table 5):\n");
  double shifts[] = {0.0, 1.0, 2.0, 2.5, 3.0};
  int num_shifts = sizeof(shifts) / sizeof(shifts[0]);

  FILE *f_shift = fopen("data/shift_results.csv", "w");
  if (f_shift) {
    fprintf(f_shift, "x0,y_x0,Pade33_Val,AbsError\n");
  }

  for (int i = 0; i < num_shifts; i++) {
    double x0 = shifts[i];
    double y_x0 = (x0 == 0.0) ? cfg->logistic_y0
                              : logistic_rk4(cfg->logistic_y0, cfg->logistic_k,
                                             cfg->logistic_r, 0.0001, x0);
    double val = evaluate_pade(3, 3, cfg->x_eval, x0, cfg->logistic_y0,
                               cfg->logistic_k, cfg->logistic_r);
    double err = fabs(val - exact_logistic);
    printf("  x0 = %.1f: y(x0) = %.6f, [3/3](%.1f) = %.6f, Abs Error = %e\n",
           x0, y_x0, cfg->x_eval, val, err);
    if (f_shift) {
      fprintf(f_shift, "%.1f,%.6f,%.6f,%e\n", x0, y_x0, val, err);
    }
  }
  if (f_shift) {
    fclose(f_shift);
  }
}

static void run_part_3_timing(const benchmark_config_t *cfg) {
  printf("\n--- Part 3: Timing Benchmarks (Average of %d runs) ---\n",
         cfg->runs);
  benchmark_timer_t timer;
  double t_newton_euler;
  double t_newton_rk4;
  double t_newton_series;
  double t_newton_alg1;
  double t_airy_euler;
  double t_airy_rk4;
  double t_airy_series;
  double t_airy_alg1;
  double t_log_euler;
  double t_log_rk4;
  double t_log_series;
  double t_log_alg1;

  double *trajectory = (double *)malloc(SPARSE_PTS * sizeof(double));
  if (!trajectory) {
    return;
  }
  double sum_val = 0.0;

  /* Newton's Cooling */
  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    sum_val += newton_euler(cfg->newton_y0, cfg->newton_tm, cfg->newton_k,
                            0.001, cfg->x_eval);
  }
  timer_end(&timer);
  t_newton_euler = timer_elapsed_us(&timer) / cfg->runs;

  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    sum_val += newton_rk4(cfg->newton_y0, cfg->newton_tm, cfg->newton_k, 0.001,
                          cfg->x_eval);
  }
  timer_end(&timer);
  t_newton_rk4 = timer_elapsed_us(&timer) / cfg->runs;

  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    newton_power_series_trajectory(cfg->series_n, SPARSE_PTS, cfg->x_eval,
                                   cfg->newton_y0, cfg->newton_tm,
                                   cfg->newton_k, trajectory);
    sum_val += trajectory[SPARSE_PTS - 1];
  }
  timer_end(&timer);
  t_newton_series = timer_elapsed_us(&timer) / cfg->runs;

  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    alg1_result_t res =
        run_algorithm_1_newton(cfg->newton_y0, cfg->newton_tm, cfg->newton_k,
                               cfg->x0, cfg->x_eval, cfg->series_n, 0);
    sum_val += res.pade_val;
  }
  timer_end(&timer);
  t_newton_alg1 = timer_elapsed_us(&timer) / cfg->runs;

  printf("Newton's Cooling:\n  Euler: %.3f μs, RK4: %.3f μs, Series: %.3f μs, "
         "Algorithm 1: %.3f μs\n",
         t_newton_euler, t_newton_rk4, t_newton_series, t_newton_alg1);

  /* Airy Equation */
  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    sum_val += airy_euler(cfg->airy_y0, cfg->airy_v0, 0.001, cfg->x_eval);
  }
  timer_end(&timer);
  t_airy_euler = timer_elapsed_us(&timer) / cfg->runs;

  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    sum_val += airy_rk4(cfg->airy_y0, cfg->airy_v0, 0.001, cfg->x_eval);
  }
  timer_end(&timer);
  t_airy_rk4 = timer_elapsed_us(&timer) / cfg->runs;

  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    airy_power_series_trajectory(cfg->series_n, SPARSE_PTS, cfg->x_eval,
                                 cfg->airy_y0, cfg->airy_v0, trajectory);
    sum_val += trajectory[SPARSE_PTS - 1];
  }
  timer_end(&timer);
  t_airy_series = timer_elapsed_us(&timer) / cfg->runs;

  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    alg1_result_t res = run_algorithm_1_airy(
        cfg->airy_y0, cfg->airy_v0, cfg->x0, cfg->x_eval, cfg->series_n, 0);
    sum_val += res.pade_val;
  }
  timer_end(&timer);
  t_airy_alg1 = timer_elapsed_us(&timer) / cfg->runs;

  printf("Airy Equation:\n  Euler: %.3f μs, RK4: %.3f μs, Series: %.3f μs, "
         "Algorithm 1: %.3f μs\n",
         t_airy_euler, t_airy_rk4, t_airy_series, t_airy_alg1);

  /* Logistic Growth */
  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    sum_val += logistic_euler(cfg->logistic_y0, cfg->logistic_k,
                              cfg->logistic_r, 0.001, cfg->x_eval);
  }
  timer_end(&timer);
  t_log_euler = timer_elapsed_us(&timer) / cfg->runs;

  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    sum_val += logistic_rk4(cfg->logistic_y0, cfg->logistic_k, cfg->logistic_r,
                            0.001, cfg->x_eval);
  }
  timer_end(&timer);
  t_log_rk4 = timer_elapsed_us(&timer) / cfg->runs;

  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    logistic_power_series_trajectory(cfg->series_n, SPARSE_PTS, cfg->x_eval,
                                     cfg->logistic_y0, cfg->logistic_k,
                                     cfg->logistic_r, trajectory);
    sum_val += trajectory[SPARSE_PTS - 1];
  }
  timer_end(&timer);
  t_log_series = timer_elapsed_us(&timer) / cfg->runs;

  timer_start(&timer);
  for (int r_idx = 0; r_idx < cfg->runs; r_idx++) {
    alg1_result_t res = run_algorithm_1_logistic(
        cfg->logistic_y0, cfg->logistic_k, cfg->logistic_r, cfg->x0,
        cfg->x_eval, cfg->tol, 0);
    sum_val += res.pade_val;
  }
  timer_end(&timer);
  t_log_alg1 = timer_elapsed_us(&timer) / cfg->runs;

  printf("Logistic Growth:\n  Euler: %.3f μs, RK4: %.3f μs, Series: %.3f μs, "
         "Algorithm 1: %.3f μs\n",
         t_log_euler, t_log_rk4, t_log_series, t_log_alg1);

  /* Save timing to CSV */
  ensure_data_dir();
  FILE *f_time = fopen("data/timing_results.csv", "w");
  if (f_time) {
    fprintf(f_time, "Model,Euler_μs,RK4_μs,PowerSeries_μs,Algorithm1_μs\n");
    fprintf(f_time, "Newton,%.6f,%.6f,%.6f,%.6f\n", t_newton_euler,
            t_newton_rk4, t_newton_series, t_newton_alg1);
    fprintf(f_time, "Airy,%.6f,%.6f,%.6f,%.6f\n", t_airy_euler, t_airy_rk4,
            t_airy_series, t_airy_alg1);
    fprintf(f_time, "Logistic,%.6f,%.6f,%.6f,%.6f\n", t_log_euler, t_log_rk4,
            t_log_series, t_log_alg1);
    fclose(f_time);
  }

  printf("  (Benchmark validation checksum: %.6f)\n", sum_val);
  free(trajectory);
}

static void run_part_4_alg1_trace(const benchmark_config_t *cfg) {
  printf("\n==================================================================="
         "=====\n");
  printf(
      "  PART 4: ALGORITHM 1 DECISION & EXECUTION TRACE (ALL 3 ODE MODELS)\n");
  printf("====================================================================="
         "===\n");

  run_algorithm_1_newton(cfg->newton_y0, cfg->newton_tm, cfg->newton_k, cfg->x0,
                         cfg->x_eval, cfg->series_n, 1);

  run_algorithm_1_airy(cfg->airy_y0, cfg->airy_v0, cfg->x0, cfg->x_eval,
                       cfg->series_n, 1);

  run_algorithm_1_logistic(cfg->logistic_y0, cfg->logistic_k, cfg->logistic_r,
                           cfg->x0, cfg->x_eval, cfg->tol, 1);
}

/* ============================================================================
 * Main Driver
 * ============================================================================
 */

int main(int argc, char **argv) {
  benchmark_config_t cfg;
  parse_cli_args(argc, argv, &cfg);

  printf("=== ODE SOLVER ACCURACY AND TIMING BENCHMARKS ===\n\n");

  if (cfg.target_part == 0 || cfg.target_part == 1) {
    run_part_1_accuracy(&cfg);
  }
  if (cfg.target_part == 0 || cfg.target_part == 2) {
    run_part_2_pade(&cfg);
  }
  if (cfg.target_part == 0 || cfg.target_part == 3) {
    run_part_3_timing(&cfg);
  }
  if (cfg.target_part == 0 || cfg.target_part == 4) {
    run_part_4_alg1_trace(&cfg);
  }

  return 0;
}
