/**
 * MONTE CARLO UNCERTAINTY PROPAGATION TEST
 * 
 * Tests sensor-based initialization with noise injection.
 * Analytical solution computed from EXACT initial conditions (not reconstructed).
 * 
 * Usage: ./mc_test [num_trials] [noise_level] [num_cells]
 *        ./mc_test 100 0.05 500
 * 
 * Compile with OpenMP:
 *   g++ -std=c++17 -O3 -fopenmp -o mc_test mc_test.cpp
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================================
// CONFIGURATION
// ============================================================================

enum class ToroTest {
  TEST1_SOD,
  TEST2_123,
  TEST3_BLAST_LEFT,
  TEST4_SLOW_SHOCK,
  TEST5_COLLISION
};

// Select test
const ToroTest CURRENT_TEST = ToroTest::TEST2_123;

struct TestParams {
  double rho_L, u_L, p_L;
  double rho_R, u_R, p_R;
  double x_diaphragm;
  double endTime;
  double gamma;
  double CFL;
  std::string name;
};

TestParams getTestParams(ToroTest test) {
  TestParams params;
  params.gamma = 1.4;
  params.CFL = 0.5;
  params.x_diaphragm = 0.5;

  switch(test) {
    case ToroTest::TEST1_SOD:
      params.rho_L = 1.0; params.u_L = 0.0; params.p_L = 1.0;
      params.rho_R = 0.125; params.u_R = 0.0; params.p_R = 0.1;
      params.endTime = 0.25;
      params.name = "Test 1: Sod Shock Tube";
      break;

    case ToroTest::TEST2_123:
      params.rho_L = 1.0; params.u_L = -2.0; params.p_L = 0.4;
      params.rho_R = 1.0; params.u_R = 2.0; params.p_R = 0.4;
      params.endTime = 0.15;
      params.name = "Test 2: 123 Problem";
      break;
    
    case ToroTest::TEST3_BLAST_LEFT:
      params.rho_L = 1.0; params.u_L = 0.0; params.p_L = 1000.0;
      params.rho_R = 1.0; params.u_R = 0.0; params.p_R = 0.01;
      params.endTime = 0.012;
      params.name = "Test 3: Left Blast Wave";
      break;

    case ToroTest::TEST4_SLOW_SHOCK:
      params.rho_L = 1.0; params.u_L = 0.0; params.p_L = 0.01;
      params.rho_R = 1.0; params.u_R = 0.0; params.p_R = 100.0;
      params.endTime = 0.035;
      params.name = "Test 4: Slow Shock";
      break;

    case ToroTest::TEST5_COLLISION:
      params.rho_L = 5.99924; params.u_L = 19.5975; params.p_L = 460.894;
      params.rho_R = 5.99242; params.u_R = -6.19633; params.p_R = 46.0950;
      params.endTime = 0.035;
      params.name = "Test 5: Collision";
      break;
  }

  return params;
}

const TestParams TEST = getTestParams(CURRENT_TEST);

// ============================================================================
// SENSOR CONFIGURATION
// ============================================================================

struct Sensor {
  double x;
  double rho;
  double u;
  double p;
};

std::vector<Sensor> getExactSensors(const TestParams& test, const std::vector<double>& positions) {
  std::vector<Sensor> sensors;
  for (double x : positions) {
    Sensor s;
    s.x = x;
    if (x < test.x_diaphragm) {
      s.rho = test.rho_L;
      s.u = test.u_L;
      s.p = test.p_L;
    } else {
      s.rho = test.rho_R;
      s.u = test.u_R;
      s.p = test.p_R;
    }
    sensors.push_back(s);
  }
  return sensors;
}

// ============================================================================
// GRID
// ============================================================================

struct Grid {
  int n;
  double length;
  double dx;
  double t;

  std::vector<double> x;
  std::vector<double> rho;
  std::vector<double> u;
  std::vector<double> p;
  std::vector<double> rho_u;
  std::vector<double> E;

  Grid(int numCells, double L) : n(numCells), length(L), dx(L / numCells), t(0.0) {
    x.resize(n);
    rho.resize(n);
    u.resize(n);
    p.resize(n);
    rho_u.resize(n);
    E.resize(n);

    for (int i = 0; i < n; ++i) {
      x[i] = (i + 0.5) * dx;
    }
  }

  void reset() {
    t = 0.0;
    std::fill(rho.begin(), rho.end(), 0.0);
    std::fill(u.begin(), u.end(), 0.0);
    std::fill(p.begin(), p.end(), 0.0);
    std::fill(rho_u.begin(), rho_u.end(), 0.0);
    std::fill(E.begin(), E.end(), 0.0);
  }

  // Initialize from exact Riemann IC (for analytical reference)
  void initializeExact(const TestParams& test) {
    for (int i = 0; i < n; ++i) {
      if (x[i] < test.x_diaphragm) {
        rho[i] = test.rho_L;
        u[i] = test.u_L;
        p[i] = test.p_L;
      } else {
        rho[i] = test.rho_R;
        u[i] = test.u_R;
        p[i] = test.p_R;
      }
    }
    primToConserved();
  }

  // Initialize from sparse sensor data (piecewise constant)
  void initializeFromSensors(const std::vector<Sensor>& sensors) {
    if (sensors.empty()) return;

    // Sort sensors by position
    std::vector<Sensor> sorted = sensors;
    std::sort(sorted.begin(), sorted.end(), [](const Sensor& a, const Sensor& b) {
      return a.x < b.x;
    });

    // Find discontinuity using ALL primitive variables (not just pressure)
    // This handles cases like Test 2 where pressure is uniform but velocity jumps
    double max_jump = 0.0;
    size_t disc_idx = 0;
    for (size_t i = 0; i < sorted.size() - 1; ++i) {
      // Density jump (normalized)
      double rho_avg = 0.5 * (sorted[i].rho + sorted[i + 1].rho);
      double rho_jump = std::abs(sorted[i + 1].rho - sorted[i].rho) / std::max(rho_avg, 1e-10);
      
      // Pressure jump (normalized)
      double p_avg = 0.5 * (sorted[i].p + sorted[i + 1].p);
      double p_jump = std::abs(sorted[i + 1].p - sorted[i].p) / std::max(p_avg, 1e-10);
      
      // Velocity jump (normalized by magnitude + offset to handle zero velocity)
      double u_scale = std::max({std::abs(sorted[i].u), std::abs(sorted[i + 1].u), 0.1});
      double u_jump = std::abs(sorted[i + 1].u - sorted[i].u) / u_scale;
      
      // Combined jump score - weight all variables
      double jump = 0.3 * rho_jump + 0.3 * p_jump + 0.4 * u_jump;
      
      if (jump > max_jump) {
        max_jump = jump;
        disc_idx = i;
      }
    }

    double x_disc = 0.5 * (sorted[disc_idx].x + sorted[disc_idx + 1].x);

    // Average states on each side
    double rho_L = 0, u_L = 0, p_L = 0;
    double rho_R = 0, u_R = 0, p_R = 0;
    int count_L = 0, count_R = 0;

    for (const auto& s : sorted) {
      if (s.x < x_disc) {
        rho_L += s.rho;
        u_L += s.u;
        p_L += s.p;
        count_L++;
      } else {
        rho_R += s.rho;
        u_R += s.u;
        p_R += s.p;
        count_R++;
      }
    }

    if (count_L > 0) { rho_L /= count_L; u_L /= count_L; p_L /= count_L; }
    if (count_R > 0) { rho_R /= count_R; u_R /= count_R; p_R /= count_R; }

    // Apply piecewise constant
    for (int i = 0; i < n; ++i) {
      if (x[i] < x_disc) {
        rho[i] = rho_L;
        u[i] = u_L;
        p[i] = p_L;
      } else {
        rho[i] = rho_R;
        u[i] = u_R;
        p[i] = p_R;
      }
    }

    primToConserved();
  }

  void primToConserved() {
    for (int i = 0; i < n; ++i) {
      rho_u[i] = rho[i] * u[i];
      double e = p[i] / ((TEST.gamma - 1.0) * rho[i]);
      E[i] = rho[i] * (e + 0.5 * u[i] * u[i]);
    }
  }

  void consToPrim() {
    for (int i = 0; i < n; ++i) {
      u[i] = rho_u[i] / rho[i];
      p[i] = (TEST.gamma - 1.0) * (E[i] - 0.5 * rho[i] * u[i] * u[i]);
      p[i] = std::max(p[i], 1e-10);
    }
  }

  std::vector<double> getEntropy() const {
    std::vector<double> s(n);
    for (int i = 0; i < n; ++i) {
      s[i] = std::log(p[i] / std::pow(rho[i], TEST.gamma));
    }
    return s;
  }
};

// ============================================================================
// NUMERICAL SOLVER (HLLC + MUSCL + RK2)
// ============================================================================

double slopeLimit(double v_minus, double v_center, double v_plus) {
  double delta_minus = v_center - v_minus;
  double delta_plus = v_plus - v_center;

  if (delta_plus * delta_minus <= 0.0) return 0.0;

  double delta_center = 0.5 * (delta_plus + delta_minus);
  double sign = (delta_center > 0.0) ? 1.0 : -1.0;
  return sign * std::min({2.0 * std::abs(delta_minus),
                          2.0 * std::abs(delta_plus),
                          std::abs(delta_center)});
}

void solveHLLC(double rho_L, double u_L, double p_L,
               double rho_R, double u_R, double p_R,
               double gamma,
               double& F_rho, double& F_rhou, double& F_E) {
  double a_L = std::sqrt(gamma * p_L / rho_L);
  double a_R = std::sqrt(gamma * p_R / rho_R);

  double e_L = p_L / ((gamma - 1.0) * rho_L);
  double e_R = p_R / ((gamma - 1.0) * rho_R);
  double E_L = rho_L * (e_L + 0.5 * u_L * u_L);
  double E_R = rho_R * (e_R + 0.5 * u_R * u_R);

  double S_L = std::min(u_L - a_L, u_R - a_R);
  double S_R = std::max(u_L + a_L, u_R + a_R);

  double S_star = (p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) /
                  (rho_L * (S_L - u_L) - rho_R * (S_R - u_R));

  double F_L_rho = rho_L * u_L;
  double F_L_rhou = rho_L * u_L * u_L + p_L;
  double F_L_E = u_L * (E_L + p_L);

  double F_R_rho = rho_R * u_R;
  double F_R_rhou = rho_R * u_R * u_R + p_R;
  double F_R_E = u_R * (E_R + p_R);

  if (S_L >= 0.0) {
    F_rho = F_L_rho;
    F_rhou = F_L_rhou;
    F_E = F_L_E;
  } else if (S_R <= 0.0) {
    F_rho = F_R_rho;
    F_rhou = F_R_rhou;
    F_E = F_R_E;
  } else if (S_star >= 0.0) {
    double factor = rho_L * (S_L - u_L) / (S_L - S_star);
    double U_star_rho = factor;
    double U_star_rhou = factor * S_star;
    double U_star_E = factor * (E_L / rho_L + (S_star - u_L) * (S_star + p_L / (rho_L * (S_L - u_L))));

    F_rho = F_L_rho + S_L * (U_star_rho - rho_L);
    F_rhou = F_L_rhou + S_L * (U_star_rhou - rho_L * u_L);
    F_E = F_L_E + S_L * (U_star_E - E_L);
  } else {
    double factor = rho_R * (S_R - u_R) / (S_R - S_star);
    double U_star_rho = factor;
    double U_star_rhou = factor * S_star;
    double U_star_E = factor * (E_R / rho_R + (S_star - u_R) * (S_star + p_R / (rho_R * (S_R - u_R))));

    F_rho = F_R_rho + S_R * (U_star_rho - rho_R);
    F_rhou = F_R_rhou + S_R * (U_star_rhou - rho_R * u_R);
    F_E = F_R_E + S_R * (U_star_E - E_R);
  }
}

void computeFluxesMUSCL(Grid& grid, std::vector<double>& F_rho,
                        std::vector<double>& F_rhou, std::vector<double>& F_E) {
  for (int i = 0; i <= grid.n; ++i) {
    double rho_L, u_L, p_L, rho_R, u_R, p_R;

    if (i == 0 || i == grid.n) {
      int idx = (i == 0) ? 0 : grid.n - 1;
      rho_L = rho_R = grid.rho[idx];
      u_L = u_R = grid.u[idx];
      p_L = p_R = grid.p[idx];
    } else {
      int i_L = i - 1;
      int i_R = i;

      // Left cell right edge
      double v_minus = (i_L > 0) ? grid.rho[i_L - 1] : grid.rho[i_L];
      double v_plus = (i_L < grid.n - 1) ? grid.rho[i_L + 1] : grid.rho[i_L];
      double slope = slopeLimit(v_minus, grid.rho[i_L], v_plus);
      rho_L = grid.rho[i_L] + 0.5 * slope;

      v_minus = (i_L > 0) ? grid.u[i_L - 1] : grid.u[i_L];
      v_plus = (i_L < grid.n - 1) ? grid.u[i_L + 1] : grid.u[i_L];
      slope = slopeLimit(v_minus, grid.u[i_L], v_plus);
      u_L = grid.u[i_L] + 0.5 * slope;

      v_minus = (i_L > 0) ? grid.p[i_L - 1] : grid.p[i_L];
      v_plus = (i_L < grid.n - 1) ? grid.p[i_L + 1] : grid.p[i_L];
      slope = slopeLimit(v_minus, grid.p[i_L], v_plus);
      p_L = grid.p[i_L] + 0.5 * slope;

      rho_L = std::max(rho_L, 1e-10);
      p_L = std::max(p_L, 1e-10);

      // Right cell left edge
      v_minus = (i_R > 0) ? grid.rho[i_R - 1] : grid.rho[i_R];
      v_plus = (i_R < grid.n - 1) ? grid.rho[i_R + 1] : grid.rho[i_R];
      slope = slopeLimit(v_minus, grid.rho[i_R], v_plus);
      rho_R = grid.rho[i_R] - 0.5 * slope;

      v_minus = (i_R > 0) ? grid.u[i_R - 1] : grid.u[i_R];
      v_plus = (i_R < grid.n - 1) ? grid.u[i_R + 1] : grid.u[i_R];
      slope = slopeLimit(v_minus, grid.u[i_R], v_plus);
      u_R = grid.u[i_R] - 0.5 * slope;

      v_minus = (i_R > 0) ? grid.p[i_R - 1] : grid.p[i_R];
      v_plus = (i_R < grid.n - 1) ? grid.p[i_R + 1] : grid.p[i_R];
      slope = slopeLimit(v_minus, grid.p[i_R], v_plus);
      p_R = grid.p[i_R] - 0.5 * slope;

      rho_R = std::max(rho_R, 1e-10);
      p_R = std::max(p_R, 1e-10);
    }

    solveHLLC(rho_L, u_L, p_L, rho_R, u_R, p_R, TEST.gamma, F_rho[i], F_rhou[i], F_E[i]);
  }
}

double computeTimestep(const Grid& grid) {
  double max_speed = 0.0;
  for (int i = 0; i < grid.n; ++i) {
    double a = std::sqrt(TEST.gamma * grid.p[i] / grid.rho[i]);
    max_speed = std::max(max_speed, std::abs(grid.u[i]) + a);
  }
  return TEST.CFL * grid.dx / max_speed;
}

void timeStepRK2(Grid& grid, double endTime) {
  double dt = computeTimestep(grid);
  if (grid.t + dt > endTime) dt = endTime - grid.t;

  std::vector<double> rho_n = grid.rho;
  std::vector<double> rho_u_n = grid.rho_u;
  std::vector<double> E_n = grid.E;

  std::vector<double> F_rho(grid.n + 1);
  std::vector<double> F_rhou(grid.n + 1);
  std::vector<double> F_E(grid.n + 1);

  // Stage 1
  computeFluxesMUSCL(grid, F_rho, F_rhou, F_E);

  std::vector<double> k1_rho(grid.n), k1_rho_u(grid.n), k1_E(grid.n);
  for (int i = 0; i < grid.n; ++i) {
    k1_rho[i] = -(1.0 / grid.dx) * (F_rho[i + 1] - F_rho[i]);
    k1_rho_u[i] = -(1.0 / grid.dx) * (F_rhou[i + 1] - F_rhou[i]);
    k1_E[i] = -(1.0 / grid.dx) * (F_E[i + 1] - F_E[i]);

    grid.rho[i] = rho_n[i] + dt * k1_rho[i];
    grid.rho_u[i] = rho_u_n[i] + dt * k1_rho_u[i];
    grid.E[i] = E_n[i] + dt * k1_E[i];

    grid.rho[i] = std::max(grid.rho[i], 1e-10);
    grid.E[i] = std::max(grid.E[i], 1e-10);
  }
  grid.consToPrim();

  // Stage 2
  computeFluxesMUSCL(grid, F_rho, F_rhou, F_E);

  for (int i = 0; i < grid.n; ++i) {
    double k2_rho = -(1.0 / grid.dx) * (F_rho[i + 1] - F_rho[i]);
    double k2_rho_u = -(1.0 / grid.dx) * (F_rhou[i + 1] - F_rhou[i]);
    double k2_E = -(1.0 / grid.dx) * (F_E[i + 1] - F_E[i]);

    grid.rho[i] = rho_n[i] + 0.5 * dt * (k1_rho[i] + k2_rho);
    grid.rho_u[i] = rho_u_n[i] + 0.5 * dt * (k1_rho_u[i] + k2_rho_u);
    grid.E[i] = E_n[i] + 0.5 * dt * (k1_E[i] + k2_E);

    grid.rho[i] = std::max(grid.rho[i], 1e-10);
    grid.E[i] = std::max(grid.E[i], 1e-10);
  }
  grid.consToPrim();
  grid.t += dt;
}

void runSimulation(Grid& grid, double endTime) {
  while (grid.t < endTime) {
    timeStepRK2(grid, endTime);
  }
}

// ============================================================================
// EXACT RIEMANN SOLVER (for analytical solution)
// ============================================================================

struct PrimitiveVars {
  double rho, u, p;
};

double solvePressureStar(double rho_L, double u_L, double p_L, double a_L,
                         double rho_R, double u_R, double p_R, double a_R,
                         double gamma) {
  const double TOL = 1e-10;
  const int MAX_ITER = 100;

  // Initial guess (PVRS)
  double p_star = std::max(TOL, 0.5 * (p_L + p_R) - 0.125 * (u_R - u_L) * (rho_L + rho_R) * (a_L + a_R));

  for (int iter = 0; iter < MAX_ITER; ++iter) {
    auto f = [gamma](double p, double p_k, double rho_k, double a_k) -> double {
      if (p > p_k) {
        double A = 2.0 / ((gamma + 1.0) * rho_k);
        double B = p_k * (gamma - 1.0) / (gamma + 1.0);
        return (p - p_k) * std::sqrt(A / (p + B));
      } else {
        return (2.0 * a_k / (gamma - 1.0)) *
               (std::pow(p / p_k, (gamma - 1.0) / (2.0 * gamma)) - 1.0);
      }
    };

    auto f_prime = [gamma](double p, double p_k, double rho_k, double a_k) -> double {
      if (p > p_k) {
        double A = 2.0 / ((gamma + 1.0) * rho_k);
        double B = p_k * (gamma - 1.0) / (gamma + 1.0);
        return std::sqrt(A / (p + B)) * (1.0 - 0.5 * (p - p_k) / (p + B));
      } else {
        return (1.0 / (rho_k * a_k)) * std::pow(p / p_k, -(gamma + 1.0) / (2.0 * gamma));
      }
    };

    double f_L = f(p_star, p_L, rho_L, a_L);
    double f_R = f(p_star, p_R, rho_R, a_R);
    double F = f_L + f_R + (u_R - u_L);
    double F_prime = f_prime(p_star, p_L, rho_L, a_L) + f_prime(p_star, p_R, rho_R, a_R);

    double p_new = std::max(TOL, p_star - F / F_prime);

    if (std::abs(p_new - p_star) / p_star < TOL) {
      return p_new;
    }
    p_star = p_new;
  }

  return p_star;
}

PrimitiveVars sampleExactSolution(double rho_L, double u_L, double p_L,
                                   double rho_R, double u_R, double p_R,
                                   double x, double t, double x0, double gamma) {
  PrimitiveVars result;
  const double TOL = 1e-10;

  double a_L = std::sqrt(gamma * p_L / rho_L);
  double a_R = std::sqrt(gamma * p_R / rho_R);

  double p_star = solvePressureStar(rho_L, u_L, p_L, a_L, rho_R, u_R, p_R, a_R, gamma);

  auto f = [gamma](double p, double p_k, double rho_k, double a_k) -> double {
    if (p > p_k) {
      double A = 2.0 / ((gamma + 1.0) * rho_k);
      double B = p_k * (gamma - 1.0) / (gamma + 1.0);
      return (p - p_k) * std::sqrt(A / (p + B));
    } else {
      return (2.0 * a_k / (gamma - 1.0)) *
             (std::pow(p / p_k, (gamma - 1.0) / (2.0 * gamma)) - 1.0);
    }
  };

  double u_star = 0.5 * (u_L + u_R) + 0.5 * (f(p_star, p_R, rho_R, a_R) - f(p_star, p_L, rho_L, a_L));

  double S = (t > TOL) ? (x - x0) / t : ((x < x0) ? -1e10 : 1e10);

  if (S <= u_star) {
    // Left of contact
    if (p_star > p_L) {
      // Left shock
      double S_L = u_L - a_L * std::sqrt((gamma + 1.0) / (2.0 * gamma) * (p_star / p_L) +
                                          (gamma - 1.0) / (2.0 * gamma));
      if (S <= S_L) {
        result.rho = rho_L;
        result.u = u_L;
        result.p = p_L;
      } else {
        result.rho = rho_L * ((p_star / p_L + (gamma - 1.0) / (gamma + 1.0)) /
                             ((gamma - 1.0) / (gamma + 1.0) * (p_star / p_L) + 1.0));
        result.u = u_star;
        result.p = p_star;
      }
    } else {
      // Left rarefaction
      double S_HL = u_L - a_L;
      if (S <= S_HL) {
        result.rho = rho_L;
        result.u = u_L;
        result.p = p_L;
      } else {
        double a_star_L = a_L * std::pow(p_star / p_L, (gamma - 1.0) / (2.0 * gamma));
        double S_TL = u_star - a_star_L;
        if (S > S_TL) {
          result.rho = rho_L * std::pow(p_star / p_L, 1.0 / gamma);
          result.u = u_star;
          result.p = p_star;
        } else {
          // Inside fan
          result.u = 2.0 / (gamma + 1.0) * (a_L + 0.5 * (gamma - 1.0) * u_L + S);
          double c = 2.0 / (gamma + 1.0) * (a_L + 0.5 * (gamma - 1.0) * (u_L - S));
          result.rho = rho_L * std::pow(c / a_L, 2.0 / (gamma - 1.0));
          result.p = p_L * std::pow(c / a_L, 2.0 * gamma / (gamma - 1.0));
        }
      }
    }
  } else {
    // Right of contact
    if (p_star > p_R) {
      // Right shock
      double S_R = u_R + a_R * std::sqrt((gamma + 1.0) / (2.0 * gamma) * (p_star / p_R) +
                                          (gamma - 1.0) / (2.0 * gamma));
      if (S >= S_R) {
        result.rho = rho_R;
        result.u = u_R;
        result.p = p_R;
      } else {
        result.rho = rho_R * ((p_star / p_R + (gamma - 1.0) / (gamma + 1.0)) /
                             ((gamma - 1.0) / (gamma + 1.0) * (p_star / p_R) + 1.0));
        result.u = u_star;
        result.p = p_star;
      }
    } else {
      // Right rarefaction
      double S_HR = u_R + a_R;
      if (S >= S_HR) {
        result.rho = rho_R;
        result.u = u_R;
        result.p = p_R;
      } else {
        double a_star_R = a_R * std::pow(p_star / p_R, (gamma - 1.0) / (2.0 * gamma));
        double S_TR = u_star + a_star_R;
        if (S < S_TR) {
          result.rho = rho_R * std::pow(p_star / p_R, 1.0 / gamma);
          result.u = u_star;
          result.p = p_star;
        } else {
          // Inside fan
          result.u = 2.0 / (gamma + 1.0) * (-a_R + 0.5 * (gamma - 1.0) * u_R + S);
          double c = 2.0 / (gamma + 1.0) * (a_R - 0.5 * (gamma - 1.0) * (u_R - S));
          result.rho = rho_R * std::pow(c / a_R, 2.0 / (gamma - 1.0));
          result.p = p_R * std::pow(c / a_R, 2.0 * gamma / (gamma - 1.0));
        }
      }
    }
  }

  return result;
}

void computeAnalyticalSolution(const Grid& grid, const TestParams& test,
                                std::vector<double>& rho_ex,
                                std::vector<double>& u_ex,
                                std::vector<double>& p_ex,
                                std::vector<double>& s_ex) {
  rho_ex.resize(grid.n);
  u_ex.resize(grid.n);
  p_ex.resize(grid.n);
  s_ex.resize(grid.n);

  for (int i = 0; i < grid.n; ++i) {
    PrimitiveVars exact = sampleExactSolution(
      test.rho_L, test.u_L, test.p_L,
      test.rho_R, test.u_R, test.p_R,
      grid.x[i], test.endTime, test.x_diaphragm, test.gamma);

    rho_ex[i] = exact.rho;
    u_ex[i] = exact.u;
    p_ex[i] = exact.p;
    s_ex[i] = std::log(exact.p / std::pow(exact.rho, test.gamma));
  }
}

// ============================================================================
// MONTE CARLO
// ============================================================================

struct MCResults {
  std::vector<double> x;

  std::vector<double> mean_rho, mean_u, mean_p, mean_s;
  std::vector<double> std_rho, std_u, std_p, std_s;
  std::vector<double> ci95_lower_rho, ci95_upper_rho;
  std::vector<double> ci95_lower_u, ci95_upper_u;
  std::vector<double> ci95_lower_p, ci95_upper_p;
  std::vector<double> ci95_lower_s, ci95_upper_s;

  int num_trials;
  double noise_level;
  double computation_time_ms;
};

std::vector<Sensor> addNoiseToSensors(const std::vector<Sensor>& base, double noise_level, std::mt19937& gen) {
  std::vector<Sensor> noisy = base;

  for (auto& s : noisy) {
    std::normal_distribution<double> noise_rho(0.0, noise_level * s.rho);
    std::normal_distribution<double> noise_u(0.0, noise_level * (std::abs(s.u) + 0.1));
    std::normal_distribution<double> noise_p(0.0, noise_level * s.p);

    s.rho = std::max(0.01, s.rho + noise_rho(gen));
    s.u = s.u + noise_u(gen);
    s.p = std::max(0.01, s.p + noise_p(gen));
  }

  return noisy;
}

MCResults runMonteCarlo(int num_cells, const std::vector<Sensor>& base_sensors,
                        int num_trials, double noise_level, const TestParams& test) {
  auto start_time = std::chrono::high_resolution_clock::now();

  MCResults results;
  results.num_trials = num_trials;
  results.noise_level = noise_level;

  const int N = num_cells;

  // Initialize
  Grid template_grid(N, 1.0);
  results.x = template_grid.x;

  std::vector<double> sum_rho(N, 0.0), sum_u(N, 0.0), sum_p(N, 0.0), sum_s(N, 0.0);
  std::vector<double> sum_sq_rho(N, 0.0), sum_sq_u(N, 0.0), sum_sq_p(N, 0.0), sum_sq_s(N, 0.0);

#ifdef _OPENMP
  int num_threads = omp_get_max_threads();
  std::cout << "\nMonte Carlo: " << num_trials << " trials, "
            << (noise_level * 100) << "% noise, " << N << " cells, "
            << num_threads << " threads\n";
#else
  std::cout << "\nMonte Carlo: " << num_trials << " trials, "
            << (noise_level * 100) << "% noise, " << N << " cells (serial)\n";
#endif

  int completed_trials = 0;

#ifdef _OPENMP
  #pragma omp parallel
  {
    // Thread-local accumulators
    std::vector<double> local_sum_rho(N, 0.0), local_sum_u(N, 0.0);
    std::vector<double> local_sum_p(N, 0.0), local_sum_s(N, 0.0);
    std::vector<double> local_sum_sq_rho(N, 0.0), local_sum_sq_u(N, 0.0);
    std::vector<double> local_sum_sq_p(N, 0.0), local_sum_sq_s(N, 0.0);

    // Thread-local RNG with unique seed
    std::mt19937 gen(std::random_device{}() + omp_get_thread_num() * 1000);

    #pragma omp for schedule(dynamic, 4)
    for (int trial = 0; trial < num_trials; ++trial) {
      // Add noise to sensors
      std::vector<Sensor> noisy_sensors = addNoiseToSensors(base_sensors, noise_level, gen);

      // Initialize grid from noisy sensors
      Grid grid(N, 1.0);
      grid.initializeFromSensors(noisy_sensors);

      // Run simulation
      runSimulation(grid, test.endTime);

      // Accumulate locally
      std::vector<double> s = grid.getEntropy();
      for (int i = 0; i < N; ++i) {
        local_sum_rho[i] += grid.rho[i];
        local_sum_u[i] += grid.u[i];
        local_sum_p[i] += grid.p[i];
        local_sum_s[i] += s[i];

        local_sum_sq_rho[i] += grid.rho[i] * grid.rho[i];
        local_sum_sq_u[i] += grid.u[i] * grid.u[i];
        local_sum_sq_p[i] += grid.p[i] * grid.p[i];
        local_sum_sq_s[i] += s[i] * s[i];
      }

      #pragma omp atomic
      completed_trials++;

      // Progress reporting (only from one thread occasionally)
      if (completed_trials % (num_trials / 10 + 1) == 0) {
        #pragma omp critical
        {
          std::cout << "  Progress: " << completed_trials << "/" << num_trials
                    << " (" << (100 * completed_trials / num_trials) << "%)\n";
        }
      }
    }

    // Reduce thread-local sums into global sums
    #pragma omp critical
    {
      for (int i = 0; i < N; ++i) {
        sum_rho[i] += local_sum_rho[i];
        sum_u[i] += local_sum_u[i];
        sum_p[i] += local_sum_p[i];
        sum_s[i] += local_sum_s[i];
        sum_sq_rho[i] += local_sum_sq_rho[i];
        sum_sq_u[i] += local_sum_sq_u[i];
        sum_sq_p[i] += local_sum_sq_p[i];
        sum_sq_s[i] += local_sum_sq_s[i];
      }
    }
  }
#else
  // Serial fallback
  std::mt19937 gen(std::random_device{}());

  for (int trial = 0; trial < num_trials; ++trial) {
    // Add noise to sensors
    std::vector<Sensor> noisy_sensors = addNoiseToSensors(base_sensors, noise_level, gen);

    // Initialize grid from noisy sensors
    Grid grid(N, 1.0);
    grid.initializeFromSensors(noisy_sensors);

    // Run simulation
    runSimulation(grid, test.endTime);

    // Accumulate
    std::vector<double> s = grid.getEntropy();
    for (int i = 0; i < N; ++i) {
      sum_rho[i] += grid.rho[i];
      sum_u[i] += grid.u[i];
      sum_p[i] += grid.p[i];
      sum_s[i] += s[i];

      sum_sq_rho[i] += grid.rho[i] * grid.rho[i];
      sum_sq_u[i] += grid.u[i] * grid.u[i];
      sum_sq_p[i] += grid.p[i] * grid.p[i];
      sum_sq_s[i] += s[i] * s[i];
    }

    completed_trials++;
    if (completed_trials % (num_trials / 10 + 1) == 0) {
      std::cout << "  Progress: " << completed_trials << "/" << num_trials
                << " (" << (100 * completed_trials / num_trials) << "%)\n";
    }
  }
#endif

  // Compute statistics
  double n = static_cast<double>(num_trials);
  const double z95 = 1.96;

  results.mean_rho.resize(N);
  results.mean_u.resize(N);
  results.mean_p.resize(N);
  results.mean_s.resize(N);
  results.std_rho.resize(N);
  results.std_u.resize(N);
  results.std_p.resize(N);
  results.std_s.resize(N);
  results.ci95_lower_rho.resize(N);
  results.ci95_upper_rho.resize(N);
  results.ci95_lower_u.resize(N);
  results.ci95_upper_u.resize(N);
  results.ci95_lower_p.resize(N);
  results.ci95_upper_p.resize(N);
  results.ci95_lower_s.resize(N);
  results.ci95_upper_s.resize(N);

  for (int i = 0; i < N; ++i) {
    results.mean_rho[i] = sum_rho[i] / n;
    results.mean_u[i] = sum_u[i] / n;
    results.mean_p[i] = sum_p[i] / n;
    results.mean_s[i] = sum_s[i] / n;

    double var_rho = std::max(0.0, (sum_sq_rho[i] / n) - results.mean_rho[i] * results.mean_rho[i]);
    double var_u = std::max(0.0, (sum_sq_u[i] / n) - results.mean_u[i] * results.mean_u[i]);
    double var_p = std::max(0.0, (sum_sq_p[i] / n) - results.mean_p[i] * results.mean_p[i]);
    double var_s = std::max(0.0, (sum_sq_s[i] / n) - results.mean_s[i] * results.mean_s[i]);

    results.std_rho[i] = std::sqrt(var_rho);
    results.std_u[i] = std::sqrt(var_u);
    results.std_p[i] = std::sqrt(var_p);
    results.std_s[i] = std::sqrt(var_s);

    results.ci95_lower_rho[i] = results.mean_rho[i] - z95 * results.std_rho[i];
    results.ci95_upper_rho[i] = results.mean_rho[i] + z95 * results.std_rho[i];
    results.ci95_lower_u[i] = results.mean_u[i] - z95 * results.std_u[i];
    results.ci95_upper_u[i] = results.mean_u[i] + z95 * results.std_u[i];
    results.ci95_lower_p[i] = results.mean_p[i] - z95 * results.std_p[i];
    results.ci95_upper_p[i] = results.mean_p[i] + z95 * results.std_p[i];
    results.ci95_lower_s[i] = results.mean_s[i] - z95 * results.std_s[i];
    results.ci95_upper_s[i] = results.mean_s[i] + z95 * results.std_s[i];
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  results.computation_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

  std::cout << "Monte Carlo complete: " << results.computation_time_ms << " ms\n";

  return results;
}

// ============================================================================
// OUTPUT
// ============================================================================

void writeDataFiles(const MCResults& mc,
                    const std::vector<double>& rho_ex,
                    const std::vector<double>& u_ex,
                    const std::vector<double>& p_ex,
                    const std::vector<double>& s_ex,
                    const std::vector<Sensor>& sensors) {
  
  // Monte Carlo results
  std::ofstream mc_file("mc_results.dat");
  mc_file << "# x  mean_rho  std_rho  ci95_lo_rho  ci95_hi_rho  mean_u  std_u  ci95_lo_u  ci95_hi_u  "
          << "mean_p  std_p  ci95_lo_p  ci95_hi_p  mean_s  std_s  ci95_lo_s  ci95_hi_s\n";
  for (size_t i = 0; i < mc.x.size(); ++i) {
    mc_file << mc.x[i] << " "
            << mc.mean_rho[i] << " " << mc.std_rho[i] << " "
            << mc.ci95_lower_rho[i] << " " << mc.ci95_upper_rho[i] << " "
            << mc.mean_u[i] << " " << mc.std_u[i] << " "
            << mc.ci95_lower_u[i] << " " << mc.ci95_upper_u[i] << " "
            << mc.mean_p[i] << " " << mc.std_p[i] << " "
            << mc.ci95_lower_p[i] << " " << mc.ci95_upper_p[i] << " "
            << mc.mean_s[i] << " " << mc.std_s[i] << " "
            << mc.ci95_lower_s[i] << " " << mc.ci95_upper_s[i] << "\n";
  }
  mc_file.close();

  // Analytical solution
  std::ofstream ana_file("analytical.dat");
  ana_file << "# x  rho  u  p  s\n";
  for (size_t i = 0; i < mc.x.size(); ++i) {
    ana_file << mc.x[i] << " " << rho_ex[i] << " " << u_ex[i] << " " << p_ex[i] << " " << s_ex[i] << "\n";
  }
  ana_file.close();

  // Sensor positions
  std::ofstream sen_file("sensors.dat");
  sen_file << "# x  rho  u  p\n";
  for (const auto& s : sensors) {
    sen_file << s.x << " " << s.rho << " " << s.u << " " << s.p << "\n";
  }
  sen_file.close();

  std::cout << "\nData files written:\n";
  std::cout << "  mc_results.dat  - Monte Carlo statistics\n";
  std::cout << "  analytical.dat  - Exact Riemann solution\n";
  std::cout << "  sensors.dat     - Sensor positions/values\n";
}

void generatePlotScript(const TestParams& test, const MCResults& mc) {
  std::ofstream gp("plot_mc.gp");

  gp << "set terminal pngcairo size 1400,1000 enhanced font 'Arial,11'\n";
  gp << "set output 'mc_uncertainty.png'\n";
  gp << "set multiplot layout 2,2 title '" << test.name << " - Monte Carlo Uncertainty ("
     << mc.num_trials << " trials, " << (mc.noise_level * 100) << "% noise)'\n";
  gp << "set grid\n\n";

  // Density
  gp << "set title 'Density'\n";
  gp << "set xlabel 'Position (m)'\n";
  gp << "set ylabel 'ρ (kg/m³)'\n";
  gp << "plot 'mc_results.dat' u 1:4:5 w filledcurves lc rgb '#aaccff' title '95% CI', \\\n";
  gp << "     'mc_results.dat' u 1:2 w l lw 2 lc rgb 'blue' title 'MC Mean', \\\n";
  gp << "     'analytical.dat' u 1:2 w l lw 2 dt 2 lc rgb 'black' title 'Analytical', \\\n";
  gp << "     'sensors.dat' u 1:2 w p pt 7 ps 1.5 lc rgb 'green' title 'Sensors'\n\n";

  // Velocity
  gp << "set title 'Velocity'\n";
  gp << "set xlabel 'Position (m)'\n";
  gp << "set ylabel 'u (m/s)'\n";
  gp << "plot 'mc_results.dat' u 1:8:9 w filledcurves lc rgb '#aaccff' title '95% CI', \\\n";
  gp << "     'mc_results.dat' u 1:6 w l lw 2 lc rgb 'blue' title 'MC Mean', \\\n";
  gp << "     'analytical.dat' u 1:3 w l lw 2 dt 2 lc rgb 'black' title 'Analytical', \\\n";
  gp << "     'sensors.dat' u 1:3 w p pt 7 ps 1.5 lc rgb 'green' title 'Sensors'\n\n";

  // Pressure
  gp << "set title 'Pressure'\n";
  gp << "set xlabel 'Position (m)'\n";
  gp << "set ylabel 'p (Pa)'\n";
  gp << "plot 'mc_results.dat' u 1:12:13 w filledcurves lc rgb '#aaccff' title '95% CI', \\\n";
  gp << "     'mc_results.dat' u 1:10 w l lw 2 lc rgb 'blue' title 'MC Mean', \\\n";
  gp << "     'analytical.dat' u 1:4 w l lw 2 dt 2 lc rgb 'black' title 'Analytical', \\\n";
  gp << "     'sensors.dat' u 1:4 w p pt 7 ps 1.5 lc rgb 'green' title 'Sensors'\n\n";

  // Entropy
  gp << "set title 'Entropy'\n";
  gp << "set xlabel 'Position (m)'\n";
  gp << "set ylabel 's (J/kg·K)'\n";
  gp << "plot 'mc_results.dat' u 1:16:17 w filledcurves lc rgb '#aaccff' title '95% CI', \\\n";
  gp << "     'mc_results.dat' u 1:14 w l lw 2 lc rgb 'blue' title 'MC Mean', \\\n";
  gp << "     'analytical.dat' u 1:5 w l lw 2 dt 2 lc rgb 'black' title 'Analytical'\n\n";

  gp << "unset multiplot\n";
  gp.close();

  int ret = system("gnuplot plot_mc.gp 2>/dev/null");
  if (ret == 0) {
    std::cout << "\nPlot generated: mc_uncertainty.png\n";
  } else {
    std::cout << "\nWARNING: gnuplot not found. Run manually: gnuplot plot_mc.gp\n";
  }
}

void printStatistics(const MCResults& mc,
                     const std::vector<double>& rho_ex,
                     const std::vector<double>& u_ex,
                     const std::vector<double>& p_ex,
                     const std::vector<double>& s_ex) {
  int N = mc.x.size();

  // Compute errors of MC mean vs analytical
  double L1_rho = 0, L1_u = 0, L1_p = 0, L1_s = 0;
  double Linf_rho = 0, Linf_u = 0, Linf_p = 0, Linf_s = 0;
  double max_std_rho = 0, max_std_u = 0, max_std_p = 0, max_std_s = 0;

  for (int i = 0; i < N; ++i) {
    double err_rho = std::abs(mc.mean_rho[i] - rho_ex[i]);
    double err_u = std::abs(mc.mean_u[i] - u_ex[i]);
    double err_p = std::abs(mc.mean_p[i] - p_ex[i]);
    double err_s = std::abs(mc.mean_s[i] - s_ex[i]);

    L1_rho += err_rho;
    L1_u += err_u;
    L1_p += err_p;
    L1_s += err_s;

    Linf_rho = std::max(Linf_rho, err_rho);
    Linf_u = std::max(Linf_u, err_u);
    Linf_p = std::max(Linf_p, err_p);
    Linf_s = std::max(Linf_s, err_s);

    max_std_rho = std::max(max_std_rho, mc.std_rho[i]);
    max_std_u = std::max(max_std_u, mc.std_u[i]);
    max_std_p = std::max(max_std_p, mc.std_p[i]);
    max_std_s = std::max(max_std_s, mc.std_s[i]);
  }

  std::cout << "\n=== MC Mean vs Analytical (Error Norms) ===\n";
  std::cout << std::scientific << std::setprecision(4);
  std::cout << "         L1          Linf\n";
  std::cout << "rho:  " << L1_rho / N << "   " << Linf_rho << "\n";
  std::cout << "u:    " << L1_u / N << "   " << Linf_u << "\n";
  std::cout << "p:    " << L1_p / N << "   " << Linf_p << "\n";
  std::cout << "s:    " << L1_s / N << "   " << Linf_s << "\n";

  std::cout << "\n=== MC Standard Deviation (Max) ===\n";
  std::cout << "rho:  " << max_std_rho << "\n";
  std::cout << "u:    " << max_std_u << "\n";
  std::cout << "p:    " << max_std_p << "\n";
  std::cout << "s:    " << max_std_s << "\n";
}

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char* argv[]) {
  // Parse arguments
  int num_trials = (argc > 1) ? std::atoi(argv[1]) : 100;
  double noise_level = (argc > 2) ? std::atof(argv[2]) : 0.05;
  int num_cells = (argc > 3) ? std::atoi(argv[3]) : 500;

  std::cout << "\n";
  std::cout << "============================================================\n";
  std::cout << "  Monte Carlo Uncertainty Propagation Test\n";
  std::cout << "============================================================\n";
  std::cout << "\n";
  std::cout << "Test Case:   " << TEST.name << "\n";
  std::cout << "Grid Cells:  " << num_cells << "\n";
  std::cout << "MC Trials:   " << num_trials << "\n";
  std::cout << "Noise Level: " << (noise_level * 100) << "%\n";
  std::cout << "End Time:    " << TEST.endTime << " s\n";
  std::cout << "\n";

  // Define sensor positions
  std::vector<double> sensor_positions = {0.1, 0.25, 0.4, 0.6, 0.75, 0.9};
  std::cout << "Sensor positions: ";
  for (double x : sensor_positions) std::cout << x << " ";
  std::cout << "\n";

  // Get exact sensor readings (from perfect IC)
  std::vector<Sensor> exact_sensors = getExactSensors(TEST, sensor_positions);

  std::cout << "\nExact sensor values (t=0):\n";
  std::cout << "  x       rho      u        p\n";
  for (const auto& s : exact_sensors) {
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  " << s.x << "   " << s.rho << "    " << s.u << "    " << s.p << "\n";
  }

  // Run Monte Carlo
  MCResults mc = runMonteCarlo(num_cells, exact_sensors, num_trials, noise_level, TEST);

  // Compute analytical solution from EXACT initial conditions
  std::cout << "\nComputing analytical solution from exact IC...\n";
  Grid ref_grid(num_cells, 1.0);
  std::vector<double> rho_ex, u_ex, p_ex, s_ex;
  computeAnalyticalSolution(ref_grid, TEST, rho_ex, u_ex, p_ex, s_ex);

  // Output
  writeDataFiles(mc, rho_ex, u_ex, p_ex, s_ex, exact_sensors);
  printStatistics(mc, rho_ex, u_ex, p_ex, s_ex);
  generatePlotScript(TEST, mc);

  std::cout << "\nDone!\n\n";

  return 0;
}