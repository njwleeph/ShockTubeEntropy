#include "solver.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

/**
 * Primitive <-> Conservative conversions
 */
ShockSolver::ConservativeVars ShockSolver::PrimitiveVars::toConservative(double gamma) const {
  ConservativeVars U;
  U.rho = rho;
  U.rho_u = rho * u;
  double e = p / ((gamma - 1.0) * rho);
  U.E = rho * (e + 0.5 * u * u);
  return U;
}

/**
 * Constructor
 */
ShockSolver::ShockSolver(const Config& cfg)
  : cfg_(cfg),
    dx_(cfg.length / cfg.numCells),
    t_(0.0),
    step_count_(0),
    flux_type_("hllc") {
  
  rho_.resize(cfg_.numCells);
  u_.resize(cfg_.numCells);
  p_.resize(cfg_.numCells);
  rho_u_.resize(cfg_.numCells);
  E_.resize(cfg_.numCells);
  x_.resize(cfg_.numCells);

  for (int i = 0; i < cfg_.numCells; i++) {
    x_[i] = (i + 0.5) * dx_;
  }

  initial_conditions_.is_valid = false;
}

/**
 * Initialization
 */
void ShockSolver::initializeShockTube(double rho_L, double u_L, double p_L,
                                      double rho_R, double u_R, double p_R,
                                      double x_diaphragm, double endTime) {

  initial_conditions_.rho_L = rho_L;
  initial_conditions_.u_L = u_L;
  initial_conditions_.p_L = p_L;
  initial_conditions_.rho_R = rho_R;
  initial_conditions_.u_R = u_R;
  initial_conditions_.p_R = p_R;
  initial_conditions_.x_diaphragm = x_diaphragm;
  initial_conditions_.endTime = endTime;
  initial_conditions_.is_valid = true;

  for (int i = 0; i < cfg_.numCells; ++i) {
    if (x_[i] < x_diaphragm) {
      rho_[i] = rho_L;
      u_[i] = u_L;
      p_[i] = p_L;
    } else {
      rho_[i] = rho_R;
      u_[i] = u_R;
      p_[i] = p_R;
    }
  }
  primitivesToConservative();
}

void ShockSolver::initializeCustom(const std::vector<PrimitiveVars>& initial) {
  if (initial.size() != static_cast<size_t>(cfg_.numCells)) {
    throw std::runtime_error("Initial conditions size mismatch");
  }

  initial_conditions_.is_valid = false;

  for (int i = 0; i < cfg_.numCells; ++i) {
    rho_[i] = initial[i].rho;
    u_[i] = initial[i].u;
    p_[i] = initial[i].p;
  }
  primitivesToConservative();
}

/**
 * Initialize verification cross-reference with Toro
 */
void ShockSolver::initializeToroTest(ToroTest test) {
  double rho_L, u_L, p_L, rho_R, u_R, p_R, x_diaphragm, endTime;
  
  switch (test) {
    case ToroTest::TEST1_SOD:
      // Test 1: Sod's problem (standard)
      rho_L = 1.0; u_L = 0.0; p_L = 1.0;
      rho_R = 0.125; u_R = 0.0; p_R = 0.1;
      x_diaphragm = 0.5;
      endTime = 0.25;
      break;
      
    case ToroTest::TEST2_123:
      // Test 2: 123 problem (strong rarefaction)
      rho_L = 1.0; u_L = -2.0; p_L = 0.4;
      rho_R = 1.0; u_R = 2.0; p_R = 0.4;
      x_diaphragm = 0.5;
      endTime = 0.15;
      break;
      
    case ToroTest::TEST3_BLAST_LEFT:
      // Test 3: Left half of blast wave problem
      rho_L = 1.0; u_L = 0.0; p_L = 1000.0;
      rho_R = 1.0; u_R = 0.0; p_R = 0.01;
      x_diaphragm = 0.5;
      endTime = 0.012;
      break;
      
    case ToroTest::TEST4_SLOW_SHOCK:
      // Test 4: Collision of two rarefaction waves
      rho_L = 1.0; u_L = 0.0; p_L = 0.01;
      rho_R = 1.0; u_R = 0.0; p_R = 100.0;
      x_diaphragm = 0.5;
      endTime = 0.035;
      break;
      
    case ToroTest::TEST5_COLLISION:
      // Test 5: Stationary contact discontinuity
      rho_L = 5.99924; u_L = 19.59745; p_L = 460.894;
      rho_R = 5.99242; u_R = -6.19633; p_R = 46.0950;
      x_diaphragm = 0.5;
      endTime = 0.035;
      break;
  }
  
  std::cout << "Initializing Toro Test " << static_cast<int>(test) + 1 << std::endl;
  initializeShockTube(rho_L, u_L, p_L, rho_R, u_R, p_R, x_diaphragm, endTime);
}

/**
 * Setter
 */
void ShockSolver::setFlux(const std::string& flux_name) {
  if (flux_name != "HLLC" && flux_name != "EntropyStable") throw std::runtime_error("Unknown flux type: " + flux_name);

  flux_type_ = flux_name;
}

/**
 * Variable Conversions
 */
void ShockSolver::primitivesToConservative() {
  for (int i = 0; i < cfg_.numCells; ++i) {
    rho_u_[i] = rho_[i] * u_[i];
    double e = p_[i] / ((cfg_.gamma - 1.0) * rho_[i]);
    E_[i] = rho_[i] * (e + 0.5 * u_[i] * u_[i]);
  }
}

void ShockSolver::conservativesToPrimitives() {
  for (int i = 0; i < cfg_.numCells; ++i) {
    u_[i] = rho_u_[i] / rho_[i];
    p_[i] = (cfg_.gamma - 1.0) * (E_[i] - 0.5 * rho_[i] * u_[i] * u_[i]);
    if (p_[i] < 0.0) {
      p_[i] = 1e-10;
    }
  }
}

/**
 * Time Integration
 */
double ShockSolver::computeTimestep() const {
  double max_speed = 0.0;
  for (int i = 0; i < cfg_.numCells; ++i) {
    double a = std::sqrt(cfg_.gamma * p_[i] / rho_[i]);
    double speed = std::abs(u_[i]) + a;
    max_speed = std::max(max_speed, speed);
  }
  
  return cfg_.CFL * dx_ / max_speed;
}

void ShockSolver::step() {
  double dt = computeTimestep();
  if (t_ + dt > cfg_.endTime) dt = cfg_.endTime - t_;

  applyBoundaryConditions();

  // Save initial state
  std::vector<double> rho_n = rho_;
  std::vector<double> rho_u_n = rho_u_;
  std::vector<double> E_n = E_;

  std::vector<double> F_rho(cfg_.numCells + 1);
  std::vector<double> F_rhou(cfg_.numCells + 1);
  std::vector<double> F_E(cfg_.numCells + 1);
  computeFluxesMUSCL(F_rho, F_rhou, F_E);

  std::vector<double> k1_rho(cfg_.numCells), k1_rho_u(cfg_.numCells), k1_E(cfg_.numCells);
  for (int i = 0; i < cfg_.numCells; ++i) {
    k1_rho[i] = -(1.0 / dx_) * (F_rho[i + 1] - F_rho[i]);
    k1_rho_u[i] = -(1.0 / dx_) * (F_rhou[i + 1] - F_rhou[i]);
    k1_E[i] = -(1.0 / dx_) * (F_E[i + 1] - F_E[i]);

    // U*
    rho_[i] = rho_n[i] + dt * k1_rho[i];
    rho_u_[i] = rho_u_n[i] + dt * k1_rho_u[i];
    E_[i] = E_n[i] + dt * k1_E[i];

    rho_[i] = std::max(rho_[i], 1e-10);
    E_[i] = std::max(E_[i], 1e-10);
  }
  conservativesToPrimitives();

  // Stage 2: k2
  computeFluxesMUSCL(F_rho, F_rhou, F_E);

  std::vector<double> k2_rho(cfg_.numCells), k2_rho_u(cfg_.numCells), k2_E(cfg_.numCells);
  for (int i = 0; i < cfg_.numCells; ++i) {
    k2_rho[i] = -(1.0 / dx_) * (F_rho[i + 1] - F_rho[i]);
    k2_rho_u[i] = -(1.0 / dx_) * (F_rhou[i + 1] - F_rhou[i]);
    k2_E[i] = -(1.0 / dx_) * (F_E[i + 1] - F_E[i]);

    rho_[i] = rho_n[i] + 0.5 * dt * (k1_rho[i] + k2_rho[i]);
    rho_u_[i] = rho_u_n[i] + 0.5 * dt * (k1_rho_u[i] + k2_rho_u[i]);
    E_[i] = E_n[i] + 0.5 * dt * (k1_E[i] + k2_E[i]);

    rho_[i] = std::max(rho_[i], 1e-10);
    E_[i] = std::max(E_[i], 1e-10);
  }
  conservativesToPrimitives();
  t_ += dt;
}

void ShockSolver::run() {
  std::cout << "Running ShockSolver...";
  std::cout << "\nGrid: " << cfg_.numCells << " cells, dx = " << dx_ << " m\n";
  std::cout << "CFL: " << cfg_.CFL 
            << ", gamma = " << cfg_.gamma << "\n\n";

  while (t_ < cfg_.endTime) {
    step();
    step_count_++;
    if (step_count_ % 100 == 0) {
      std::cout << "Step " << step_count_ 
                << ", t = " << t_ << " s"
                << ", dt = " << computeTimestep() << " s\n";
    }
  }

  std::cout << "\nSimulation complete: " << step_count_ << " steps\n";
  std::cout << "Final time: " << t_ << " s\n";
}

/**
 * Flux Computation
 */
void ShockSolver::solveEntropyConservative(double rho_L, double u_L, double p_L,
                                           double rho_R, double u_R, double p_R,
                                           double& F_rho, double& F_rhou, double& F_E) {
  // Define log Mean locally
  auto logMean = [](double a, double b) {
    if (std::abs(a - b) < 1e-10) return a;

    a = std::max(a, 1e-10);
    b = std::max(b, 1e-10);

    double xi = b / a;
    double f = (xi - 1.0) / (xi + 1.0);
    double u = f * f;

    // Edge case if a = b (use taylor expand)
    if (u < 1e-2) {
      double F = 1.0 + u / 3.0 + u * u / 5.0 + u * u * u / 7.0;
      return 0.5 * (a + b) / F;
    }

    return (b - a) / std::log(xi);
  };

  // Safety
  rho_L = std::max(rho_L, 1e-10);
  rho_R = std::max(rho_R, 1e-10);
  p_L = std::max(p_L, 1e-10);
  p_R = std::max(p_R, 1e-10);

  // Log means
  double rho_log = logMean(rho_L, rho_R);
  double p_log = logMean(p_L, p_R);

  // Arithmetic averages
  double u_avg = 0.5 * (u_L + u_R);

  // Specific internal energy
  double e_L = p_L / ((cfg_.gamma - 1.0) * rho_L);
  double e_R = p_R / ((cfg_.gamma - 1.0) * rho_R);
  double e_avg = 0.5 * (e_L + e_R);

  // Kinetic energy
  double ke_avg = 0.5 * u_avg * u_avg;

  // Entropy Conservative Flux
  F_rho = rho_log * u_avg;
  F_rhou = rho_log * u_avg * u_avg + p_log;
  F_E = rho_log * u_avg * (e_avg + ke_avg) + p_log * u_avg;
}

/**
 * Compute dissipation (Ismael-Roe 2009)
 */
void ShockSolver::computeDissipation(double rho_L, double u_L, double p_L,
                                double rho_R, double u_R, double p_R,
                                double& D_rho, double& D_rhou, double& D_E) {
  // Average States
  double rho_avg = 0.5 * (rho_L + rho_R);
  double u_avg = 0.5 * (u_L + u_R);
  double p_avg = 0.5 * (p_L + p_R);

  // Average sound speed
  double a_avg = std::sqrt(cfg_.gamma * p_avg / rho_avg);

  // Jumps in prim variables
  double drho = rho_R - rho_L;
  double du = u_R - u_L;
  double dp = p_R - p_L;

  // Compute specific enthalpy
  double H_L = (cfg_.gamma / (cfg_.gamma - 1.0)) * (p_L / rho_L) + 0.5 * u_L * u_L;
  double H_R = (cfg_.gamma / (cfg_.gamma - 1.0)) * (p_R / rho_R) + 0.5 * u_R * u_R;
  double H_avg = 0.5 * (H_L + H_R);

  // Eigenvalues of Euler system:
  double lambda1 = std::abs(u_avg - a_avg);   // Left-going wave
  double lambda2 = std::abs(u_avg);           // Contact
  double lambda3 = std::abs(u_avg + a_avg);   // Right-going wave

  // Characteristic jumps (project primitive jumps onto eigenvectors)
  double w1 = 0.5 * (dp / (a_avg * a_avg) - rho_avg * du / a_avg);          // Left wave amplitude
  double w2 = drho - dp / (a_avg * a_avg);                                  // Contact wave amplitude
  double w3 = 0.5 * (dp / (a_avg * a_avg) + rho_avg * du / a_avg);          // Right wave amplitude

  // Dissipation
  D_rho = lambda1 * w1 * 1.0 + lambda2 * w2 * 1.0 + lambda3 * w3 * 1.0;
  D_rhou = lambda1 * w1 * (u_avg - a_avg) + lambda2 * w2 * u_avg + lambda3 * w3 * (u_avg + a_avg);
  D_E = lambda1 * w1 * (H_avg - u_avg * a_avg) + lambda2 * w2 * (0.5 * u_avg * u_avg) + lambda3 * w3 * (H_avg + u_avg * a_avg);
}

/**
 * Full entropy-stable flux
 */
void ShockSolver::solveEntropyStable(double rho_L, double u_L, double p_L,
                                     double rho_R, double u_R, double p_R,
                                     double& F_rho, double& F_rhou, double& F_E) {
  double F_EC_rho, F_EC_rhou, F_EC_E;
  solveEntropyConservative(rho_L, u_L, p_L, rho_R, u_R, p_R, F_EC_rho, F_EC_rhou, F_EC_E);

  double D_rho, D_rhou, D_E;
  computeDissipation(rho_L, u_L, p_L, rho_R, u_R, p_R, D_rho, D_rhou, D_E);

  F_rho = F_EC_rho - 0.5 * D_rho;
  F_rhou = F_EC_rhou - 0.5 * D_rhou;
  F_E = F_EC_E - 0.5 * D_E;
}

/**
 * MUSCL - RK2 Spatial Reconstruction
 */
void ShockSolver::computeFluxesMUSCL(std::vector<double>& F_rho, std::vector<double>& F_rhou, std::vector<double>& F_E) {
  for (int i = 0; i <= cfg_.numCells; ++i) {
    double rho_L, u_L, p_L, rho_R, u_R, p_R;

    if (i == 0 || i == cfg_.numCells) {
      // Boundary: use cell center values
      int idx = (i == 0) ? 0 : cfg_.numCells - 1;
      rho_L = rho_R = rho_[idx];
      u_L = u_R = u_[idx];
      p_L = p_R = p_[idx];
    } else {
      int i_L = i - 1;
      int i_R = i;

      // Detect high-pressure regions to avoid pressure oscillations created from MUSCL Scheme
      double p_threshold = 1000.0;    // Pa
      bool high_pressure = (p_[i_L] > p_threshold || p_[i_R] > p_threshold);

      if (high_pressure) {
        rho_L = rho_[i_L];
        u_L = u_[i_L];
        p_L = p_[i_L];
        rho_R = rho_[i_R];
        u_R = u_[i_R];
        p_R = p_[i_R];
      } else {
        // Reconstruct Left cell - right edge
        double v_minus = (i_L > 0) ? rho_[i_L - 1] : rho_[i_L];
        double v_plus = (i_L < cfg_.numCells - 1) ? rho_[i_L + 1] : rho_[i_L];
        double slope_rho = slopeLimit(v_minus, rho_[i_L], v_plus);
        rho_L = rho_[i_L] + 0.5 * slope_rho;

        v_minus = (i_L > 0) ? u_[i_L - 1] : u_[i_L];
        v_plus = (i_L < cfg_.numCells - 1) ? u_[i_L + 1] : u_[i_L];
        double slope_u = slopeLimit(v_minus, u_[i_L], v_plus);
        u_L = u_[i_L] + 0.5 * slope_u;

        v_minus = (i_L > 0) ? p_[i_L - 1] : p_[i_L];
        v_plus = (i_L < cfg_.numCells - 1) ? p_[i_L + 1] : p_[i_L];
        double slope_p = slopeLimit(v_minus, p_[i_L], v_plus);
        p_L = p_[i_L] + 0.5 * slope_p;

        rho_L = std::max(rho_L, 1e-10);
        p_L = std::max(p_L, 1e-10);

        // Reconstruct Right cell - left edge
        v_minus = (i_R > 0) ? rho_[i_R - 1] : rho_[i_R];
        v_plus = (i_R < cfg_.numCells - 1) ? rho_[i_R + 1] : rho_[i_R];
        slope_rho = slopeLimit(v_minus, rho_[i_R], v_plus);
        rho_R = rho_[i_R] + 0.5 * slope_rho;

        v_minus = (i_R > 0) ? u_[i_R - 1] : u_[i_R];
        v_plus = (i_R < cfg_.numCells - 1) ? u_[i_R + 1] : u_[i_R];
        slope_u = slopeLimit(v_minus, u_[i_R], v_plus);
        u_R = u_[i_R] + 0.5 * slope_u;

        v_minus = (i_R > 0) ? p_[i_R - 1] : p_[i_R];
        v_plus = (i_R < cfg_.numCells - 1) ? p_[i_R + 1] : p_[i_R];
        slope_p = slopeLimit(v_minus, p_[i_R], v_plus);
        p_R = p_[i_R] + 0.5 * slope_p;

        rho_R = std::max(rho_R, 1e-10);
        p_R = std::max(p_R, 1e-10); 
      }
    }

    if (flux_type_ == "HLLC") {
      solveHLLC(rho_L, u_L, p_L, rho_R, u_R, p_R, F_rho[i], F_rhou[i], F_E[i]);
    } else {
      solveEntropyStable(rho_L, u_L, p_L, rho_R, u_R, p_R, F_rho[i], F_rhou[i], F_E[i]);
    } 
  }
}

/**
 * HLLC Riemann Solver
 */
void ShockSolver::solveHLLC(double rho_L, double u_L, double p_L,
                            double rho_R, double u_R, double p_R,
                            double& F_rho, double& F_rhou, double& F_E) {
  double a_L = std::sqrt(cfg_.gamma * p_L / rho_L);
  double a_R = std::sqrt(cfg_.gamma * p_R / rho_R);

  double p_star = 0.5 * (p_L + p_R) - 0.125 * (u_R - u_L) * (rho_L + rho_R) * (a_L + a_R);
  p_star = std::max(0.0, p_star);

  double q_L = (p_star > p_L) ? 
    std::sqrt(1.0 + (cfg_.gamma + 1.0) / (2.0 * cfg_.gamma) * (p_star / p_L - 1.0)) : 1.0;
  double q_R = (p_star > p_R) ?
    std::sqrt(1.0 + (cfg_.gamma + 1.0) / (2.0 * cfg_.gamma) * (p_star / p_R - 1.0)) : 1.0;
  
  double S_L = u_L - a_L * q_L;
  double S_R = u_R + a_R * q_R;

  double S_star = (p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) /
                  (rho_L * (S_L - u_L) - rho_R * (S_R - u_R));
  
  if (S_L >= 0.0) {
    F_rho = rho_L * u_L;
    F_rhou = rho_L * u_L * u_L + p_L;
    double e_L = p_L / ((cfg_.gamma - 1.0) * rho_L);
    double E_L = rho_L * (e_L + 0.5 * u_L * u_L);
    F_E = u_L * (E_L + p_L);
  } else if (S_R <= 0.0) {
    F_rho = rho_R * u_R;
    F_rhou = rho_R * u_R * u_R + p_R;
    double e_R = p_R / ((cfg_.gamma - 1.0) * rho_R);
    double E_R = rho_R * (e_R + 0.5 * u_R * u_R);
    F_E = u_R * (E_R + p_R);
  } else if (S_star >= 0.0) {
    double rho_star_L = rho_L * (S_L - u_L) / (S_L - S_star);
    double e_L = p_L / ((cfg_.gamma - 1.0) * rho_L);
    double E_L = rho_L * (e_L + 0.5 * u_L * u_L);
    double E_star_L = rho_star_L * (E_L / rho_L + (S_star - u_L) * (S_star + p_L / (rho_L * (S_L - u_L))));

    double F_L_rho = rho_L * u_L;
    double F_L_rhou = rho_L * u_L * u_L + p_L;
    double F_L_E = u_L * (E_L + p_L);

    F_rho = F_L_rho + S_L * (rho_star_L - rho_L);
    F_rhou = F_L_rhou + S_L * (rho_star_L * S_star - rho_L * u_L);
    F_E = F_L_E + S_L * (E_star_L - E_L);
  } else {
    double rho_star_R = rho_R * (S_R - u_R) / (S_R - S_star);
    double e_R = p_R / ((cfg_.gamma - 1.0) * rho_R);
    double E_R = rho_R * (e_R + 0.5 * u_R * u_R);
    double E_star_R = rho_star_R * (E_R / rho_R + (S_star - u_R) * (S_star + p_R / (rho_R * (S_R - u_R))));

    double F_R_rho = rho_R * u_R;
    double F_R_rhou = rho_R * u_R * u_R + p_R;
    double F_R_E = u_R * (E_R + p_R);

    F_rho = F_R_rho + S_R * (rho_star_R - rho_R);
    F_rhou = F_R_rhou + S_R * (rho_star_R * S_star - rho_R * u_R);
    F_E = F_R_E + S_R * (E_star_R - E_R);
  }
}

/**
 * Pressure star helper
 */
double ShockSolver::solvePressureStar(double rho_L, double u_L, double p_L, double a_L,
                                     double rho_R, double u_R, double p_R, double a_R) const {
  double p_star = std::pow((a_L + a_R - 0.5 * (cfg_.gamma - 1.0) * (u_R - u_L)) / 
                           (a_L / std::pow(p_L, (cfg_.gamma - 1.0) / (2.0 * cfg_.gamma)) +
                            a_R / std::pow(p_R, (cfg_.gamma - 1.0) / (2.0 * cfg_.gamma))),
                            2.0 * cfg_.gamma / (cfg_.gamma - 1.0));
  p_star = std::max(TOL, p_star);

  for (int iter = 0; iter < MAX_ITER; ++iter) {
    auto f = [this](double p, double p_k, double rho_k, double a_k) -> double {
      if (p > p_k) {
        double A = 2.0 / ((cfg_.gamma + 1.0) * rho_k);
        double B = p_k * (cfg_.gamma - 1.0) / (cfg_.gamma + 1.0);
        return (p - p_k) * std::sqrt(A / (p + B));
      } else {
        return (2.0 * a_k / (cfg_.gamma - 1.0)) *
               (std::pow(p / p_k, (cfg_.gamma - 1.0) / (2.0 * cfg_.gamma)) - 1.0);
      }
    };

    auto f_prime = [this](double p, double p_k, double rho_k, double a_k) -> double {
      if (p > p_k) {
        double A = 2.0 / ((cfg_.gamma + 1.0) * rho_k);
        double B = p_k * (cfg_.gamma - 1.0) / (cfg_.gamma + 1.0);
        return std::sqrt(A / (p + B)) * (1.0 - 0.5 * (p - p_k) / (p + B));
      } else {
        return (1.0 / (rho_k * a_k)) * std::pow(p / p_k, -(cfg_.gamma + 1.0) / (2.0 * cfg_.gamma));
      }
    };

    double f_L = f(p_star, p_L, rho_L, a_L);
    double f_R = f(p_star, p_R, rho_R, a_R);
    double f_L_prime = f_prime(p_star, p_L, rho_L, a_L);
    double f_R_prime = f_prime(p_star, p_R, rho_R, a_R);

    double F = f_L + f_R + (u_R - u_L);
    double F_prime = f_L_prime + f_R_prime;

    double p_new = p_star - F / F_prime;
    p_new = std::max(TOL, p_new);

    if (std::abs(p_new - p_star) / p_star < TOL) {
      return p_new;
    }

    p_star = p_new;
  }

  return p_star;
}

/**
 * Slope Limiter (MC)
 */
double ShockSolver::slopeLimit(double v_minus, double v_center, double v_plus) const {
  double delta_minus = v_center - v_minus;
  double delta_plus = v_plus - v_center;

  // If slopes have opposite signs, we're at an extremum -> set slope to 0
  if (delta_plus * delta_minus <= 0.0) return 0.0;
  double delta_center = 0.5 * (delta_plus + delta_minus);
  double sign = (delta_center > 0.0) ? 1.0 : -1.0;
  
  return sign * std::min({2.0 * std::abs(delta_minus),
                          std::abs(delta_center),
                          2.0 * std::abs(delta_plus)});
}

/**
 * Boundary Conditions
 */
void ShockSolver::applyBoundaryConditions() {
}

/**
 * Diagnostics
 */
double ShockSolver::getTotalMass() const {
  double total = 0.0;
  for (int i = 0; i < cfg_.numCells; ++i) {
    total += rho_[i] * dx_;
  }
  return total;
}

double ShockSolver::getTotalEnergy() const {
  double total = 0.0;
  for (int i = 0; i < cfg_.numCells; ++i) {
    double e = p_[i] / ((cfg_.gamma - 1.0) * rho_[i]);
    double E = rho_[i] * (e + 0.5 * u_[i] * u_[i]);
    total += E * dx_;
  }
  return total;
}

/**
 * Entropy calculation
 */
std::vector<double> ShockSolver::getEntropy() const {
  std::vector<double> entropy(cfg_.numCells);
  
  for (int i = 0; i < cfg_.numCells; ++i) {
    double rho_safe = std::max(rho_[i], TOL);
    double p_safe = std::max(p_[i], TOL);
    
    entropy[i] = std::log(p_safe / std::pow(rho_safe, cfg_.gamma));
  }
  
  return entropy;
}

double ShockSolver::getEntropyGeneration() const {
  auto entropy = getEntropy();
  
  double maxEntropy = *std::max_element(entropy.begin(), entropy.end());
  double minEntropy = *std::min_element(entropy.begin(), entropy.end());
  
  return maxEntropy - minEntropy;
}

void ShockSolver::placeSensors(const std::vector<double>& positions) {
  sensor_positions_ = positions;
  sensor_indices_.clear();

  for (double pos : positions) {
    // Find closest cell to sensor position
    double min_dist = 1e10;
    int closest_idx = 0;

    for (int i = 0; i < cfg_.numCells; ++i) {
      double dist = std::abs(x_[i] - pos);
      if (dist < min_dist) {
        min_dist = dist;
        closest_idx = i;
      }
    }

    sensor_indices_.push_back(closest_idx);
    std::cout << "Sensor at x= " << pos << " cell " << closest_idx
              << " (actual x= " << x_[closest_idx] << ")" << std::endl;
  }
}

std::vector<ShockSolver::SensorReading> ShockSolver::getSensorReadings() const {
  std::vector<SensorReading> readings;
  auto entropy = getEntropy();

  for (size_t i = 0; i < sensor_indices_.size(); ++i) {
    int idx = sensor_indices_[i];
    
    SensorReading reading;
    reading.position = x_[idx];
    reading.density = rho_[idx];
    reading.velocity = u_[idx];
    reading.pressure = p_[idx];
    reading.entropy = entropy[idx];
    reading.time = t_;

    readings.push_back(reading);
  }

  return readings;
}

/**
 * Sparse Initialization - Reconstruct field from limited sensor data
 */

/**
 * Reference corresponding Toro test at t = 0
 */
std::vector<ShockSolver::SparseDataPoint> ShockSolver::getSensorDataForToroTest(
  ToroTest test,
  const std::vector<double>& sensor_positions) {
  
  double rho_L, u_L, p_L, rho_R, u_R, p_R, x_diaphragm;

  switch (test) {
    case ToroTest::TEST1_SOD:
      rho_L = 1.0; u_L = 0.0; p_L = 1.0;
      rho_R = 0.125; u_R = 0.0; p_R = 0.1;
      x_diaphragm = 0.5;
      break;
      
    case ToroTest::TEST2_123:
      rho_L = 1.0; u_L = -2.0; p_L = 0.4;
      rho_R = 1.0; u_R = 2.0; p_R = 0.4;
      x_diaphragm = 0.5;
      break;
      
    case ToroTest::TEST3_BLAST_LEFT:
      rho_L = 1.0; u_L = 0.0; p_L = 1000.0;
      rho_R = 1.0; u_R = 0.0; p_R = 0.01;
      x_diaphragm = 0.5;
      break;
      
    case ToroTest::TEST4_SLOW_SHOCK:
      rho_L = 1.0; u_L = 0.0; p_L = 0.01;
      rho_R = 1.0; u_R = 0.0; p_R = 100.0;
      x_diaphragm = 0.5;
      break;
      
    case ToroTest::TEST5_COLLISION:
      rho_L = 5.99924; u_L = 19.5975; p_L = 460.894;
      rho_R = 5.99242; u_R = -6.19633; p_R = 46.0950;
      x_diaphragm = 0.5;
      break;
      
    default:
      rho_L = 1.0; u_L = 0.0; p_L = 1.0;
      rho_R = 0.125; u_R = 0.0; p_R = 0.1;
      x_diaphragm = 0.5;
  }

  std::vector<SparseDataPoint> sensor_data;

  for (double x : sensor_positions) {
    SparseDataPoint point;
    point.x = x;

    if (x < x_diaphragm) {
      point.rho = rho_L;
      point.u = u_L;
      point.p = p_L;
    } else {
      point.rho = rho_R;
      point.u = u_R;
      point.p = p_R;
    }

    sensor_data.push_back(point);
  }

  return sensor_data;
}

/**
 * Initialize from sparse sensor data using Interp methods
 */
void ShockSolver::initializeFromSparseData(
  const std::vector<SparseDataPoint>& sparse_data,
  const std::string& interpolation_method) {
  
  if (sparse_data.empty()) {
    throw std::runtime_error("No sparse data provided");
  } 

  // Sort by position
  std::vector<SparseDataPoint> sorted_data = sparse_data;
  std::sort(sorted_data.begin(), sorted_data.end(),
            [](const SparseDataPoint& a, const SparseDataPoint& b) {
              return a.x < b.x;
            });
  
  // Store for validation - estimate diaphragm position from data
  // Find where biggest density jump is
  double max_jump = 0.0;
  double estimated_diaphragm = 0.5;
  for (size_t i = 0; i < sorted_data.size() - 1; ++i) {
    // Normalize jumps by magnitude
    double rho_avg = 0.5 * (sorted_data[i].rho + sorted_data[i + 1].rho);
    double u_avg = 0.5 * (sorted_data[i].u + sorted_data[i + 1].u);
    double p_avg = 0.5 * (sorted_data[i].p + sorted_data[i + 1].p);

    double rho_jump = std::abs(sorted_data[i + 1].rho - sorted_data[i].rho) / (rho_avg + 1e-10);
    double u_jump = std::abs(sorted_data[i + 1].u - sorted_data[i].u) / (u_avg + 1e-10);
    double p_jump = std::abs(sorted_data[i + 1].p - sorted_data[i].p) / (p_avg + 1e-10);

    // Use maximum normalized jump across all variables
    double total_jump = std::max({rho_jump, u_jump, p_jump});

    if (total_jump > max_jump) {
      max_jump = total_jump;
      estimated_diaphragm = 0.5 * (sorted_data[i].x + sorted_data[i + 1].x);
    }
  }

  estimated_diaphragm = 0.5;

  // Store initial conditions for analytical solution comparison
  // Use first and last sensor as left/right states
  initial_conditions_.rho_L = sorted_data.front().rho;
  initial_conditions_.rho_R = sorted_data.back().rho;
  initial_conditions_.u_L = sorted_data.front().u;
  initial_conditions_.u_R = sorted_data.back().u;
  initial_conditions_.p_L = sorted_data.front().p;
  initial_conditions_.p_R = sorted_data.back().p;
  initial_conditions_.x_diaphragm = estimated_diaphragm;
  initial_conditions_.is_valid = true;

  std::cout << "Initializing from " << sparse_data.size() << " sparse data points" << std::endl;
  std::cout << "Interpolation Method: " << interpolation_method << std::endl;
  std::cout << "Estimated diaphragm position: " << estimated_diaphragm << std::endl;

  // Interpolate to fill cells
  for (int i = 0; i < cfg_.numCells; ++i) {
    double x = x_[i];

    // Find bracketing points
    if (x <= sorted_data.front().x) {
      // Extrapolate from first point
      rho_[i] = sorted_data.front().rho;
      u_[i] = sorted_data.front().u;
      p_[i] = sorted_data.front().p;
    } else if (x >= sorted_data.back().x) {
      // Extrapolate from last point
      rho_[i] = sorted_data.back().rho;
      u_[i] = sorted_data.back().u;
      p_[i] = sorted_data.back().p;
    } else {
      // Find bracketing sensors
      size_t left_idx = 0;
      for (size_t j = 0; j < sorted_data.size() - 1; ++j) {
        if (x >= sorted_data[j].x && x < sorted_data[j + 1].x) {
          left_idx = j;
          break;
        }
      }

      const SparseDataPoint& left = sorted_data[left_idx];
      const SparseDataPoint& right = sorted_data[left_idx + 1];

      if (interpolation_method == "piecewise_constant") {
        // Use nearest neighbor (left point for tie)
        double mid = 0.5 * (left.x + right.x);
        if (x < mid) {
          rho_[i] = left.rho;
          u_[i] = left.u;
          p_[i] = left.p; 
        } else {
          rho_[i] = right.rho;
          u_[i] = right.u;
          p_[i] = right.p;
        }
      } else {
        // Linear interpolation (default)
        double t = (x - left.x) / (right.x - left.x);
        rho_[i] = left.rho + t * (right.rho - left.rho);
        u_[i] = left.u + t * (right.u - left.u);
        p_[i] = left.p + t * (right.p - left.p);
      }
    }
  }

  primitivesToConservative();

  std::cout << "Sparse initialization complete" << std::endl;
}



/**
 * VALIDATION FRAMEWORK - Analytical Solution Comparison
 * 
 * Based on Toro's "Riemann Solvers and Numerical Methods for Fluid Dynamics"
 * Chapter 4: The Riemann Problem for the Euler Equations
 */

void ShockSolver::computeAnalyticalSolution(std::vector<double>& rho_exact,
                                             std::vector<double>& u_exact,
                                             std::vector<double>& p_exact) const {
  if (!initial_conditions_.is_valid) {
    throw std::runtime_error("Cannot compute analytical solution: initial conditions not set");
  }

  // Resize output vectors
  rho_exact.resize(cfg_.numCells);
  u_exact.resize(cfg_.numCells);
  p_exact.resize(cfg_.numCells);

  // Extract initial conditions
  double rho_L = initial_conditions_.rho_L;
  double u_L = initial_conditions_.u_L;
  double p_L = initial_conditions_.p_L;
  double rho_R = initial_conditions_.rho_R;
  double u_R = initial_conditions_.u_R;
  double p_R = initial_conditions_.p_R;
  double x_diaphragm = initial_conditions_.x_diaphragm;

  // Compute sound speeds
  double a_L = std::sqrt(cfg_.gamma * p_L / rho_L);
  double a_R = std::sqrt(cfg_.gamma * p_R / rho_R);

  // Solve for star region (p*, u*)
  double p_star = solvePressureStar(rho_L, u_L, p_L, a_L, rho_R, u_R, p_R, a_R);

  // Compute u_star using f functions
  auto f = [this](double p, double p_k, double rho_k, double a_k) -> double {
    if (p > p_k) {
      double A = 2.0 / ((cfg_.gamma + 1.0) * rho_k);
      double B = p_k * (cfg_.gamma - 1.0) / (cfg_.gamma + 1.0);
      return (p - p_k) * std::sqrt(A / (p + B));
    } else {
      return (2.0 * a_k / (cfg_.gamma - 1.0)) * 
             (std::pow(p / p_k, (cfg_.gamma - 1.0) / (2.0 * cfg_.gamma)) - 1.0);
    }
  };

  double f_L = f(p_star, p_L, rho_L, a_L);
  double f_R = f(p_star, p_R, rho_R, a_R);
  double u_star = 0.5 * (u_L + u_R) + 0.5 * (f_R - f_L);

  // Sample exact solution at each grid point
  for (int i = 0; i < cfg_.numCells; ++i) {
    double S;
    if (t_ > TOL) {
      S = (x_[i] - x_diaphragm) / t_;
    } else {
      S = (x_[i] < x_diaphragm) ? -1e10 : 1e10;
    }

    double rho_sample, u_sample, p_sample;

    if (S <= u_star) {
      if (p_star > p_L) {
        // Left shock
        double S_L = u_L - a_L * std::sqrt((cfg_.gamma + 1.0) / (2.0 * cfg_.gamma) * (p_star / p_L) +
                                            (cfg_.gamma - 1.0) / (2.0 * cfg_.gamma));
        if (S <= S_L) {
          rho_sample = rho_L;
          u_sample = u_L;
          p_sample = p_L;
        } else {
          rho_sample = rho_L * ((p_star / p_L + (cfg_.gamma - 1.0) / (cfg_.gamma + 1.0)) /
                               ((cfg_.gamma - 1.0) / (cfg_.gamma + 1.0) * (p_star / p_L) + 1.0));
          u_sample = u_star;
          p_sample = p_star;                     
        }                                     
      } else {
        // Left rarefaction
        double S_HL = u_L - a_L;
        if (S <= S_HL) {
          rho_sample = rho_L;
          u_sample = u_L;
          p_sample = p_L;
        } else {
          double a_star_L = a_L * std::pow(p_star / p_L, (cfg_.gamma - 1.0) / (2.0 * cfg_.gamma));
          double S_TL = u_star - a_star_L;
          if (S > S_TL) {
            rho_sample = rho_L * std::pow(p_star / p_L, 1.0 / cfg_.gamma);
            u_sample = u_star;
            p_sample = p_star;
          } else {
            // Inside rarefaction fan
            u_sample = 2.0 / (cfg_.gamma + 1.0) * (a_L + 0.5 * (cfg_.gamma - 1.0) * u_L + S);
            double c = 2.0 / (cfg_.gamma + 1.0) * (a_L + 0.5 * (cfg_.gamma - 1.0) * (u_L - S));
            rho_sample = rho_L * std::pow(c / a_L, 2.0 / (cfg_.gamma - 1.0));
            p_sample = p_L * std::pow(c / a_L, 2.0 * cfg_.gamma / (cfg_.gamma - 1.0));
          }
        }
      }
    } else {
      // Right side
      if (p_star > p_R) {
        // Right shock
        double S_R = u_R + a_R * std::sqrt((cfg_.gamma + 1.0) / (2.0 * cfg_.gamma) * (p_star / p_R) + 
                                            (cfg_.gamma - 1.0) / (2.0 * cfg_.gamma));
        if (S >= S_R) {
          rho_sample = rho_R;
          u_sample = u_R;
          p_sample = p_R;
        } else {
          rho_sample = rho_R * ((p_star / p_R + (cfg_.gamma - 1.0) / (cfg_.gamma + 1.0))/
                               ((cfg_.gamma - 1.0) / (cfg_.gamma + 1.0) * (p_star / p_R) + 1.0));
          u_sample = u_star;
          p_sample = p_star;
        }
      } else {
        // Right rarefaction
        double S_HR = u_R + a_R;
        if (S >= S_HR) {
          rho_sample = rho_R;
          u_sample = u_R;
          p_sample = p_R;
        } else {
          double a_star_R = a_R * std::pow(p_star / p_R, (cfg_.gamma - 1.0) / (2.0 * cfg_.gamma));
          double S_TR = u_star + a_star_R;
          if (S < S_TR) {
            rho_sample = rho_R * std::pow(p_star / p_R, 1.0 / cfg_.gamma);
            u_sample = u_star;
            p_sample = p_star;
          } else {
            // Inside rarefaction fan
            u_sample = 2.0 / (cfg_.gamma + 1.0) * (-a_R + 0.5 * (cfg_.gamma - 1.0) * u_R + S);
            double c = 2.0 / (cfg_.gamma + 1.0) * (a_R - 0.5 * (cfg_.gamma - 1.0) * (u_R - S));
            rho_sample = rho_R * std::pow(c / a_R, 2.0 / (cfg_.gamma - 1.0));
            p_sample = p_R * std::pow(c / a_R, 2.0 * cfg_.gamma / (cfg_.gamma - 1.0));
          }
        }
      }
    }

    rho_exact[i] = rho_sample;
    u_exact[i] = u_sample;
    p_exact[i] = p_sample;
  }
}

/**
 * Validate numerical solution against analytical solution
 * 
 * Compute L1, L2, and Linf error norms:
 * L1 (mean abs error): average error magnitude
 * L2 (RMS error): root-mean-square error
 * Linf (max error): worst-case pointwise error
 * 
 */
ShockSolver::ValidationMetrics ShockSolver::validateAgainstAnalytical() const {
  ValidationMetrics metrics;

  // Compute analytical solution
  std::vector<double> rho_exact, u_exact, p_exact;
  computeAnalyticalSolution(rho_exact, u_exact, p_exact);

  // Initialize error accumulators
  double sum_rho_L1 = 0.0, sum_rho_L2 = 0.0, max_rho = 0.0;
  double sum_u_L1 = 0.0, sum_u_L2 = 0.0, max_u = 0.0;
  double sum_p_L1 = 0.0, sum_p_L2 = 0.0, max_p = 0.0;
  double sum_s_L1 = 0.0, sum_s_L2 = 0.0, max_s = 0.0;

  // Compute errors for each grid point
  for (int i = 0; i < cfg_.numCells; ++i) {
    // Absolute errors
    double err_rho = std::abs(rho_[i] - rho_exact[i]);
    double err_u = std::abs(u_[i] - u_exact[i]);
    double err_p = std::abs(p_[i] - p_exact[i]);
    
    // Entropy: s = ln(p / rho^gamma)
    double s_num = std::log(p_[i] / std::pow(rho_[i], cfg_.gamma));
    double s_exact = std::log(p_exact[i] / std::pow(rho_exact[i], cfg_.gamma));
    double err_s = std::abs(s_num - s_exact);

    // L1 norm (sum of abs errors)
    sum_rho_L1 += err_rho;
    sum_u_L1 += err_u;
    sum_p_L1 += err_p;
    sum_s_L1 += err_s;

    // L2 norm (sum of squared errors)
    sum_rho_L2 += err_rho * err_rho;
    sum_u_L2 += err_u * err_u;
    sum_p_L2 += err_p * err_p;
    sum_s_L2 += err_s * err_s;

    // Linf norm (max absolute error)
    if (err_rho > max_rho) { max_rho = err_rho; }
    if (err_u > max_u) { max_u = err_u; }
    if (err_p > max_p) { max_p = err_p; }
    if (err_s > max_s) { max_s = err_s; }
  }

  // Normalize by number of cells 
  int N = cfg_.numCells;
  metrics.L1_error_density = sum_rho_L1 / N;
  metrics.L1_error_velocity = sum_u_L1 / N;
  metrics.L1_error_pressure = sum_p_L1 / N;
  metrics.L1_error_entropy = sum_s_L1 / N;

  metrics.L2_error_density = std::sqrt(sum_rho_L2 / N);
  metrics.L2_error_velocity = std::sqrt(sum_u_L2 / N);
  metrics.L2_error_pressure = std::sqrt(sum_p_L2 / N);
  metrics.L2_error_entropy = std::sqrt(sum_s_L2 / N);

  metrics.Linf_error_density = max_rho;
  metrics.Linf_error_velocity = max_u;
  metrics.Linf_error_pressure = max_p;
  metrics.Linf_error_entropy = max_s;

  metrics.num_points_validated = N;

  return metrics;
}