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
    scheme_type_("godunov"),
    solver_type_("hllc") {
  
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
                                      double x_diaphragm) {

  initial_conditions_.rho_L = rho_L;
  initial_conditions_.u_L = u_L;
  initial_conditions_.p_L = p_L;
  initial_conditions_.rho_R = rho_R;
  initial_conditions_.u_R = u_R;
  initial_conditions_.p_R = p_R;
  initial_conditions_.x_diaphragm = x_diaphragm;
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
  double rho_L, u_L, p_L, rho_R, u_R, p_R, x_diaphragm;
  
  switch (test) {
    case ToroTest::TEST1_SOD:
      // Test 1: Sod's problem (standard)
      rho_L = 1.0; u_L = 0.0; p_L = 1.0;
      rho_R = 0.125; u_R = 0.0; p_R = 0.1;
      x_diaphragm = 0.5;
      break;
      
    case ToroTest::TEST2_123:
      // Test 2: 123 problem (strong rarefaction)
      rho_L = 1.0; u_L = -2.0; p_L = 0.4;
      rho_R = 1.0; u_R = 2.0; p_R = 0.4;
      x_diaphragm = 0.5;
      break;
      
    case ToroTest::TEST3_BLAST_LEFT:
      // Test 3: Left half of blast wave problem
      rho_L = 1.0; u_L = 0.0; p_L = 1000.0;
      rho_R = 1.0; u_R = 0.0; p_R = 0.01;
      x_diaphragm = 0.5;
      break;
      
    case ToroTest::TEST4_COLLISION:
      // Test 4: Collision of two rarefaction waves
      rho_L = 5.99924; u_L = 19.5975; p_L = 460.894;
      rho_R = 5.99242; u_R = -6.19633; p_R = 46.0950;
      x_diaphragm = 0.4;
      break;
      
    case ToroTest::TEST5_STATIONARY:
      // Test 5: Stationary contact discontinuity
      rho_L = 1.0; u_L = -19.59745; p_L = 1000.0;
      rho_R = 1.0; u_R = -19.59745; p_R = 0.01;
      x_diaphragm = 0.8;
      break;
  }
  
  std::cout << "Initializing Toro Test " << static_cast<int>(test) + 1 << std::endl;
  initializeShockTube(rho_L, u_L, p_L, rho_R, u_R, p_R, x_diaphragm);
}


/**
 * Setters
 */
void ShockSolver::setScheme(const std::string& scheme_name) {
  if (scheme_name != "godunov" && scheme_name != "muscl") {
    throw std::runtime_error("Unknown scheme: " + scheme_name);
  }
  scheme_type_ = scheme_name;
}

void ShockSolver::setRiemannSolver(const std::string& solver_name) {
  if (solver_name != "exact" && solver_name != "hllc") {
    throw std::runtime_error("Unknown Riemann solver: " + solver_name);
  }
  solver_type_ = solver_name;
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
  
  // MUSCL needs lower CFL for stability (second-order extrapolation)
  double effective_cfl = cfg_.CFL;
  if (scheme_type_ == "muscl") {
    effective_cfl = std::min(cfg_.CFL, 0.5);
  }
  
  return effective_cfl * dx_ / max_speed;
}

void ShockSolver::step() {
  double dt = computeTimestep();
  if (t_ + dt > cfg_.endTime) {
    dt = cfg_.endTime - t_;
  }

  applyBoundaryConditions();

  std::vector<double> F_rho(cfg_.numCells + 1);
  std::vector<double> F_rhou(cfg_.numCells + 1);
  std::vector<double> F_E(cfg_.numCells + 1);

  computeFluxes(F_rho, F_rhou, F_E);
  updateSolution(F_rho, F_rhou, F_E, dt);
  conservativesToPrimitives();

  t_ += dt;
  step_count_++;
}

void ShockSolver::run() {
  std::cout << "Running ShockSolver with " << scheme_type_
            << " scheme and " << solver_type_ << " Riemann solver\n";
  std::cout << "Grid: " << cfg_.numCells << " cells, dx = " << dx_ << " m\n";
  
  double effective_cfl = cfg_.CFL;
  if (scheme_type_ == "muscl") {
    effective_cfl = std::min(cfg_.CFL, 0.5);
  }
  std::cout << "CFL = " << effective_cfl;
  if (scheme_type_ == "muscl" && cfg_.CFL > 0.5) {
    std::cout << " (capped from " << cfg_.CFL << " for MUSCL stability)";
  }
  std::cout << ", gamma = " << cfg_.gamma << "\n\n";

  while (t_ < cfg_.endTime) {
    step();
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
void ShockSolver::computeFluxes(std::vector<double>& F_rho,
                                std::vector<double>& F_rhou,
                                std::vector<double>& F_E) {
  if (scheme_type_ == "godunov") {
    for (int i = 0; i <= cfg_.numCells; ++i) {
      double rho_L, u_L, p_L, rho_R, u_R, p_R;
      
      if (i == 0) {
        rho_L = rho_[0];
        u_L = u_[0];
        p_L = p_[0];
        rho_R = rho_[0];
        u_R = u_[0];
        p_R = p_[0];
      } else if (i == cfg_.numCells) {
        rho_L = rho_[cfg_.numCells - 1];
        u_L = u_[cfg_.numCells - 1];
        p_L = p_[cfg_.numCells - 1];
        rho_R = rho_[cfg_.numCells - 1];
        u_R = u_[cfg_.numCells - 1];
        p_R = p_[cfg_.numCells - 1];
      } else {
        rho_L = rho_[i - 1];
        u_L = u_[i - 1];
        p_L = p_[i - 1];
        rho_R = rho_[i];
        u_R = u_[i];
        p_R = p_[i];
      }
      
      solveRiemann(rho_L, u_L, p_L, rho_R, u_R, p_R,
                  F_rho[i], F_rhou[i], F_E[i]);
    }
  } else if (scheme_type_ == "muscl") {
    // MUSCL-Hancock: reconstruct, then evolve by dt/2, then solve Riemann
    double dt = computeTimestep();
    if (t_ + dt > cfg_.endTime) {
      dt = cfg_.endTime - t_;
    }
    double dtdx = dt / dx_;
    
    for (int i = 0; i <= cfg_.numCells; ++i) {
      double rho_L, u_L, p_L, rho_R, u_R, p_R;

      if (i == 0) {
        rho_L = rho_[0];
        u_L = u_[0];
        p_L = p_[0];
        rho_R = rho_[0];
        u_R = u_[0];
        p_R = p_[0];
      } else if (i == cfg_.numCells) {
        rho_L = rho_[cfg_.numCells - 1];
        u_L = u_[cfg_.numCells - 1];
        p_L = p_[cfg_.numCells - 1];
        rho_R = rho_[cfg_.numCells - 1];
        u_R = u_[cfg_.numCells - 1];
        p_R = p_[cfg_.numCells - 1];
      } else {
        // Interface i sits between cell i-1 and cell i
        // Left state = RIGHT edge of cell i-1
        // Right state = LEFT edge of cell i
        double rho_L_left, u_L_left, p_L_left;   // left edge of cell i-1
        double rho_R_right, u_R_right, p_R_right; // right edge of cell i
        
        // Get both edges of cell i-1
        reconstructMUSCL(i - 1, rho_L_left, u_L_left, p_L_left, rho_L, u_L, p_L);
        
        // Get both edges of cell i
        reconstructMUSCL(i, rho_R, u_R, p_R, rho_R_right, u_R_right, p_R_right);
        
        // MUSCL-Hancock predictor: evolve edge states by dt/2
        // Using primitive variable formulation - compute all changes first, then apply
        
        // For cell i-1: evolve the RIGHT edge state
        {
          double drho = rho_L - rho_L_left;  // slope across cell i-1
          double du = u_L - u_L_left;
          double dp = p_L - p_L_left;
          
          // Compute changes (using original values)
          double d_rho = -0.5 * dtdx * (u_L * drho + rho_L * du);
          double d_u = -0.5 * dtdx * (u_L * du + dp / rho_L);
          double d_p = -0.5 * dtdx * (u_L * dp + cfg_.gamma * p_L * du);
          
          // Apply all at once
          rho_L += d_rho;
          u_L += d_u;
          p_L += d_p;
          
          // Ensure positivity
          rho_L = std::max(rho_L, 1e-10);
          p_L = std::max(p_L, 1e-10);
        }
        
        // For cell i: evolve the LEFT edge state
        {
          double drho = rho_R_right - rho_R;  // slope across cell i
          double du = u_R_right - u_R;
          double dp = p_R_right - p_R;
          
          // Compute changes (using original values)
          double d_rho = -0.5 * dtdx * (u_R * drho + rho_R * du);
          double d_u = -0.5 * dtdx * (u_R * du + dp / rho_R);
          double d_p = -0.5 * dtdx * (u_R * dp + cfg_.gamma * p_R * du);
          
          // Apply all at once
          rho_R += d_rho;
          u_R += d_u;
          p_R += d_p;
          
          // Ensure positivity
          rho_R = std::max(rho_R, 1e-10);
          p_R = std::max(p_R, 1e-10);
        }
      }

      solveRiemann(rho_L, u_L, p_L, rho_R, u_R, p_R,
                  F_rho[i], F_rhou[i], F_E[i]);
    }
  }
}

/**
 * Riemann Solver Dispatch
 */
void ShockSolver::solveRiemann(double rho_L, double u_L, double p_L,
                               double rho_R, double u_R, double p_R,
                               double& F_rho, double& F_rhou, double& F_E) {
  if (solver_type_ == "hllc") {
    solveHLLC(rho_L, u_L, p_L, rho_R, u_R, p_R, F_rho, F_rhou, F_E);
  } else if (solver_type_ == "exact") {
    solveExact(rho_L, u_L, p_L, rho_R, u_R, p_R, F_rho, F_rhou, F_E);
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
 * Exact Riemann Solver
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

void ShockSolver::sampleExactSolution(double rho_L, double u_L, double p_L, double a_L,
                                      double rho_R, double u_R, double p_R, double a_R,
                                      double p_star, double u_star,
                                      double& rho_sample, double& u_sample, double& p_sample) const {
  double S = 0.0;

  if (S <= u_star) {
    if (p_star > p_L) {
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
          u_sample = 2.0 / (cfg_.gamma + 1.0) * (a_L + 0.5 * (cfg_.gamma - 1.0) * u_L + S);
          double c = 2.0 / (cfg_.gamma + 1.0) * (a_L + 0.5 * (cfg_.gamma - 1.0) * (u_L - S));
          rho_sample = rho_L * std::pow(c / a_L, 2.0 / (cfg_.gamma - 1.0));
          p_sample = p_L * std::pow(c / a_L, 2.0 * cfg_.gamma / (cfg_.gamma - 1.0));
        }
      }
    }
  } else {
    if (p_star > p_R) {
      double S_R = u_R + a_R * std::sqrt((cfg_.gamma + 1.0) / (2.0 * cfg_.gamma) * (p_star / p_R) + 
                                         (cfg_.gamma - 1.0) / (2.0 * cfg_.gamma));
      if (S >= S_R) {
        rho_sample = rho_R;
        u_sample = u_R;
        p_sample = p_R;
      } else {
        rho_sample = rho_R * ((p_star / p_R + (cfg_.gamma - 1.0) / (cfg_.gamma + 1.0)) /
                             ((cfg_.gamma - 1.0) / (cfg_.gamma + 1.0) * (p_star / p_R) + 1.0));
        u_sample = u_star;
        p_sample = p_star;
      }
    } else {
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
          u_sample = 2.0 / (cfg_.gamma + 1.0) * (-a_R + 0.5 * (cfg_.gamma - 1.0) * u_R + S);
          double c = 2.0 / (cfg_.gamma + 1.0) * (a_R - 0.5 * (cfg_.gamma - 1.0) * (u_R - S));
          rho_sample = rho_R * std::pow(c / a_R, 2.0 / (cfg_.gamma - 1.0));
          p_sample = p_R * std::pow(c / a_R, 2.0 * cfg_.gamma / (cfg_.gamma - 1.0));
        }
      }
    }
  }
}

void ShockSolver::solveExact(double rho_L, double u_L, double p_L,
                             double rho_R, double u_R, double p_R,
                             double& F_rho, double& F_rhou, double& F_E) {
  double a_L = std::sqrt(cfg_.gamma * p_L / rho_L);
  double a_R = std::sqrt(cfg_.gamma * p_R / rho_R);

  double p_star = solvePressureStar(rho_L, u_L, p_L, a_L, rho_R, u_R, p_R, a_R);

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

  double rho_sample, u_sample, p_sample;
  sampleExactSolution(rho_L, u_L, p_L, a_L, rho_R, u_R, p_R, a_R,
                      p_star, u_star, rho_sample, u_sample, p_sample);
  
  F_rho = rho_sample * u_sample;
  F_rhou = rho_sample * u_sample * u_sample + p_sample;
  double e_sample = p_sample / ((cfg_.gamma - 1.0) * rho_sample);
  double E_sample = rho_sample * (e_sample + 0.5 * u_sample * u_sample);
  F_E = u_sample * (E_sample + p_sample);
}

/**
 * MUSCL Reconstruction
 */
void ShockSolver::reconstructMUSCL(int i,
                                    double& rho_L, double& u_L, double& p_L,
                                    double& rho_R, double& u_R, double& p_R) {
  double slope_rho = slopeLimit(
    i > 0 ? rho_[i - 1] : rho_[i],
    rho_[i],
    i < cfg_.numCells - 1 ? rho_[i + 1] : rho_[i]
  );
  
  double slope_u = slopeLimit(
    i > 0 ? u_[i - 1] : u_[i],
    u_[i],
    i < cfg_.numCells - 1 ? u_[i + 1] : u_[i]
  );

  double slope_p = slopeLimit(
    i > 0 ? p_[i - 1] : p_[i],
    p_[i],
    i < cfg_.numCells - 1 ? p_[i + 1] : p_[i]
  );

  rho_L = rho_[i] - 0.5 * slope_rho;
  u_L = u_[i] - 0.5 * slope_u;
  p_L = p_[i] - 0.5 * slope_p;

  rho_R = rho_[i] + 0.5 * slope_rho;
  u_R = u_[i] + 0.5 * slope_u;
  p_R = p_[i] + 0.5 * slope_p;

  rho_L = std::max(rho_L, TOL);
  rho_R = std::max(rho_R, TOL);
  p_L = std::max(p_L, TOL);
  p_R = std::max(p_R, TOL);
}

double ShockSolver::slopeLimit(double v_minus, double v_center, double v_plus) const {
  double delta_minus = v_center - v_minus;
  double delta_plus = v_plus - v_center;

  if (delta_plus * delta_minus <= 0.0) {
    return 0.0;
  }

  double delta_center = 0.5 * (v_plus - v_minus);

  double sign = (delta_center > 0.0) ? 1.0 : -1.0;
  return sign * std::min({std::abs(delta_center), 
                          2.0 * std::abs(delta_minus),
                          2.0 * std::abs(delta_plus)});
}

/**
 * Boundary Conditions
 */
void ShockSolver::applyBoundaryConditions() {
}

/**
 * Update Solution
 */
void ShockSolver::updateSolution(const std::vector<double>& F_rho,
                                  const std::vector<double>& F_rhou,
                                  const std::vector<double>& F_E,
                                  double dt) {
  for (int i = 0; i < cfg_.numCells; ++i) {
    rho_[i] = rho_[i] - (dt / dx_) * (F_rho[i + 1] - F_rho[i]);
    rho_u_[i] = rho_u_[i] - (dt / dx_) * (F_rhou[i + 1] - F_rhou[i]);
    E_[i] = E_[i] - (dt / dx_) * (F_E[i + 1] - F_E[i]);

    rho_[i] = std::max(rho_[i], TOL);
    E_[i] = std::max(E_[i], TOL);
  }
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
      
    case ToroTest::TEST4_COLLISION:
      rho_L = 5.99924; u_L = 19.5975; p_L = 460.894;
      rho_R = 5.99242; u_R = -6.19633; p_R = 46.0950;
      x_diaphragm = 0.4;
      break;
      
    case ToroTest::TEST5_STATIONARY:
      rho_L = 1.0; u_L = -19.59745; p_L = 1000.0;
      rho_R = 1.0; u_R = -19.59745; p_R = 0.01;
      x_diaphragm = 0.8;
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
    double jump = std::abs(sorted_data[i + 1].rho - sorted_data[i].rho);
    if (jump > max_jump) {
      max_jump = jump;
      estimated_diaphragm = 0.5 * (sorted_data[i].x + sorted_data[i + 1].x);
    }
  }

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

  // Compute errors for each grid point
  for (int i = 0; i < cfg_.numCells; ++i) {
    // Absolute errors
    double err_rho = std::abs(rho_[i] - rho_exact[i]);
    double err_u = std::abs(u_[i] - u_exact[i]);
    double err_p = std::abs(p_[i] - p_exact[i]);

    // L1 norm (sum of abs errors)
    sum_rho_L1 += err_rho;
    sum_u_L1 += err_u;
    sum_p_L1 += err_p;

    // L2 norm (sum of squared errors)
    sum_rho_L2 += err_rho * err_rho;
    sum_u_L2 += err_u * err_u;
    sum_p_L2 += err_p * err_p;

    // Linf norm (max absolute error)
    if (err_rho > max_rho) { max_rho = err_rho; }
    if (err_u > max_u) { max_u = err_u; }
    if (err_p > max_p) { max_p = err_p; }
  }

  // Normalize by number of cells 
  int N = cfg_.numCells;
  metrics.L1_error_density = sum_rho_L1 / N;
  metrics.L1_error_velocity = sum_u_L1 / N;
  metrics.L1_error_pressure = sum_p_L1 / N;

  metrics.L2_error_density = std::sqrt(sum_rho_L2 / N);
  metrics.L2_error_velocity = std::sqrt(sum_u_L2 / N);
  metrics.L2_error_pressure = std::sqrt(sum_p_L2 / N);

  metrics.Linf_error_density = max_rho;
  metrics.Linf_error_velocity = max_u;
  metrics.Linf_error_pressure = max_p;

  metrics.num_points_validated = N;

  return metrics;
}