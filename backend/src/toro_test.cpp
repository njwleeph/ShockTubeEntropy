/**
 * TORO TEST 2
 * 
 * NEW VERSION:
 * - Stripped old MUSCL schemes, HLLC, and RK2 methods. 
 * - Pull from repo to access old methods or reference old solver.cpp
 * - Clean implemented entropy-conservative method
 */

#include "solver.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>

/**
 * Configuration (Toro Test 2: 123 Problem)
 */
const double rho_L_config = 1.0, u_L_config = -2.0, p_L_config = 0.4;
const double rho_R_config = 1.0, u_R_config = 2.0, p_R_config = 0.4;
const double x_diapgragm_config = 0.5;
const double endTime = 0.2;
const double gamma_config = 1.4;
const double CFL = 0.5;

/**
 * Grid Structure
 */
struct Grid {
  int n;              // Grid cells
  double length;      // Grid length
  double dx;          // Grid spacing
  double t;           // Time [s]

  std::vector<double> x;        // Position
  std::vector<double> rho;      // Density
  std::vector<double> u;        // Velocity
  std::vector<double> p;        // Pressure
  std::vector<double> rho_u;    // Conserved momentum
  std::vector<double> E;        // Conserved Energy

  Grid(int numCells, double L) : n(numCells), length(L), dx(L / numCells), t(0.0) {
    x.resize(n);
    rho.resize(n);
    u.resize(n);
    p.resize(n);
    rho_u.resize(n);
    E.resize(n);

    for (int i = 0; i < n; ++i) { x[i] = (i + 0.5) * dx; }
  }

  void initialize() {
    for (int i = 0; i < n; ++i) {
      if (x[i] < x_diapgragm_config) {
        rho[i] = rho_L_config;
        u[i] = u_L_config;
        p[i] = p_L_config;
      } else {
        rho[i] = rho_R_config;
        u[i] = u_R_config;
        p[i] = p_R_config;
      }
    }
    primToConserved();
  }

  void primToConserved() {
    for (int i = 0; i < n; ++i) {
      rho_u[i] = rho[i] * u[i];
      double e = p[i] / ((gamma_config - 1.0) * rho[i]);
      E[i] = rho[i] * (e + 0.5 * u[i] * u[i]);
    }
  }

  void consToPrim() {
    for (int i = 0; i < n; ++i) {
      u[i] = rho_u[i] / rho[i];
      p[i] = (gamma_config - 1.0) * (E[i] - 0.5 * rho[i] * u[i] * u[i]);
      p[i] = std::max(p[i], 1e-10);
    }
  }
};

/**
 * MC Slope Limiter
 */
double slopeLimit(double v_minus, double v_center, double v_plus) {
  double delta_minus = v_center - v_minus;
  double delta_plus = v_plus - v_center;

  if (delta_plus * delta_minus <= 0.0) return 0.0;  // At extremum

  // MC (Monotized Central) Limiter
  double delta_center = 0.5 * (delta_plus + delta_minus);
  double sign = (delta_center > 0.0) ? 1.0 : -1.0;
  return sign * std::min({2.0 * std::abs(delta_minus),
                          2.0 * std::abs(delta_plus),
                          std::abs(delta_center)});
}

/**
 * Entropy Conservative Flux (Chandrashekar 2013)
 */
double logMean(double a, double b) {
  double xi = b / a;
  double f = (xi - 1.0) / (xi + 1.0);
  double u = f * f;

  // Handle edge case a = b w/ Taylor series
  if (u < 1e-2) {
    double F = 1.0 * u / 3.0 + u * u / 5.0 + u * u * u / 7.0;
    return 0.5 * (a + b) / F;
  }

  return (b - a) / std::log(xi);
}

void solveEntropyConservative(double rho_L, double u_L, double p_L,
                              double rho_R, double u_R, double p_R,
                              double& F_rho, double& F_rhou, double& F_E) {
  // Logarithmic means (preserve Entropy)
  double rho_log = logMean(rho_L, rho_R);
  double p_log = logMean(p_L, p_R);

  // Arithmetic averages
  double u_avg = 0.5 * (u_L + u_R);

  // Specific internal energy
  double e_L = p_L / ((gamma_config - 1.0) * rho_L);
  double e_R = p_R / ((gamma_config - 1.0) * rho_R);
  double e_avg = 0.5 * (e_L + e_R);

  // Kinetic energy
  double ke_avg = 0.5 * u_avg * u_avg;

  // Entropy-conservative flux
  F_rho = rho_log * u_avg;
  F_rhou = rho_log * u_avg * u_avg + p_log;
  F_E = rho_log * u_avg * (e_avg + ke_avg) + p_log * u_avg;
}

/**
 * MUSCL - RK2 Spatial Reconstruction
 */
void computeFluxesMUSCL(Grid& grid, std::vector<double>& F_rho,
                        std::vector<double>& F_rhou, std::vector<double>& F_E) {
  
  for (int i = 0; i <= grid.n; ++i) {
    double rho_L, u_L, p_L, rho_R, u_R, p_R;
    
    if (i == 0 || i == grid.n) {
      // Boundary: use cell center values
      int idx = (i == 0) ? 0 : grid.n - 1;
      rho_L = rho_R = grid.rho[idx];
      u_L = u_R = grid.u[idx];
      p_L = p_R = grid.p[idx];
    } else {
      // Interior: MUSCL reconstruction
      int i_L = i - 1;
      int i_R = i;
      
      // Reconstruct LEFT cell - right edge
      {
        double v_minus = (i_L > 0) ? grid.rho[i_L - 1] : grid.rho[i_L];
        double v_plus = (i_L < grid.n - 1) ? grid.rho[i_L + 1] : grid.rho[i_L];
        double slope_rho = slopeLimit(v_minus, grid.rho[i_L], v_plus);
        rho_L = grid.rho[i_L] + 0.5 * slope_rho;
        
        v_minus = (i_L > 0) ? grid.u[i_L - 1] : grid.u[i_L];
        v_plus = (i_L < grid.n - 1) ? grid.u[i_L + 1] : grid.u[i_L];
        double slope_u = slopeLimit(v_minus, grid.u[i_L], v_plus);
        u_L = grid.u[i_L] + 0.5 * slope_u;
        
        v_minus = (i_L > 0) ? grid.p[i_L - 1] : grid.p[i_L];
        v_plus = (i_L < grid.n - 1) ? grid.p[i_L + 1] : grid.p[i_L];
        double slope_p = slopeLimit(v_minus, grid.p[i_L], v_plus);
        p_L = grid.p[i_L] + 0.5 * slope_p;
        
        rho_L = std::max(rho_L, 1e-10);
        p_L = std::max(p_L, 1e-10);
      }
      
      // Reconstruct RIGHT cell - left edge
      {
        double v_minus = (i_R > 0) ? grid.rho[i_R - 1] : grid.rho[i_R];
        double v_plus = (i_R < grid.n - 1) ? grid.rho[i_R + 1] : grid.rho[i_R];
        double slope_rho = slopeLimit(v_minus, grid.rho[i_R], v_plus);
        rho_R = grid.rho[i_R] - 0.5 * slope_rho;
        
        v_minus = (i_R > 0) ? grid.u[i_R - 1] : grid.u[i_R];
        v_plus = (i_R < grid.n - 1) ? grid.u[i_R + 1] : grid.u[i_R];
        double slope_u = slopeLimit(v_minus, grid.u[i_R], v_plus);
        u_R = grid.u[i_R] - 0.5 * slope_u;
        
        v_minus = (i_R > 0) ? grid.p[i_R - 1] : grid.p[i_R];
        v_plus = (i_R < grid.n - 1) ? grid.p[i_R + 1] : grid.p[i_R];
        double slope_p = slopeLimit(v_minus, grid.p[i_R], v_plus);
        p_R = grid.p[i_R] - 0.5 * slope_p;
        
        rho_R = std::max(rho_R, 1e-10);
        p_R = std::max(p_R, 1e-10);
      }
    }
    
    solveEntropyConservative(rho_L, u_L, p_L, rho_R, u_R, p_R, F_rho[i], F_rhou[i], F_E[i]);
  }
}

/**
 * RK2 Time Integration
 */
double computeTimestep(const Grid& grid) {
  double max_speed = 0.0;
  for (int i = 0; i < grid.n; ++i) {
    double a = std::sqrt(gamma_config * grid.p[i] / grid.rho[i]);
    max_speed = std::max(max_speed, std::abs(grid.u[i]) + a);
  }
  
  return CFL * grid.dx / max_speed;
}

void timeStep(Grid& grid) {
  double dt = computeTimestep(grid);
  if (grid.t + dt > endTime) dt = endTime - grid.t;

  // Save initial state
  std::vector<double> rho_n = grid.rho;
  std::vector<double> rho_u_n = grid.rho_u;
  std::vector<double> E_n = grid.E;

  // State 1: k1
  std::vector<double> F_rho(grid.n + 1);
  std::vector<double> F_rhou(grid.n + 1);
  std::vector<double> F_E(grid.n + 1);
  computeFluxesMUSCL(grid, F_rho, F_rhou, F_E);

  std::vector<double> k1_rho(grid.n), k1_rho_u(grid.n), k1_E(grid.n);
  for (int i = 0; i < grid.n; ++i) {
    k1_rho[i] = -(1.0 / grid.dx) * (F_rho[i + 1] - F_rho[i]);
    k1_rho_u[i] = -(1.0 / grid.dx) * (F_rhou[i + 1] - F_rhou[i]);
    k1_E[i] = -(1.0 / grid.dx) * (F_E[i + 1] - F_E[i]);

    // U*
    grid.rho[i] = rho_n[i] + dt * k1_rho[i];
    grid.rho_u[i] = rho_u_n[i] + dt * k1_rho_u[i];
    grid.E[i] = E_n[i] + dt * k1_E[i];

    grid.rho[i] = std::max(grid.rho[i], 1e-10);
    grid.E[i] = std::max(grid.E[i], 1e-10);
  }
  grid.consToPrim();

  // Stage 2: k2
  computeFluxesMUSCL(grid, F_rho, F_rhou, F_E);

  std::vector<double> k2_rho(grid.n), k2_rho_u(grid.n), k2_E(grid.n);
  for (int i = 0; i < grid.n; ++i) {
    k2_rho[i] = -(1.0 / grid.dx) * (F_rho[i + 1] - F_rho[i]);
    k2_rho_u[i] = -(1.0 / grid.dx) * (F_rhou[i + 1] - F_rhou[i]);
    k2_E[i] = -(1.0 / grid.dx) * (F_E[i + 1] - F_E[i]);

    // U^(n + 1)
    grid.rho[i] = rho_n[i] + 0.5 * dt * (k1_rho[i] + k2_rho[i]);
    grid.rho_u[i] = rho_u_n[i] + 0.5 * dt * (k1_rho_u[i] + k2_rho_u[i]);
    grid.E[i] = E_n[i] + 0.5 * dt * (k1_E[i] + k2_E[i]);

    grid.rho[i] = std::max(grid.rho[i], 1e-10);
    grid.E[i] = std::max(grid.E[i], 1e-10);
  }
  grid.consToPrim();
  grid.t += dt;
}

void runSimulation(Grid& grid) {
  int step = 0;
  while (grid.t < endTime) {
    timeStep(grid);
    step++;
  }
  std::cout << "Simulation complete: " << step << " steps, t = " << grid.t << " s\n\n";
}

/**
 * Analytical Solution
 */
void computeAnalytical(const Grid& grid, 
                       std::vector<double>& rho_exact, std::vector<double>& u_exact, std::vector<double>& p_exact) {
  ShockSolver::Config cfg;
  cfg.numCells = grid.n;
  cfg.length = grid.length;
  cfg.endTime = endTime;
  cfg.CFL = CFL;
  cfg.gamma = gamma_config;

  ShockSolver solver(cfg);
  solver.initializeShockTube(rho_L_config, u_L_config, p_L_config, rho_R_config, u_R_config, p_R_config, x_diapgragm_config);
  solver.run();

  solver.computeAnalyticalSolution(rho_exact, u_exact, p_exact);
}

/**
 * Error Analysis
 */
void computeErrors(const Grid& grid, 
                   const std::vector<double>& rho_exact, const std::vector<double>& u_exact, const std::vector<double>& p_exact) {
  double sum_rho_L1 = 0.0, sum_u_L1 = 0.0, sum_p_L1 = 0.0, sum_s_L1 = 0.0;
  double max_rho = 0.0, max_u = 0.0, max_p = 0.0, max_s = 0.0;
  int max_s_idx = 0;

  for (int i = 0; i < grid.n; ++i) {
    double err_rho = std::abs(grid.rho[i] - rho_exact[i]);
    double err_u = std::abs(grid.u[i] - u_exact[i]);
    double err_p = std::abs(grid.p[i] - p_exact[i]);

    double s_num = std::log(grid.p[i] / std::pow(grid.rho[i], gamma_config));
    double s_exact = std::log(p_exact[i] / std::pow(rho_exact[i], gamma_config));
    double err_s = std::abs(s_num - s_exact);

    max_rho = std::max(max_rho, err_rho);
    max_u = std::max(max_u, err_u);
    max_p = std::max(max_p, err_p);

    sum_rho_L1 += err_rho;
    sum_u_L1 += err_u;
    sum_p_L1 += err_p;
    sum_s_L1 += err_s;

    if (err_s > max_s) {
      max_s = err_s;
      max_s_idx = i;
    }
  }

  std::cout << "=== Error Norms ===" << std::endl;
  std::cout << "        L1        Linf" << std::endl;
  std::cout << "rho: " << std::scientific << std::setprecision(4)
            << sum_rho_L1 / grid.n << "  " << max_rho << std::endl;
  std::cout << "u: " << sum_u_L1 / grid.n << "  " << max_u << std::endl;
  std::cout << "p: " << sum_p_L1 / grid.n << "  " << max_p << std::endl;
  std::cout << "s: " << sum_s_L1 / grid.n << "  " << max_s << "\n\n";

  std::cout << "=== Max Entropy Error ===" << std::endl;
  std::cout << "Location: cell " << max_s_idx << " (x = " << grid.x[max_s_idx] << ")" << std::endl;
  double s_num_max = std::log(grid.p[max_s_idx] / std::pow(grid.rho[max_s_idx], gamma_config));
  double s_exact_max = std::log(p_exact[max_s_idx] / std::pow(rho_exact[max_s_idx], gamma_config));
  std::cout << "    s_numerical = " << s_num_max << std::endl;
  std::cout << "    s_exact = " << s_exact_max << std::endl;
  std::cout << "    error = " << s_num_max - s_exact_max << std::endl;
}

/**
 * Plotting
 */
void plotResults(const Grid& grid, const std::vector<double>& rho_ex,
                 const std::vector<double>& u_ex, const std::vector<double>& p_ex) {
  
  std::ofstream data_num("plot_numerical.dat");
  std::ofstream data_ex("plot_analytical.dat");
  
  for (int i = 0; i < grid.n; ++i) {
    double s_num = std::log(grid.p[i] / std::pow(grid.rho[i], gamma_config));
    double s_ex = std::log(p_ex[i] / std::pow(rho_ex[i], gamma_config));
    
    data_num << grid.x[i] << " " << grid.rho[i] << " " << grid.u[i] << " " 
             << grid.p[i] << " " << s_num << "\n";
    data_ex << grid.x[i] << " " << rho_ex[i] << " " << u_ex[i] << " " 
            << p_ex[i] << " " << s_ex << "\n";
  }
  
  data_num.close();
  data_ex.close();
  
  std::ofstream gp("plot_script.gp");
  gp << "set terminal pngcairo size 1200,900 enhanced font 'Arial,12'\n";
  gp << "set output 'toro_test2_results.png'\n";
  gp << "set multiplot layout 2,2 title 'Toro Test 2: MUSCL-RK2 + Entropy-Conservative Flux'\n";
  gp << "set grid\n\n";
  
  gp << "set title 'Density'\n";
  gp << "set xlabel 'Position (m)'\n";
  gp << "set ylabel 'Density (kg/m^3)'\n";
  gp << "plot 'plot_analytical.dat' u 1:2 w l lw 2 lt 2 lc rgb 'gray' title 'Analytical', \\\n";
  gp << "     'plot_numerical.dat' u 1:2 w l lw 2 lc rgb 'blue' title 'Numerical'\n\n";
  
  gp << "set title 'Velocity'\n";
  gp << "set xlabel 'Position (m)'\n";
  gp << "set ylabel 'Velocity (m/s)'\n";
  gp << "plot 'plot_analytical.dat' u 1:3 w l lw 2 lt 2 lc rgb 'gray' title 'Analytical', \\\n";
  gp << "     'plot_numerical.dat' u 1:3 w l lw 2 lc rgb 'green' title 'Numerical'\n\n";
  
  gp << "set title 'Pressure'\n";
  gp << "set xlabel 'Position (m)'\n";
  gp << "set ylabel 'Pressure (Pa)'\n";
  gp << "plot 'plot_analytical.dat' u 1:4 w l lw 2 lt 2 lc rgb 'gray' title 'Analytical', \\\n";
  gp << "     'plot_numerical.dat' u 1:4 w l lw 2 lc rgb 'red' title 'Numerical'\n\n";
  
  gp << "set title 'Entropy (Chandrashekar 2013 Entropy-Conservative Flux)'\n";
  gp << "set xlabel 'Position (m)'\n";
  gp << "set ylabel 'Entropy (J/kg·K)'\n";
  gp << "plot 'plot_analytical.dat' u 1:5 w l lw 2 lt 2 lc rgb 'gray' title 'Analytical', \\\n";
  gp << "     'plot_numerical.dat' u 1:5 w l lw 2 lc rgb 'purple' title 'Numerical'\n\n";
  
  gp << "unset multiplot\n";
  gp.close();
  
  int ret = system("gnuplot plot_script.gp 2>/dev/null");
  
  if (ret == 0) {
    std::cout << "=== PLOT GENERATED ===\n";
    std::cout << "File: toro_test2_results.png\n\n";
  } else {
    std::cout << "WARNING: gnuplot not found\n";
    std::cout << "Data saved to: plot_numerical.dat, plot_analytical.dat\n\n";
  }
}

/**
 * MAIN
 */
int main(int argc, char* argv[]) {
  int num_cells = (argc > 1) ? std::atoi(argv[1]) : 100;

  std::cout << std::string(70, '=') << std::endl;
  std::cout << "TORO TEST 2 USING ENTROPY-CONSERVATIVE SHOCK TUBE SIMULATOR" << std::endl;
  std::cout << std::string(70, '=') << std::endl;
  std::cout << "Method: MUSCL-RK2 + Chandrashekar (2013) Entropy-Conservative Flux Methods" << std::endl;
  std::cout << "Limiter: MC (Monotized Central)" << std::endl;
  std::cout << "Cells: " << num_cells << std::endl;
  std::cout << "CFL: " << CFL << std::endl;

  // Initalize
  Grid grid(num_cells, 1.0);
  grid.initialize();

  std::cout << "Running simulation...\n";
  runSimulation(grid);

  // Analyze
  std::cout << "Computing analytical solution...\n";
  std::vector<double> rho_exact, u_exact, p_exact;
  computeAnalytical(grid, rho_exact, u_exact, p_exact);

  computeErrors(grid, rho_exact, u_exact, p_exact);
  plotResults(grid, rho_exact, u_exact, p_exact);

  std::cout << "DONE\n\n";

  return 0;

}