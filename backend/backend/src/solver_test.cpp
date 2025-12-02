#include "solver.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

/**
 * MUSCL Solver Test
 */

void printHeader(const std::string& title) {
  std::cout << "\n" << std::string(60, '=') << std::endl;
  std::cout << "  " << title << std::endl;
  std::cout << std::string(60, '=') << std::endl;
}

void testSolverMUSCL() {
  printHeader("Testing Godunov vs MUSCL");
  
  ShockSolver::Config cfg;
  cfg.length = 1.0;
  cfg.numCells = 100;
  cfg.gamma = 1.4;
  cfg.CFL = 0.9;
  cfg.endTime = 0.2;
  
  double rho_L = 1.0, u_L = 0.0, p_L = 1.0;
  double rho_R = 0.125, u_R = 0.0, p_R = 0.1;
  double x_diaphragm = 0.5;
  
  std::cout << "\nSod shock tube, 100 cells, t=0.2" << std::endl;
  
  // Run Godunov
  std::cout << "\n--- Godunov (CFL=0.9) ---" << std::endl;
  ShockSolver godunov(cfg);
  godunov.setScheme("godunov");
  godunov.setRiemannSolver("hllc");
  godunov.initializeShockTube(rho_L, u_L, p_L, rho_R, u_R, p_R, x_diaphragm);
  godunov.run();
  
  // Run MUSCL with CFL=0.9
  std::cout << "\n--- MUSCL (CFL=0.9) ---" << std::endl;
  ShockSolver muscl(cfg);
  muscl.setScheme("muscl");
  muscl.setRiemannSolver("hllc");
  muscl.initializeShockTube(rho_L, u_L, p_L, rho_R, u_R, p_R, x_diaphragm);
  muscl.run();
  
  // Run MUSCL with CFL=0.5
  std::cout << "\n--- MUSCL (CFL=0.5) ---" << std::endl;
  cfg.CFL = 0.5;
  ShockSolver muscl_low(cfg);
  muscl_low.setScheme("muscl");
  muscl_low.setRiemannSolver("hllc");
  muscl_low.initializeShockTube(rho_L, u_L, p_L, rho_R, u_R, p_R, x_diaphragm);
  muscl_low.run();
  
  // Get analytical solution
  std::vector<double> rho_exact, u_exact, p_exact;
  godunov.computeAnalyticalSolution(rho_exact, u_exact, p_exact);
  
  const auto& rho_g = godunov.getDensity();
  const auto& rho_m = muscl.getDensity();
  const auto& rho_ml = muscl_low.getDensity();
  
  // Find max errors
  double max_g = 0, max_m = 0, max_ml = 0;
  int max_g_i = 0, max_m_i = 0, max_ml_i = 0;
  for (int i = 0; i < 100; ++i) {
    double eg = std::abs(rho_g[i] - rho_exact[i]);
    double em = std::abs(rho_m[i] - rho_exact[i]);
    double eml = std::abs(rho_ml[i] - rho_exact[i]);
    if (eg > max_g) { max_g = eg; max_g_i = i; }
    if (em > max_m) { max_m = em; max_m_i = i; }
    if (eml > max_ml) { max_ml = eml; max_ml_i = i; }
  }
  
  const auto& x = godunov.getPositions();
  std::cout << "\n--- Max errors ---" << std::endl;
  std::cout << "Godunov:        " << max_g << " at cell " << max_g_i << " (x=" << x[max_g_i] << ")" << std::endl;
  std::cout << "MUSCL CFL=0.9:  " << max_m << " at cell " << max_m_i << " (x=" << x[max_m_i] << ")" << std::endl;
  std::cout << "MUSCL CFL=0.5:  " << max_ml << " at cell " << max_ml_i << " (x=" << x[max_ml_i] << ")" << std::endl;
  
  // Print around max MUSCL error
  std::cout << "\n--- Around MUSCL CFL=0.9 max error (cell " << max_m_i << ") ---" << std::endl;
  std::cout << std::setw(6) << "Cell" << std::setw(8) << "x" 
            << std::setw(10) << "Exact" 
            << std::setw(10) << "Godunov" 
            << std::setw(10) << "MUSCL.9"
            << std::setw(10) << "MUSCL.5" << std::endl;
  
  int start = std::max(0, max_m_i - 5);
  int end = std::min(99, max_m_i + 5);
  for (int i = start; i <= end; ++i) {
    std::cout << std::fixed << std::setprecision(4)
              << std::setw(6) << i 
              << std::setw(8) << x[i]
              << std::setw(10) << rho_exact[i]
              << std::setw(10) << rho_g[i] 
              << std::setw(10) << rho_m[i]
              << std::setw(10) << rho_ml[i] << std::endl;
  }
  
  // Error norms
  auto g_metrics = godunov.validateAgainstAnalytical();
  auto m_metrics = muscl.validateAgainstAnalytical();
  auto ml_metrics = muscl_low.validateAgainstAnalytical();
  
  std::cout << "\n--- L1 Error (density) ---" << std::endl;
  std::cout << "Godunov:        " << std::scientific << g_metrics.L1_error_density << std::endl;
  std::cout << "MUSCL CFL=0.9:  " << m_metrics.L1_error_density << std::endl;
  std::cout << "MUSCL CFL=0.5:  " << ml_metrics.L1_error_density << std::endl;
  
  std::cout << "\n--- Ratios (vs Godunov) ---" << std::endl;
  std::cout << "MUSCL CFL=0.9:  " << std::fixed << std::setprecision(3) 
            << m_metrics.L1_error_density / g_metrics.L1_error_density << std::endl;
  std::cout << "MUSCL CFL=0.5:  " 
            << ml_metrics.L1_error_density / g_metrics.L1_error_density << std::endl;
  std::cout << "(Should be < 1.0 if MUSCL is working)" << std::endl;
}

int main() {
  std::cout << R"(
  ╔═══════════════════════════════════════════╗
  ║       MUSCL Solver Debug Test Suite       ║
  ╚═══════════════════════════════════════════╝
  )" << std::endl;

  testSolverMUSCL();

  std::cout << "\nDone." << std::endl;
  
  return 0;
}