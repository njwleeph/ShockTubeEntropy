#pragma once
#include <vector>
#include <array>
#include <memory>
#include <string>

/**
 * 1D Euler equation solver for compressible flow with shock discont.
 * 
 * Solver the conservative form of Euler equations:
 *   ∂U/∂t = ∂F(U)/∂x = 0
 * 
 * where U = [ρ, ρu, E]^T and F = [ρu, ρu² + p, u(E + p)]^T
 * 
 * Features:
 * - Multiple numerical schemes (Godunov, MUSCL-Hancock)
 * - Multiple Riemann solvers (Exact, HLLC)
 * - Accurate shock caputring with TVD slop limiters
 * - Validated against analytical solutions (Sod's shock tube)
 */

class ShockSolver {
public:
  struct Config {
    double length;      // Domain length [m]
    int numCells;       // grid cells
    double gamma;       // specific heat ratio (1.4 for air 1.22 for combustion)
    double CFL;         // CFL number for timestep stability 
    double endTime;     // Simulation end time [s]
  };

/**
 * Conservative variables (rho, momentum, energy)
 */
  struct ConservativeVars {
    double rho;         // Density
    double rho_u;       // Momentum
    double E;           // Total energy

    /**
     * Compute velocity from momentum
     */
    double velocity() const { return rho_u / rho;}

    /**
     * Compute pressure using ideal gas equation of state
     * gamma Specific heat ratio
     */
    double pressure(double gamma) const {
      double u = velocity();
      return (gamma - 1.0) * (E - 0.5 * rho * u * u);
    }

    /**
     * Compute sound speed
     * gamma Specific heat ratio
     */
    double soundSpeed(double gamma) const {
      return std::sqrt(gamma * pressure(gamma)/ rho);
    }
  };

  /**
   * Primitive variables
   */
  struct PrimitiveVars {
    double rho;     // Density [kg/m^3]
    double u;       // Velocity [m/s]
    double p;       // Pressure [Pa]

    /**
     * Convert to conservative form
     * gamma Specific ehat ratio
     */
    ConservativeVars toConservative(double gamma) const;
  };

  /**
   * Flux vector (F_rho, F_rhou, F_E)
   */
  struct FluxVector {
    double F_rho;       // Mass flux
    double F_rhou;      // Momentum flux
    double F_E;         // Energy flux

    FluxVector() : F_rho(0), F_rhou(0), F_E(0) {}
    FluxVector(double f1, double f2, double f3) : F_rho(f1), F_rhou(f2), F_E(f3) {}
  };

  /**
   * Constructor
   * Configuration parameters
   */
  ShockSolver(const Config& cfg);

  /**
   * Initialize with classic shock tube problem (Sod, 1978)
   * rho_L Left state density
   * u_L Left state velocity
   * p_L Left state pressure
   * rho_R Right state density
   * u_R Right state velocity
   * p_R Right state pressure
   * x_diaphragm Position of initial discontinuity
   */
  void initializeShockTube(double rho_L, double u_L, double p_L,
                          double rho_R, double u_R, double p_R,
                          double x_diaphragm);

  /**
   * Initialize with custom initial conditions
   * Initial Vector of primitive states (one per cell)
   */
  void initializeCustom(const std::vector<PrimitiveVars>& initial);

  /**
   * Set numerical scheme
   * "godunov" (first-order) or "muscl" (second-order)
   */
  void setScheme(const std::string& scheme_name);

  /**
   * Set Riemann solver
   * solver_name "exact" or "hllc"
   */
  void setRiemannSolver(const std::string& solver_name);

  /**
   * Run simulation until end time
   */
  void run();

  /**
   * Advance solution by one timestep
   */
  void step();

  /**
   * Getters
   */
  const std::vector<double>& getDensity() const { return rho_; }
  const std::vector<double>& getVelocity() const { return u_; }
  const std::vector<double>& getPressure() const { return p_; }
  const std::vector<double>& getPositions() const { return x_; }
  double getCurrentTime() const { return t_; }
  int getStepCount() const { return step_count_; }
  double getTotalMass() const;
  double getTotalEnergy() const;
  double getEntropyGeneration() const;
  std::vector<double> getEntropy() const;

  /**
   * To cross-reference solutions use Toro's solutions
   */
  enum class ToroTest {
    TEST1_SOD,            // Standard Sod Problem
    TEST2_123,            // 123 problem (strong rarefraction)
    TEST3_BLAST_LEFT,     // Left half of blast wave
    TEST4_COLLISION,      // Two rarefraction waves
    TEST5_STATIONARY      // Stationary contact
  };

  /**
   * Initialize with Toro's standard test problems
   */
  void initializeToroTest(ToroTest test);

  /**
   * Sensor data structure
   */
  struct SensorReading{
    double position;
    double density;
    double velocity;
    double pressure;
    double entropy;
    double time;
  };

  /**
   * Place sensors
   */
  void placeSensors(const std::vector<double>& positions);

  /**
   * Get current sensor readings
   */
  std::vector<SensorReading> getSensorReadings() const;

  /**
   * Get sensor positions
   */
  const std::vector<double>& getSensorPositions() const { return sensor_positions_; }

  /**
   * Validation against analytical solution
   */
  struct ValidationMetrics {
    double L1_error_density;
    double L2_error_density;
    double Linf_error_density;
    double L1_error_velocity;
    double L2_error_velocity;
    double Linf_error_velocity;
    double L1_error_pressure;
    double L2_error_pressure;
    double Linf_error_pressure;
    double max_relative_error_density;
    double max_relative_error_velocity;
    double max_relative_error_pressure;
    bool passes_3percent_threshold;
    int num_points_validated;
  };

  /**
   * Compute analytical solution at current time
   * Uses exact Riemann solver to generate reference solution
   * Based on Toro's analytical solution methodology
   */
  void computeAnalyticalSolution(std::vector<double>& rho_exact,
                                 std::vector<double>& u_exact,
                                 std::vector<double>& p_exact) const;
  
  /**
   * Validate numerical solution against analytical solution
   * Returns error metrics (L1, L2, Linf norms)
   * Checks max relative error
   */
  ValidationMetrics validateAgainstAnalytical() const;

  /**
   * Get initial conditions (needed for analytical solution)
   */
  struct InitialConditions {
    double rho_L, u_L, p_L;
    double rho_R, u_R, p_R;
    double x_diaphragm;
    bool is_valid;
  };

  const InitialConditions& getInitialConditions() const { return initial_conditions_; }


private:
  Config cfg_;             
  double dx_;               // Cell spacing [m]
  double t_;                // Current time [s]
  int step_count_;          // Number of timesteps taken

  // Solution arrays
  std::vector<double> rho_;     // Density field
  std::vector<double> u_;       // Velocity field
  std::vector<double> p_;       // Pressure field
  std::vector<double> x_;       // Position field

  // Conservative variables (internal)
  std::vector<double> rho_u_;    // Momentum field
  std::vector<double> E_;       // Total energy field

  // Sensor variables
  std::vector<double> sensor_positions_;
  std::vector<double> sensor_indices_;

  // Numerical method selection
  std::string scheme_type_;     // "godunov" or "muscl"
  std::string solver_type_;     // "exact" or "hllc"

  // Store initial conditions for analytical solution
  InitialConditions initial_conditions_;

  /**
   * Helpers - inline for performance
   */
  inline int idx(int i) const { return i; }

  /**
   * Convert primitive to conservative variables
   */
  void primitivesToConservative();

  /**
   * Convert conservative to primitive variables
   */
  void conservativesToPrimitives();

  /**
   * Compute timestep from CFL condition
   * Stable timestep [s]
   */
  double computeTimestep() const;

  /**
   * Compute fluxes at all cell interfaces
   * F_rho Mass flux array (size numCells + 1)
   * F_rhou Momentum flux array (size numCells + 1)
   * F_E Energy flux array (size numCells + 1)
   */
  void computeFluxes(std::vector<double>& F_rho,
                     std::vector<double>& F_rhou,
                     std::vector<double>& F_E);

  /**
   * Solve Riemann problem at single interface
   * rho_L Left state density
   * u_L Left state velocity
   * p_L Left state pressure
   * rho_R Right state density
   * u_R Right state velocity
   * p_R Right state pressure
   * F_rho Output: mass flux
   * F_rhou Output: momentum flux
   * F_E Output: energy flux
   */
  void solveRiemann(double rho_L, double u_L, double p_L,
                   double rho_R, double u_R, double p_R,
                   double& F_rho, double& F_rhou, double& F_E);

  /**
   * Exact Riemann solver (iterative)
   */
  void solveExact(double rho_l, double u_L, double p_L,
                 double rho_R, double u_R, double p_R,
                 double& F_rho, double& F_rhou, double& F_E);
  
  /**
   * HLLC approximate Riemann solver
   */
  void solveHLLC(double rho_L, double u_L, double p_L,
                double rho_R, double u_R, double p_R,
                double& F_rho, double& F_rhou, double& F_E);
  
  /**
   * MUSCL reconstruction with slope limiters
   * i Cell index
   * rho_L Output: reconstructed left density at i + 1 / 2
   * u_L Output: reconstructed left velocity at i + 1 / 2
   * p_L Output: reconstructed left pressure at i + 1 / 2
   * rho_R Output: reconstructed right density at i + 1 / 2
   * u_R Output: reconstructed right velocity at i + 1 / 2
   * p_R Output: reconstructed right pressure at i + 1 / 2
   */
  void reconstructMUSCL(int i,
                        double& rho_L, double& u_L, double& p_L,
                        double& rho_R, double& u_R, double& p_R);

  /**
   * Update conservative variables with computed fluxes
   * F_rho Mass flux array
   * F_rhou Momentum flux array
   * F_E Energy flux array
   * dt Timestep
   */
  void updateSolution(const std::vector<double>&  F_rho,
                     const std::vector<double>& F_rhou,
                     const std::vector<double>& F_E,
                     double dt);
  /**
   * Solver for p_star in exact Riemann Solver (Newton-Raphson)
   */
  double solvePressureStar(double rho_L, double u_L, double p_L, double a_L,
                          double rho_R, double u_R, double p_R, double a_R) const;

  /**
   * Sample exact Riemann solution at x/t = 0
   */
  void sampleExactSolution(double rho_L, double u_L, double p_L, double a_L,
                          double rho_R, double u_R, double p_R, double a_R,
                          double p_star, double u_star,
                          double& rho_sample, double& u_sample, double& p_sample) const;
  
  /**
   * Slope Limiter Function
   */
  double slopeLimit(double v_minus, double v_center, double v_plus) const;

  /**
   * Apply boundary conditions (transmissive)
   */
  void applyBoundaryConditions();
  // Constants
  static constexpr double TOL = 1e-6;                   // Tolerance for solvers
  static constexpr int MAX_ITER = 20;            // Max iterations
};