#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include <string>
#include <omp.h>

class ShockSolver {
public:
  /**
   * Configuration
   */
  struct Config {
    double length = 1.0;
    int numCells = 1000;
    double gamma = 1.4;
    double CFL = 0.5;
    double endTime = 0.2;
  };

  /**
   * Conservative variables
   */
  struct ConservativeVars {
    double rho;
    double rho_u;
    double E;
  };

  /**
   * Primitive variables
   */
  struct PrimitiveVars {
    double rho;
    double u;
    double p;

    ConservativeVars toConservative(double gamma) const;
  };

  /**
   * Monte Carlo results 
   */
  struct MonteCarloResults {
    std::vector<double> x;    // Grid position

    // Mean Values
    std::vector<double> mean_rho;
    std::vector<double> mean_u;
    std::vector<double> mean_p;
    std::vector<double> mean_entropy;

    // Standard deviations
    std::vector<double> std_rho;
    std::vector<double> std_u;
    std::vector<double> std_p;
    std::vector<double> std_entropy;

    // Confidence intervals 
    std::vector<double> ci95_lower_rho, ci95_upper_rho;
    std::vector<double> ci95_lower_u, ci95_upper_u;
    std::vector<double> ci95_lower_p, ci95_upper_p;
    std::vector<double> ci95_lower_entropy, ci95_upper_entropy;

    std::vector<double> analytical_rho, analytical_u, analytical_p, analytical_entropy;
    std::vector<double> sensor_x;

    int num_trials;     // Number of Monte Carlo samples
    double noise_level; // Gaussian noise std dev as fraction
    double computation_time_ms;
    bool success;      
  };

  /**
   * Sensor reading
   */
  struct SensorReading {
    double position;
    double density;
    double velocity;
    double pressure;
    double entropy;
    double time;
  };

  /**
   * Sparse data point for reconstruction
   */
  struct SparseDataPoint {
    double x;
    double rho;
    double u;
    double p;
  };

  /**
   * Validation metrics
   */
  struct ValidationMetrics {
    double L1_error_density;
    double L1_error_velocity;
    double L1_error_pressure;
    double L1_error_entropy;

    double L2_error_density;
    double L2_error_velocity;
    double L2_error_pressure;
    double L2_error_entropy;

    double Linf_error_density;
    double Linf_error_velocity;
    double Linf_error_pressure;
    double Linf_error_entropy;

    int num_points_validated;
  };

  /**
   * Toro test cases
   */
  enum class ToroTest {
    TEST1_SOD,
    TEST2_123,
    TEST3_BLAST_LEFT,
    TEST4_SLOW_SHOCK,
    TEST5_COLLISION
  };

  /**
   * Time Integration Flux Computation Types
   */
  enum class TimeIntegration {
    RK2,      // New Runge-Kutta 2nd order (stable, robust)
    HANCOCK
  };

  /**
   * Constructor
   */
  explicit ShockSolver(const Config& cfg);

  /**
   * Initialization
   */
  void initializeShockTube(double rho_L, double u_L, double p_L,
                           double rho_R, double u_R, double p_R,
                           double x_diaphragm = 0.5, double endTime = 0.25);
  void initializeCustom(const std::vector<PrimitiveVars>& initial);
  void initializeToroTest(ToroTest test);
  void initializeFromSparseData(const std::vector<SparseDataPoint>& sparse_data,
                                 const std::string& interpolation_method = "piecewise_constant",
                                 double diaphragm_x = -1.0);

  /**
   * Static helper for Toro test sensor data
   */
  static std::vector<SparseDataPoint> getSensorDataForToroTest(
    ToroTest test,
    const std::vector<double>& sensor_positions);

  /**
   * Setter
   */
  void setFlux(const std::string& flux_name);
  void setTimeIntegration(const std::string& method);

  /**
   * Simulation control
   */
  void step();
  void run();

  /**
   * Sensor operations
   */
  void placeSensors(const std::vector<double>& positions);
  std::vector<SensorReading> getSensorReadings() const;
  std::vector<SensorReading> getAnalyticalAtSensors() const;

  /**
   * Validation
   */
  void computeAnalyticalSolution(std::vector<double>& rho_exact,
                                  std::vector<double>& u_exact,
                                  std::vector<double>& p_exact,
                                  std::vector<double>& entropy_exact) const;
  ValidationMetrics validateAgainstAnalytical() const;

  /**
   * Anomaly Detection System
   * 
   * Compares predicted sensor readings from physics simulation against
   * actual measurements to detect deviations indicating:
   *   - Sensor faults (drift, failure, noise)
   *   - Model errors (physics assumptions violated)
   */
  /**
   * Monte Carlo Uncertainty Propogation Method
   */
  MonteCarloResults runMonteCarloUncertainty(
    const std::vector<SparseDataPoint>& base_sensors,   // Nominal sensor readings without noise
    double noise_level,                                 // Guassian noise std dev as fraction
    int num_trials,                                     // Number of Monte Carlo samples
    const std::string& interpolation_method,            // linear or piecewise constant methods
    const std::string& time_integration_method         // Either RK2 for RK2 or Hancock for MUSCL-Hancock
  );

  /**
   * Getters
   */
  const std::vector<double>& getDensity() const { return rho_; }
  const std::vector<double>& getVelocity() const { return u_; }
  const std::vector<double>& getPressure() const { return p_; }
  const std::vector<double>& getPositions() const { return x_; }
  std::vector<double> getEntropy() const;

  double getCurrentTime() const { return t_; }
  int getStepCount() const { return step_count_; }
  double getTotalMass() const;
  double getTotalEnergy() const;
  double getEntropyGeneration() const;

private:
  /**
   * Internal methods
   */
  void primitivesToConservative();
  void conservativesToPrimitives();
  double computeTimestep() const;
  void applyBoundaryConditions();

  void solveHLLC(double rho_L, double u_L, double p_L,
                 double rho_R, double u_R, double p_R,
                 double& F_rho, double& F_rhou, double& F_E);

  double solvePressureStar(double rho_L, double u_L, double p_L, double a_L,
                           double rho_R, double u_R, double p_R, double a_R) const;

  double slopeLimit(double v_minus, double v_center, double v_plus) const;

  void solveEntropyConservative(double rho_L, double u_L, double p_L,
                                double rho_R, double u_R, double p_R,
                                double& F_rho, double& F_rhou, double& F_E);

  void solveEntropyStable(double rho_L, double u_L, double p_L,
                          double rho_R, double u_R, double p_R,
                          double& F_rho, double& F_rhou, double& F_E);

  void computeDissipation(double rho_L, double u_L, double p_L,
                          double rho_R, double u_R, double p_R,
                          double& D_rho, double& D_rhou, double& D_E);

  void computeFluxesMUSCL(std::vector<double>& F_rho, std::vector<double>& F_rhou, std::vector<double>& F_E);
                          

  /**
   * Member variables
   */
  Config cfg_;
  double dx_;
  double t_;
  int step_count_;

  std::string flux_type_;
  TimeIntegration time_integration_;

  void stepRK2();
  void stepHancock();

  // Primitive variables
  std::vector<double> rho_;
  std::vector<double> u_;
  std::vector<double> p_;

  // Conservative variables
  std::vector<double> rho_u_;
  std::vector<double> E_;

  // Grid
  std::vector<double> x_;

  // Sensors
  std::vector<double> sensor_positions_;
  std::vector<int> sensor_indices_;

  // Initial conditions (for anaalytical solution)
  struct InitialConditions {
    double rho_L, u_L, p_L;
    double rho_R, u_R, p_R;
    double x_diaphragm;
    double endTime;
    bool is_valid;
  } initial_conditions_;

  // Constants
  static constexpr double TOL = 1e-10;
  static constexpr int MAX_ITER = 100;
};

#endif // SOLVER_HPP