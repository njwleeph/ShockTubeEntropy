#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include <string>

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
   * Anomaly detection configuration
   */
  struct AnomalyConfig {
    double density_threshold = 0.05;    // 5% deviation
    double velocity_threshold = 0.10;   // 10% deviation (velocity can be near zero)
    double pressure_threshold = 0.05;   // 5% deviation
    double entropy_threshold = 0.05;    // 5% deviation
    double score_threshold = 0.1;       // Combined metric threshold
    
    // Weights for composite score
    double weight_density = 1.0;
    double weight_velocity = 0.5;
    double weight_pressure = 1.5;       // Pressure typically most reliable
    double weight_entropy = 1.0;        // Entropy useful for detecting non-isentropic events
  };

  /**
   * Anomaly detection result for a single sensor
   */
  struct AnomalyResult {
    double time;
    int sensor_index;
    double position;
    
    // Predicted values from simulation
    double predicted_density;
    double predicted_velocity;
    double predicted_pressure;
    double predicted_entropy;
    
    // Actual values from sensor
    double actual_density;
    double actual_velocity;
    double actual_pressure;
    double actual_entropy;
    
    // Raw residuals (predicted - actual)
    double residual_density;
    double residual_velocity;
    double residual_pressure;
    double residual_entropy;
    
    // Normalized residuals (absolute % error)
    double normalized_density;
    double normalized_velocity;
    double normalized_pressure;
    double normalized_entropy;
    
    // Overall anomaly score (weighted RMS of normalized errors)
    double anomaly_score;
    bool is_anomalous;
    
    // Which variable triggered the flag
    std::string primary_deviation;
  };

  /**
   * Summary of anomaly check across all sensors
   */
  struct AnomalySummary {
    double time;
    int total_sensors;
    int anomalous_sensors;
    double max_anomaly_score;
    int max_score_sensor_index;
    bool system_anomaly;           // True if multiple sensors flagged
    std::string severity;          // "nominal", "warning", "critical"
    std::vector<AnomalyResult> sensor_results;
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
                                  std::vector<double>& p_exact) const;
  ValidationMetrics validateAgainstAnalytical() const;

  /**
   * Anomaly Detection System
   * 
   * Compares predicted sensor readings from physics simulation against
   * actual measurements to detect deviations indicating:
   *   - Sensor faults (drift, failure, noise)
   *   - Model errors (physics assumptions violated)
   *   - Physical anomalies (unexpected flow behavior)
   */
  void setAnomalyConfig(const AnomalyConfig& config);
  AnomalyConfig getAnomalyConfig() const { return anomaly_config_; }
  
  AnomalyResult compareSensorReading(
      int sensor_idx,
      const SensorReading& predicted,
      const SensorReading& actual) const;
  
  std::vector<AnomalyResult> checkForAnomalies(
      const std::vector<SensorReading>& actual_readings) const;
  
  AnomalySummary analyzeAnomalies(
      const std::vector<SensorReading>& actual_readings) const;
  
  AnomalySummary analyzeAnomaliesWithReference(
      const std::vector<SensorReading>& actual_readings,
      const std::vector<SensorReading>& reference_readings) const;

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

  // Anomaly detection configuration
  AnomalyConfig anomaly_config_;

  // Initial conditions (for analytical solution)
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