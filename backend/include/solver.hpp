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
    double CFL = 0.9;
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

    double L2_error_density;
    double L2_error_velocity;
    double L2_error_pressure;

    double Linf_error_density;
    double Linf_error_velocity;
    double Linf_error_pressure;

    int num_points_validated;
  };

  /**
   * Toro test cases
   */
  enum class ToroTest {
    TEST1_SOD,
    TEST2_123,
    TEST3_BLAST_LEFT,
    TEST4_COLLISION,
    TEST5_STATIONARY
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
                           double x_diaphragm = 0.5);
  void initializeCustom(const std::vector<PrimitiveVars>& initial);
  void initializeToroTest(ToroTest test);
  void initializeFromSparseData(const std::vector<SparseDataPoint>& sparse_data,
                                 const std::string& interpolation_method = "linear");

  /**
   * Static helper for Toro test sensor data
   */
  static std::vector<SparseDataPoint> getSensorDataForToroTest(
    ToroTest test,
    const std::vector<double>& sensor_positions);

  /**
   * Configuration setters
   */
  void setScheme(const std::string& scheme_name);
  void setRiemannSolver(const std::string& solver_name);

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

  /**
   * Validation
   */
  void computeAnalyticalSolution(std::vector<double>& rho_exact,
                                  std::vector<double>& u_exact,
                                  std::vector<double>& p_exact) const;
  ValidationMetrics validateAgainstAnalytical() const;

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

  void computeFluxes(std::vector<double>& F_rho,
                     std::vector<double>& F_rhou,
                     std::vector<double>& F_E);

  void solveRiemann(double rho_L, double u_L, double p_L,
                    double rho_R, double u_R, double p_R,
                    double& F_rho, double& F_rhou, double& F_E);

  void solveHLLC(double rho_L, double u_L, double p_L,
                 double rho_R, double u_R, double p_R,
                 double& F_rho, double& F_rhou, double& F_E);

  void solveExact(double rho_L, double u_L, double p_L,
                  double rho_R, double u_R, double p_R,
                  double& F_rho, double& F_rhou, double& F_E);

  double solvePressureStar(double rho_L, double u_L, double p_L, double a_L,
                           double rho_R, double u_R, double p_R, double a_R) const;

  void sampleExactSolution(double rho_L, double u_L, double p_L, double a_L,
                           double rho_R, double u_R, double p_R, double a_R,
                           double p_star, double u_star,
                           double& rho_sample, double& u_sample, double& p_sample) const;

  void reconstructMUSCL(int i,
                        double& rho_L, double& u_L, double& p_L,
                        double& rho_R, double& u_R, double& p_R);

  double slopeLimit(double v_minus, double v_center, double v_plus) const;

  void updateSolution(const std::vector<double>& F_rho,
                      const std::vector<double>& F_rhou,
                      const std::vector<double>& F_E,
                      double dt);

  /**
   * Member variables
   */
  Config cfg_;
  double dx_;
  double t_;
  int step_count_;

  std::string scheme_type_;
  std::string solver_type_;

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

  // Initial conditions (for analytical solution)
  struct InitialConditions {
    double rho_L, u_L, p_L;
    double rho_R, u_R, p_R;
    double x_diaphragm;
    bool is_valid;
  } initial_conditions_;

  // Constants
  static constexpr double TOL = 1e-10;
  static constexpr int MAX_ITER = 100;
};

#endif // SOLVER_HPP