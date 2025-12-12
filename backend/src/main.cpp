#include "solver.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include "cpp-httplib/httplib.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

/**
 * Global States
 */
std::unique_ptr<ShockSolver> g_solver = nullptr;

/**
 * Helpers
 */

// Convert data to JSON
json solutionToJson(const ShockSolver& solver) {
  json result;

  const auto& rho = solver.getDensity();
  const auto& u = solver.getVelocity();
  const auto& p = solver.getPressure();
  const auto& x = solver.getPositions();
  auto entropy = solver.getEntropy();

  result["time"] = solver.getCurrentTime();
  result["step"] = solver.getStepCount();
  result["numCells"] = rho.size();

  result["x"] = json::array();
  result["density"] = json::array();
  result["velocity"] = json::array();
  result["pressure"] = json::array();
  result["entropy"] = json::array();
  
  for (size_t i = 0; i < rho.size(); ++i) {
    result["x"].push_back(x[i]);
    result["density"].push_back(rho[i]);
    result["velocity"].push_back(u[i]);
    result["pressure"].push_back(p[i]);
    result["entropy"].push_back(entropy[i]);
  }

  return result;
}

// Create error response JSON
json errorResponse(const std::string& message) {
  json response;
  response["success"] = false;
  response["error"] = message;
  return response;
}

// Create success response JSON
json successResponse(const std::string& message) {
  json response;
  response["success"] = true;
  response["message"] = message;
  return response;
}

/**
 * API Endpoints
 */

/**
 * POST /api/simulation/create
 */
void handleCreateSimulation(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/simulation/create" << std::endl;
    
    auto body = json::parse(req.body);
    std::cout << "Request body: " << body.dump(2) << std::endl;

    ShockSolver::Config cfg;
    cfg.length = body.value("length", 1.0);
    cfg.numCells = body.value("numCells", 1000);
    cfg.gamma = body.value("gamma", 1.4);
    cfg.CFL = body.value("CFL", 0.5);
    cfg.endTime = body.value("endTime", 0.2);

    std::cout << "Creating solver with config:" << std::endl;
    std::cout << "  length: " << cfg.length << std::endl;
    std::cout << "  numCells: " << cfg.numCells << std::endl;
    std::cout << "  gamma: " << cfg.gamma << std::endl;
    std::cout << "  CFL: " << cfg.CFL << std::endl;
    std::cout << "  endTime: " << cfg.endTime << std::endl;

    g_solver = std::make_unique<ShockSolver>(cfg);
    std::cout << "Solver created successfully" << std::endl;

    json response = successResponse("Simulation created successfully");
    response["config"] = {
      {"length", cfg.length},
      {"numCells", cfg.numCells},
      {"gamma", cfg.gamma},
      {"CFL", cfg.CFL},
      {"endTime", cfg.endTime}
    };

    res.set_content(response.dump(), "application/json");
    std::cout << "Response sent" << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleCreateSimulation: " << e.what() << std::endl;
    res.status = 400;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/simulation/init/shocktube
 */
void handleInitShockTube(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/simulation/init/shocktube" << std::endl;
    
    if (!g_solver) {
      std::cerr << "ERROR: Solver not created" << std::endl;
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    std::cout << "Request body: " << body.dump(2) << std::endl;

    double rho_L = body.value("rho_L", 1.0);
    double u_L = body.value("u_L", 0.0);
    double p_L = body.value("p_L", 1.0);
    double rho_R = body.value("rho_R", 0.125);
    double u_R = body.value("u_R", 0.0);
    double p_R = body.value("p_R", 0.1);
    double x_diaphragm = body.value("x_diaphragm", 0.5);
    double endTime = body.value("endTime", 0.25);

    std::cout << "Initializing shock tube:" << std::endl;
    std::cout << "  Left:  rho=" << rho_L << ", u=" << u_L << ", p=" << p_L << std::endl;
    std::cout << "  Right: rho=" << rho_R << ", u=" << u_R << ", p=" << p_R << std::endl;
    std::cout << "  Diaphragm: " << x_diaphragm << std::endl;

    g_solver->initializeShockTube(rho_L, u_L, p_L, rho_R, u_R, p_R, x_diaphragm, endTime);
    std::cout << "Shock tube initialized" << std::endl;

    json response = successResponse("Shock tube initialized");
    response["initial_conditions"] = {
      {"left", {{"rho", rho_L}, {"u", u_L}, {"p", p_L}}},
      {"right", {{"rho", rho_R}, {"u", u_R}, {"p", p_R}}},
      {"diaphragm", x_diaphragm},
      {"end time", endTime}
    };

    res.set_content(response.dump(), "application/json");
    std::cout << "Response sent" << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleInitShockTube: " << e.what() << std::endl;
    res.status = 400;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/simulation/init/toro
 */
void handleInitToro(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/simulation/init/toro" << std::endl;

    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    int test_number = body.value("test", 1);

    ShockSolver::ToroTest test;
    switch(test_number) {
      case 1: test = ShockSolver::ToroTest::TEST1_SOD; break;
      case 2: test = ShockSolver::ToroTest::TEST2_123; break;
      case 3: test = ShockSolver::ToroTest::TEST3_BLAST_LEFT; break;
      case 4: test = ShockSolver::ToroTest::TEST4_SLOW_SHOCK; break;
      case 5: test = ShockSolver::ToroTest::TEST5_COLLISION; break;
      default: test = ShockSolver::ToroTest::TEST1_SOD;
    }

    g_solver->initializeToroTest(test);

    json response = successResponse("Toro test initialized");
    response["test_number"] = test_number;
    res.set_content(response.dump(), "application/json");
  } catch (const std::exception& e) {
    res.status = 400;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/simulation/sensors/place
 */
void handlePlaceSensors(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/simulation/sensors/place" << std::endl;

    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    std::vector<double> positions = body["positions"];

    g_solver->placeSensors(positions);

    json response = successResponse("Sensors placed");
    response["sensor_count"] = positions.size();
    response["positions"] = positions;
    res.set_content(response.dump(), "application/json");
  } catch (const std::exception& e) {
    res.status = 400;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/simulation/configure
 */
void handleConfigure(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/simulation/configure" << std::endl;
    
    if (!g_solver) {
      std::cerr << "ERROR: Solver not created" << std::endl;
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    std::cout << "Request body: " << body.dump(1) << std::endl;

    std::string flux = body.value("flux", "HLLC");
    std::string time_integration = body.value("timeIntegration", "RK2");

    std::cout << "Configuring: " << std::endl;
    std::cout << "  Flux: " << flux << std::endl;
    std::cout << "  Time Integration: " << time_integration << std::endl;

    g_solver->setFlux(flux);
    g_solver->setTimeIntegration(time_integration);
    std::cout << "Configuration updated" << std::endl;

    json response = successResponse("Configuration updated");
    response["fluxType"] = flux;
    response["time_integration"] = time_integration;

    res.set_content(response.dump(), "application/json");
    std::cout << "Response sent" << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleConfigure: " << e.what() << std::endl;
    res.status = 400;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/simulation/run
 */
void handleRun(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/simulation/run <<<" << std::endl;
    std::cout << "======================================" << std::endl;
    
    if (!g_solver) {
      std::cerr << "ERROR: No solver exists!" << std::endl;
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    std::cout << "Solver exists, calling run()..." << std::endl;
    std::cout << "Current time before run: " << g_solver->getCurrentTime() << std::endl;
    std::cout << "Current step before run: " << g_solver->getStepCount() << std::endl;
    
    std::cout << "\n*** CALLING g_solver->run() ***" << std::endl;
    g_solver->run();
    std::cout << "*** g_solver->run() RETURNED ***\n" << std::endl;

    std::cout << "Current time after run: " << g_solver->getCurrentTime() << std::endl;
    std::cout << "Current step after run: " << g_solver->getStepCount() << std::endl;

    double final_time = g_solver->getCurrentTime();
    int total_steps = g_solver->getStepCount();
    double total_mass = g_solver->getTotalMass();
    double total_energy = g_solver->getTotalEnergy();

    std::cout << "\nSimulation Results:" << std::endl;
    std::cout << "  Final time: " << final_time << std::endl;
    std::cout << "  Total steps: " << total_steps << std::endl;
    std::cout << "  Total mass: " << total_mass << std::endl;
    std::cout << "  Total energy: " << total_energy << std::endl;

    json response = successResponse("Simulation completed");
    response["final_time"] = final_time;
    response["total_steps"] = total_steps;
    response["total_mass"] = total_mass;
    response["total_energy"] = total_energy;

    res.set_content(response.dump(), "application/json");
    std::cout << "Response sent" << std::endl;
    std::cout << "======================================\n" << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "\nEXCEPTION IN handleRun: " << e.what() << "\n" << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/simulation/step
 */
void handleStep(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/simulation/step" << std::endl;
    
    if (!g_solver) {
      std::cerr << "ERROR: Solver not created" << std::endl;
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    int num_steps = body.value("num_steps", 1);
    std::cout << "Stepping " << num_steps << " times" << std::endl;

    for (int i = 0; i < num_steps; ++i) {
      g_solver->step();
    }

    json response = successResponse("Steps completed");
    response["current_time"] = g_solver->getCurrentTime();
    response["step_count"] = g_solver->getStepCount();

    res.set_content(response.dump(), "application/json");
    std::cout << "Steps completed, t=" << g_solver->getCurrentTime() << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleStep: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * GET /api/simulation/status
 */
void handleGetStatus(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> GET /api/simulation/status" << std::endl;
    
    if (!g_solver) {
      std::cerr << "ERROR: Solver not created" << std::endl;
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    json response;
    response["success"] = true;
    response["time"] = g_solver->getCurrentTime();
    response["step"] = g_solver->getStepCount();
    response["total_mass"] = g_solver->getTotalMass();
    response["total_energy"] = g_solver->getTotalEnergy();

    res.set_content(response.dump(), "application/json");
    std::cout << "Status sent" << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleGetStatus: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * GET /api/simulation/data
 */
void handleGetData(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> GET /api/simulation/data" << std::endl;
    
    if (!g_solver) {
      std::cerr << "ERROR: Solver not created" << std::endl;
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    std::cout << "Converting solution to JSON..." << std::endl;
    json response = solutionToJson(*g_solver);
    response["success"] = true;

    std::cout << "Solution contains " << response["numCells"] << " cells" << std::endl;
    std::cout << "Time: " << response["time"] << ", Steps: " << response["step"] << std::endl;

    res.set_content(response.dump(), "application/json");
    std::cout << "Data sent" << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleGetData: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * GET /api/simulation/sensors/readings
 */
void handleGetSensorReadings(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> GET /api/simulation/sensors/readings" << std::endl;

    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto readings = g_solver->getSensorReadings();

    json response;
    response["success"] = true;
    response["time"] = g_solver->getCurrentTime();
    response["readings"] = json::array();

    for (const auto& reading : readings) {
      json r;
      r["position"] = reading.position;
      r["density"] = reading.density;
      r["velocity"] = reading.velocity;
      r["pressure"] = reading.pressure;
      r["entropy"] = reading.entropy;
      response["readings"].push_back(r);
    }

    res.set_content(response.dump(), "application/json");
  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * GET /api/simulation/sensors/analytical
 * Get analytical (exact Riemann) solution at sensor positions
 * This is the "ground truth" for anomaly detection comparison
 */
void handleGetAnalyticalAtSensors(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> GET /api/simulation/sensors/analytical" << std::endl;

    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto readings = g_solver->getAnalyticalAtSensors();

    json response;
    response["success"] = true;
    response["time"] = g_solver->getCurrentTime();
    response["readings"] = json::array();

    for (const auto& reading : readings) {
      json r;
      r["position"] = reading.position;
      r["density"] = reading.density;
      r["velocity"] = reading.velocity;
      r["pressure"] = reading.pressure;
      r["entropy"] = reading.entropy;
      response["readings"].push_back(r);
    }

    res.set_content(response.dump(), "application/json");
    std::cout << "Returned analytical solution at " << readings.size() << " sensor positions" << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleGetAnalyticalAtSensors: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * GET /api/simulation/validate
 * Validate numerical solution against analytical solution
 */
void handleValidate(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> GET /api/simulation/validate" << std::endl;
    
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto metrics = g_solver->validateAgainstAnalytical();

    json response;
    response["success"] = true;
    
    // L1 errors (mean absolute error)
    response["L1_density"] = metrics.L1_error_density;
    response["L1_velocity"] = metrics.L1_error_velocity;
    response["L1_pressure"] = metrics.L1_error_pressure;
    response["L1_entropy"] = metrics.L1_error_entropy;
    
    // L2 errors (RMS error)
    response["L2_density"] = metrics.L2_error_density;
    response["L2_velocity"] = metrics.L2_error_velocity;
    response["L2_pressure"] = metrics.L2_error_pressure;
    response["L2_entropy"] = metrics.L2_error_entropy;
    
    // Linf errors (max absolute error)
    response["Linf_density"] = metrics.Linf_error_density;
    response["Linf_velocity"] = metrics.Linf_error_velocity;
    response["Linf_pressure"] = metrics.Linf_error_pressure;
    response["Linf_entropy"] = metrics.Linf_error_entropy;
    
    // Legacy field for backwards compatibility
    response["L1_error"] = metrics.L1_error_density;
    
    std::cout << "Validation Results:" << std::endl;
    std::cout << "  L1 (rho):   " << metrics.L1_error_density << std::endl;
    std::cout << "  L2 (rho):   " << metrics.L2_error_density << std::endl;
    std::cout << "  Linf (rho): " << metrics.Linf_error_density << std::endl;

    res.set_content(response.dump(), "application/json");
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleValidate: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * GET /api/simulation/analytical
 * Get analytical solution for plotting overlay
 */
void handleGetAnalytical(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> GET /api/simulation/analytical" << std::endl;
    
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    std::vector<double> rho_exact, u_exact, p_exact;
    g_solver->computeAnalyticalSolution(rho_exact, u_exact, p_exact);

    const auto& x = g_solver->getPositions();

    json response;
    response["success"] = true;
    response["time"] = g_solver->getCurrentTime();
    response["numCells"] = rho_exact.size();
    
    response["x"] = json::array();
    response["density"] = json::array();
    response["velocity"] = json::array();
    response["pressure"] = json::array();
    
    for (size_t i = 0; i < rho_exact.size(); ++i) {
      response["x"].push_back(x[i]);
      response["density"].push_back(rho_exact[i]);
      response["velocity"].push_back(u_exact[i]);
      response["pressure"].push_back(p_exact[i]);
    }

    res.set_content(response.dump(), "application/json");
    std::cout << "Analytical solution sent" << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleGetAnalytical: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/toro/sensor-data
 * Get what sensors would read for corresponding Toro test at t = 0
 */
void handleGetToroSensorData(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/toro/sensor-data" << std::endl;

    auto body = json::parse(req.body);
    int test_number = body.value("test", 1);
    std::vector<double> positions = body["positions"];

    ShockSolver::ToroTest test;
    switch(test_number) {
      case 1: test = ShockSolver::ToroTest::TEST1_SOD; break;
      case 2: test = ShockSolver::ToroTest::TEST2_123; break;
      case 3: test = ShockSolver::ToroTest::TEST3_BLAST_LEFT; break;
      case 4: test = ShockSolver::ToroTest::TEST4_SLOW_SHOCK; break;
      case 5: test = ShockSolver::ToroTest::TEST5_COLLISION; break;
      default: test = ShockSolver::ToroTest::TEST1_SOD; 
    }

    auto sensor_data = ShockSolver::getSensorDataForToroTest(test, positions);

    json response;
    response["success"] = true;
    response["test"] = test_number;
    response["sensors"] = json::array();
    for (const auto& s : sensor_data) {
      json sensor;
      sensor["x"] = s.x;
      sensor["rho"] = s.rho;
      sensor["u"] = s.u;
      sensor["p"] = s.p;
      response["sensors"].push_back(sensor);
    }

    res.set_content(response.dump(), "application/json");
    std::cout << "Toro sensor data sent for " << positions.size() << " sensors" << std::endl;

  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleGetToroSensorData: " << e.what() << std::endl;
    res.status = 400;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/simulation/init/sparse
 * Initialize from sparse sensor data
 * 
 * Body: {
 *   "sensors": [{"x": 0.1, "rho": 1.0, "u": 0.0, "p": 1.0}, ...],
 *   "interpolation": "piecewise_constant",  // or "linear", "characteristic"
 *   "diaphragm_x": 0.5  // optional, -1 or omit for auto-detect
 * }
 */
void handleInitSparse(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/simulation/init/sparse" << std::endl;

    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    std::string interpolation = body.value("interpolation", "piecewise_constant");
    double diaphragm_x = body.value("diaphragm_x", -1.0);

    std::vector<ShockSolver::SparseDataPoint> sparse_data;
    for (const auto& sensor : body["sensors"]) {
      ShockSolver::SparseDataPoint point;
      point.x = sensor["x"];
      point.rho = sensor["rho"];
      point.u = sensor["u"];
      point.p = sensor["p"];
      sparse_data.push_back(point);
    }

    std::cout << "Initializing from " << sparse_data.size() << " sensors" << std::endl;
    std::cout << "Interpolation: " << interpolation << std::endl;
    std::cout << "Diaphragm x: " << (diaphragm_x < 0 ? "auto-detect" : std::to_string(diaphragm_x)) << std::endl;

    g_solver->initializeFromSparseData(sparse_data, interpolation, diaphragm_x);

    json response = successResponse("Initialized from sparse data");
    response["num_sensors"] = sparse_data.size();
    response["interpolation"] = interpolation;
    response["diaphragm_x"] = diaphragm_x;
    
    res.set_content(response.dump(), "application/json");
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleInitSparse: " << e.what() << std::endl;
    res.status = 400;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/simulation/reset
 */
void handleReset(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/simulation/reset" << std::endl;
    
    g_solver.reset();

    json response = successResponse("Simulation reset");
    res.set_content(response.dump(), "application/json");
    std::cout << "Simulation reset" << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleReset: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * ============================================================================
 * ANOMALY DETECTION ENDPOINTS
 * ============================================================================
 */

/**
 * POST /api/anomaly/configure
 * Configure anomaly detection thresholds and weights
 * 
 * Body: {
 *   "density_threshold": 0.05,
 *   "velocity_threshold": 0.10,
 *   "pressure_threshold": 0.05,
 *   "score_threshold": 0.1,
 *   "weight_density": 1.0,
 *   "weight_velocity": 0.5,
 *   "weight_pressure": 1.5
 * }
 */
void handleAnomalyConfig(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/anomaly/configure" << std::endl;
    
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    
    ShockSolver::AnomalyConfig config;
    config.density_threshold = body.value("density_threshold", 0.05);
    config.velocity_threshold = body.value("velocity_threshold", 0.10);
    config.pressure_threshold = body.value("pressure_threshold", 0.05);
    config.entropy_threshold = body.value("entropy_threshold", 0.05);
    config.score_threshold = body.value("score_threshold", 0.1);
    config.weight_density = body.value("weight_density", 1.0);
    config.weight_velocity = body.value("weight_velocity", 0.5);
    config.weight_pressure = body.value("weight_pressure", 1.5);
    config.weight_entropy = body.value("weight_entropy", 1.0);

    g_solver->setAnomalyConfig(config);

    json response = successResponse("Anomaly detection configured");
    response["config"] = {
      {"density_threshold", config.density_threshold},
      {"velocity_threshold", config.velocity_threshold},
      {"pressure_threshold", config.pressure_threshold},
      {"entropy_threshold", config.entropy_threshold},
      {"score_threshold", config.score_threshold},
      {"weight_density", config.weight_density},
      {"weight_velocity", config.weight_velocity},
      {"weight_pressure", config.weight_pressure},
      {"weight_entropy", config.weight_entropy}
    };

    res.set_content(response.dump(), "application/json");
    std::cout << "Anomaly config updated" << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleAnomalyConfig: " << e.what() << std::endl;
    res.status = 400;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * GET /api/anomaly/config
 * Get current anomaly detection configuration
 */
void handleGetAnomalyConfig(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> GET /api/anomaly/config" << std::endl;
    
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto config = g_solver->getAnomalyConfig();

    json response;
    response["success"] = true;
    response["config"] = {
      {"density_threshold", config.density_threshold},
      {"velocity_threshold", config.velocity_threshold},
      {"pressure_threshold", config.pressure_threshold},
      {"entropy_threshold", config.entropy_threshold},
      {"score_threshold", config.score_threshold},
      {"weight_density", config.weight_density},
      {"weight_velocity", config.weight_velocity},
      {"weight_pressure", config.weight_pressure},
      {"weight_entropy", config.weight_entropy}
    };

    res.set_content(response.dump(), "application/json");
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleGetAnomalyConfig: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/anomaly/check
 * Compare simulation predictions against provided sensor readings
 * 
 * Body: {
 *   "readings": [
 *     {"position": 0.2, "density": 1.0, "velocity": 0.0, "pressure": 1.0},
 *     ...
 *   ]
 * }
 * 
 * Returns detailed per-sensor anomaly analysis
 */
void handleAnomalyCheck(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/anomaly/check" << std::endl;
    
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    
    // Parse actual sensor readings from request
    std::vector<ShockSolver::SensorReading> actual_readings;
    for (const auto& r : body["readings"]) {
      ShockSolver::SensorReading reading;
      reading.position = r["position"];
      reading.density = r["density"];
      reading.velocity = r["velocity"];
      reading.pressure = r["pressure"];
      reading.entropy = r.value("entropy", 0.0);  // Optional, default to 0
      reading.time = g_solver->getCurrentTime();
      actual_readings.push_back(reading);
    }

    std::cout << "Checking " << actual_readings.size() << " sensor readings for anomalies" << std::endl;

    // Run anomaly detection
    auto results = g_solver->checkForAnomalies(actual_readings);
    
    // Build response
    json response;
    response["success"] = true;
    response["time"] = g_solver->getCurrentTime();
    response["anomaly_detected"] = false;
    response["results"] = json::array();
    
    int anomaly_count = 0;
    double max_score = 0.0;
    int max_score_idx = -1;
    
    for (const auto& r : results) {
      json result;
      result["sensor_index"] = r.sensor_index;
      result["position"] = r.position;
      
      // Predicted vs actual
      result["predicted"] = {
        {"density", r.predicted_density},
        {"velocity", r.predicted_velocity},
        {"pressure", r.predicted_pressure},
        {"entropy", r.predicted_entropy}
      };
      result["actual"] = {
        {"density", r.actual_density},
        {"velocity", r.actual_velocity},
        {"pressure", r.actual_pressure},
        {"entropy", r.actual_entropy}
      };
      
      // Residuals
      result["residual_density"] = r.residual_density;
      result["residual_velocity"] = r.residual_velocity;
      result["residual_pressure"] = r.residual_pressure;
      result["residual_entropy"] = r.residual_entropy;
      
      // Normalized errors (percentage)
      result["normalized_density"] = r.normalized_density;
      result["normalized_velocity"] = r.normalized_velocity;
      result["normalized_pressure"] = r.normalized_pressure;
      result["normalized_entropy"] = r.normalized_entropy;
      
      // Anomaly status
      result["anomaly_score"] = r.anomaly_score;
      result["is_anomalous"] = r.is_anomalous;
      result["primary_deviation"] = r.primary_deviation;
      
      if (r.is_anomalous) {
        response["anomaly_detected"] = true;
        anomaly_count++;
      }
      if (r.anomaly_score > max_score) {
        max_score = r.anomaly_score;
        max_score_idx = r.sensor_index;
      }
      
      response["results"].push_back(result);
    }
    
    // Summary statistics
    response["summary"] = {
      {"total_sensors", results.size()},
      {"anomalous_sensors", anomaly_count},
      {"max_anomaly_score", max_score},
      {"max_score_sensor_index", max_score_idx}
    };

    res.set_content(response.dump(), "application/json");
    
    std::cout << "Anomaly check complete: " << anomaly_count << "/" << results.size() 
              << " sensors flagged, max score = " << max_score << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleAnomalyCheck: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/anomaly/analyze
 * Full anomaly analysis with severity classification and system-level assessment
 * 
 * Body: {
 *   "readings": [...],       // Values to check (numerical solution)
 *   "reference": [...]       // Optional: ground truth to compare against (analytical solution)
 *                            // If not provided, compares against current simulation state
 * }
 * 
 * Returns summary with severity level
 */
void handleAnomalyAnalyze(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/anomaly/analyze" << std::endl;
    
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    
    // Parse actual sensor readings (numerical solution)
    std::vector<ShockSolver::SensorReading> actual_readings;
    for (const auto& r : body["readings"]) {
      ShockSolver::SensorReading reading;
      reading.position = r["position"];
      reading.density = r["density"];
      reading.velocity = r["velocity"];
      reading.pressure = r["pressure"];
      reading.entropy = r.value("entropy", 0.0);
      reading.time = g_solver->getCurrentTime();
      actual_readings.push_back(reading);
    }

    // Check if reference readings are provided (analytical solution)
    std::vector<ShockSolver::SensorReading> reference_readings;
    bool has_reference = body.contains("reference") && body["reference"].is_array();
    
    if (has_reference) {
      for (const auto& r : body["reference"]) {
        ShockSolver::SensorReading reading;
        reading.position = r["position"];
        reading.density = r["density"];
        reading.velocity = r["velocity"];
        reading.pressure = r["pressure"];
        reading.entropy = r.value("entropy", 0.0);
        reading.time = g_solver->getCurrentTime();
        reference_readings.push_back(reading);
      }
      std::cout << "Using provided reference (analytical) for comparison" << std::endl;
    } else {
      std::cout << "Using simulation state for comparison" << std::endl;
    }

    // Run analysis - either with reference or against simulation state
    ShockSolver::AnomalySummary summary;
    if (has_reference && reference_readings.size() == actual_readings.size()) {
      // Direct comparison: actual (numerical) vs reference (analytical)
      summary = g_solver->analyzeAnomaliesWithReference(actual_readings, reference_readings);
    } else {
      // Original behavior: compare against simulation state
      summary = g_solver->analyzeAnomalies(actual_readings);
    }
    
    // Build response
    json response;
    response["success"] = true;
    response["time"] = summary.time;
    response["comparison_mode"] = has_reference ? "numerical_vs_analytical" : "readings_vs_simulation";
    
    // Summary
    response["summary"] = {
      {"total_sensors", summary.total_sensors},
      {"anomalous_sensors", summary.anomalous_sensors},
      {"max_anomaly_score", summary.max_anomaly_score},
      {"max_score_sensor_index", summary.max_score_sensor_index},
      {"system_anomaly", summary.system_anomaly},
      {"severity", summary.severity}
    };
    
    // Detailed results
    response["sensor_results"] = json::array();
    for (const auto& r : summary.sensor_results) {
      json result;
      result["sensor_index"] = r.sensor_index;
      result["position"] = r.position;
      result["anomaly_score"] = r.anomaly_score;
      result["is_anomalous"] = r.is_anomalous;
      result["primary_deviation"] = r.primary_deviation;
      result["residuals"] = {
        {"density", r.residual_density},
        {"velocity", r.residual_velocity},
        {"pressure", r.residual_pressure},
        {"entropy", r.residual_entropy}
      };
      result["normalized_errors"] = {
        {"density", r.normalized_density},
        {"velocity", r.normalized_velocity},
        {"pressure", r.normalized_pressure},
        {"entropy", r.normalized_entropy}
      };
      response["sensor_results"].push_back(result);
    }

    res.set_content(response.dump(), "application/json");
    
    std::cout << "Anomaly analysis complete: severity = " << summary.severity 
              << ", " << summary.anomalous_sensors << "/" << summary.total_sensors 
              << " sensors flagged" << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleAnomalyAnalyze: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/anomaly/inject
 * Inject a simulated fault into sensor readings for testing
 * Useful for demonstrating anomaly detection capabilities
 * 
 * Body: {
 *   "sensor_index": 2,
 *   "fault_type": "pressure_drift",  // "pressure_drift", "density_spike", "velocity_bias", "random_noise"
 *   "magnitude": 0.2                 // Fractional magnitude of fault
 * }
 */
void handleAnomalyInject(const httplib::Request& req, httplib::Response& res) {
  try {
    std::cout << "\n>>> POST /api/anomaly/inject" << std::endl;
    
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    
    int sensor_index = body.value("sensor_index", 0);
    std::string fault_type = body.value("fault_type", "pressure_drift");
    double magnitude = body.value("magnitude", 0.2);

    // Get current sensor readings (truth)
    auto readings = g_solver->getSensorReadings();
    
    if (sensor_index < 0 || sensor_index >= static_cast<int>(readings.size())) {
      res.status = 400;
      res.set_content(errorResponse("Invalid sensor index").dump(), "application/json");
      return;
    }

    // Create faulted readings
    json response;
    response["success"] = true;
    response["fault_injected"] = {
      {"sensor_index", sensor_index},
      {"fault_type", fault_type},
      {"magnitude", magnitude}
    };
    
    response["original_reading"] = {
      {"density", readings[sensor_index].density},
      {"velocity", readings[sensor_index].velocity},
      {"pressure", readings[sensor_index].pressure},
      {"entropy", readings[sensor_index].entropy}
    };

    // Apply fault
    if (fault_type == "pressure_drift") {
      readings[sensor_index].pressure *= (1.0 + magnitude);
    } else if (fault_type == "density_spike") {
      readings[sensor_index].density *= (1.0 + magnitude);
    } else if (fault_type == "velocity_bias") {
      readings[sensor_index].velocity += magnitude * std::max(std::abs(readings[sensor_index].velocity), 1.0);
    } else if (fault_type == "entropy_spike") {
      // Entropy can be negative, so we add magnitude directly scaled by abs value
      readings[sensor_index].entropy += magnitude * std::max(std::abs(readings[sensor_index].entropy), 1.0);
    } else if (fault_type == "random_noise") {
      // Add noise to all variables
      readings[sensor_index].density *= (1.0 + magnitude * (2.0 * (rand() / (double)RAND_MAX) - 1.0));
      readings[sensor_index].velocity += magnitude * (2.0 * (rand() / (double)RAND_MAX) - 1.0);
      readings[sensor_index].pressure *= (1.0 + magnitude * (2.0 * (rand() / (double)RAND_MAX) - 1.0));
      readings[sensor_index].entropy += magnitude * (2.0 * (rand() / (double)RAND_MAX) - 1.0);
    }

    response["faulted_reading"] = {
      {"density", readings[sensor_index].density},
      {"velocity", readings[sensor_index].velocity},
      {"pressure", readings[sensor_index].pressure},
      {"entropy", readings[sensor_index].entropy}
    };

    // Return all readings with fault applied
    response["readings"] = json::array();
    for (const auto& r : readings) {
      json reading;
      reading["position"] = r.position;
      reading["density"] = r.density;
      reading["velocity"] = r.velocity;
      reading["pressure"] = r.pressure;
      reading["entropy"] = r.entropy;
      response["readings"].push_back(reading);
    }

    res.set_content(response.dump(), "application/json");
    
    std::cout << "Fault injected: " << fault_type << " with magnitude " << magnitude 
              << " at sensor " << sensor_index << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "ERROR in handleAnomalyInject: " << e.what() << std::endl;
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * POST /api/simulation/montecarlo
 */
void handleMonteCarlo(const httplib::Request& req, httplib::Response& res) {
  try {
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);

    // Parse parameters
    double noise_level = body.value("noiseLevel", 0.05);
    int num_trials = body.value("numTrials", 100);
    std::string interpolation = body.value("interpolation", "piecewise_constant");
    std::string time_integration_method = body.value("timeIntegration", "RK2");
    bool include_analytical = body.value("includeAnalytical", true);

    // Parse sensor data
    std::vector<ShockSolver::SparseDataPoint> sensors;
    if (body.contains("sensors")) {
      for (const auto& s : body["sensors"]) {
        ShockSolver::SparseDataPoint point;
        point.x = s["x"].get<double>();
        point.rho = s["rho"].get<double>();
        point.u = s["u"].get<double>();
        point.p = s["p"].get<double>();
        sensors.push_back(point);
      }
    }

    // Run Monte Carlo
    auto results = g_solver->runMonteCarloUncertainty(
      sensors, noise_level, num_trials, interpolation, time_integration_method, include_analytical
    );
    // Build response
    json response = successResponse("Monte Carlo completed");
    response["x"] = results.x;
    response["mean_rho"] = results.mean_rho;
    response["mean_u"] = results.mean_u;
    response["mean_p"] = results.mean_p;
    response["mean_entropy"] = results.mean_entropy;
    response["std_rho"] = results.std_rho;
    response["std_u"] = results.std_u;
    response["std_p"] = results.std_p;
    response["std_entropy"] = results.std_entropy;
    response["ci95_lower_rho"] = results.ci95_lower_rho;
    response["ci95_upper_rho"] = results.ci95_upper_rho;
    response["ci95_lower_u"] = results.ci95_lower_u;
    response["ci95_upper_u"] = results.ci95_upper_u;
    response["ci95_lower_p"] = results.ci95_lower_p;
    response["ci95_upper_p"] = results.ci95_upper_p;
    response["ci95_lower_s"] = results.ci95_lower_s;
    response["ci95_upper_s"] = results.ci95_upper_s;
    response["computation_time_ms"] = results.computation_time_ms;
    response["success"] = results.success;
    response["num_trials"] = results.num_trials;
    response["noise_level"] = results.noise_level;

    if (!results.analytical_rho.empty()) {
      response["analytical_rho"] = results.analytical_rho;
      response["analytical_u"] = results.analytical_u;
      response["analytical_p"] = results.analytical_p;
      response["analytical_entropy"] = results.analytical_entropy;
    }

    response["sensor_x"] = results.sensor_x;

    res.set_content(response.dump(), "application/json");
  } catch (const std::exception& e) {
    res.status = 400;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

/**
 * GET /api/health
 */
void handleHealth(const httplib::Request& req, httplib::Response& res) {
  json response;
  response["status"] = "healthy";
  response["solver_initialized"] = (g_solver != nullptr);
  res.set_content(response.dump(), "application/json");
}

/**
 * Main Server
 */
int main() {
  httplib::Server svr;

  svr.set_default_headers({
    {"Access-Control-Allow-Origin", "*"},
    {"Access-Control-Allow-Methods", "GET, POST, OPTIONS"},
    {"Access-Control-Allow-Headers", "Content-Type"}
  });

  svr.Options("/(.*)", [](const httplib::Request& req, httplib::Response& res) {
    res.status = 200;
  });

  std::cout << R"(
    ╔═══════════════════════════════════════════════════════════╗
    ║   Shock Tube Simulation API Server                        ║
    ║   Physics Engine: 1D Euler Equations                      ║
    ║   Features: HLLC/Entropy-Stable, MUSCL, Anomaly Detection ║
    ║   Port: 8080                                              ║
    ╚═══════════════════════════════════════════════════════════╝
  )" << std::endl;

  // Simulation endpoints
  svr.Post("/api/simulation/create", handleCreateSimulation);
  svr.Post("/api/simulation/init/shocktube", handleInitShockTube);
  svr.Post("/api/simulation/init/toro", handleInitToro);
  svr.Post("/api/simulation/init/sparse", handleInitSparse);
  svr.Post("/api/simulation/configure", handleConfigure);
  svr.Post("/api/simulation/sensors/place", handlePlaceSensors);
  svr.Post("/api/simulation/run", handleRun);
  svr.Post("/api/simulation/step", handleStep);
  svr.Post("/api/simulation/reset", handleReset);
  svr.Get("/api/simulation/status", handleGetStatus);
  svr.Get("/api/simulation/data", handleGetData);
  svr.Get("/api/simulation/sensors/readings", handleGetSensorReadings);
  svr.Get("/api/simulation/sensors/analytical", handleGetAnalyticalAtSensors);
  svr.Get("/api/simulation/validate", handleValidate);
  svr.Get("/api/simulation/analytical", handleGetAnalytical);
  svr.Post("/api/toro/sensor-data", handleGetToroSensorData);
  
  // Anomaly detection endpoints
  svr.Post("/api/anomaly/configure", handleAnomalyConfig);
  svr.Get("/api/anomaly/config", handleGetAnomalyConfig);
  svr.Post("/api/anomaly/check", handleAnomalyCheck);
  svr.Post("/api/anomaly/analyze", handleAnomalyAnalyze);
  svr.Post("/api/anomaly/inject", handleAnomalyInject);

  // Monte Carlo
  svr.Post("/api/simulation/montecarlo", handleMonteCarlo);
  
  // Health check
  svr.Get("/api/health", handleHealth);

  std::cout << "\nAPI Endpoints:" << std::endl;
  std::cout << "\n  === Simulation ===" << std::endl;
  std::cout << "  POST   /api/simulation/create         - Create new simulation" << std::endl;
  std::cout << "  POST   /api/simulation/init/shocktube - Initialize shock tube" << std::endl;
  std::cout << "  POST   /api/simulation/init/toro      - Initialize Toro tests" << std::endl;
  std::cout << "  POST   /api/simulation/init/sparse    - Initialize from sparse sensor data" << std::endl;
  std::cout << "  POST   /api/simulation/configure      - Configure flux type" << std::endl;
  std::cout << "  POST   /api/simulation/sensors/place  - Place sensors" << std::endl;
  std::cout << "  POST   /api/simulation/run            - Run full simulation" << std::endl;
  std::cout << "  POST   /api/simulation/step           - Step N timesteps" << std::endl;
  std::cout << "  POST   /api/simulation/reset          - Reset simulation" << std::endl;
  std::cout << "  GET    /api/simulation/status         - Get status" << std::endl;
  std::cout << "  GET    /api/simulation/data           - Get full solution" << std::endl;
  std::cout << "  GET    /api/simulation/sensors/readings - Get sensor readings" << std::endl;
  std::cout << "  GET    /api/simulation/sensors/analytical - Get analytical solution at sensors" << std::endl;
  std::cout << "  GET    /api/simulation/validate       - Validate vs analytical" << std::endl;
  std::cout << "  GET    /api/simulation/analytical     - Get analytical solution" << std::endl;
  std::cout << "  POST   /api/toro/sensor-data          - Get sensor data for Toro test" << std::endl;
  
  std::cout << "\n  === Anomaly Detection ===" << std::endl;
  std::cout << "  POST   /api/anomaly/configure         - Configure thresholds" << std::endl;
  std::cout << "  GET    /api/anomaly/config            - Get current config" << std::endl;
  std::cout << "  POST   /api/anomaly/check             - Check readings for anomalies" << std::endl;
  std::cout << "  POST   /api/anomaly/analyze           - Full analysis with severity" << std::endl;
  std::cout << "  POST   /api/anomaly/inject            - Inject fault for testing" << std::endl;
  
  std::cout << "\n  === System ===" << std::endl;
  std::cout << "  GET    /api/health                    - Health check" << std::endl;

  std::cout << "\n  === Monte Carlo Uncertainty Propogation ===" << std::endl;
  std::cout << "  POST    /api/simulation/montecarlo    - Monte Carlo Uncertainty Propogation in sensor errors" << std::endl;
  
  std::cout << "\nServer starting on http://localhost:8080" << std::endl;
  std::cout << "Waiting for requests...\n" << std::endl;

  if (!svr.listen("0.0.0.0", 8080)) {
    std::cerr << "Failed to start server" << std::endl;
    return 1;
  }

  return 0;
}