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
    cfg.CFL = body.value("CFL", 0.9);
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

    std::cout << "Initializing shock tube:" << std::endl;
    std::cout << "  Left:  rho=" << rho_L << ", u=" << u_L << ", p=" << p_L << std::endl;
    std::cout << "  Right: rho=" << rho_R << ", u=" << u_R << ", p=" << p_R << std::endl;
    std::cout << "  Diaphragm: " << x_diaphragm << std::endl;

    g_solver->initializeShockTube(rho_L, u_L, p_L, rho_R, u_R, p_R, x_diaphragm);
    std::cout << "Shock tube initialized" << std::endl;

    json response = successResponse("Shock tube initialized");
    response["initial_conditions"] = {
      {"left", {{"rho", rho_L}, {"u", u_L}, {"p", p_L}}},
      {"right", {{"rho", rho_R}, {"u", u_R}, {"p", p_R}}},
      {"diaphragm", x_diaphragm}
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
      case 4: test = ShockSolver::ToroTest::TEST4_COLLISION; break;
      case 5: test = ShockSolver::ToroTest::TEST5_STATIONARY; break;
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
    std::cout << "Request body: " << body.dump(2) << std::endl;

    std::string scheme = body.value("scheme", "godunov");
    std::string riemann_solver = body.value("riemann_solver", "hllc");

    std::cout << "Configuring:" << std::endl;
    std::cout << "  Scheme: " << scheme << std::endl;
    std::cout << "  Riemann solver: " << riemann_solver << std::endl;

    g_solver->setScheme(scheme);
    g_solver->setRiemannSolver(riemann_solver);
    std::cout << "Configuration updated" << std::endl;

    json response = successResponse("Configuration updated");
    response["scheme"] = scheme;
    response["riemann_solver"] = riemann_solver;

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
    
    // L2 errors (RMS error)
    response["L2_density"] = metrics.L2_error_density;
    response["L2_velocity"] = metrics.L2_error_velocity;
    response["L2_pressure"] = metrics.L2_error_pressure;
    
    // Linf errors (max absolute error)
    response["Linf_density"] = metrics.Linf_error_density;
    response["Linf_velocity"] = metrics.Linf_error_velocity;
    response["Linf_pressure"] = metrics.Linf_error_pressure;
    
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
      case 4: test = ShockSolver::ToroTest::TEST4_COLLISION; break;
      case 5: test = ShockSolver::ToroTest::TEST5_STATIONARY; break;
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
    std::string interpolation = body.value("interpolation", "linear");

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

    g_solver->initializeFromSparseData(sparse_data, interpolation);

    json response = successResponse("Initialized from sparse data");
    response["num_sensors"] = sparse_data.size();
    response["interpolation"] = interpolation;
    
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
    ╔═══════════════════════════════════════════╗
    ║   Shock Tube Simulation API Server        ║
    ║   Physics Engine: 1D Euler Equations      ║
    ║   Port: 8080                              ║
    ╚═══════════════════════════════════════════╝
  )" << std::endl;

  svr.Post("/api/simulation/create", handleCreateSimulation);
  svr.Post("/api/simulation/init/shocktube", handleInitShockTube);
  svr.Post("/api/simulation/init/toro", handleInitToro);
  svr.Post("/api/simulation/configure", handleConfigure);
  svr.Post("/api/simulation/run", handleRun);
  svr.Post("/api/simulation/step", handleStep);
  svr.Post("/api/simulation/sensors/place", handlePlaceSensors);
  svr.Get("/api/simulation/sensors/readings", handleGetSensorReadings);
  svr.Post("/api/simulation/reset", handleReset);
  svr.Get("/api/simulation/status", handleGetStatus);
  svr.Get("/api/simulation/data", handleGetData);
  svr.Get("/api/simulation/validate", handleValidate);
  svr.Get("/api/simulation/analytical", handleGetAnalytical);
  svr.Post("/api/toro/sensor-data", handleGetToroSensorData);
  svr.Post("/api/simulation/init/sparse", handleInitSparse);
  svr.Get("/api/health", handleHealth);

  std::cout << "\nAPI Endpoints:" << std::endl;
  std::cout << "  POST   /api/simulation/create         - Create new simulation" << std::endl;
  std::cout << "  POST   /api/simulation/init/shocktube - Initialize shock tube" << std::endl;
  std::cout << "  POST   /api/simulation/init/toro      - Initialize Toro tests" << std::endl;
  std::cout << "  POST   /api/simulation/init/sparse    - Initialize from sparse sensor data" << std::endl;
  std::cout << "  POST   /api/simulation/configure      - Set scheme/solver" << std::endl;
  std::cout << "  POST   /api/simulation/sensors/place  - Place sensors" << std::endl;
  std::cout << "  POST   /api/simulation/run            - Run full simulation" << std::endl;
  std::cout << "  POST   /api/simulation/step           - Step N timesteps" << std::endl;
  std::cout << "  GET    /api/simulation/sensors/readings - Get sensor readings" << std::endl;
  std::cout << "  POST   /api/simulation/reset          - Reset simulation" << std::endl;
  std::cout << "  GET    /api/simulation/status         - Get status" << std::endl;
  std::cout << "  GET    /api/simulation/data           - Get full solution" << std::endl;
  std::cout << "  GET    /api/simulation/validate       - Validate vs analytical" << std::endl;
  std::cout << "  GET    /api/simulation/analytical     - Get analytical solution" << std::endl;
  std::cout << "  POST   /api/toro/sensor-data          - Get sensor data according to Toro test" << std::endl;
  std::cout << "  GET    /api/health                    - Health check" << std::endl;
  std::cout << "\nServer starting on http://localhost:8080" << std::endl;
  std::cout << "Waiting for requests...\n" << std::endl;

  if (!svr.listen("0.0.0.0", 8080)) {
    std::cerr << "Failed to start server" << std::endl;
    return 1;
  }

  return 0;
}