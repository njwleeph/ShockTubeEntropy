#include "solver.hpp"
#include "cpp-httplib/httplib.h"
#include <nlohmann/json.hpp>
#include <iostream>
#include <memory>

using json = nlohmann::json;

std::unique_ptr<ShockSolver> g_solver = nullptr;

json errorResponse(const std::string& message) {
  return {{"success", false}, {"error", message}};
}

json successResponse(const std::string& message) {
  return {{"success", true}, {"message", message}};
}

ShockSolver::ToroTest parseToroTest(int test) {
  switch(test) {
    case 1: return ShockSolver::ToroTest::TEST1_SOD;
    case 2: return ShockSolver::ToroTest::TEST2_123;
    case 3: return ShockSolver::ToroTest::TEST3_BLAST_LEFT;
    case 4: return ShockSolver::ToroTest::TEST4_SLOW_SHOCK;
    case 5: return ShockSolver::ToroTest::TEST5_COLLISION;
    default: return ShockSolver::ToroTest::TEST1_SOD;
  }
}

/**
 * API Endpoints
 */

// POST /api/simulation/create
void handleCreate(const httplib::Request& req, httplib::Response& res) {
  try {
    auto body = json::parse(req.body);

    ShockSolver::Config cfg;
    cfg.length = body.value("length", 1.0);
    cfg.numCells = body.value("numCells", 500);
    cfg.gamma = body.value("gamma", 1.4);
    cfg.CFL = body.value("CFL", 0.5);
    cfg.endTime = body.value("endTime", 0.25);

    g_solver = std::make_unique<ShockSolver>(cfg);

    res.set_content(successResponse("Solver created").dump(), "application/json");
  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

// POST /api/simulation/reset
void handleReset(const httplib::Request& req, httplib::Response& res) {
  try {
    if (g_solver) g_solver.reset();
    res.set_content(successResponse("Simulation reset").dump(), "application/json");

  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

// POST /api/simulation/configure
void handleConfigure(const httplib::Request& req, httplib::Response& res) {
  try {
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);

    if (body.contains("flux")) g_solver->setFlux(body["flux"].get<std::string>());
    if (body.contains("timeIntegration")) g_solver->setTimeIntegration(body["timeIntegration"].get<std::string>());

    res.set_content(successResponse("Configuration updated").dump(), "application/json");
    
  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

// POST /api/simulation/initialize/toro
void handleInitToro(const httplib::Request& req, httplib::Response& res) {
  try {
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Solver not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);
    int test = body.value("test", 1);

    g_solver->initializeToroTest(parseToroTest(test));

    res.set_content(successResponse("Toro test initialized").dump(), "application/json");

  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

// POST /api/simulation/sparse-init
void handleSparseInit(const httplib::Request& req, httplib::Response& res) {
  try {
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);

    std::vector<ShockSolver::SparseDataPoint> sensors;
    for (const auto& s : body["sensors"]) {
      ShockSolver::SparseDataPoint pt;
      pt.x = s["x"].get<double>();
      pt.rho = s["rho"].get<double>();
      pt.u = s["u"].get<double>();
      pt.p = s["p"].get<double>();
      sensors.push_back(pt);
    }

    std::string interpolation = body.value("interpolation", "linear");

    g_solver->initializeFromSparseData(sensors, interpolation);

    res.set_content(successResponse("Sparse initialization complete").dump(), "application/json");

  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

// POST /api/simulation/run
void handleRun(const httplib::Request& req, httplib::Response& res) {
  try {
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Solver not created").dump(), "application/json");
      return;
    }

    g_solver->run();

    json response = successResponse("Simulation complete");
    response["time"] = g_solver->getCurrentTime();
    response["steps"] = g_solver->getStepCount();

    res.set_content(response.dump(), "application/json");

  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

// GET /api/simulation/data
void handleGetData(const httplib::Request& req, httplib::Response& res) {
  try {
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Solver not created").dump(), "application/json");
      return;
    }

    const auto& rho = g_solver->getDensity();
    const auto& u = g_solver->getVelocity();
    const auto& p = g_solver->getPressure();
    const auto& x = g_solver->getPositions();
    auto entropy = g_solver->getEntropy();

    json response;
    response["success"] = true;
    response["time"] = g_solver->getCurrentTime();
    response["x"] = x;
    response["density"] = rho;
    response["velocity"] = u;
    response["pressure"] = p;
    response["entropy"] = entropy;

    res.set_content(response.dump(), "application/json");

  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

// GET /api/simulation/analytical
void handleGetAnalytical(const httplib::Request& req, httplib::Response& res) {
  try {
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    } 

    std::vector<double> rho_exact, u_exact, p_exact, entropy_exact;
    g_solver->computeAnalyticalSolution(rho_exact, u_exact, p_exact, entropy_exact);

    const auto& x = g_solver->getPositions();

    json response;
    response["success"] = true;
    response["x"] = x;
    response["density"] = rho_exact;
    response["velocity"] = u_exact;
    response["pressure"] = p_exact;
    response["entropy"] = entropy_exact;

    res.set_content(response.dump(), "application/json");

  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()), "application/json");
  }
}

// GET /api/simulation/validate
void handleValidate(const httplib::Request& req, httplib::Response& res) {
  try {
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto metrics = g_solver->validateAgainstAnalytical();

    json response;
    response["success"] = true;
    response["L1_error_density"] = metrics.L1_error_density;
    response["L1_error_velocity"] = metrics.L1_error_velocity;
    response["L1_error_pressure"] = metrics.L1_error_pressure;
    response["L1_error_entropy"] = metrics.L1_error_entropy;

    response["L2_error_density"] = metrics.L2_error_density;
    response["L2_error_velocity"] = metrics.L2_error_velocity;
    response["L2_error_pressure"] = metrics.L2_error_pressure;
    response["L2_error_entropy"] = metrics.L2_error_entropy;

    response["Linf_error_density"] = metrics.Linf_error_density;
    response["Linf_error_velocity"] = metrics.Linf_error_velocity;
    response["Linf_error_pressure"] = metrics.Linf_error_pressure;
    response["Linf_error_entropy"] = metrics.Linf_error_entropy;

    response["num_points_validated"] = metrics.num_points_validated;

    res.set_content(response.dump(), "application/json");

  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

// POST /api/simulation/montecarlo
void handleMonteCarlo(const httplib::Request& req, httplib::Response& res) {
  try {
    if (!g_solver) {
      res.status = 400;
      res.set_content(errorResponse("Simulation not created").dump(), "application/json");
      return;
    }

    auto body = json::parse(req.body);

    std::vector<ShockSolver::SparseDataPoint> sensors;
    for (const auto& s : body["sensors"]) {
      ShockSolver::SparseDataPoint point;
      point.x = s["x"].get<double>();
      point.rho = s["rho"].get<double>();
      point.u = s["u"].get<double>();
      point.p = s["p"].get<double>();
      sensors.push_back(point);
    }

    double noise_level = body.value("noiseLevel", 0.05);
    int num_trials = body.value("numTrials", 100);
    std::string interpolation = body.value("interpolation", "piecewise_constant");
    std::string time_integration_method = body.value("timeIntegration", "RK2");

    auto results = g_solver->runMonteCarloUncertainty(
      sensors, noise_level, num_trials, interpolation, time_integration_method
    );

    json response;
    response["success"] = results.success;
    response["noise_level"] = results.noise_level;
    response["num_trials"] = results.num_trials;
    response["computation_time_ms"] = results.computation_time_ms;

    response["x"] = results.x;
    response["sensor_x"] = results.sensor_x;

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
    response["ci95_lower_entropy"] = results.ci95_lower_entropy;
    response["ci95_upper_entropy"] = results.ci95_upper_entropy;

    response["analytical_rho"] = results.analytical_rho;
    response["analytical_u"] = results.analytical_u;
    response["analytical_p"] = results.analytical_p;
    response["analytical_entropy"] = results.analytical_entropy;

    res.set_content(response.dump(), "application/json");

  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

// POST /api/toro/sensor_data
void handleGetToroSensorData(const httplib::Request& req, httplib::Response& res) {
  try {
    auto body = json::parse(req.body);

    int test = body.value("test", 1);

    std::vector<double> positions;
    for (const auto& p : body["positions"]) positions.push_back(p.get<double>());

    auto sensors = ShockSolver::getSensorDataForToroTest(parseToroTest(test), positions);

    json response;
    response["success"] = true;
    response["sensors"] = json::array();

    for (const auto& s : sensors) {
      response["sensors"].push_back({
        {"x", s.x},
        {"rho", s.rho},
        {"u", s.u},
        {"p", s.p}
      });
    }

    res.set_content(response.dump(), "application/json");

  } catch (const std::exception& e) {
    res.status = 500;
    res.set_content(errorResponse(e.what()).dump(), "application/json");
  }
}

// GET /api/health
void handleHealth(const httplib::Request& req, httplib::Response& res) {
  json response;
  response["status"] = "ok";
  response["solver_initialized"] = (g_solver != nullptr);
  res.set_content(response.dump(), "application/json");
}

/**
 * MAIN
 */
int main() {
  httplib::Server svr;

  // CORS
  svr.set_pre_routing_handler([](const httplib::Request& req, httplib::Response& res) {
    res.set_header("Access-Control-Allow-Origin", "*");
    res.set_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
    res.set_header("Access-Control-Allow-Headers", "Content-Type");

    if (req.method == "OPTIONS") {
      res.status = 204;
      return httplib::Server::HandlerResponse::Handled;
    }
    return httplib::Server::HandlerResponse::Unhandled;
  });

  // Register endpoints
  svr.Post("/api/simulation/create", handleCreate);
  svr.Post("/api/simulation/reset", handleReset);
  svr.Post("/api/simulation/configure", handleConfigure);
  svr.Post("/api/simulation/initialize/toro", handleInitToro);
  svr.Post("/api/simulation/sparse-init", handleSparseInit);
  svr.Post("/api/simulation/run", handleRun);
  svr.Get("/api/simulation/data", handleGetData);
  svr.Get("/api/simulation/analytical", handleGetAnalytical);
  svr.Get("/api/simulation/validate", handleValidate);
  svr.Post("/api/simulation/montecarlo", handleMonteCarlo);
  svr.Post("/api/toro/sensor-data", handleGetToroSensorData);
  svr.Get("/api/health", handleHealth);

  std::cout << "\n";
  std::cout << "========================================\n";
  std::cout << "  Shock Tube Entropy Simulator\n";
  std::cout << "========================================\n";
  std::cout << "\n";
  std::cout << "API Endpoints:\n";
  std::cout << "  POST /api/simulation/create        - Create solver\n";
  std::cout << "  POST /api/simulation/reset         - Reset simulation\n";
  std::cout << "  POST /api/simulation/configure     - Set flux/time integration\n";
  std::cout << "  POST /api/simulation/initialize/toro - Init Toro test\n";
  std::cout << "  POST /api/simulation/sparse-init   - Init from sensors\n";
  std::cout << "  POST /api/simulation/run           - Run simulation\n";
  std::cout << "  GET  /api/simulation/data          - Get numerical solution\n";
  std::cout << "  GET  /api/simulation/analytical    - Get exact solution\n";
  std::cout << "  GET  /api/simulation/validate      - Validate vs analytical\n";
  std::cout << "  POST /api/simulation/montecarlo    - Run Monte Carlo\n";
  std::cout << "  POST /api/toro/sensor-data         - Get sensor data for Toro test\n";
  std::cout << "  GET  /api/health                   - Health check\n";
  std::cout << "\n";
  std::cout << "Server starting on http://localhost:8080\n";
  std::cout << "\n";

  svr.listen("0.0.0.0", 8080);

  return 0;
}