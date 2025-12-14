import React, { useState, useRef, useEffect } from 'react';

const API_BASE = 'http://localhost:8080';

// Toro test configurations
const TORO_TESTS = {
  1: {
    name: 'Test 1: Sod Shock Tube',
    endTime: 0.25,
    left: { rho: 1.0, u: 0.0, p: 1.0 },
    right: { rho: 0.125, u: 0.0, p: 0.1 },
    x0: 0.5,
  },
  2: {
    name: 'Test 2: 123 Problem',
    endTime: 0.15,
    left: { rho: 1.0, u: -2.0, p: 0.4 },
    right: { rho: 1.0, u: 2.0, p: 0.4 },
    x0: 0.5,
  },
  3: {
    name: 'Test 3: Left Blast Wave',
    endTime: 0.012,
    left: { rho: 1.0, u: 0.0, p: 1000.0 },
    right: { rho: 1.0, u: 0.0, p: 0.01 },
    x0: 0.5,
  },
  4: {
    name: 'Test 4: Right Blast Wave',
    endTime: 0.035,
    left: { rho: 1.0, u: 0.0, p: 0.01 },
    right: { rho: 1.0, u: 0.0, p: 100.0 },
    x0: 0.5,
  },
  5: {
    name: 'Test 5: Collision',
    endTime: 0.035,
    left: { rho: 5.99924, u: 19.5975, p: 460.894 },
    right: { rho: 5.99242, u: -6.19633, p: 46.0950 },
    x0: 0.5,
  },
};

const DEFAULT_SENSOR_POSITIONS = [0.1, 0.25, 0.4, 0.6, 0.75, 0.9];

function generateSensorsFromTest(testId) {
  const test = TORO_TESTS[testId];
  return DEFAULT_SENSOR_POSITIONS.map(x => ({
    x,
    rho: x < test.x0 ? test.left.rho : test.right.rho,
    u: x < test.x0 ? test.left.u : test.right.u,
    p: x < test.x0 ? test.left.p : test.right.p,
  }));
}

export default function App() {
  // Test selection
  const [selectedTest, setSelectedTest] = useState(1);
  
  // Solver config
  const [config, setConfig] = useState({
    numCells: 500,
    gamma: 1.4,
    CFL: 0.5,
    flux: 'HLLC',
    timeIntegration: 'RK2',
    interpolation: 'piecewise_constant',
  });

  // Monte Carlo settings
  const [mcSettings, setMcSettings] = useState({
    numTrials: 100,
    noiseLevel: 0.05,
  });

  // Sensors
  const [sensors, setSensors] = useState(generateSensorsFromTest(1));

  // Results
  const [mcData, setMcData] = useState(null);
  const [numericalData, setNumericalData] = useState(null);
  const [analyticalData, setAnalyticalData] = useState(null);
  const [errors, setErrors] = useState(null);

  // UI state
  const [isRunning, setIsRunning] = useState(false);
  const [statusMessage, setStatusMessage] = useState('Ready');

  // Canvas refs
  const canvasRefs = {
    density: useRef(null),
    velocity: useRef(null),
    pressure: useRef(null),
    entropy: useRef(null),
  };

  // Load sensors when test changes
  const loadTestSensors = async () => {
    try {
      // Use backend API to get sensor data for selected test
      const positions = sensors.map(s => s.x);
      const response = await fetch(`${API_BASE}/api/toro/sensor-data`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ test: selectedTest, positions }),
      });
      const data = await response.json();
      if (data.success) {
        setSensors(data.sensors);
        setStatusMessage(`Loaded sensors for ${TORO_TESTS[selectedTest].name}`);
      } else {
        // Fallback to local generation
        setSensors(generateSensorsFromTest(selectedTest));
        setStatusMessage(`Loaded sensors locally for ${TORO_TESTS[selectedTest].name}`);
      }
    } catch (err) {
      // Fallback to local generation
      setSensors(generateSensorsFromTest(selectedTest));
      setStatusMessage(`Loaded sensors locally for ${TORO_TESTS[selectedTest].name}`);
    }
  };

  // Update sensor
  const updateSensor = (index, field, value) => {
    const newSensors = [...sensors];
    newSensors[index] = { ...newSensors[index], [field]: parseFloat(value) || 0 };
    setSensors(newSensors);
  };

  // Add/remove sensors
  const addSensor = () => {
    const lastX = sensors[sensors.length - 1]?.x || 0;
    const test = TORO_TESTS[selectedTest];
    const newX = Math.min(lastX + 0.1, 0.95);
    setSensors([...sensors, {
      x: newX,
      rho: newX < test.x0 ? test.left.rho : test.right.rho,
      u: newX < test.x0 ? test.left.u : test.right.u,
      p: newX < test.x0 ? test.left.p : test.right.p,
    }]);
  };

  const removeSensor = (index) => {
    if (sensors.length > 2) {
      setSensors(sensors.filter((_, i) => i !== index));
    }
  };

  // Run Monte Carlo simulation
  const runSimulation = async () => {
    setIsRunning(true);
    setStatusMessage('Creating solver...');

    try {
      const endTime = TORO_TESTS[selectedTest].endTime;

      // 1. Create solver
      const createRes = await fetch(`${API_BASE}/api/simulation/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          length: 1.0,
          numCells: config.numCells,
          gamma: config.gamma,
          CFL: config.CFL,
          endTime,
        }),
      });
      if (!createRes.ok) throw new Error('Failed to create solver');

      // 2. Configure flux and time integration
      setStatusMessage('Configuring solver...');
      const configRes = await fetch(`${API_BASE}/api/simulation/configure`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          flux: config.flux,
          timeIntegration: config.timeIntegration,
        }),
      });
      if (!configRes.ok) throw new Error('Failed to configure solver');

      // 3. Initialize from sensors (sets initial conditions for analytical)
      setStatusMessage('Initializing from sensors...');
      const initRes = await fetch(`${API_BASE}/api/simulation/sparse-init`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sensors,
          interpolation: config.interpolation,
        }),
      });
      if (!initRes.ok) throw new Error('Failed to initialize from sensors');

      // 4. Run Monte Carlo
      setStatusMessage(`Running Monte Carlo (${mcSettings.numTrials} trials)...`);
      const mcResponse = await fetch(`${API_BASE}/api/simulation/montecarlo`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sensors,
          noiseLevel: mcSettings.noiseLevel,
          numTrials: mcSettings.numTrials,
          interpolation: config.interpolation,
          timeIntegration: config.timeIntegration,
        }),
      });
      const mcResult = await mcResponse.json();
      console.log('MC Result:', mcResult);
      setMcData(mcResult);

      // 5. Re-initialize and run single numerical simulation for comparison
      setStatusMessage('Running numerical simulation...');
      
      // Reset solver
      await fetch(`${API_BASE}/api/simulation/reset`, { method: 'POST' });
      
      // Recreate
      await fetch(`${API_BASE}/api/simulation/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          length: 1.0,
          numCells: config.numCells,
          gamma: config.gamma,
          CFL: config.CFL,
          endTime,
        }),
      });

      await fetch(`${API_BASE}/api/simulation/configure`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          flux: config.flux,
          timeIntegration: config.timeIntegration,
        }),
      });

      await fetch(`${API_BASE}/api/simulation/sparse-init`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sensors,
          interpolation: config.interpolation,
        }),
      });

      await fetch(`${API_BASE}/api/simulation/run`, { method: 'POST' });

      const numResponse = await fetch(`${API_BASE}/api/simulation/data`);
      const numResult = await numResponse.json();
      console.log('Numerical Result:', numResult);
      setNumericalData(numResult);

      // 6. Set analytical from MC response
      if (mcResult.analytical_rho) {
        setAnalyticalData({
          x: mcResult.x,
          density: mcResult.analytical_rho,
          velocity: mcResult.analytical_u,
          pressure: mcResult.analytical_p,
          entropy: mcResult.analytical_entropy,
        });
      }

      // 7. Get validation errors
      const valResponse = await fetch(`${API_BASE}/api/simulation/validate`);
      const valResult = await valResponse.json();
      console.log('Validation Result:', valResult);
      setErrors(valResult);

      setStatusMessage(`Complete! MC: ${mcResult.computation_time_ms?.toFixed(0)}ms`);

    } catch (err) {
      console.error('Simulation error:', err);
      setStatusMessage(`Error: ${err.message}`);
    } finally {
      setIsRunning(false);
    }
  };

  // Drawing function
  const drawPlot = (canvasRef, title, yLabel, data) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    const W = canvas.width;
    const H = canvas.height;

    // Layout
    const margin = { left: 65, right: 110, top: 35, bottom: 45 };
    const plotW = W - margin.left - margin.right;
    const plotH = H - margin.top - margin.bottom;

    // Clear
    ctx.fillStyle = '#fff';
    ctx.fillRect(0, 0, W, H);

    // Compute Y range from all data
    let yMin = Infinity, yMax = -Infinity;
    
    const updateRange = (arr) => {
      if (!arr) return;
      arr.forEach(v => {
        if (isFinite(v)) {
          yMin = Math.min(yMin, v);
          yMax = Math.max(yMax, v);
        }
      });
    };

    updateRange(data.mcMean);
    updateRange(data.mcLower);
    updateRange(data.mcUpper);
    updateRange(data.numerical);
    updateRange(data.analytical);

    if (!isFinite(yMin) || !isFinite(yMax)) {
      yMin = 0; yMax = 1;
    }

    const yPad = (yMax - yMin) * 0.1 || 0.1;
    yMin -= yPad;
    yMax += yPad;

    // Transform
    const toX = (x) => margin.left + x * plotW;
    const toY = (y) => margin.top + (1 - (y - yMin) / (yMax - yMin)) * plotH;

    // Grid
    ctx.strokeStyle = '#e5e7eb';
    ctx.lineWidth = 1;
    for (let i = 0; i <= 5; i++) {
      const x = toX(i / 5);
      ctx.beginPath();
      ctx.moveTo(x, margin.top);
      ctx.lineTo(x, margin.top + plotH);
      ctx.stroke();

      const y = margin.top + (i / 5) * plotH;
      ctx.beginPath();
      ctx.moveTo(margin.left, y);
      ctx.lineTo(margin.left + plotW, y);
      ctx.stroke();
    }

    // Sensor positions
    if (data.sensorX) {
      ctx.strokeStyle = '#22c55e';
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 4]);
      data.sensorX.forEach(sx => {
        const x = toX(sx);
        ctx.beginPath();
        ctx.moveTo(x, margin.top);
        ctx.lineTo(x, margin.top + plotH);
        ctx.stroke();
      });
      ctx.setLineDash([]);
    }

    // 95% CI band
    if (data.x && data.mcLower && data.mcUpper && data.mcLower.length > 0) {
      ctx.fillStyle = 'rgba(59, 130, 246, 0.2)';
      ctx.beginPath();
      ctx.moveTo(toX(data.x[0]), toY(data.mcLower[0]));
      for (let i = 1; i < data.x.length; i++) {
        ctx.lineTo(toX(data.x[i]), toY(data.mcLower[i]));
      }
      for (let i = data.x.length - 1; i >= 0; i--) {
        ctx.lineTo(toX(data.x[i]), toY(data.mcUpper[i]));
      }
      ctx.closePath();
      ctx.fill();
    }

    // Helper to draw line
    const drawLine = (xArr, yArr, color, dash = [], lineWidth = 2) => {
      if (!xArr || !yArr || xArr.length === 0) return;
      ctx.strokeStyle = color;
      ctx.lineWidth = lineWidth;
      ctx.setLineDash(dash);
      ctx.beginPath();
      let started = false;
      for (let i = 0; i < xArr.length; i++) {
        if (isFinite(yArr[i])) {
          const px = toX(xArr[i]);
          const py = toY(yArr[i]);
          if (!started) {
            ctx.moveTo(px, py);
            started = true;
          } else {
            ctx.lineTo(px, py);
          }
        }
      }
      ctx.stroke();
      ctx.setLineDash([]);
    };

    // MC Mean (blue)
    if (data.mcMean) drawLine(data.x, data.mcMean, '#3b82f6', [], 2);

    // Numerical (red)
    if (data.numerical) drawLine(data.numX || data.x, data.numerical, '#ef4444', [], 2);

    // Analytical (gray dashed)
    if (data.analytical) drawLine(data.x, data.analytical, '#6b7280', [6, 4], 2);

    // Axes
    ctx.strokeStyle = '#374151';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(margin.left, margin.top);
    ctx.lineTo(margin.left, margin.top + plotH);
    ctx.lineTo(margin.left + plotW, margin.top + plotH);
    ctx.stroke();

    // X-axis labels
    ctx.fillStyle = '#374151';
    ctx.font = '11px system-ui';
    ctx.textAlign = 'center';
    for (let i = 0; i <= 5; i++) {
      const val = i / 5;
      ctx.fillText(val.toFixed(1), toX(val), margin.top + plotH + 15);
    }
    ctx.font = '12px system-ui';
    ctx.fillText('Position (m)', margin.left + plotW / 2, H - 8);

    // Y-axis labels
    ctx.textAlign = 'right';
    ctx.font = '11px system-ui';
    for (let i = 0; i <= 5; i++) {
      const val = yMin + (yMax - yMin) * (1 - i / 5);
      const y = margin.top + (i / 5) * plotH;
      ctx.fillText(val.toPrecision(3), margin.left - 8, y + 4);
    }

    // Y-axis label
    ctx.save();
    ctx.translate(14, margin.top + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.font = '12px system-ui';
    ctx.fillText(yLabel, 0, 0);
    ctx.restore();

    // Title
    ctx.fillStyle = '#111827';
    ctx.font = 'bold 14px system-ui';
    ctx.textAlign = 'center';
    ctx.fillText(title, margin.left + plotW / 2, 20);

    // Legend
    const legendX = margin.left + plotW + 8;
    let legendY = margin.top + 5;
    ctx.font = '10px system-ui';
    ctx.textAlign = 'left';

    const legendItem = (color, dash, label, hasBand = false) => {
      if (hasBand) {
        ctx.fillStyle = 'rgba(59, 130, 246, 0.2)';
        ctx.fillRect(legendX, legendY - 5, 20, 10);
      }
      ctx.strokeStyle = color;
      ctx.lineWidth = 2;
      ctx.setLineDash(dash);
      ctx.beginPath();
      ctx.moveTo(legendX, legendY);
      ctx.lineTo(legendX + 20, legendY);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = '#374151';
      ctx.fillText(label, legendX + 25, legendY + 3);
      legendY += 16;
    };

    if (data.mcMean) legendItem('#3b82f6', [], 'MC Mean', true);
    if (data.numerical) legendItem('#ef4444', [], 'Numerical');
    if (data.analytical) legendItem('#6b7280', [6, 4], 'Analytical');
    
    if (data.sensorX) {
      ctx.strokeStyle = '#22c55e';
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(legendX, legendY);
      ctx.lineTo(legendX + 20, legendY);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = '#374151';
      ctx.fillText('Sensors', legendX + 25, legendY + 3);
    }
  };

  // Redraw plots when data changes
  useEffect(() => {
    const sensorX = sensors.map(s => s.x);

    const buildPlotData = (mcMeanKey, mcLowerKey, mcUpperKey, numKey, anaKey) => ({
      x: mcData?.x || numericalData?.x || analyticalData?.x,
      numX: numericalData?.x,
      mcMean: mcData?.[mcMeanKey],
      mcLower: mcData?.[mcLowerKey],
      mcUpper: mcData?.[mcUpperKey],
      numerical: numericalData?.[numKey],
      analytical: analyticalData?.[anaKey],
      sensorX,
    });

    drawPlot(canvasRefs.density, 'Density', 'ρ (kg/m³)',
      buildPlotData('mean_rho', 'ci95_lower_rho', 'ci95_upper_rho', 'density', 'density'));

    drawPlot(canvasRefs.velocity, 'Velocity', 'u (m/s)',
      buildPlotData('mean_u', 'ci95_lower_u', 'ci95_upper_u', 'velocity', 'velocity'));

    drawPlot(canvasRefs.pressure, 'Pressure', 'p (Pa)',
      buildPlotData('mean_p', 'ci95_lower_p', 'ci95_upper_p', 'pressure', 'pressure'));

    drawPlot(canvasRefs.entropy, 'Entropy', 's (J/kg·K)',
      buildPlotData('mean_entropy', 'ci95_lower_entropy', 'ci95_upper_entropy', 'entropy', 'entropy'));

  }, [mcData, numericalData, analyticalData, sensors]);

  return (
    <div style={{ display: 'flex', height: '100vh', fontFamily: 'system-ui, sans-serif', background: '#f8fafc' }}>
      
      {/* Sidebar */}
      <div style={{ width: '340px', background: '#fff', borderRight: '1px solid #e2e8f0', overflowY: 'auto', padding: '16px' }}>
        
        <h1 style={{ fontSize: '18px', fontWeight: 'bold', marginBottom: '16px', color: '#1e293b' }}>
          Shock Tube Simulator
        </h1>

        {/* Toro Test Selection */}
        <Section title="Toro Test">
          <select
            value={selectedTest}
            onChange={(e) => setSelectedTest(parseInt(e.target.value))}
            style={selectStyle}
          >
            {Object.entries(TORO_TESTS).map(([id, test]) => (
              <option key={id} value={id}>{test.name}</option>
            ))}
          </select>
          <button onClick={loadTestSensors} style={{ ...buttonStyle, marginTop: '8px', width: '100%' }}>
            Load Test Sensors
          </button>
        </Section>

        {/* Solver Config */}
        <Section title="Solver Configuration">
          <Label>Grid Cells</Label>
          <input
            type="number"
            value={config.numCells}
            onChange={(e) => setConfig({ ...config, numCells: parseInt(e.target.value) || 100 })}
            style={inputStyle}
          />

          <Label>Flux Scheme</Label>
          <select
            value={config.flux}
            onChange={(e) => setConfig({ ...config, flux: e.target.value })}
            style={selectStyle}
          >
            <option value="HLLC">HLLC</option>
            <option value="EntropyStable">Entropy Stable</option>
          </select>

          <Label>Time Integration</Label>
          <select
            value={config.timeIntegration}
            onChange={(e) => setConfig({ ...config, timeIntegration: e.target.value })}
            style={selectStyle}
          >
            <option value="RK2">RK2 (Heun)</option>
            <option value="Hancock">MUSCL-Hancock</option>
          </select>

          <Label>Interpolation</Label>
          <select
            value={config.interpolation}
            onChange={(e) => setConfig({ ...config, interpolation: e.target.value })}
            style={selectStyle}
          >
            <option value="piecewise_constant">Piecewise Constant</option>
            <option value="linear">Linear</option>
          </select>
        </Section>

        {/* Monte Carlo Settings */}
        <Section title="Monte Carlo">
          <Label>Number of Trials</Label>
          <input
            type="number"
            value={mcSettings.numTrials}
            onChange={(e) => setMcSettings({ ...mcSettings, numTrials: parseInt(e.target.value) || 10 })}
            style={inputStyle}
          />

          <Label>Noise Level (%)</Label>
          <input
            type="number"
            value={mcSettings.noiseLevel * 100}
            onChange={(e) => setMcSettings({ ...mcSettings, noiseLevel: (parseFloat(e.target.value) || 0) / 100 })}
            style={inputStyle}
            step="0.5"
          />
        </Section>

        {/* Sensor Table */}
        <Section title="Sensors">
          <table style={{ width: '100%', fontSize: '11px', borderCollapse: 'collapse' }}>
            <thead>
              <tr style={{ background: '#f1f5f9' }}>
                <th style={thStyle}>x</th>
                <th style={thStyle}>ρ</th>
                <th style={thStyle}>u</th>
                <th style={thStyle}>p</th>
                <th style={thStyle}></th>
              </tr>
            </thead>
            <tbody>
              {sensors.map((s, i) => (
                <tr key={i}>
                  <td style={tdStyle}>
                    <input type="number" value={s.x} step="0.01"
                      onChange={(e) => updateSensor(i, 'x', e.target.value)}
                      style={cellInputStyle} />
                  </td>
                  <td style={tdStyle}>
                    <input type="number" value={s.rho} step="0.1"
                      onChange={(e) => updateSensor(i, 'rho', e.target.value)}
                      style={cellInputStyle} />
                  </td>
                  <td style={tdStyle}>
                    <input type="number" value={s.u} step="0.1"
                      onChange={(e) => updateSensor(i, 'u', e.target.value)}
                      style={cellInputStyle} />
                  </td>
                  <td style={tdStyle}>
                    <input type="number" value={s.p} step="0.1"
                      onChange={(e) => updateSensor(i, 'p', e.target.value)}
                      style={cellInputStyle} />
                  </td>
                  <td style={tdStyle}>
                    <button onClick={() => removeSensor(i)} style={removeButtonStyle}>×</button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
          <button onClick={addSensor} style={{ ...buttonStyle, marginTop: '8px', width: '100%' }}>
            + Add Sensor
          </button>
        </Section>

        {/* Run Button */}
        <button
          onClick={runSimulation}
          disabled={isRunning}
          style={{
            ...buttonStyle,
            width: '100%',
            padding: '12px',
            fontSize: '14px',
            fontWeight: 'bold',
            background: isRunning ? '#94a3b8' : '#3b82f6',
            marginTop: '8px',
          }}
        >
          {isRunning ? 'Running...' : 'Run Monte Carlo'}
        </button>

        <p style={{ fontSize: '12px', color: '#64748b', marginTop: '8px', textAlign: 'center' }}>
          {statusMessage}
        </p>

        {/* Error Metrics - Numerical vs Analytical */}
        {errors && errors.success && numericalData && analyticalData && (
          <Section title="Numerical vs Analytical (% Error)">
            <table style={{ width: '100%', fontSize: '10px', borderCollapse: 'collapse' }}>
              <thead>
                <tr style={{ background: '#f1f5f9' }}>
                  <th style={thStyle}></th>
                  <th style={thStyle}>L1</th>
                  <th style={thStyle}>L2</th>
                  <th style={thStyle}>L∞</th>
                </tr>
              </thead>
              <tbody>
                <ErrorRow label="ρ" {...computeErrors(numericalData.density, analyticalData.density)} />
                <ErrorRow label="u" {...computeErrors(numericalData.velocity, analyticalData.velocity)} />
                <ErrorRow label="p" {...computeErrors(numericalData.pressure, analyticalData.pressure)} />
                <ErrorRow label="s" {...computeErrors(numericalData.entropy, analyticalData.entropy)} />
              </tbody>
            </table>
          </Section>
        )}

        {/* MC Mean vs Analytical Errors */}
        {mcData && analyticalData && mcData.mean_rho && analyticalData.density && (
          <Section title="MC Mean vs Analytical (% Error)">
            <table style={{ width: '100%', fontSize: '10px', borderCollapse: 'collapse' }}>
              <thead>
                <tr style={{ background: '#f1f5f9' }}>
                  <th style={thStyle}></th>
                  <th style={thStyle}>L1</th>
                  <th style={thStyle}>L2</th>
                  <th style={thStyle}>L∞</th>
                </tr>
              </thead>
              <tbody>
                <ErrorRow label="ρ" {...computeErrors(mcData.mean_rho, analyticalData.density)} />
                <ErrorRow label="u" {...computeErrors(mcData.mean_u, analyticalData.velocity)} />
                <ErrorRow label="p" {...computeErrors(mcData.mean_p, analyticalData.pressure)} />
                <ErrorRow label="s" {...computeErrors(mcData.mean_entropy, analyticalData.entropy)} />
              </tbody>
            </table>
          </Section>
        )}
      </div>

      {/* Main Plot Area */}
      <div style={{ flex: 1, padding: '16px', overflowY: 'auto' }}>
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '16px' }}>
          <canvas ref={canvasRefs.density} width={560} height={380} style={canvasStyle} />
          <canvas ref={canvasRefs.velocity} width={560} height={380} style={canvasStyle} />
          <canvas ref={canvasRefs.pressure} width={560} height={380} style={canvasStyle} />
          <canvas ref={canvasRefs.entropy} width={560} height={380} style={canvasStyle} />
        </div>
      </div>
    </div>
  );
}

// Compute L1, L2, Linf errors between two arrays (both absolute and relative %)
function computeErrors(computed, exact) {
  if (!computed || !exact || computed.length !== exact.length) {
    return { l1: null, l2: null, linf: null };
  }
  
  const n = computed.length;
  let sumAbs = 0, sumSq = 0, maxErr = 0;
  let sumRelAbs = 0, sumRelSq = 0, maxRelErr = 0;
  
  // Compute range for normalization
  const maxVal = Math.max(...exact.map(Math.abs));
  const scale = maxVal > 1e-10 ? maxVal : 1;
  
  for (let i = 0; i < n; i++) {
    const err = Math.abs(computed[i] - exact[i]);
    sumAbs += err;
    sumSq += err * err;
    if (err > maxErr) maxErr = err;
    
    // Relative error (normalized by range)
    const relErr = err / scale;
    sumRelAbs += relErr;
    sumRelSq += relErr * relErr;
    if (relErr > maxRelErr) maxRelErr = relErr;
  }
  
  return {
    l1: sumAbs / n,
    l2: Math.sqrt(sumSq / n),
    linf: maxErr,
    l1_pct: (sumRelAbs / n) * 100,
    l2_pct: Math.sqrt(sumRelSq / n) * 100,
    linf_pct: maxRelErr * 100
  };
}

// Helper Components
function Section({ title, children }) {
  return (
    <div style={{ marginBottom: '16px', padding: '12px', background: '#f8fafc', borderRadius: '8px', border: '1px solid #e2e8f0' }}>
      <h3 style={{ fontSize: '13px', fontWeight: '600', marginBottom: '10px', color: '#475569' }}>{title}</h3>
      {children}
    </div>
  );
}

function Label({ children }) {
  return <label style={{ display: 'block', fontSize: '11px', color: '#64748b', marginBottom: '4px', marginTop: '8px' }}>{children}</label>;
}

function ErrorRow({ label, l1, l2, linf, l1_pct, l2_pct, linf_pct }) {
  const fmt = (v) => v != null ? v.toExponential(2) : '-';
  const fmtPct = (v) => v != null ? v.toFixed(2) + '%' : '-';
  
  // If we have percentage values, show those; otherwise show absolute
  if (l1_pct != null) {
    return (
      <tr>
        <td style={{ ...tdStyle, fontWeight: '600' }}>{label}</td>
        <td style={tdStyle}>{fmtPct(l1_pct)}</td>
        <td style={tdStyle}>{fmtPct(l2_pct)}</td>
        <td style={tdStyle}>{fmtPct(linf_pct)}</td>
      </tr>
    );
  }
  
  return (
    <tr>
      <td style={{ ...tdStyle, fontWeight: '600' }}>{label}</td>
      <td style={tdStyle}>{fmt(l1)}</td>
      <td style={tdStyle}>{fmt(l2)}</td>
      <td style={tdStyle}>{fmt(linf)}</td>
    </tr>
  );
}

// Styles
const inputStyle = {
  width: '100%',
  padding: '6px 8px',
  fontSize: '12px',
  border: '1px solid #d1d5db',
  borderRadius: '4px',
  boxSizing: 'border-box',
};

const selectStyle = {
  ...inputStyle,
  background: '#fff',
};

const buttonStyle = {
  padding: '8px 12px',
  fontSize: '12px',
  background: '#3b82f6',
  color: '#fff',
  border: 'none',
  borderRadius: '6px',
  cursor: 'pointer',
};

const thStyle = {
  padding: '6px 4px',
  textAlign: 'center',
  fontWeight: '600',
  color: '#475569',
};

const tdStyle = {
  padding: '4px 2px',
  textAlign: 'center',
};

const cellInputStyle = {
  width: '100%',
  padding: '4px',
  fontSize: '11px',
  border: '1px solid #e2e8f0',
  borderRadius: '3px',
  textAlign: 'center',
  boxSizing: 'border-box',
};

const removeButtonStyle = {
  padding: '2px 6px',
  fontSize: '12px',
  background: '#fee2e2',
  color: '#dc2626',
  border: 'none',
  borderRadius: '3px',
  cursor: 'pointer',
};

const canvasStyle = {
  background: '#fff',
  borderRadius: '8px',
  boxShadow: '0 1px 3px rgba(0,0,0,0.1)',
};