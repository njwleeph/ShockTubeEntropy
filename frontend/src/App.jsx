import React, { useState, useRef, useEffect } from 'react';

const API_BASE = 'http://localhost:8080';

const TORO_TESTS = [
  { id: 1, name: 'Test 1: Sod Shock Tube', endTime: 0.25 },
  { id: 2, name: 'Test 2: 123 Problem', endTime: 0.15 },
  { id: 3, name: 'Test 3: Left Blast Wave', endTime: 0.012 },
  { id: 4, name: 'Test 4: Collision', endTime: 0.035 },
  { id: 5, name: 'Test 5: Stationary Contact', endTime: 0.035 },
];

const DEFAULT_SENSORS = [
  { x: 0.1, rho: 1.0, u: 0.0, p: 1.0 },
  { x: 0.3, rho: 1.0, u: 0.0, p: 1.0 },
  { x: 0.45, rho: 1.0, u: 0.0, p: 1.0 },
  { x: 0.55, rho: 0.125, u: 0.0, p: 0.1 },
  { x: 0.7, rho: 0.125, u: 0.0, p: 0.1 },
  { x: 0.9, rho: 0.125, u: 0.0, p: 0.1 },
];

export default function App() {
  // Configuration
  const [config, setConfig] = useState({
    toroTest: 1,
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
  const [sensors, setSensors] = useState(DEFAULT_SENSORS);

  // Results
  const [numericalData, setNumericalData] = useState(null);
  const [analyticalData, setAnalyticalData] = useState(null);
  const [monteCarloData, setMonteCarloData] = useState(null);
  const [validationResults, setValidationResults] = useState(null);
  const [simInfo, setSimInfo] = useState(null);

  // UI State
  const [isRunning, setIsRunning] = useState(false);
  const [isRunningMC, setIsRunningMC] = useState(false);
  const [error, setError] = useState(null);

  // Canvas refs
  const densityCanvasRef = useRef(null);
  const velocityCanvasRef = useRef(null);
  const pressureCanvasRef = useRef(null);
  const entropyCanvasRef = useRef(null);

  // Get end time for selected test
  const getEndTime = () => {
    const test = TORO_TESTS.find(t => t.id === config.toroTest);
    return test ? test.endTime : 0.25;
  };

  // Load sensor data from Toro test
  const loadToroSensorData = async () => {
    try {
      const positions = sensors.map(s => s.x);
      const response = await fetch(`${API_BASE}/api/toro/sensor-data`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ test: config.toroTest, positions }),
      });
      const data = await response.json();
      if (data.success) {
        setSensors(data.sensors);
      }
    } catch (err) {
      setError('Failed to load Toro sensor data: ' + err.message);
    }
  };

  // Run simulation
  const runSimulation = async () => {
    setIsRunning(true);
    setError(null);

    try {
      // 1. Create solver
      await fetch(`${API_BASE}/api/simulation/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          length: 1.0,
          numCells: config.numCells,
          gamma: config.gamma,
          CFL: config.CFL,
          endTime: getEndTime(),
        }),
      });

      // 2. Configure flux and time integration
      await fetch(`${API_BASE}/api/simulation/configure`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          flux: config.flux,
          timeIntegration: config.timeIntegration,
        }),
      });

      // 3. Initialize from sparse sensor data
      await fetch(`${API_BASE}/api/simulation/sparse-init`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sensors: sensors,
          interpolation: config.interpolation,
        }),
      });

      // 4. Run simulation
      const runResponse = await fetch(`${API_BASE}/api/simulation/run`, {
        method: 'POST',
      });
      const runData = await runResponse.json();
      setSimInfo({ time: runData.time, steps: runData.steps });

      // 5. Get numerical data
      const numResponse = await fetch(`${API_BASE}/api/simulation/data`);
      const numData = await numResponse.json();
      setNumericalData(numData);

      // 6. Get analytical data
      const anaResponse = await fetch(`${API_BASE}/api/simulation/analytical`);
      const anaData = await anaResponse.json();
      setAnalyticalData(anaData);

      // 7. Get validation metrics
      const valResponse = await fetch(`${API_BASE}/api/simulation/validate`);
      const valData = await valResponse.json();
      setValidationResults(valData);

      // Clear MC data when running regular simulation
      setMonteCarloData(null);

    } catch (err) {
      setError('Simulation failed: ' + err.message);
    } finally {
      setIsRunning(false);
    }
  };

  // Run Monte Carlo
  const runMonteCarlo = async () => {
    setIsRunningMC(true);
    setError(null);

    try {
      // 1. Create solver
      await fetch(`${API_BASE}/api/simulation/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          length: 1.0,
          numCells: config.numCells,
          gamma: config.gamma,
          CFL: config.CFL,
          endTime: getEndTime(),
        }),
      });

      // 2. Configure
      await fetch(`${API_BASE}/api/simulation/configure`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          flux: config.flux,
          timeIntegration: config.timeIntegration,
        }),
      });

      // 3. Initialize (needed to set initial conditions for analytical)
      await fetch(`${API_BASE}/api/simulation/sparse-init`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sensors: sensors,
          interpolation: config.interpolation,
        }),
      });

      // 4. Run Monte Carlo
      const mcResponse = await fetch(`${API_BASE}/api/simulation/montecarlo`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sensors: sensors,
          noiseLevel: mcSettings.noiseLevel,
          numTrials: mcSettings.numTrials,
          interpolation: config.interpolation,
          timeIntegration: config.timeIntegration,
        }),
      });
      const mcData = await mcResponse.json();
      setMonteCarloData(mcData);

      // Also run single simulation for numerical overlay
      await fetch(`${API_BASE}/api/simulation/sparse-init`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sensors: sensors,
          interpolation: config.interpolation,
        }),
      });

      const runResponse = await fetch(`${API_BASE}/api/simulation/run`, {
        method: 'POST',
      });
      const runData = await runResponse.json();
      setSimInfo({ time: runData.time, steps: runData.steps, mcTime: mcData.computation_time_ms });

      const numResponse = await fetch(`${API_BASE}/api/simulation/data`);
      const numData = await numResponse.json();
      setNumericalData(numData);

      // Analytical comes from MC response
      setAnalyticalData({
        x: mcData.x,
        density: mcData.analytical_rho,
        velocity: mcData.analytical_u,
        pressure: mcData.analytical_p,
        entropy: mcData.analytical_entropy,
      });

      // Get validation
      const valResponse = await fetch(`${API_BASE}/api/simulation/validate`);
      const valData = await valResponse.json();
      setValidationResults(valData);

    } catch (err) {
      setError('Monte Carlo failed: ' + err.message);
    } finally {
      setIsRunningMC(false);
    }
  };

  // Sensor table handlers
  const updateSensor = (index, field, value) => {
    const newSensors = [...sensors];
    newSensors[index] = { ...newSensors[index], [field]: parseFloat(value) || 0 };
    setSensors(newSensors);
  };

  const addSensor = () => {
    const lastX = sensors.length > 0 ? sensors[sensors.length - 1].x : 0;
    setSensors([...sensors, { x: Math.min(lastX + 0.1, 0.99), rho: 1.0, u: 0.0, p: 1.0 }]);
  };

  const removeSensor = (index) => {
    if (sensors.length > 2) {
      setSensors(sensors.filter((_, i) => i !== index));
    }
  };

  // Drawing function
  const drawPlot = (canvasRef, title, yLabel, datasets, sensorPositions) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;

    // Layout
    const padding = { left: 70, right: 120, top: 40, bottom: 50 };
    const plotWidth = width - padding.left - padding.right;
    const plotHeight = height - padding.top - padding.bottom;

    // Clear
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, width, height);

    // Find data range
    let xMin = 0, xMax = 1;
    let yMin = Infinity, yMax = -Infinity;

    datasets.forEach(ds => {
      if (ds.y) {
        ds.y.forEach(v => {
          if (isFinite(v)) {
            yMin = Math.min(yMin, v);
            yMax = Math.max(yMax, v);
          }
        });
      }
      if (ds.yLower && ds.yUpper) {
        ds.yLower.forEach(v => { if (isFinite(v)) yMin = Math.min(yMin, v); });
        ds.yUpper.forEach(v => { if (isFinite(v)) yMax = Math.max(yMax, v); });
      }
    });

    // Handle empty or invalid range
    if (!isFinite(yMin) || !isFinite(yMax)) {
      yMin = 0;
      yMax = 1;
    }

    // Add margin
    const yRange = yMax - yMin || 1;
    yMin -= yRange * 0.1;
    yMax += yRange * 0.1;

    // Transform functions
    const toCanvasX = (x) => padding.left + (x - xMin) / (xMax - xMin) * plotWidth;
    const toCanvasY = (y) => padding.top + (1 - (y - yMin) / (yMax - yMin)) * plotHeight;

    // Draw grid
    ctx.strokeStyle = '#e5e7eb';
    ctx.lineWidth = 1;
    
    // Vertical grid lines
    for (let i = 0; i <= 5; i++) {
      const x = toCanvasX(xMin + (xMax - xMin) * i / 5);
      ctx.beginPath();
      ctx.moveTo(x, padding.top);
      ctx.lineTo(x, padding.top + plotHeight);
      ctx.stroke();
    }

    // Horizontal grid lines
    for (let i = 0; i <= 5; i++) {
      const y = toCanvasY(yMin + (yMax - yMin) * i / 5);
      ctx.beginPath();
      ctx.moveTo(padding.left, y);
      ctx.lineTo(padding.left + plotWidth, y);
      ctx.stroke();
    }

    // Draw sensor positions
    if (sensorPositions && sensorPositions.length > 0) {
      ctx.strokeStyle = '#22c55e';
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 4]);
      sensorPositions.forEach(sx => {
        const x = toCanvasX(sx);
        ctx.beginPath();
        ctx.moveTo(x, padding.top);
        ctx.lineTo(x, padding.top + plotHeight);
        ctx.stroke();
      });
      ctx.setLineDash([]);
    }

    // Draw datasets
    datasets.forEach(ds => {
      if (!ds.x || ds.x.length === 0) return;

      // Draw CI band if present
      if (ds.yLower && ds.yUpper) {
        ctx.fillStyle = ds.bandColor || 'rgba(59, 130, 246, 0.2)';
        ctx.beginPath();
        ctx.moveTo(toCanvasX(ds.x[0]), toCanvasY(ds.yLower[0]));
        for (let i = 1; i < ds.x.length; i++) {
          ctx.lineTo(toCanvasX(ds.x[i]), toCanvasY(ds.yLower[i]));
        }
        for (let i = ds.x.length - 1; i >= 0; i--) {
          ctx.lineTo(toCanvasX(ds.x[i]), toCanvasY(ds.yUpper[i]));
        }
        ctx.closePath();
        ctx.fill();
      }

      // Draw line
      if (ds.y) {
        ctx.strokeStyle = ds.color || '#000';
        ctx.lineWidth = ds.lineWidth || 2;
        ctx.setLineDash(ds.dash || []);
        ctx.beginPath();
        let started = false;
        for (let i = 0; i < ds.x.length; i++) {
          if (isFinite(ds.y[i])) {
            const cx = toCanvasX(ds.x[i]);
            const cy = toCanvasY(ds.y[i]);
            if (!started) {
              ctx.moveTo(cx, cy);
              started = true;
            } else {
              ctx.lineTo(cx, cy);
            }
          }
        }
        ctx.stroke();
        ctx.setLineDash([]);
      }
    });

    // Draw axes
    ctx.strokeStyle = '#374151';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(padding.left, padding.top);
    ctx.lineTo(padding.left, padding.top + plotHeight);
    ctx.lineTo(padding.left + plotWidth, padding.top + plotHeight);
    ctx.stroke();

    // X-axis labels
    ctx.fillStyle = '#374151';
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'center';
    for (let i = 0; i <= 5; i++) {
      const val = xMin + (xMax - xMin) * i / 5;
      const x = toCanvasX(val);
      ctx.fillText(val.toFixed(2), x, padding.top + plotHeight + 15);
    }
    ctx.font = '12px sans-serif';
    ctx.fillText('Position (m)', padding.left + plotWidth / 2, height - 8);

    // Y-axis labels
    ctx.textAlign = 'right';
    ctx.font = '11px sans-serif';
    for (let i = 0; i <= 5; i++) {
      const val = yMin + (yMax - yMin) * i / 5;
      const y = toCanvasY(val);
      ctx.fillText(val.toPrecision(3), padding.left - 8, y + 4);
    }

    // Y-axis label (rotated)
    ctx.save();
    ctx.translate(15, padding.top + plotHeight / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.font = '12px sans-serif';
    ctx.fillText(yLabel, 0, 0);
    ctx.restore();

    // Title
    ctx.fillStyle = '#111827';
    ctx.font = 'bold 14px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText(title, padding.left + plotWidth / 2, 20);

    // Legend
    const legendX = padding.left + plotWidth + 10;
    let legendY = padding.top + 10;
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'left';

    datasets.forEach(ds => {
      if (ds.label) {
        // Draw band sample if present
        if (ds.yLower && ds.yUpper) {
          ctx.fillStyle = ds.bandColor || 'rgba(59, 130, 246, 0.2)';
          ctx.fillRect(legendX, legendY - 5, 20, 10);
        }

        // Draw line sample
        ctx.strokeStyle = ds.color || '#000';
        ctx.lineWidth = 2;
        ctx.setLineDash(ds.dash || []);
        ctx.beginPath();
        ctx.moveTo(legendX, legendY);
        ctx.lineTo(legendX + 20, legendY);
        ctx.stroke();
        ctx.setLineDash([]);

        // Label
        ctx.fillStyle = '#374151';
        ctx.fillText(ds.label, legendX + 25, legendY + 4);
        legendY += 18;
      }
    });

    // Sensor legend
    if (sensorPositions && sensorPositions.length > 0) {
      ctx.strokeStyle = '#22c55e';
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(legendX, legendY);
      ctx.lineTo(legendX + 20, legendY);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = '#374151';
      ctx.fillText('Sensors', legendX + 25, legendY + 4);
    }
  };

  // Effect to redraw plots
  useEffect(() => {
    const sensorX = sensors.map(s => s.x);

    // Build datasets for each plot
    const buildDatasets = (numKey, anaKey, mcMeanKey, mcLowerKey, mcUpperKey) => {
      const datasets = [];

      // Monte Carlo CI band (draw first so it's behind)
      if (monteCarloData && monteCarloData[mcMeanKey]) {
        datasets.push({
          x: monteCarloData.x,
          y: monteCarloData[mcMeanKey],
          yLower: monteCarloData[mcLowerKey],
          yUpper: monteCarloData[mcUpperKey],
          color: '#3b82f6',
          bandColor: 'rgba(59, 130, 246, 0.2)',
          label: 'MC Mean ± 95% CI',
          lineWidth: 2,
        });
      }

      // Numerical solution
      if (numericalData && numericalData[numKey]) {
        datasets.push({
          x: numericalData.x,
          y: numericalData[numKey],
          color: '#ef4444',
          label: 'Numerical',
          lineWidth: 2,
        });
      }

      // Analytical solution
      if (analyticalData && analyticalData[anaKey]) {
        datasets.push({
          x: analyticalData.x,
          y: analyticalData[anaKey],
          color: '#6b7280',
          dash: [6, 4],
          label: 'Analytical',
          lineWidth: 2,
        });
      }

      return datasets;
    };

    drawPlot(
      densityCanvasRef,
      'Density',
      'ρ (kg/m³)',
      buildDatasets('density', 'density', 'mean_rho', 'ci95_lower_rho', 'ci95_upper_rho'),
      sensorX
    );

    drawPlot(
      velocityCanvasRef,
      'Velocity',
      'u (m/s)',
      buildDatasets('velocity', 'velocity', 'mean_u', 'ci95_lower_u', 'ci95_upper_u'),
      sensorX
    );

    drawPlot(
      pressureCanvasRef,
      'Pressure',
      'p (Pa)',
      buildDatasets('pressure', 'pressure', 'mean_p', 'ci95_lower_p', 'ci95_upper_p'),
      sensorX
    );

    drawPlot(
      entropyCanvasRef,
      'Entropy',
      's (J/kg·K)',
      buildDatasets('entropy', 'entropy', 'mean_entropy', 'ci95_lower_entropy', 'ci95_upper_entropy'),
      sensorX
    );
  }, [numericalData, analyticalData, monteCarloData, sensors]);

  // Styles
  const styles = {
    container: {
      display: 'flex',
      height: '100vh',
      fontFamily: 'system-ui, -apple-system, sans-serif',
      backgroundColor: '#f8fafc',
    },
    sidebar: {
      width: '320px',
      backgroundColor: '#ffffff',
      borderRight: '1px solid #e2e8f0',
      overflowY: 'auto',
      padding: '16px',
    },
    main: {
      flex: 1,
      padding: '16px',
      overflowY: 'auto',
    },
    section: {
      marginBottom: '20px',
      padding: '12px',
      backgroundColor: '#f8fafc',
      borderRadius: '8px',
      border: '1px solid #e2e8f0',
    },
    sectionTitle: {
      fontSize: '13px',
      fontWeight: '600',
      color: '#1e293b',
      marginBottom: '12px',
      paddingBottom: '8px',
      borderBottom: '1px solid #e2e8f0',
    },
    label: {
      display: 'block',
      fontSize: '11px',
      fontWeight: '500',
      color: '#64748b',
      marginBottom: '4px',
    },
    select: {
      width: '100%',
      padding: '8px',
      fontSize: '12px',
      border: '1px solid #d1d5db',
      borderRadius: '6px',
      marginBottom: '10px',
      backgroundColor: '#fff',
    },
    input: {
      width: '100%',
      padding: '8px',
      fontSize: '12px',
      border: '1px solid #d1d5db',
      borderRadius: '6px',
      marginBottom: '10px',
      boxSizing: 'border-box',
    },
    button: {
      width: '100%',
      padding: '10px',
      fontSize: '13px',
      fontWeight: '600',
      color: '#fff',
      backgroundColor: '#3b82f6',
      border: 'none',
      borderRadius: '6px',
      cursor: 'pointer',
      marginBottom: '8px',
    },
    buttonSecondary: {
      width: '100%',
      padding: '10px',
      fontSize: '13px',
      fontWeight: '600',
      color: '#374151',
      backgroundColor: '#e5e7eb',
      border: 'none',
      borderRadius: '6px',
      cursor: 'pointer',
      marginBottom: '8px',
    },
    buttonMC: {
      width: '100%',
      padding: '10px',
      fontSize: '13px',
      fontWeight: '600',
      color: '#fff',
      backgroundColor: '#8b5cf6',
      border: 'none',
      borderRadius: '6px',
      cursor: 'pointer',
      marginBottom: '8px',
    },
    plotGrid: {
      display: 'grid',
      gridTemplateColumns: 'repeat(2, 1fr)',
      gap: '16px',
    },
    plotContainer: {
      backgroundColor: '#fff',
      borderRadius: '8px',
      padding: '12px',
      border: '1px solid #e2e8f0',
    },
    canvas: {
      width: '100%',
      height: 'auto',
    },
    table: {
      width: '100%',
      fontSize: '11px',
      borderCollapse: 'collapse',
    },
    th: {
      padding: '6px 4px',
      textAlign: 'left',
      borderBottom: '1px solid #e2e8f0',
      color: '#64748b',
      fontWeight: '600',
    },
    td: {
      padding: '4px',
      borderBottom: '1px solid #f1f5f9',
    },
    tdInput: {
      width: '100%',
      padding: '4px',
      fontSize: '11px',
      border: '1px solid #e2e8f0',
      borderRadius: '4px',
      boxSizing: 'border-box',
    },
    removeBtn: {
      padding: '2px 6px',
      fontSize: '10px',
      color: '#ef4444',
      backgroundColor: '#fef2f2',
      border: '1px solid #fecaca',
      borderRadius: '4px',
      cursor: 'pointer',
    },
    addBtn: {
      padding: '6px 12px',
      fontSize: '11px',
      color: '#3b82f6',
      backgroundColor: '#eff6ff',
      border: '1px solid #bfdbfe',
      borderRadius: '4px',
      cursor: 'pointer',
      marginTop: '8px',
    },
    error: {
      padding: '10px',
      backgroundColor: '#fef2f2',
      color: '#dc2626',
      borderRadius: '6px',
      fontSize: '12px',
      marginBottom: '12px',
    },
    infoRow: {
      display: 'flex',
      justifyContent: 'space-between',
      fontSize: '11px',
      padding: '4px 0',
      borderBottom: '1px solid #f1f5f9',
    },
    infoLabel: {
      color: '#64748b',
    },
    infoValue: {
      color: '#1e293b',
      fontWeight: '500',
    },
    errorSection: {
      marginBottom: '8px',
    },
    errorTitle: {
      fontSize: '10px',
      fontWeight: '600',
      color: '#475569',
      marginBottom: '4px',
    },
  };

  return (
    <div style={styles.container}>
      {/* Sidebar */}
      <div style={styles.sidebar}>
        <h2 style={{ fontSize: '16px', fontWeight: '700', marginBottom: '16px', color: '#1e293b' }}>
          Shock Tube Simulator
        </h2>

        {error && <div style={styles.error}>{error}</div>}

        {/* Configuration */}
        <div style={styles.section}>
          <div style={styles.sectionTitle}>Configuration</div>

          <label style={styles.label}>Toro Test Case</label>
          <select
            style={styles.select}
            value={config.toroTest}
            onChange={(e) => setConfig({ ...config, toroTest: parseInt(e.target.value) })}
          >
            {TORO_TESTS.map(t => (
              <option key={t.id} value={t.id}>{t.name}</option>
            ))}
          </select>

          <label style={styles.label}>Grid Cells</label>
          <input
            type="number"
            style={styles.input}
            value={config.numCells}
            onChange={(e) => setConfig({ ...config, numCells: parseInt(e.target.value) || 500 })}
          />

          <label style={styles.label}>Flux Scheme</label>
          <select
            style={styles.select}
            value={config.flux}
            onChange={(e) => setConfig({ ...config, flux: e.target.value })}
          >
            <option value="HLLC">HLLC</option>
            <option value="EntropyStable">Entropy Stable</option>
          </select>

          <label style={styles.label}>Time Integration</label>
          <select
            style={styles.select}
            value={config.timeIntegration}
            onChange={(e) => setConfig({ ...config, timeIntegration: e.target.value })}
          >
            <option value="RK2">RK2</option>
            <option value="Hancock">MUSCL-Hancock</option>
          </select>

          <label style={styles.label}>Interpolation</label>
          <select
            style={styles.select}
            value={config.interpolation}
            onChange={(e) => setConfig({ ...config, interpolation: e.target.value })}
          >
            <option value="piecewise_constant">Piecewise Constant</option>
            <option value="linear">Linear</option>
          </select>
        </div>

        {/* Monte Carlo Settings */}
        <div style={styles.section}>
          <div style={styles.sectionTitle}>Monte Carlo Settings</div>

          <label style={styles.label}>Number of Trials</label>
          <input
            type="number"
            style={styles.input}
            value={mcSettings.numTrials}
            onChange={(e) => setMcSettings({ ...mcSettings, numTrials: parseInt(e.target.value) || 100 })}
          />

          <label style={styles.label}>Noise Level (%)</label>
          <input
            type="number"
            step="0.01"
            style={styles.input}
            value={mcSettings.noiseLevel * 100}
            onChange={(e) => setMcSettings({ ...mcSettings, noiseLevel: (parseFloat(e.target.value) || 5) / 100 })}
          />
        </div>

        {/* Sensors */}
        <div style={styles.section}>
          <div style={styles.sectionTitle}>Sensors ({sensors.length})</div>

          <button style={styles.buttonSecondary} onClick={loadToroSensorData}>
            Load Toro Test Data
          </button>

          <div style={{ maxHeight: '200px', overflowY: 'auto', marginTop: '8px' }}>
            <table style={styles.table}>
              <thead>
                <tr>
                  <th style={styles.th}>x</th>
                  <th style={styles.th}>ρ</th>
                  <th style={styles.th}>u</th>
                  <th style={styles.th}>p</th>
                  <th style={styles.th}></th>
                </tr>
              </thead>
              <tbody>
                {sensors.map((s, i) => (
                  <tr key={i}>
                    <td style={styles.td}>
                      <input
                        type="number"
                        step="0.01"
                        style={styles.tdInput}
                        value={s.x}
                        onChange={(e) => updateSensor(i, 'x', e.target.value)}
                      />
                    </td>
                    <td style={styles.td}>
                      <input
                        type="number"
                        step="0.01"
                        style={styles.tdInput}
                        value={s.rho}
                        onChange={(e) => updateSensor(i, 'rho', e.target.value)}
                      />
                    </td>
                    <td style={styles.td}>
                      <input
                        type="number"
                        step="0.1"
                        style={styles.tdInput}
                        value={s.u}
                        onChange={(e) => updateSensor(i, 'u', e.target.value)}
                      />
                    </td>
                    <td style={styles.td}>
                      <input
                        type="number"
                        step="0.1"
                        style={styles.tdInput}
                        value={s.p}
                        onChange={(e) => updateSensor(i, 'p', e.target.value)}
                      />
                    </td>
                    <td style={styles.td}>
                      <button style={styles.removeBtn} onClick={() => removeSensor(i)}>×</button>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <button style={styles.addBtn} onClick={addSensor}>+ Add Sensor</button>
        </div>

        {/* Run Buttons */}
        <div style={styles.section}>
          <button
            style={{ ...styles.button, opacity: isRunning ? 0.6 : 1 }}
            onClick={runSimulation}
            disabled={isRunning || isRunningMC}
          >
            {isRunning ? 'Running...' : 'Run Simulation'}
          </button>

          <button
            style={{ ...styles.buttonMC, opacity: isRunningMC ? 0.6 : 1 }}
            onClick={runMonteCarlo}
            disabled={isRunning || isRunningMC}
          >
            {isRunningMC ? 'Running MC...' : 'Run Monte Carlo'}
          </button>
        </div>

        {/* Simulation Info */}
        {simInfo && (
          <div style={styles.section}>
            <div style={styles.sectionTitle}>Simulation Info</div>
            <div style={styles.infoRow}>
              <span style={styles.infoLabel}>Final Time</span>
              <span style={styles.infoValue}>{simInfo.time?.toFixed(6)} s</span>
            </div>
            <div style={styles.infoRow}>
              <span style={styles.infoLabel}>Time Steps</span>
              <span style={styles.infoValue}>{simInfo.steps}</span>
            </div>
            {simInfo.mcTime && (
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>MC Compute Time</span>
                <span style={styles.infoValue}>{simInfo.mcTime.toFixed(1)} ms</span>
              </div>
            )}
          </div>
        )}

        {/* Monte Carlo Stats */}
        {monteCarloData && (
          <div style={styles.section}>
            <div style={styles.sectionTitle}>Monte Carlo Statistics</div>
            <div style={styles.infoRow}>
              <span style={styles.infoLabel}>Trials</span>
              <span style={styles.infoValue}>{monteCarloData.num_trials}</span>
            </div>
            <div style={styles.infoRow}>
              <span style={styles.infoLabel}>Noise Level</span>
              <span style={styles.infoValue}>{(monteCarloData.noise_level * 100).toFixed(1)}%</span>
            </div>
            <div style={styles.infoRow}>
              <span style={styles.infoLabel}>Compute Time</span>
              <span style={styles.infoValue}>{monteCarloData.computation_time_ms?.toFixed(1)} ms</span>
            </div>
          </div>
        )}

        {/* Validation Results */}
        {validationResults && (
          <div style={styles.section}>
            <div style={styles.sectionTitle}>Validation Errors</div>

            <div style={styles.errorSection}>
              <div style={styles.errorTitle}>L1 Error (Mean Absolute)</div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Density</span>
                <span style={styles.infoValue}>{validationResults.L1_error_density?.toExponential(3)}</span>
              </div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Velocity</span>
                <span style={styles.infoValue}>{validationResults.L1_error_velocity?.toExponential(3)}</span>
              </div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Pressure</span>
                <span style={styles.infoValue}>{validationResults.L1_error_pressure?.toExponential(3)}</span>
              </div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Entropy</span>
                <span style={styles.infoValue}>{validationResults.L1_error_entropy?.toExponential(3)}</span>
              </div>
            </div>

            <div style={styles.errorSection}>
              <div style={styles.errorTitle}>L2 Error (RMS)</div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Density</span>
                <span style={styles.infoValue}>{validationResults.L2_error_density?.toExponential(3)}</span>
              </div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Velocity</span>
                <span style={styles.infoValue}>{validationResults.L2_error_velocity?.toExponential(3)}</span>
              </div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Pressure</span>
                <span style={styles.infoValue}>{validationResults.L2_error_pressure?.toExponential(3)}</span>
              </div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Entropy</span>
                <span style={styles.infoValue}>{validationResults.L2_error_entropy?.toExponential(3)}</span>
              </div>
            </div>

            <div style={styles.errorSection}>
              <div style={styles.errorTitle}>L∞ Error (Maximum)</div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Density</span>
                <span style={styles.infoValue}>{validationResults.Linf_error_density?.toExponential(3)}</span>
              </div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Velocity</span>
                <span style={styles.infoValue}>{validationResults.Linf_error_velocity?.toExponential(3)}</span>
              </div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Pressure</span>
                <span style={styles.infoValue}>{validationResults.Linf_error_pressure?.toExponential(3)}</span>
              </div>
              <div style={styles.infoRow}>
                <span style={styles.infoLabel}>Entropy</span>
                <span style={styles.infoValue}>{validationResults.Linf_error_entropy?.toExponential(3)}</span>
              </div>
            </div>
          </div>
        )}
      </div>

      {/* Main Plot Area */}
      <div style={styles.main}>
        <div style={styles.plotGrid}>
          <div style={styles.plotContainer}>
            <canvas ref={densityCanvasRef} width={550} height={350} style={styles.canvas} />
          </div>
          <div style={styles.plotContainer}>
            <canvas ref={velocityCanvasRef} width={550} height={350} style={styles.canvas} />
          </div>
          <div style={styles.plotContainer}>
            <canvas ref={pressureCanvasRef} width={550} height={350} style={styles.canvas} />
          </div>
          <div style={styles.plotContainer}>
            <canvas ref={entropyCanvasRef} width={550} height={350} style={styles.canvas} />
          </div>
        </div>
      </div>
    </div>
  );
}