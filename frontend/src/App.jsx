import React, { useState, useEffect, useRef } from 'react';

const ShockTubeSimulator = () => {
  const API_URL = 'http://localhost:8080';

  // Simulation config
  const [config, setConfig] = useState({
    numCells: 500,
    gamma: 1.4,
    CFL: 0.5,
    endTime: 0.25,
    fluxType: 'HLLC',
    timeIntegration: 'RK2',
    interpolation: 'linear'
  });

  // Toro test selection
  const [selectedTest, setSelectedTest] = useState(1);

  // Toro test end times
  const TORO_TEST_END_TIMES = {
    1: 0.25,
    2: 0.15,
    3: 0.012,
    4: 0.035,
    5: 0.035
  };

  // Sensor data - editable
  const [sensors, setSensors] = useState([
    { x: 0.10, rho: 1.0, u: 0.0, p: 1.0 },
    { x: 0.25, rho: 1.0, u: 0.0, p: 1.0 },
    { x: 0.40, rho: 1.0, u: 0.0, p: 1.0 },
    { x: 0.60, rho: 0.125, u: 0.0, p: 0.1 },
    { x: 0.75, rho: 0.125, u: 0.0, p: 0.1 },
    { x: 0.90, rho: 0.125, u: 0.0, p: 0.1 }
  ]);

  // Results
  const [numericalData, setNumericalData] = useState(null);
  const [analyticalData, setAnalyticalData] = useState(null);
  const [validationResults, setValidationResults] = useState(null);
  const [status, setStatus] = useState({ ready: false, running: false });

  // Monte Carlo config and results
  const [monteCarloConfig, setMonteCarloConfig] = useState({
    noiseLevel: 0.05,
    numTrials: 100
  });
  const [monteCarloResults, setMonteCarloResults] = useState(null);

  // Canvas refs
  const densityCanvasRef = useRef(null);
  const velocityCanvasRef = useRef(null);
  const pressureCanvasRef = useRef(null);
  const entropyCanvasRef = useRef(null);

  // Add sensor
  const addSensor = () => {
    if (sensors.length >= 10) return;
    const lastX = sensors.length > 0 ? sensors[sensors.length - 1].x : 0.5;
    const newX = Math.min(0.95, lastX + 0.1);
    setSensors([...sensors, { x: newX, rho: 0.5, u: 0.0, p: 0.5 }]);
  };

  // Remove sensor
  const removeSensor = (index) => {
    if (sensors.length <= 2) return;
    setSensors(sensors.filter((_, i) => i !== index));
  };

  // Update sensor value
  const updateSensor = (index, field, value) => {
    const updated = [...sensors];
    updated[index] = { ...updated[index], [field]: parseFloat(value) || 0 };
    setSensors(updated);
  };

  // Load sensor data from Toro test
  const loadToroTest = async () => {
    try {
      const positions = sensors.map(s => s.x);
      const response = await fetch(`${API_URL}/api/toro/sensor-data`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ test: selectedTest, positions })
      });
      const data = await response.json();

      if (data.success) {
        setSensors(data.sensors);
        setConfig(prev => ({
          ...prev,
          endTime: TORO_TEST_END_TIMES[selectedTest] || 0.25
        }));
      }
    } catch (error) {
      console.error('Error loading Toro test:', error);
    }
  };

  // Run standard Toro test simulation (exact initial conditions)
  const runSimulation = async () => {
    try {
      setStatus({ ready: false, running: true });
      setValidationResults(null);
      setMonteCarloResults(null);

      await fetch(`${API_URL}/api/simulation/reset`, { method: 'POST' });

      await fetch(`${API_URL}/api/simulation/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          length: 1.0,
          numCells: config.numCells,
          gamma: config.gamma,
          CFL: config.CFL,
          endTime: config.endTime
        })
      });

      await fetch(`${API_URL}/api/simulation/configure`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          flux: config.fluxType,
          timeIntegration: config.timeIntegration
        })
      });

      await fetch(`${API_URL}/api/simulation/initialize/toro`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ test: selectedTest })
      });

      const runResponse = await fetch(`${API_URL}/api/simulation/run`, { method: 'POST' });
      const runData = await runResponse.json();

      if (runData.success) {
        const dataResponse = await fetch(`${API_URL}/api/simulation/data`);
        const simData = await dataResponse.json();
        if (simData.success) {
          setNumericalData({
            x: simData.x,
            density: simData.density,
            velocity: simData.velocity,
            pressure: simData.pressure,
            entropy: simData.entropy
          });
        }

        const anaResponse = await fetch(`${API_URL}/api/simulation/analytical`);
        const anaData = await anaResponse.json();
        if (anaData.success) {
          setAnalyticalData({
            x: anaData.x,
            density: anaData.density,
            velocity: anaData.velocity,
            pressure: anaData.pressure
          });
        }

        const valResponse = await fetch(`${API_URL}/api/simulation/validate`);
        const valData = await valResponse.json();
        if (valData.success) {
          setValidationResults(valData);
        }
      }

      setStatus({ ready: true, running: false });

    } catch (error) {
      console.error('Error running simulation:', error);
      setStatus({ ready: false, running: false });
    }
  };

  // Run simulation from sparse sensor data
  const runFromSensors = async () => {
    try {
      setStatus({ ready: false, running: true });
      setValidationResults(null);
      setMonteCarloResults(null);

      await fetch(`${API_URL}/api/simulation/reset`, { method: 'POST' });

      await fetch(`${API_URL}/api/simulation/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          length: 1.0,
          numCells: config.numCells,
          gamma: config.gamma,
          CFL: config.CFL,
          endTime: config.endTime
        })
      });

      await fetch(`${API_URL}/api/simulation/configure`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          flux: config.fluxType,
          timeIntegration: config.timeIntegration
        })
      });

      // Initialize Toro test for analytical reference
      await fetch(`${API_URL}/api/simulation/initialize/toro`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ test: selectedTest })
      });

      // Initialize from sparse sensor data
      await fetch(`${API_URL}/api/simulation/sparse-init`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sensors: sensors,
          interpolation: config.interpolation
        })
      });

      const runResponse = await fetch(`${API_URL}/api/simulation/run`, { method: 'POST' });
      const runData = await runResponse.json();

      if (runData.success) {
        const dataResponse = await fetch(`${API_URL}/api/simulation/data`);
        const simData = await dataResponse.json();
        if (simData.success) {
          setNumericalData({
            x: simData.x,
            density: simData.density,
            velocity: simData.velocity,
            pressure: simData.pressure,
            entropy: simData.entropy
          });
        }

        const anaResponse = await fetch(`${API_URL}/api/simulation/analytical`);
        const anaData = await anaResponse.json();
        if (anaData.success) {
          setAnalyticalData({
            x: anaData.x,
            density: anaData.density,
            velocity: anaData.velocity,
            pressure: anaData.pressure
          });
        }

        const valResponse = await fetch(`${API_URL}/api/simulation/validate`);
        const valData = await valResponse.json();
        if (valData.success) {
          setValidationResults(valData);
        }
      }

      setStatus({ ready: true, running: false });

    } catch (error) {
      console.error('Error running from sensors:', error);
      setStatus({ ready: false, running: false });
    }
  };

  // Run Monte Carlo analysis
  const runMonteCarlo = async () => {
    try {
      setStatus({ ready: false, running: true });
      setNumericalData(null);
      setAnalyticalData(null);
      setValidationResults(null);

      await fetch(`${API_URL}/api/simulation/reset`, { method: 'POST' });

      await fetch(`${API_URL}/api/simulation/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          length: 1.0,
          numCells: config.numCells,
          gamma: config.gamma,
          CFL: config.CFL,
          endTime: config.endTime
        })
      });

      await fetch(`${API_URL}/api/simulation/configure`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          flux: config.fluxType,
          timeIntegration: config.timeIntegration
        })
      });

      await fetch(`${API_URL}/api/simulation/initialize/toro`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ test: selectedTest })
      });

      const mcResponse = await fetch(`${API_URL}/api/simulation/montecarlo`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sensors: sensors,
          noiseLevel: monteCarloConfig.noiseLevel,
          numTrials: monteCarloConfig.numTrials,
          interpolation: config.interpolation,
          includeAnalytical: true
        })
      });

      const mcData = await mcResponse.json();

      if (mcData.success) {
        setMonteCarloResults(mcData);
      } else {
        console.error('Monte Carlo failed:', mcData.error);
      }

      setStatus({ ready: true, running: false });

    } catch (error) {
      console.error('Error running Monte Carlo:', error);
      setStatus({ ready: false, running: false });
    }
  };

  // Draw comparison plot (for regular simulation)
  const drawComparisonPlot = (canvasRef, xNum, yNum, xAna, yAna, ylabel, colorNum, colorAna) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    canvas.width = 500;
    canvas.height = 300;

    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    const padding = { left: 60, right: 100, top: 30, bottom: 45 };
    const plotWidth = canvas.width - padding.left - padding.right;
    const plotHeight = canvas.height - padding.top - padding.bottom;

    const allY = [...(yNum || []), ...(yAna || [])];
    const allX = [...(xNum || []), ...(xAna || [])];

    if (allY.length === 0) return;

    const xMin = Math.min(...allX);
    const xMax = Math.max(...allX);
    const yMin = Math.min(...allY);
    const yMax = Math.max(...allY);
    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;
    const yPadding = yRange * 0.1;

    const xScale = (x) => padding.left + ((x - xMin) / xRange) * plotWidth;
    const yScale = (y) => canvas.height - padding.bottom - ((y - (yMin - yPadding)) / (yRange + 2 * yPadding)) * plotHeight;

    // Grid
    ctx.strokeStyle = '#e0e0e0';
    ctx.lineWidth = 1;
    for (let i = 0; i <= 5; i++) {
      const y = padding.top + (plotHeight / 5) * i;
      ctx.beginPath();
      ctx.moveTo(padding.left, y);
      ctx.lineTo(canvas.width - padding.right, y);
      ctx.stroke();
    }

    // Draw sensor positions
    ctx.strokeStyle = '#10b981';
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    sensors.forEach(sensor => {
      const sx = xScale(sensor.x);
      ctx.beginPath();
      ctx.moveTo(sx, padding.top);
      ctx.lineTo(sx, canvas.height - padding.bottom);
      ctx.stroke();
    });
    ctx.setLineDash([]);

    // Analytical (dashed)
    if (xAna && yAna && yAna.length > 0) {
      ctx.strokeStyle = colorAna;
      ctx.lineWidth = 2;
      ctx.setLineDash([8, 4]);
      ctx.beginPath();
      ctx.moveTo(xScale(xAna[0]), yScale(yAna[0]));
      for (let i = 1; i < xAna.length; i++) {
        ctx.lineTo(xScale(xAna[i]), yScale(yAna[i]));
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // Numerical (solid)
    if (xNum && yNum && yNum.length > 0) {
      ctx.strokeStyle = colorNum;
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.moveTo(xScale(xNum[0]), yScale(yNum[0]));
      for (let i = 1; i < xNum.length; i++) {
        ctx.lineTo(xScale(xNum[i]), yScale(yNum[i]));
      }
      ctx.stroke();
    }

    // Axes
    ctx.strokeStyle = '#333';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(padding.left, padding.top);
    ctx.lineTo(padding.left, canvas.height - padding.bottom);
    ctx.lineTo(canvas.width - padding.right, canvas.height - padding.bottom);
    ctx.stroke();

    // Labels
    ctx.fillStyle = '#333';
    ctx.font = 'bold 12px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText(ylabel, canvas.width / 2, 18);
    ctx.fillText('Position (m)', canvas.width / 2, canvas.height - 5);

    // Legend
    const legendX = canvas.width - padding.right + 10;
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'left';

    ctx.strokeStyle = colorNum;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(legendX, 50);
    ctx.lineTo(legendX + 20, 50);
    ctx.stroke();
    ctx.fillText('Numerical', legendX + 25, 54);

    ctx.strokeStyle = colorAna;
    ctx.setLineDash([8, 4]);
    ctx.beginPath();
    ctx.moveTo(legendX, 70);
    ctx.lineTo(legendX + 20, 70);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillText('Analytical', legendX + 25, 74);

    ctx.strokeStyle = '#10b981';
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(legendX, 90);
    ctx.lineTo(legendX + 20, 90);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillText('Sensors', legendX + 25, 94);
  };

  // Draw Monte Carlo plot with uncertainty bands
  const drawMonteCarloPlot = (canvasRef, title, x, mean, lower, upper, analytical) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    canvas.width = 500;
    canvas.height = 300;

    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    const padding = { left: 60, right: 100, top: 30, bottom: 45 };
    const plotWidth = canvas.width - padding.left - padding.right;
    const plotHeight = canvas.height - padding.top - padding.bottom;

    let allY = [...mean];
    if (lower && upper) allY = [...allY, ...lower, ...upper];
    if (analytical) allY = [...allY, ...analytical];

    const xMin = Math.min(...x);
    const xMax = Math.max(...x);
    const yMin = Math.min(...allY) * 0.95;
    const yMax = Math.max(...allY) * 1.05;

    const xScale = (val) => padding.left + ((val - xMin) / (xMax - xMin)) * plotWidth;
    const yScale = (val) => canvas.height - padding.bottom - ((val - yMin) / (yMax - yMin)) * plotHeight;

    // Grid
    ctx.strokeStyle = '#e0e0e0';
    ctx.lineWidth = 1;
    for (let i = 0; i <= 5; i++) {
      const y = padding.top + (plotHeight / 5) * i;
      ctx.beginPath();
      ctx.moveTo(padding.left, y);
      ctx.lineTo(canvas.width - padding.right, y);
      ctx.stroke();
    }

    // Uncertainty band
    if (lower && upper) {
      ctx.fillStyle = 'rgba(59, 130, 246, 0.2)';
      ctx.beginPath();
      ctx.moveTo(xScale(x[0]), yScale(upper[0]));
      for (let i = 1; i < x.length; i++) {
        ctx.lineTo(xScale(x[i]), yScale(upper[i]));
      }
      for (let i = x.length - 1; i >= 0; i--) {
        ctx.lineTo(xScale(x[i]), yScale(lower[i]));
      }
      ctx.closePath();
      ctx.fill();
    }

    // Sensor markers
    ctx.strokeStyle = '#10b981';
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    sensors.forEach(sensor => {
      const sx = xScale(sensor.x);
      ctx.beginPath();
      ctx.moveTo(sx, padding.top);
      ctx.lineTo(sx, canvas.height - padding.bottom);
      ctx.stroke();
    });
    ctx.setLineDash([]);

    // Analytical (dashed red)
    if (analytical && analytical.length > 0) {
      ctx.strokeStyle = '#dc2626';
      ctx.lineWidth = 2;
      ctx.setLineDash([8, 4]);
      ctx.beginPath();
      ctx.moveTo(xScale(x[0]), yScale(analytical[0]));
      for (let i = 1; i < x.length; i++) {
        ctx.lineTo(xScale(x[i]), yScale(analytical[i]));
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // Mean (solid blue)
    ctx.strokeStyle = '#2563eb';
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    ctx.moveTo(xScale(x[0]), yScale(mean[0]));
    for (let i = 1; i < x.length; i++) {
      ctx.lineTo(xScale(x[i]), yScale(mean[i]));
    }
    ctx.stroke();

    // Axes
    ctx.strokeStyle = '#333';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(padding.left, padding.top);
    ctx.lineTo(padding.left, canvas.height - padding.bottom);
    ctx.lineTo(canvas.width - padding.right, canvas.height - padding.bottom);
    ctx.stroke();

    // Labels
    ctx.fillStyle = '#333';
    ctx.font = 'bold 12px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText(title, canvas.width / 2, 18);

    // Legend
    const legendX = canvas.width - padding.right + 10;
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'left';

    ctx.strokeStyle = '#2563eb';
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    ctx.moveTo(legendX, 50);
    ctx.lineTo(legendX + 20, 50);
    ctx.stroke();
    ctx.fillText('Mean', legendX + 25, 54);

    ctx.strokeStyle = '#dc2626';
    ctx.setLineDash([8, 4]);
    ctx.beginPath();
    ctx.moveTo(legendX, 70);
    ctx.lineTo(legendX + 20, 70);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillText('Analytical', legendX + 25, 74);

    ctx.fillStyle = 'rgba(59, 130, 246, 0.3)';
    ctx.fillRect(legendX, 85, 20, 10);
    ctx.fillStyle = '#333';
    ctx.fillText('95% CI', legendX + 25, 94);
  };

  // Compute entropy from p and rho
  const computeEntropy = (pressure, density, gamma = 1.4) => {
    if (!pressure || !density) return null;
    return pressure.map((p, i) => {
      const rho = density[i];
      if (rho <= 0 || p <= 0) return 0;
      return p / Math.pow(rho, gamma);
    });
  };

  // Redraw plots when data changes
  useEffect(() => {
    if (monteCarloResults) {
      const { x, mean_rho, mean_u, mean_p, mean_entropy,
              ci95_lower_rho, ci95_upper_rho,
              ci95_lower_u, ci95_upper_u,
              ci95_lower_p, ci95_upper_p,
              ci95_lower_entropy, ci95_upper_entropy,
              analytical_rho, analytical_u, analytical_p, analytical_entropy } = monteCarloResults;

      drawMonteCarloPlot(densityCanvasRef, 'Density (kg/m3)', x, mean_rho, ci95_lower_rho, ci95_upper_rho, analytical_rho);
      drawMonteCarloPlot(velocityCanvasRef, 'Velocity (m/s)', x, mean_u, ci95_lower_u, ci95_upper_u, analytical_u);
      drawMonteCarloPlot(pressureCanvasRef, 'Pressure (Pa)', x, mean_p, ci95_lower_p, ci95_upper_p, analytical_p);
      drawMonteCarloPlot(entropyCanvasRef, 'Entropy (J/kg-K)', x, mean_entropy, ci95_lower_entropy, ci95_upper_entropy, analytical_entropy);

    } else if (numericalData || analyticalData) {
      const analyticalEntropy = computeEntropy(analyticalData?.pressure, analyticalData?.density, config.gamma);

      drawComparisonPlot(densityCanvasRef, numericalData?.x, numericalData?.density, analyticalData?.x, analyticalData?.density, 'Density (kg/m3)', '#2563eb', '#94a3b8');
      drawComparisonPlot(velocityCanvasRef, numericalData?.x, numericalData?.velocity, analyticalData?.x, analyticalData?.velocity, 'Velocity (m/s)', '#16a34a', '#94a3b8');
      drawComparisonPlot(pressureCanvasRef, numericalData?.x, numericalData?.pressure, analyticalData?.x, analyticalData?.pressure, 'Pressure (Pa)', '#dc2626', '#94a3b8');
      drawComparisonPlot(entropyCanvasRef, numericalData?.x, numericalData?.entropy, analyticalData?.x, analyticalEntropy, 'Entropy (J/kg-K)', '#9333ea', '#94a3b8');
    }
  }, [monteCarloResults, numericalData, analyticalData, sensors, config.gamma]);

  const styles = {
    container: {
      display: 'flex',
      minHeight: '100vh',
      fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif',
      backgroundColor: '#f8fafc'
    },
    main: {
      flex: 1,
      padding: '24px',
      overflowY: 'auto'
    },
    sidebar: {
      width: '320px',
      backgroundColor: '#fff',
      borderLeft: '1px solid #e2e8f0',
      padding: '24px',
      overflowY: 'auto'
    },
    card: {
      backgroundColor: '#fff',
      borderRadius: '8px',
      border: '1px solid #e2e8f0',
      padding: '20px',
      marginBottom: '20px'
    },
    cardTitle: {
      fontSize: '16px',
      fontWeight: '600',
      color: '#334155',
      marginTop: 0,
      marginBottom: '16px'
    },
    row: {
      display: 'flex',
      gap: '12px',
      marginBottom: '12px',
      alignItems: 'flex-end'
    },
    label: {
      fontSize: '12px',
      fontWeight: '500',
      color: '#64748b',
      marginBottom: '4px'
    },
    select: {
      padding: '8px 12px',
      fontSize: '13px',
      border: '1px solid #cbd5e1',
      borderRadius: '4px',
      backgroundColor: '#fff'
    },
    input: {
      padding: '8px 10px',
      fontSize: '13px',
      border: '1px solid #cbd5e1',
      borderRadius: '4px',
      width: '80px'
    },
    button: {
      padding: '10px 20px',
      fontSize: '14px',
      fontWeight: '500',
      border: 'none',
      borderRadius: '6px',
      cursor: 'pointer'
    },
    buttonPrimary: {
      backgroundColor: '#2563eb',
      color: '#fff'
    },
    buttonSecondary: {
      backgroundColor: '#e2e8f0',
      color: '#334155'
    },
    buttonSmall: {
      padding: '6px 12px',
      fontSize: '12px'
    },
    sensorRow: {
      display: 'grid',
      gridTemplateColumns: '1fr 1fr 1fr 1fr 30px',
      gap: '8px',
      marginBottom: '8px',
      alignItems: 'center'
    },
    sensorInput: {
      padding: '6px 8px',
      fontSize: '12px',
      border: '1px solid #cbd5e1',
      borderRadius: '4px',
      width: '100%',
      boxSizing: 'border-box'
    },
    removeBtn: {
      padding: '4px 8px',
      fontSize: '12px',
      border: 'none',
      borderRadius: '4px',
      backgroundColor: '#fee2e2',
      color: '#991b1b',
      cursor: 'pointer'
    },
    canvasWrapper: {
      marginBottom: '16px'
    },
    sidebarSection: {
      marginBottom: '24px'
    },
    sectionTitle: {
      fontSize: '12px',
      fontWeight: '600',
      color: '#64748b',
      textTransform: 'uppercase',
      marginBottom: '12px'
    },
    valueRow: {
      display: 'flex',
      justifyContent: 'space-between',
      marginBottom: '6px',
      fontSize: '13px'
    },
    valueLabel: {
      color: '#64748b'
    },
    valueData: {
      color: '#1e293b',
      fontWeight: '500'
    },
    plotsGrid: {
      display: 'grid',
      gridTemplateColumns: 'repeat(2, 1fr)',
      gap: '16px'
    }
  };

  return (
    <div style={styles.container}>
      <div style={styles.main}>
        <h1 style={{ marginTop: 0, marginBottom: '24px', color: '#1e293b' }}>
          Shock Tube Entropy Simulator
        </h1>

        {/* Test Selection */}
        <div style={styles.card}>
          <h3 style={styles.cardTitle}>Test Configuration</h3>
          <div style={styles.row}>
            <div>
              <div style={styles.label}>Toro Test Case</div>
              <select
                style={styles.select}
                value={selectedTest}
                onChange={(e) => setSelectedTest(parseInt(e.target.value))}
              >
                <option value={1}>Test 1 - Sod</option>
                <option value={2}>Test 2 - 123 Problem</option>
                <option value={3}>Test 3 - Left Blast</option>
                <option value={4}>Test 4 - Slow Shock</option>
                <option value={5}>Test 5 - Collision</option>
              </select>
            </div>
            <div>
              <div style={styles.label}>Grid Cells</div>
              <select
                style={styles.select}
                value={config.numCells}
                onChange={(e) => setConfig({...config, numCells: parseInt(e.target.value)})}
              >
                <option value={100}>100</option>
                <option value={200}>200</option>
                <option value={500}>500</option>
                <option value={1000}>1000</option>
              </select>
            </div>
            <div>
              <div style={styles.label}>Flux Scheme</div>
              <select
                style={styles.select}
                value={config.fluxType}
                onChange={(e) => setConfig({...config, fluxType: e.target.value})}
              >
                <option value="HLLC">HLLC</option>
                <option value="EntropyStable">Entropy-Stable</option>
              </select>
            </div>
            <div>
              <div style={styles.label}>Time Integration</div>
              <select
                style={styles.select}
                value={config.timeIntegration}
                onChange={(e) => setConfig({...config, timeIntegration: e.target.value})}
              >
                <option value="RK2">RK2</option>
                <option value="Hancock">MUSCL-Hancock</option>
              </select>
            </div>
          </div>
          <div style={{ marginTop: '12px' }}>
            <button
              style={{
                ...styles.button,
                ...styles.buttonPrimary,
                opacity: status.running ? 0.6 : 1
              }}
              onClick={runSimulation}
              disabled={status.running}
            >
              {status.running ? 'Running...' : 'Run Toro Test'}
            </button>
          </div>
        </div>

        {/* Sparse Sensor Data */}
        <div style={styles.card}>
          <h3 style={styles.cardTitle}>Sparse Sensor Data</h3>

          <div style={{ marginBottom: '12px' }}>
            <button
              style={{ ...styles.button, ...styles.buttonSecondary, ...styles.buttonSmall }}
              onClick={loadToroTest}
            >
              Load Toro Test Data
            </button>
          </div>

          <div style={{ marginBottom: '12px' }}>
            <div style={{ ...styles.sensorRow, fontWeight: '600', color: '#64748b', fontSize: '11px' }}>
              <span>Position</span>
              <span>Density</span>
              <span>Velocity</span>
              <span>Pressure</span>
              <span></span>
            </div>
            {sensors.map((sensor, idx) => (
              <div key={idx} style={styles.sensorRow}>
                <input
                  type="number"
                  style={styles.sensorInput}
                  value={sensor.x}
                  onChange={(e) => updateSensor(idx, 'x', e.target.value)}
                  step="0.05"
                />
                <input
                  type="number"
                  style={styles.sensorInput}
                  value={sensor.rho}
                  onChange={(e) => updateSensor(idx, 'rho', e.target.value)}
                  step="0.1"
                />
                <input
                  type="number"
                  style={styles.sensorInput}
                  value={sensor.u}
                  onChange={(e) => updateSensor(idx, 'u', e.target.value)}
                  step="0.1"
                />
                <input
                  type="number"
                  style={styles.sensorInput}
                  value={sensor.p}
                  onChange={(e) => updateSensor(idx, 'p', e.target.value)}
                  step="0.1"
                />
                <button
                  style={styles.removeBtn}
                  onClick={() => removeSensor(idx)}
                  disabled={sensors.length <= 2}
                >
                  x
                </button>
              </div>
            ))}
          </div>

          <div style={{ display: 'flex', gap: '12px', alignItems: 'center' }}>
            <button
              style={{ ...styles.button, ...styles.buttonSecondary, ...styles.buttonSmall }}
              onClick={addSensor}
              disabled={sensors.length >= 10}
            >
              + Add Sensor
            </button>
            <div>
              <select
                style={styles.select}
                value={config.interpolation}
                onChange={(e) => setConfig({...config, interpolation: e.target.value})}
              >
                <option value="linear">Linear Interpolation</option>
                <option value="piecewise_constant">Piecewise Constant</option>
              </select>
            </div>
          </div>

          <div style={{ marginTop: '12px' }}>
            <button
              style={{
                ...styles.button,
                ...styles.buttonPrimary,
                backgroundColor: '#059669',
                opacity: status.running ? 0.6 : 1
              }}
              onClick={runFromSensors}
              disabled={status.running}
            >
              {status.running ? 'Running...' : 'Run from Sensors'}
            </button>
          </div>
        </div>

        {/* Monte Carlo */}
        <div style={styles.card}>
          <h3 style={styles.cardTitle}>Monte Carlo Uncertainty Analysis</h3>
          <div style={styles.row}>
            <div>
              <div style={styles.label}>Sensor Noise (%)</div>
              <input
                type="number"
                style={styles.input}
                value={monteCarloConfig.noiseLevel * 100}
                onChange={(e) => setMonteCarloConfig({
                  ...monteCarloConfig,
                  noiseLevel: parseFloat(e.target.value) / 100
                })}
                step="1"
                min="0"
                max="20"
              />
            </div>
            <div>
              <div style={styles.label}>Number of Trials</div>
              <select
                style={styles.select}
                value={monteCarloConfig.numTrials}
                onChange={(e) => setMonteCarloConfig({
                  ...monteCarloConfig,
                  numTrials: parseInt(e.target.value)
                })}
              >
                <option value={50}>50</option>
                <option value={100}>100</option>
                <option value={200}>200</option>
                <option value={500}>500</option>
              </select>
            </div>
          </div>
          <div style={{ marginTop: '12px' }}>
            <button
              style={{
                ...styles.button,
                ...styles.buttonPrimary,
                backgroundColor: '#7c3aed',
                opacity: status.running ? 0.6 : 1
              }}
              onClick={runMonteCarlo}
              disabled={status.running}
            >
              {status.running ? 'Running...' : 'Run Monte Carlo'}
            </button>
          </div>
        </div>

        {/* Results */}
        <div style={styles.card}>
          <h3 style={styles.cardTitle}>Results</h3>
          {(numericalData || monteCarloResults) ? (
            <div style={styles.plotsGrid}>
              <div style={styles.canvasWrapper}>
                <canvas ref={densityCanvasRef} />
              </div>
              <div style={styles.canvasWrapper}>
                <canvas ref={velocityCanvasRef} />
              </div>
              <div style={styles.canvasWrapper}>
                <canvas ref={pressureCanvasRef} />
              </div>
              <div style={styles.canvasWrapper}>
                <canvas ref={entropyCanvasRef} />
              </div>
            </div>
          ) : (
            <div style={{ padding: '60px', textAlign: 'center', color: '#94a3b8' }}>
              Run a simulation to view results
            </div>
          )}
        </div>
      </div>

      {/* Sidebar */}
      <div style={styles.sidebar}>
        <h2 style={{ fontSize: '18px', fontWeight: '600', marginTop: 0, marginBottom: '20px', color: '#1e293b' }}>
          Analysis
        </h2>

        {/* Validation Results */}
        {validationResults && (
          <div style={styles.sidebarSection}>
            <div style={styles.sectionTitle}>Validation vs Analytical</div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Density Error</span>
              <span style={{
                ...styles.valueData,
                color: validationResults.density_error_percent < 3 ? '#166534' : '#991b1b'
              }}>
                {validationResults.density_error_percent?.toFixed(2)}%
              </span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Velocity Error</span>
              <span style={{
                ...styles.valueData,
                color: validationResults.velocity_error_percent < 3 ? '#166534' : '#991b1b'
              }}>
                {validationResults.velocity_error_percent?.toFixed(2)}%
              </span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Pressure Error</span>
              <span style={{
                ...styles.valueData,
                color: validationResults.pressure_error_percent < 3 ? '#166534' : '#991b1b'
              }}>
                {validationResults.pressure_error_percent?.toFixed(2)}%
              </span>
            </div>
            <div style={{ ...styles.valueRow, marginTop: '8px', paddingTop: '8px', borderTop: '1px solid #e2e8f0' }}>
              <span style={styles.valueLabel}>L1 Error</span>
              <span style={styles.valueData}>{validationResults.L1_error?.toFixed(4)}</span>
            </div>
          </div>
        )}

        {/* Monte Carlo Results */}
        {monteCarloResults && (
          <div style={styles.sidebarSection}>
            <div style={styles.sectionTitle}>Monte Carlo Results</div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Trials</span>
              <span style={styles.valueData}>{monteCarloResults.num_trials}</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Noise Level</span>
              <span style={styles.valueData}>{(monteCarloResults.noise_level * 100).toFixed(1)}%</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Compute Time</span>
              <span style={styles.valueData}>{monteCarloResults.computation_time_ms?.toFixed(0)} ms</span>
            </div>
            <div style={{ marginTop: '12px', fontSize: '12px', color: '#64748b' }}>
              <strong>Max Std Dev:</strong>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Density</span>
              <span style={styles.valueData}>{Math.max(...monteCarloResults.std_rho).toFixed(4)}</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Velocity</span>
              <span style={styles.valueData}>{Math.max(...monteCarloResults.std_u).toFixed(4)}</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Pressure</span>
              <span style={styles.valueData}>{Math.max(...monteCarloResults.std_p).toFixed(2)}</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Entropy</span>
              <span style={styles.valueData}>{Math.max(...monteCarloResults.std_entropy).toFixed(4)}</span>
            </div>
          </div>
        )}

        {/* Configuration */}
        <div style={styles.sidebarSection}>
          <div style={styles.sectionTitle}>Configuration</div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Test Case</span>
            <span style={styles.valueData}>Toro {selectedTest}</span>
          </div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Grid Cells</span>
            <span style={styles.valueData}>{config.numCells}</span>
          </div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Flux</span>
            <span style={styles.valueData}>{config.fluxType}</span>
          </div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Time Integration</span>
            <span style={styles.valueData}>{config.timeIntegration}</span>
          </div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>CFL</span>
            <span style={styles.valueData}>{config.CFL}</span>
          </div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>End Time</span>
            <span style={styles.valueData}>{config.endTime} s</span>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ShockTubeSimulator;