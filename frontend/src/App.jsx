import React, { useState, useEffect, useRef } from 'react';

const ShockTubeSimulator = () => {
  const API_URL = 'http://localhost:8080';

  // Simulation config
  const [config, setConfig] = useState({
    numCells: 1000,
    gamma: 1.4,
    CFL: 0.5,
    endTime: 0.25,
    fluxType: 'HLLC',
    interpolation: 'piecewise_constant'
  });

  // Toro test selection
  const [selectedTest, setSelectedTest] = useState(1);

  const TORO_TEST_END_TIMES = {
    1: 0.25,    // SOD
    2: 0.15,    // 123
    3: 0.012,   // Blast Left
    4: 0.035,   // Slow Shock
    5: 0.035    // Collision
  };

  // Sensor data
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

  // Canvas refs
  const densityCanvasRef = useRef(null);
  const velocityCanvasRef = useRef(null);
  const pressureCanvasRef = useRef(null);
  const entropyCanvasRef = useRef(null);

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

        const newEndTime = TORO_TEST_END_TIMES[selectedTest];
        setConfig(prev => ({ ...prev, endTime: newEndTime }));

        console.log(`Loaded Test ${selectedTest} with endTime = ${newEndTime}s`);
      }
    } catch (error) {
      console.error('Error loading Toro test:', error);
    }
  };

  // Run complete simulation
  const runSimulation = async () => {
    try {
      setStatus({ ready: false, running: true });
      setValidationResults(null);

      // Reset
      await fetch(`${API_URL}/api/simulation/reset`, { method: 'POST' });

      // Create
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

      // Configure
      await fetch(`${API_URL}/api/simulation/configure`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ flux: config.fluxType })
      });

      // Initialize from sparse data
      await fetch(`${API_URL}/api/simulation/init/sparse`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sensors: sensors,
          interpolation: config.interpolation
        })
      });

      // Run
      await fetch(`${API_URL}/api/simulation/run`, { method: 'POST' });

      // Get numerical solution
      const numResponse = await fetch(`${API_URL}/api/simulation/data`);
      const numData = await numResponse.json();
      if (numData.success) {
        setNumericalData(numData);
      }

      // Get analytical solution
      const anaResponse = await fetch(`${API_URL}/api/simulation/analytical`);
      const anaData = await anaResponse.json();
      if (anaData.success) {
        setAnalyticalData(anaData);
      }

      // Get validation
      const valResponse = await fetch(`${API_URL}/api/simulation/validate`);
      const valData = await valResponse.json();
      if (valData.success) {
        setValidationResults(valData);
      }

      setStatus({ ready: true, running: false });

    } catch (error) {
      console.error('Error running simulation:', error);
      setStatus({ ready: false, running: false });
    }
  };

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

  // Draw comparison plot
  const drawComparisonPlot = (canvasRef, xNum, yNum, xAna, yAna, ylabel, colorNum, colorAna) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    canvas.width = 900;
    canvas.height = 280;

    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    const padding = { left: 70, right: 120, top: 30, bottom: 50 };
    const plotWidth = canvas.width - padding.left - padding.right;
    const plotHeight = canvas.height - padding.top - padding.bottom;

    // Combine data for axis scaling
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

    // Grid lines
    ctx.strokeStyle = '#e0e0e0';
    ctx.lineWidth = 1;
    for (let i = 0; i <= 5; i++) {
      const y = padding.top + (plotHeight / 5) * i;
      ctx.beginPath();
      ctx.moveTo(padding.left, y);
      ctx.lineTo(canvas.width - padding.right, y);
      ctx.stroke();
    }

    // Axes
    ctx.strokeStyle = '#333';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(padding.left, canvas.height - padding.bottom);
    ctx.lineTo(canvas.width - padding.right, canvas.height - padding.bottom);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(padding.left, padding.top);
    ctx.lineTo(padding.left, canvas.height - padding.bottom);
    ctx.stroke();

    // Draw analytical solution first (dashed, behind)
    if (xAna && yAna && xAna.length > 0) {
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

    // Draw numerical solution (solid, on top)
    if (xNum && yNum && xNum.length > 0) {
      ctx.strokeStyle = colorNum;
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      ctx.moveTo(xScale(xNum[0]), yScale(yNum[0]));
      for (let i = 1; i < xNum.length; i++) {
        ctx.lineTo(xScale(xNum[i]), yScale(yNum[i]));
      }
      ctx.stroke();
    }

    // Draw sensor positions
    ctx.fillStyle = '#e74c3c';
    for (const sensor of sensors) {
      const sx = xScale(sensor.x);
      ctx.beginPath();
      ctx.moveTo(sx, padding.top);
      ctx.lineTo(sx, canvas.height - padding.bottom);
      ctx.strokeStyle = '#e74c3c';
      ctx.lineWidth = 1;
      ctx.setLineDash([3, 3]);
      ctx.stroke();
      ctx.setLineDash([]);
      
      // Sensor marker
      ctx.beginPath();
      ctx.arc(sx, padding.top + 8, 5, 0, 2 * Math.PI);
      ctx.fill();
    }

    // Axis labels
    ctx.fillStyle = '#333';
    ctx.font = 'bold 13px Arial';
    ctx.textAlign = 'center';
    ctx.fillText('Position (m)', padding.left + plotWidth / 2, canvas.height - 10);

    ctx.save();
    ctx.translate(18, canvas.height / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(ylabel, 0, 0);
    ctx.restore();

    // Y-axis values
    ctx.textAlign = 'right';
    ctx.font = '11px Arial';
    for (let i = 0; i <= 5; i++) {
      const yVal = (yMin - yPadding) + ((yRange + 2 * yPadding) / 5) * i;
      const y = canvas.height - padding.bottom - (plotHeight / 5) * i;
      ctx.fillText(yVal.toFixed(3), padding.left - 8, y + 4);
    }

    // X-axis values
    ctx.textAlign = 'center';
    for (let i = 0; i <= 5; i++) {
      const xVal = xMin + (xRange / 5) * i;
      const x = padding.left + (plotWidth / 5) * i;
      ctx.fillText(xVal.toFixed(2), x, canvas.height - padding.bottom + 20);
    }

    // Legend
    const legendX = canvas.width - padding.right + 15;
    const legendY = padding.top + 20;
    
    ctx.font = 'bold 11px Arial';
    ctx.textAlign = 'left';
    
    // Numerical
    ctx.strokeStyle = colorNum;
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    ctx.moveTo(legendX, legendY);
    ctx.lineTo(legendX + 25, legendY);
    ctx.stroke();
    ctx.fillStyle = '#333';
    ctx.fillText('Numerical', legendX + 30, legendY + 4);
    
    // Analytical
    ctx.strokeStyle = colorAna;
    ctx.lineWidth = 2;
    ctx.setLineDash([8, 4]);
    ctx.beginPath();
    ctx.moveTo(legendX, legendY + 25);
    ctx.lineTo(legendX + 25, legendY + 25);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillText('Analytical', legendX + 30, legendY + 29);
    
    // Sensors
    ctx.fillStyle = '#e74c3c';
    ctx.beginPath();
    ctx.arc(legendX + 12, legendY + 50, 5, 0, 2 * Math.PI);
    ctx.fill();
    ctx.fillStyle = '#333';
    ctx.fillText('Sensors', legendX + 30, legendY + 54);
  };

  // Compute entropy from primitive variables: s = ln(p / rho^gamma)
  const computeEntropy = (pressure, density, gamma = 1.4) => {
    if (!pressure || !density) return null;
    return pressure.map((p, i) => {
      const rho = density[i];
      if (rho <= 0 || p <= 0) return 0;
      return Math.log(p / Math.pow(rho, gamma));
    });
  };

  // Redraw plots when data changes
  useEffect(() => {
    if (numericalData || analyticalData) {
      drawComparisonPlot(
        densityCanvasRef,
        numericalData?.x, numericalData?.density,
        analyticalData?.x, analyticalData?.density,
        'Density (kg/m³)', '#2563eb', '#94a3b8'
      );
      drawComparisonPlot(
        velocityCanvasRef,
        numericalData?.x, numericalData?.velocity,
        analyticalData?.x, analyticalData?.velocity,
        'Velocity (m/s)', '#16a34a', '#94a3b8'
      );
      drawComparisonPlot(
        pressureCanvasRef,
        numericalData?.x, numericalData?.pressure,
        analyticalData?.x, analyticalData?.pressure,
        'Pressure (Pa)', '#dc2626', '#94a3b8'
      );
      
      // Compute analytical entropy from analytical p and rho
      const analyticalEntropy = computeEntropy(
        analyticalData?.pressure, 
        analyticalData?.density, 
        config.gamma
      );
      
      drawComparisonPlot(
        entropyCanvasRef,
        numericalData?.x, numericalData?.entropy,
        analyticalData?.x, analyticalEntropy,
        'Entropy (J/kg·K)', '#9333ea', '#94a3b8'
      );
    }
  }, [numericalData, analyticalData, sensors, config.gamma]);

  // Helper to format error values
  const formatError = (value) => {
    if (value === undefined || value === null) return '-';
    if (Math.abs(value) < 0.0001) return value.toExponential(2);
    return value.toFixed(4);
  };

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
      width: '360px',
      backgroundColor: '#fff',
      borderLeft: '1px solid #e2e8f0',
      padding: '24px',
      overflowY: 'auto'
    },
    header: {
      marginBottom: '24px'
    },
    title: {
      fontSize: '24px',
      fontWeight: '700',
      color: '#1e293b',
      margin: '0 0 4px 0'
    },
    subtitle: {
      fontSize: '14px',
      color: '#64748b',
      margin: 0
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
      alignItems: 'center'
    },
    label: {
      fontSize: '13px',
      fontWeight: '500',
      color: '#475569',
      marginBottom: '4px'
    },
    input: {
      width: '70px',
      padding: '6px 8px',
      fontSize: '13px',
      border: '1px solid #cbd5e1',
      borderRadius: '4px'
    },
    select: {
      padding: '8px 12px',
      fontSize: '13px',
      border: '1px solid #cbd5e1',
      borderRadius: '4px',
      backgroundColor: '#fff'
    },
    button: {
      padding: '10px 20px',
      fontSize: '14px',
      fontWeight: '500',
      border: 'none',
      borderRadius: '6px',
      cursor: 'pointer',
      transition: 'background-color 0.2s'
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
      padding: '4px 10px',
      fontSize: '12px'
    },
    sensorRow: {
      display: 'grid',
      gridTemplateColumns: '60px 80px 80px 80px 30px',
      gap: '8px',
      alignItems: 'center',
      marginBottom: '8px',
      fontSize: '13px'
    },
    sensorInput: {
      width: '100%',
      padding: '5px 6px',
      fontSize: '12px',
      border: '1px solid #cbd5e1',
      borderRadius: '3px',
      textAlign: 'right'
    },
    removeBtn: {
      background: 'none',
      border: 'none',
      color: '#ef4444',
      cursor: 'pointer',
      fontSize: '16px',
      padding: '2px 6px'
    },
    canvasWrapper: {
      border: '1px solid #e2e8f0',
      borderRadius: '6px',
      overflow: 'hidden',
      marginBottom: '16px',
      backgroundColor: '#fff'
    },
    sidebarSection: {
      marginBottom: '24px',
      paddingBottom: '20px',
      borderBottom: '1px solid #e2e8f0'
    },
    sectionTitle: {
      fontSize: '14px',
      fontWeight: '600',
      color: '#334155',
      marginBottom: '12px'
    },
    valueRow: {
      display: 'flex',
      justifyContent: 'space-between',
      padding: '6px 0',
      fontSize: '13px'
    },
    valueLabel: {
      color: '#64748b'
    },
    valueData: {
      fontFamily: 'monospace',
      fontWeight: '500'
    },
    errorTable: {
      width: '100%',
      borderCollapse: 'collapse',
      fontSize: '12px',
      marginTop: '8px'
    },
    errorTableHeader: {
      textAlign: 'left',
      padding: '6px 4px',
      borderBottom: '2px solid #e2e8f0',
      color: '#64748b',
      fontWeight: '600'
    },
    errorTableCell: {
      padding: '6px 4px',
      borderBottom: '1px solid #f1f5f9',
      fontFamily: 'monospace',
      textAlign: 'right'
    },
    errorTableCellLabel: {
      padding: '6px 4px',
      borderBottom: '1px solid #f1f5f9',
      color: '#475569',
      textAlign: 'left'
    }
  };

  return (
    <div style={styles.container}>
      <div style={styles.main}>
        <div style={styles.header}>
          <h1 style={styles.title}>1D Shock Tube Simulator</h1>
          <p style={styles.subtitle}>Compressible Flow Reconstruction from Sparse Telemetry</p>
        </div>

        {/* Sensor Configuration */}
        <div style={styles.card}>
          <h3 style={styles.cardTitle}>Sensor Configuration</h3>
          
          <div style={{ ...styles.row, marginBottom: '16px' }}>
            <select
              style={styles.select}
              value={selectedTest}
              onChange={(e) => setSelectedTest(Number(e.target.value))}
            >
              <option value={1}>Test 1: Sod Problem</option>
              <option value={2}>Test 2: 123 Problem</option>
              <option value={3}>Test 3: Blast Wave</option>
              <option value={4}>Test 4: Slow Shock</option>
              <option value={5}>Test 5: Collision</option>
            </select>
            <button
              style={{ ...styles.button, ...styles.buttonSecondary }}
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

          <button
            style={{ ...styles.button, ...styles.buttonSecondary, ...styles.buttonSmall }}
            onClick={addSensor}
            disabled={sensors.length >= 10}
          >
            + Add Sensor
          </button>
        </div>

        {/* Simulation Settings */}
        <div style={styles.card}>
          <h3 style={styles.cardTitle}>Simulation Settings</h3>
          
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: '16px' }}>
            <div>
              <div style={styles.label}>Grid Cells</div>
              <select
                style={{ ...styles.select, width: '100%' }}
                value={config.numCells}
                onChange={(e) => setConfig({ ...config, numCells: Number(e.target.value) })}
              >
                <option value={200}>200</option>
                <option value={500}>500</option>
                <option value={1000}>1000</option>
                <option value={2000}>2000</option>
              </select>
            </div>
            <div>
              <div style={styles.label}>Flux Type</div>
              <select
                style={{ ...styles.select, width: '100%' }}
                value={config.fluxType}
                onChange={(e) => setConfig({ ...config, fluxType: e.target.value })}
              >
                <option value="HLLC">HLLC</option>
                <option value="EntropyStable">Entropy Stable (Ismael-Roe)</option>
              </select>
            </div>
            <div>
              <div style={styles.label}>Interpolation</div>
              <select
                style={{ ...styles.select, width: '100%' }}
                value={config.interpolation}
                onChange={(e) => setConfig({ ...config, interpolation: e.target.value })}
              >
                <option value="piecewise_constant">Piecewise Constant</option>
                <option value="linear">Linear</option>
              </select>
            </div>
          </div>

          <div style={{ marginTop: '20px' }}>
            <button
              style={{
                ...styles.button,
                ...styles.buttonPrimary,
                opacity: status.running ? 0.6 : 1
              }}
              onClick={runSimulation}
              disabled={status.running}
            >
              {status.running ? 'Running...' : 'Run Simulation'}
            </button>
          </div>
        </div>

        {/* Results */}
        <div style={styles.card}>
          <h3 style={styles.cardTitle}>Results</h3>
          
          {numericalData ? (
            <>
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
            </>
          ) : (
            <div style={{ padding: '60px', textAlign: 'center', color: '#94a3b8' }}>
              Configure sensors and run simulation to view results
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
            
            {/* Error Norms Table */}
            <div style={{ fontSize: '12px', fontWeight: '600', color: '#64748b', marginBottom: '8px' }}>
              Error Norms
            </div>
            <table style={styles.errorTable}>
              <thead>
                <tr>
                  <th style={styles.errorTableHeader}></th>
                  <th style={{ ...styles.errorTableHeader, textAlign: 'right' }}>L₁</th>
                  <th style={{ ...styles.errorTableHeader, textAlign: 'right' }}>L₂</th>
                  <th style={{ ...styles.errorTableHeader, textAlign: 'right' }}>L∞</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td style={styles.errorTableCellLabel}>Density</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.L1_density)}</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.L2_density)}</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.Linf_density)}</td>
                </tr>
                <tr>
                  <td style={styles.errorTableCellLabel}>Velocity</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.L1_velocity)}</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.L2_velocity)}</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.Linf_velocity)}</td>
                </tr>
                <tr>
                  <td style={styles.errorTableCellLabel}>Pressure</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.L1_pressure)}</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.L2_pressure)}</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.Linf_pressure)}</td>
                </tr>
                <tr>
                  <td style={styles.errorTableCellLabel}>Entropy</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.L1_entropy)}</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.L2_entropy)}</td>
                  <td style={styles.errorTableCell}>{formatError(validationResults.Linf_entropy)}</td>
                </tr>
              </tbody>
            </table>
          </div>
        )}

        {/* Simulation Info */}
        {numericalData && (
          <div style={styles.sidebarSection}>
            <div style={styles.sectionTitle}>Simulation Info</div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Time</span>
              <span style={styles.valueData}>{numericalData.time.toFixed(4)} s</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Steps</span>
              <span style={styles.valueData}>{numericalData.step}</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Grid Cells</span>
              <span style={styles.valueData}>{numericalData.numCells}</span>
            </div>
          </div>
        )}

        {/* Sensor Summary */}
        <div style={styles.sidebarSection}>
          <div style={styles.sectionTitle}>Sensor Summary</div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Active Sensors</span>
            <span style={styles.valueData}>{sensors.length}</span>
          </div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Coverage</span>
            <span style={styles.valueData}>
              {sensors.length > 0 
                ? `${(Math.min(...sensors.map(s => s.x))).toFixed(2)} - ${(Math.max(...sensors.map(s => s.x))).toFixed(2)} m`
                : '-'
              }
            </span>
          </div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Interpolation</span>
            <span style={styles.valueData}>{config.interpolation === 'linear' ? 'Linear' : 'Piecewise'}</span>
          </div>
        </div>

        {/* Configuration */}
        <div style={styles.sidebarSection}>
          <div style={styles.sectionTitle}>Configuration</div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Flux</span>
            <span style={styles.valueData}>{config.fluxType.toUpperCase()}</span>
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