import React, { useState, useEffect, useRef } from 'react';

const ShockTubeViewer = () => {
  const [data, setData] = useState(null);
  const [status, setStatus] = useState({ initialized: false, running: false });
  const [config, setConfig] = useState({
    length: 1.0,
    numCells: 1000,
    gamma: 1.4,
    CFL: 0.9,
    endTime: 0.2,
    scheme: 'godunov',
    riemannSolver: 'hllc',
    rho_L: 1.0,
    u_L: 0.0,
    p_L: 1.0,
    rho_R: 0.125,
    u_R: 0.0,
    p_R: 0.1,
    x_diaphragm: 0.5
  });
  const [simulationInfo, setSimulationInfo] = useState(null);
  const [sensorPositions, setSensorPositions] = useState([0.1, 0.3, 0.5, 0.7, 0.9]);
  const [sensorReadings, setSensorReadings] = useState(null);
  const [selectedToroTest, setSelectedToroTest] = useState(1);
  const [validationResults, setValidationResults] = useState(null);

  const densityCanvasRef = useRef(null);
  const velocityCanvasRef = useRef(null);
  const pressureCanvasRef = useRef(null);
  const entropyCanvasRef = useRef(null);

  const API_URL = 'http://localhost:8080';

  const fetchData = async () => {
    try {
      const response = await fetch(`${API_URL}/api/simulation/data`);
      const jsonData = await response.json();
      if (jsonData.success) {
        setData(jsonData);
      }
    } catch (error) {
      console.error('Error fetching data:', error);
    }
  };

  const fetchStatus = async () => {
    try {
      const response = await fetch(`${API_URL}/api/simulation/status`);
      const statusData = await response.json();
      if (statusData.success) {
        setSimulationInfo(statusData);
      }
    } catch (error) {
      console.error('Error fetching status:', error);
    }
  };

  const fetchValidation = async () => {
    try {
      const response = await fetch(`${API_URL}/api/simulation/validate`);
      const data = await response.json();
      if (data.success) {
        setValidationResults(data);
      }
    } catch (error) {
      console.error('Error fetching validation:', error);
    }
  };

  const initializeSimulation = async () => {
    try {
      await fetch(`${API_URL}/api/simulation/reset`, { method: 'POST' });

      await fetch(`${API_URL}/api/simulation/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          length: config.length,
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
          scheme: config.scheme,
          riemann_solver: config.riemannSolver
        })
      });

      await fetch(`${API_URL}/api/simulation/init/shocktube`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          rho_L: config.rho_L,
          u_L: config.u_L,
          p_L: config.p_L,
          rho_R: config.rho_R,
          u_R: config.u_R,
          p_R: config.p_R,
          x_diaphragm: config.x_diaphragm
        })
      });

      setStatus({ initialized: true, running: false });
      setValidationResults(null);
      await fetchData();
    } catch (error) {
      console.error('Error initializing:', error);
    }
  };

  const placeSensors = async () => {
    try {
      await fetch(`${API_URL}/api/simulation/sensors/place`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ positions: sensorPositions })
      });
      console.log('Sensors placed');
    } catch (error) {
      console.error('Error placing sensors:', error);
    }
  };

  const fetchSensorReadings = async () => {
    try {
      const response = await fetch(`${API_URL}/api/simulation/sensors/readings`);
      const data = await response.json();
      if (data.success) {
        setSensorReadings(data);
      }
    } catch (error) {
      console.error('Error fetching sensor readings:', error);
    }
  };

  const initializeToroTest = async () => {
    try {
      await fetch(`${API_URL}/api/simulation/reset`, { method: 'POST' });

      await fetch(`${API_URL}/api/simulation/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          length: config.length,
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
          scheme: config.scheme,
          riemann_solver: config.riemannSolver
        })
      });

      await fetch(`${API_URL}/api/simulation/init/toro`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ test: selectedToroTest })
      });

      await placeSensors();

      setStatus({ initialized: true, running: false });
      setValidationResults(null);
      await fetchData();
    } catch (error) {
      console.error('Error initializing Toro Test:', error);
    }
  };

  const runSimulation = async () => {
    try {
      setStatus(prev => ({ ...prev, running: true }));

      const response = await fetch(`${API_URL}/api/simulation/run`, { method: 'POST' });
      const result = await response.json();

      if (result.success) {
        setSimulationInfo(result);
        await fetchData();
        await fetchSensorReadings();
        await fetchValidation();
      }

      setStatus(prev => ({ ...prev, running: false }));
    } catch (error) {
      console.error('Error running simulation:', error);
      setStatus(prev => ({ ...prev, running: false }));
    }
  };

  useEffect(() => {
    const checkHealth = async () => {
      try {
        await fetch(`${API_URL}/api/health`);
      } catch (error) {
        console.error('Cannot connect to API server');
      }
    };
    checkHealth();
  }, []);

  const drawPlot = (canvasRef, xData, yData, ylabel, color) => {
    const canvas = canvasRef.current;
    if (!canvas || !xData || !yData) return;

    const ctx = canvas.getContext('2d');
    canvas.width = 1000;
    canvas.height = 250;

    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    const padding = { left: 60, right: 20, top: 20, bottom: 40 };
    const plotWidth = canvas.width - padding.left - padding.right;
    const plotHeight = canvas.height - padding.top - padding.bottom;

    const xMin = Math.min(...xData);
    const xMax = Math.max(...xData);
    const yMin = Math.min(...yData);
    const yMax = Math.max(...yData);

    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;

    const xScale = (x) => padding.left + ((x - xMin) / xRange) * plotWidth;
    const yScale = (y) => canvas.height - padding.bottom - ((y - yMin) / yRange) * plotHeight;

    ctx.strokeStyle = '#ddd';
    ctx.lineWidth = 1;
    for (let i = 0; i <= 5; i++) {
      const y = padding.top + (plotHeight / 5) * i;
      ctx.beginPath();
      ctx.moveTo(padding.left, y);
      ctx.lineTo(canvas.width - padding.right, y);
      ctx.stroke();
    }

    ctx.strokeStyle = '#000';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(padding.left, canvas.height - padding.bottom);
    ctx.lineTo(canvas.width - padding.right, canvas.height - padding.bottom);
    ctx.stroke();

    ctx.beginPath();
    ctx.moveTo(padding.left, padding.top);
    ctx.lineTo(padding.left, canvas.height - padding.bottom);
    ctx.stroke();

    ctx.strokeStyle = color;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(xScale(xData[0]), yScale(yData[0]));
    for (let i = 1; i < xData.length; i++) {
      ctx.lineTo(xScale(xData[i]), yScale(yData[i]));
    }
    ctx.stroke();

    ctx.fillStyle = '#000';
    ctx.font = 'bold 12px Arial';
    ctx.textAlign = 'center';
    ctx.fillText('Position (m)', canvas.width / 2, canvas.height - 5);

    ctx.save();
    ctx.translate(15, canvas.height / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.fillText(ylabel, 0, 0);
    ctx.restore();

    ctx.textAlign = 'right';
    ctx.font = '11px Arial';
    for (let i = 0; i <= 5; i++) {
      const yVal = yMin + (yRange / 5) * i;
      const y = canvas.height - padding.bottom - (plotHeight / 5) * i;
      ctx.fillText(yVal.toFixed(3), padding.left - 5, y + 4);
    }

    ctx.textAlign = 'center';
    for (let i = 0; i <= 5; i++) {
      const xVal = xMin + (xRange / 5) * i;
      const x = padding.left + (plotWidth / 5) * i;
      ctx.fillText(xVal.toFixed(2), x, canvas.height - padding.bottom + 20);
    }

    if (sensorPositions && sensorPositions.length > 0) {
      ctx.fillStyle = '#ff0000';
      ctx.strokeStyle = '#ff0000';
      ctx.lineWidth = 2;

      for (const sensorPos of sensorPositions) {
        const x = xScale(sensorPos);

        ctx.beginPath();
        ctx.moveTo(x, padding.top);
        ctx.lineTo(x, canvas.height - padding.bottom);
        ctx.setLineDash([5, 5]);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.beginPath();
        ctx.arc(x, padding.top + 10, 5, 0, 2 * Math.PI);
        ctx.fill();
      }
    }
  };

  useEffect(() => {
    if (!data) return;

    drawPlot(densityCanvasRef, data.x, data.density, 'Density (kg/m³)', '#2563eb');
    drawPlot(velocityCanvasRef, data.x, data.velocity, 'Velocity (m/s)', '#16a34a');
    drawPlot(pressureCanvasRef, data.x, data.pressure, 'Pressure (Pa)', '#dc2626');
    
    if (data.entropy) {
      drawPlot(entropyCanvasRef, data.x, data.entropy, 'Specific Entropy (J/kg·K)', '#9333ea');
    }
  }, [data, sensorPositions]);

  const diagnostics = React.useMemo(() => {
    if (!data) return null;

    const maxDensity = Math.max(...data.density);
    const minDensity = Math.min(...data.density);
    const maxVelocity = Math.max(...data.velocity);
    const maxPressure = Math.max(...data.pressure);

    let entropyGeneration = null;
    if (data.entropy) {
      const maxEntropy = Math.max(...data.entropy);
      const minEntropy = Math.min(...data.entropy);
      entropyGeneration = (maxEntropy - minEntropy).toExponential(3);
    }

    return {
      maxDensity: maxDensity.toFixed(6),
      minDensity: minDensity.toFixed(6),
      maxVelocity: maxVelocity.toFixed(6),
      maxPressure: maxPressure.toFixed(6),
      entropyGeneration
    };
  }, [data]);

  const styles = {
    container: {
      display: 'flex',
      height: '100vh',
      fontFamily: 'Arial, sans-serif',
      backgroundColor: '#f5f5f5'
    },
    mainContent: {
      flex: 1,
      padding: '30px',
      overflowY: 'auto'
    },
    sidebar: {
      width: '300px',
      backgroundColor: '#ffffff',
      borderLeft: '1px solid #ddd',
      padding: '20px',
      boxShadow: '-2px 0 5px rgba(0,0,0,0.1)',
      overflowY: 'auto'
    },
    header: {
      marginBottom: '20px',
      paddingBottom: '15px',
      borderBottom: '2px solid #333'
    },
    title: {
      fontSize: '24px',
      fontWeight: 'bold',
      margin: '0 0 10px 0',
      color: '#333'
    },
    statusBar: {
      backgroundColor: '#ffffff',
      padding: '15px',
      borderRadius: '5px',
      marginBottom: '20px',
      boxShadow: '0 2px 4px rgba(0,0,0,0.1)'
    },
    controlPanel: {
      backgroundColor: '#ffffff',
      padding: '20px',
      borderRadius: '5px',
      marginBottom: '20px',
      boxShadow: '0 2px 4px rgba(0,0,0,0.1)'
    },
    button: {
      padding: '10px 20px',
      marginRight: '10px',
      marginBottom: '10px',
      border: 'none',
      borderRadius: '4px',
      cursor: 'pointer',
      fontSize: '14px',
      fontWeight: '500',
      transition: 'background-color 0.2s'
    },
    buttonPrimary: {
      backgroundColor: '#4CAF50',
      color: 'white'
    },
    buttonDanger: {
      backgroundColor: '#f44336',
      color: 'white'
    },
    buttonSecondary: {
      backgroundColor: '#2196F3',
      color: 'white'
    },
    buttonDisabled: {
      backgroundColor: '#ccc',
      cursor: 'not-allowed'
    },
    visualizationContainer: {
      backgroundColor: '#ffffff',
      padding: '20px',
      borderRadius: '5px',
      boxShadow: '0 2px 4px rgba(0,0,0,0.1)',
      marginBottom: '20px'
    },
    canvasWrapper: {
      border: '2px solid #333',
      backgroundColor: '#fff',
      marginBottom: '20px'
    },
    sidebarSection: {
      marginBottom: '25px',
      paddingBottom: '20px',
      borderBottom: '1px solid #eee'
    },
    sectionTitle: {
      fontSize: '16px',
      fontWeight: 'bold',
      marginBottom: '12px',
      color: '#555'
    },
    valueRow: {
      display: 'flex',
      justifyContent: 'space-between',
      padding: '8px 0',
      fontSize: '14px'
    },
    valueLabel: {
      color: '#666',
      fontWeight: '500'
    },
    valueData: {
      color: '#333',
      fontFamily: 'monospace'
    },
    formGroup: {
      marginBottom: '15px'
    },
    label: {
      display: 'block',
      marginBottom: '5px',
      fontWeight: 'bold',
      fontSize: '14px'
    },
    input: {
      width: '100%',
      padding: '8px',
      fontSize: '14px',
      borderRadius: '4px',
      border: '1px solid #ccc'
    },
    select: {
      width: '100%',
      padding: '8px',
      fontSize: '14px',
      borderRadius: '4px',
      border: '1px solid #ccc'
    }
  };

  return (
    <div style={styles.container}>
      <div style={styles.mainContent}>
        <div style={styles.header}>
          <h1 style={styles.title}>1D Shock Tube Simulation</h1>
          <p style={{ margin: 0, color: '#666' }}>Compressible Euler Equations</p>
        </div>

        <div style={styles.statusBar}>
          <div style={{ display: 'flex', justifyContent: 'space-between' }}>
            <div>
              <strong>Status:</strong> {status.initialized ? (status.running ? 'Running' : 'Ready') : 'Not Initialized'}
            </div>
            <div>
              <strong>Scheme:</strong> {config.scheme.toUpperCase()} | <strong>Solver:</strong> {config.riemannSolver.toUpperCase()}
            </div>
            <div>
              <strong>Grid:</strong> {config.numCells} cells
            </div>
          </div>
        </div>

        <div style={styles.controlPanel}>
          <h3 style={{ marginTop: 0, marginBottom: '15px' }}>Test Selection</h3>

          <div style={styles.formGroup}>
            <label style={styles.label}>Toro Test Case</label>
            <select
              style={styles.select}
              value={selectedToroTest}
              onChange={(e) => setSelectedToroTest(Number(e.target.value))}
            >
              <option value="1">Test 1: Sod Problem</option>
              <option value="2">Test 2: 123 Problem (Strong Rarefraction)</option>
              <option value="3">Test 3: Blast Wave (Left Half)</option>
              <option value="4">Test 4: Collision</option>
              <option value="5">Test 5: Stationary Contact</option>
            </select>
          </div>
          
          <button
            onClick={initializeToroTest}
            style={{
              ...styles.button,
              backgroundColor: '#9C27B0',
              color: 'white'
            }}
          >
            Initialize Toro Test
          </button>

          <button
            onClick={initializeSimulation}
            style={{
              ...styles.button,
              ...styles.buttonSecondary
            }}
          >
            Initialize Custom
          </button>
        </div>

        <div style={styles.controlPanel}>
          <h3 style={{ marginTop: 0, marginBottom: '15px' }}>Sensor Configuration</h3>

          <div style={styles.formGroup}>
            <label style={styles.label}>Sensor Positions (comma-separated, 0-1)</label>
            <input
              style={styles.input}
              type="text"
              value={sensorPositions.join(', ')}
              onChange={(e) => {
                const positions = e.target.value.split(',').map(s => parseFloat(s.trim())).filter(n => !isNaN(n));
                setSensorPositions(positions);
              }}
            />
          </div>

          <button
            onClick={placeSensors}
            disabled={!status.initialized}
            style={{
              ...styles.button,
              ...(status.initialized ? styles.buttonSecondary : styles.buttonDisabled)
            }}
          >
            Update Sensors
          </button>
        </div>

        <div style={styles.controlPanel}>
          <h3 style={{ marginTop: 0, marginBottom: '15px' }}>Simulation Parameters</h3>

          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '15px' }}>
            <div style={styles.formGroup}>
              <label style={styles.label}>Grid Cells</label>
              <select
                style={styles.select}
                value={config.numCells}
                onChange={(e) => setConfig({ ...config, numCells: Number(e.target.value) })}
              >
                <option value="100">100</option>
                <option value="500">500</option>
                <option value="1000">1000</option>
                <option value="2000">2000</option>
              </select>
            </div>

            <div style={styles.formGroup}>
              <label style={styles.label}>Scheme</label>
              <select
                style={styles.select}
                value={config.scheme}
                onChange={(e) => setConfig({ ...config, scheme: e.target.value })}
              >
                <option value="godunov">Godunov (1st Order)</option>
                <option value="muscl">MUSCL (2nd Order)</option>
              </select>
            </div>

            <div style={styles.formGroup}>
              <label style={styles.label}>Riemann Solver</label>
              <select
                style={styles.select}
                value={config.riemannSolver}
                onChange={(e) => setConfig({ ...config, riemannSolver: e.target.value })}
              >
                <option value="hllc">HLLC</option>
                <option value="exact">Exact</option>
              </select>
            </div>

            <div style={styles.formGroup}>
              <label style={styles.label}>CFL Number</label>
              <input
                style={styles.input}
                type="number"
                value={config.CFL}
                onChange={(e) => setConfig({ ...config, CFL: parseFloat(e.target.value) })}
                step="0.1"
              />
            </div>

            <div style={styles.formGroup}>
              <label style={styles.label}>End Time (s)</label>
              <input
                style={styles.input}
                type="number"
                value={config.endTime}
                onChange={(e) => setConfig({ ...config, endTime: parseFloat(e.target.value) })}
                step="0.01"
              />
            </div>

            <div style={styles.formGroup}>
              <label style={styles.label}>Gamma</label>
              <input
                style={styles.input}
                type="number"
                value={config.gamma}
                onChange={(e) => setConfig({ ...config, gamma: parseFloat(e.target.value) })}
                step="0.01"
              />
            </div>
          </div>

          <div style={{ marginTop: '20px' }}>
            <h4 style={{ marginBottom: '10px' }}>Initial Conditions</h4>
            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '20px' }}>
              <div>
                <h5 style={{ marginBottom: '10px', color: '#666' }}>Left State</h5>
                <div style={styles.formGroup}>
                  <label style={styles.label}>Density (kg/m³)</label>
                  <input
                    style={styles.input}
                    type="number"
                    value={config.rho_L}
                    onChange={(e) => setConfig({ ...config, rho_L: parseFloat(e.target.value) })}
                    step="0.1"
                  />
                </div>
                <div style={styles.formGroup}>
                  <label style={styles.label}>Velocity (m/s)</label>
                  <input
                    style={styles.input}
                    type="number"
                    value={config.u_L}
                    onChange={(e) => setConfig({ ...config, u_L: parseFloat(e.target.value) })}
                    step="0.1"
                  />
                </div>
                <div style={styles.formGroup}>
                  <label style={styles.label}>Pressure (Pa)</label>
                  <input
                    style={styles.input}
                    type="number"
                    value={config.p_L}
                    onChange={(e) => setConfig({ ...config, p_L: parseFloat(e.target.value) })}
                    step="0.1"
                  />
                </div>
              </div>

              <div>
                <h5 style={{ marginBottom: '10px', color: '#666' }}>Right State</h5>
                <div style={styles.formGroup}>
                  <label style={styles.label}>Density (kg/m³)</label>
                  <input
                    style={styles.input}
                    type="number"
                    value={config.rho_R}
                    onChange={(e) => setConfig({ ...config, rho_R: parseFloat(e.target.value) })}
                    step="0.01"
                  />
                </div>
                <div style={styles.formGroup}>
                  <label style={styles.label}>Velocity (m/s)</label>
                  <input
                    style={styles.input}
                    type="number"
                    value={config.u_R}
                    onChange={(e) => setConfig({ ...config, u_R: parseFloat(e.target.value) })}
                    step="0.1"
                  />
                </div>
                <div style={styles.formGroup}>
                  <label style={styles.label}>Pressure (Pa)</label>
                  <input
                    style={styles.input}
                    type="number"
                    value={config.p_R}
                    onChange={(e) => setConfig({ ...config, p_R: parseFloat(e.target.value) })}
                    step="0.01"
                  />
                </div>
              </div>
            </div>
          </div>

          <div style={{ marginTop: '20px' }}>
            <button
              onClick={initializeSimulation}
              style={{
                ...styles.button,
                ...styles.buttonSecondary
              }}
            >
              Initialize
            </button>

            <button
              onClick={runSimulation}
              disabled={!status.initialized || status.running}
              style={{
                ...styles.button,
                ...(status.initialized && !status.running ? styles.buttonPrimary : styles.buttonDisabled)
              }}
            >
              {status.running ? 'Running...' : 'Run Simulation'}
            </button>

            <button
              onClick={async () => {
                await fetch(`${API_URL}/api/simulation/reset`, { method: 'POST' });
                setStatus({ initialized: false, running: false });
                setData(null);
                setSimulationInfo(null);
                setValidationResults(null);
              }}
              style={{
                ...styles.button,
                backgroundColor: '#FF9800',
                color: 'white'
              }}
            >
              Reset
            </button>
          </div>
        </div>

        <div style={styles.visualizationContainer}>
          <h3 style={{ marginTop: 0, marginBottom: '15px' }}>Results</h3>
          {data ? (
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
              {data.entropy && (
                <div style={styles.canvasWrapper}>
                  <canvas ref={entropyCanvasRef} />
                </div>
              )}
            </>
          ) : (
            <div style={{ padding: '40px', textAlign: 'center', color: '#999' }}>
              Initialize and run simulation to view results
            </div>
          )}
        </div>
      </div>

      <div style={styles.sidebar}>
        <h2 style={{ marginTop: 0, marginBottom: '20px', fontSize: '18px' }}>Diagnostics</h2>

        {validationResults && (
          <div style={styles.sidebarSection}>
            <div style={styles.sectionTitle}>Validation (vs Analytical)</div>
            
            <div style={{
              padding: '12px',
              borderRadius: '6px',
              marginBottom: '15px',
              textAlign: 'center',
              fontWeight: 'bold',
              fontSize: '16px',
              backgroundColor: validationResults.passes ? '#d4edda' : '#f8d7da',
              color: validationResults.passes ? '#155724' : '#721c24',
              border: `2px solid ${validationResults.passes ? '#c3e6cb' : '#f5c6cb'}`
            }}>
              {validationResults.passes ? 'PASSES' : 'FAILS'} 3% Threshold
            </div>
            
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Density Error:</span>
              <span style={{
                ...styles.valueData,
                color: validationResults.density_error_percent < 3 ? '#155724' : '#721c24',
                fontWeight: 'bold'
              }}>
                {validationResults.density_error_percent.toFixed(2)}%
              </span>
            </div>
            
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Velocity Error:</span>
              <span style={{
                ...styles.valueData,
                color: validationResults.velocity_error_percent < 3 ? '#155724' : '#721c24',
                fontWeight: 'bold'
              }}>
                {validationResults.velocity_error_percent.toFixed(2)}%
              </span>
            </div>
            
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Pressure Error:</span>
              <span style={{
                ...styles.valueData,
                color: validationResults.pressure_error_percent < 3 ? '#155724' : '#721c24',
                fontWeight: 'bold'
              }}>
                {validationResults.pressure_error_percent.toFixed(2)}%
              </span>
            </div>
            
            <div style={{
              marginTop: '10px',
              padding: '8px',
              backgroundColor: '#f8f9fa',
              borderRadius: '4px',
              fontSize: '12px',
              color: '#666'
            }}>
              L1 Error: {validationResults.L1_error.toFixed(4)}
            </div>
          </div>
        )}

        {simulationInfo && (
          <div style={styles.sidebarSection}>
            <div style={styles.sectionTitle}>Simulation Info</div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Final Time:</span>
              <span style={styles.valueData}>{simulationInfo.time?.toFixed(4)} s</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Total Steps:</span>
              <span style={styles.valueData}>{simulationInfo.step}</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Total Mass:</span>
              <span style={styles.valueData}>{simulationInfo.total_mass?.toFixed(6)}</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Total Energy:</span>
              <span style={styles.valueData}>{simulationInfo.total_energy?.toFixed(6)}</span>
            </div>
          </div>
        )}

        {diagnostics && (
          <div style={styles.sidebarSection}>
            <div style={styles.sectionTitle}>Flow Properties</div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Max Density:</span>
              <span style={styles.valueData}>{diagnostics.maxDensity}</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Min Density:</span>
              <span style={styles.valueData}>{diagnostics.minDensity}</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Max Velocity:</span>
              <span style={styles.valueData}>{diagnostics.maxVelocity}</span>
            </div>
            <div style={styles.valueRow}>
              <span style={styles.valueLabel}>Max Pressure:</span>
              <span style={styles.valueData}>{diagnostics.maxPressure}</span>
            </div>
            {diagnostics.entropyGeneration && (
              <div style={styles.valueRow}>
                <span style={styles.valueLabel}>Entropy Generation:</span>
                <span style={styles.valueData}>{diagnostics.entropyGeneration}</span>
              </div>
            )}
          </div>
        )}

        {sensorReadings && (
          <div style={styles.sidebarSection}>
            <div style={styles.sectionTitle}>Sensor Readings</div>
            {sensorReadings.readings.map((reading, idx) => (
              <div key={idx} style={{ 
                marginBottom: '15px', 
                padding: '10px', 
                backgroundColor: '#f8f9fa', 
                borderRadius: '4px',
                border: '1px solid #e0e0e0'
              }}>
                <div style={{ fontWeight: 'bold', marginBottom: '5px', color: '#d32f2f' }}>
                  Sensor {idx + 1} @ x={reading.position.toFixed(3)}m
                </div>
                <div style={styles.valueRow}>
                  <span style={styles.valueLabel}>Density:</span>
                  <span style={styles.valueData}>{reading.density.toFixed(4)}</span>
                </div>
                <div style={styles.valueRow}>
                  <span style={styles.valueLabel}>Velocity:</span>
                  <span style={styles.valueData}>{reading.velocity.toFixed(4)}</span>
                </div>
                <div style={styles.valueRow}>
                  <span style={styles.valueLabel}>Pressure:</span>
                  <span style={styles.valueData}>{reading.pressure.toFixed(4)}</span>
                </div>
                <div style={styles.valueRow}>
                  <span style={styles.valueLabel}>Entropy:</span>
                  <span style={styles.valueData}>{reading.entropy.toFixed(4)}</span>
                </div>
              </div>
            ))}
          </div>
        )}

        <div style={styles.sidebarSection}>
          <div style={styles.sectionTitle}>Configuration</div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Domain Length:</span>
            <span style={styles.valueData}>{config.length} m</span>
          </div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Grid Cells:</span>
            <span style={styles.valueData}>{config.numCells}</span>
          </div>
          <div style={styles.valueRow}>
            <span style={styles.valueLabel}>Cell Width:</span>
            <span style={styles.valueData}>{(config.length / config.numCells).toFixed(6)} m</span>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ShockTubeViewer;