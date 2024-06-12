# GR Particle Path Tool

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [GUI Features](#gui-features)
- [Metrics](#metrics)
- [Contribution](#contribution)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Introduction
The GR Particle Path Research Tool is designed for researchers and enthusiasts in the field of general relativity and black hole physics. This tool provides a comprehensive suite for simulating and visualizing the paths of particles around black holes, allowing for detailed analysis and exploration of various metrics.

## Features
- **Metric Selection**: Choose from various predefined metrics such as Schwarzschild, Kerr, Reissner-Nordström, and more.
- **Parameter Input**: Set initial conditions and parameters for simulations.
- **Event Horizon Detection**: Automatically detect and highlight the event horizon.
- **Data Export and Analysis**: Export simulation data in various formats and integrate with analysis tools.
- **Real-Time Parameter Adjustment**: Adjust parameters in real-time and see immediate effects.
- **Comparative Analysis**: Compare different simulations side-by-side.
- **Gravitational Wave Simulation**: Visualize gravitational wave signals emitted by black holes.
- **User-Defined Metrics**: Add custom metrics for simulations.
- **Advanced Visualization**: Utilize heat maps, vector fields, and 4D representations.
- **Batch Processing**: Run multiple simulations with different parameters.
- **Integration with Astronomy Databases**: Use real-world black hole data in simulations.
- **Machine Learning Integration**: Predict black hole behaviors using machine learning models.
- **Simulation Playback Control**: Pause, rewind, and fast-forward simulations.
- **User Customizable Interface**: Customize GUI layout and save configurations.
- **Detailed Documentation and Tutorials**: Access comprehensive documentation and interactive tutorials.

## Installation
To install the GR Particle Path Research Tool, follow these steps:

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/your-username/GR_Particle_Path_Research_Tool.git
   cd GR_Particle_Path_Research_Tool
   ```

2. **Create a Virtual Environment:**
   ```bash
   python3 -m venv venv
   source venv/bin/activate   # On Windows: venv\Scripts\activate
   ```

3. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

4. **Install the Tool:**
   ```bash
   pip install -e .
   ```

## Usage
### Running Simulations
To run a simulation, use the `run_simulation.py` script:
```bash
python run_simulation.py
```

### Launching the GUI
To launch the graphical user interface:
```bash
python gui.py
```

## GUI Features
### Metric Selection
- Select from various predefined metrics: Schwarzschild, Kerr, Reissner-Nordström, FLRW, etc.

### Parameters and Initial Conditions
- Set initial conditions for mass, spin, charge, initial positions, and velocities.

### Additional Features
- **Event Horizon Detection**: Automatically detect and highlight the event horizon.
- **Energy and Angular Momentum**: Calculate and display energy and angular momentum of the particle.
- **Parameter Sweep**: Automatically vary parameters over specified ranges and observe the outcomes.
- **Save and Load Simulations**: Save the state of a simulation to disk and reload it later.
- **Export Data**: Export simulation data in formats like CSV and JSON for further analysis.
- **Real-Time Parameter Adjustment**: Adjust parameters while the simulation is running and see the results immediately.
- **Comparative Analysis**: Compare multiple simulations side-by-side.

## Metrics
### Available Metrics
- **Schwarzschild Metric**
- **Kerr Metric**
- **Reissner-Nordström Metric**
- **Kerr-Newman Metric**
- **FLRW Metric**
- **Anti-de Sitter Metric**
- **de Sitter Metric**
- **Minkowski Metric**
- **Vaidya Metric**
- **Bertotti-Robinson Metric**
- **Göttingen Metric**

## Contribution
We welcome contributions! 

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
