import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from gr_particle_path_tool.metrics import schwarzschild_metric
from gr_particle_path_tool.utils import christoffel_symbols, geodesic_equations, solve_geodesic

def run_simulation():
    try:
        M = float(mass_entry.get())
        r0 = float(r_entry.get())
        theta0 = float(theta_entry.get())
        phi0 = float(phi_entry.get())
        dr_dtau0 = float(dr_dtau_entry.get())
        dtheta_dtau0 = float(dtheta_dtau_entry.get())
        
        # Calculate the angular velocity for a stable circular orbit
        dphi_dtau0 = np.sqrt(M / r0**3)
        
        # Define Schwarzschild metric
        metric = schwarzschild_metric().subs(sp.symbols('M'), M)
        coords = sp.symbols('t r theta phi')
        parameters = {sp.symbols('M'): M}
        christoffel = christoffel_symbols(metric, coords, parameters)
        
        # Define initial conditions
        initial_conditions = [0, r0, theta0, phi0, 0, dr_dtau0, dtheta_dtau0, dphi_dtau0]
        t_span = (0, 100)
        t_eval = np.linspace(0, 100, 1000)
        
        # Geodesic equations and solve
        geodesic_func = geodesic_equations(christoffel, coords, parameters)
        solution = solve_geodesic(geodesic_func, initial_conditions, t_span, t_eval)
        
        # Extract solution
        t = solution.t
        r = solution.y[1]
        theta = solution.y[2]
        phi = solution.y[3]
        
        # Convert to Cartesian coordinates for plotting
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)
        
        # Plot the geodesic path
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x, y, z, label='Geodesic Path')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Particle Path around Schwarzschild Black Hole')
        ax.legend()
        
        plt.show()
        
    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numerical values.")

# GUI setup
root = tk.Tk()
root.title("Geodesic Path Simulator")

# Input fields
ttk.Label(root, text="Mass (M):").grid(column=0, row=0, padx=10, pady=5)
mass_entry = ttk.Entry(root)
mass_entry.grid(column=1, row=0, padx=10, pady=5)
mass_entry.insert(0, "1.0")

ttk.Label(root, text="Initial r:").grid(column=0, row=1, padx=10, pady=5)
r_entry = ttk.Entry(root)
r_entry.grid(column=1, row=1, padx=10, pady=5)
r_entry.insert(0, "10.0")

ttk.Label(root, text="Initial θ (radians):").grid(column=0, row=2, padx=10, pady=5)
theta_entry = ttk.Entry(root)
theta_entry.grid(column=1, row=2, padx=10, pady=5)
theta_entry.insert(0, f"{np.pi / 2}")

ttk.Label(root, text="Initial φ (radians):").grid(column=0, row=3, padx=10, pady=5)
phi_entry = ttk.Entry(root)
phi_entry.grid(column=1, row=3, padx=10, pady=5)
phi_entry.insert(0, "0.0")

ttk.Label(root, text="Initial dr/dτ:").grid(column=0, row=4, padx=10, pady=5)
dr_dtau_entry = ttk.Entry(root)
dr_dtau_entry.grid(column=1, row=4, padx=10, pady=5)
dr_dtau_entry.insert(0, "0.0")

ttk.Label(root, text="Initial dθ/dτ:").grid(column=0, row=5, padx=10, pady=5)
dtheta_dtau_entry = ttk.Entry(root)
dtheta_dtau_entry.grid(column=1, row=5, padx=10, pady=5)
dtheta_dtau_entry.insert(0, "0.0")

run_button = ttk.Button(root, text="Run Simulation", command=run_simulation)
run_button.grid(column=0, row=6, columnspan=2, pady=10)

root.mainloop()
