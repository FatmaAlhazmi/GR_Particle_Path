import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import sympy as sp
from sympy import symbols, Derivative
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from gr_particle_path_tool.metrics import (
    schwarzschild_metric, kerr_metric, reissner_nordstrom_metric, kerr_newman_metric,
    flrw_metric, de_sitter_metric, anti_de_sitter_metric, minkowski_metric,
    vaidya_metric, gottingen_metric, bertotti_robinson_metric
)
from gr_particle_path_tool.utils import christoffel_symbols, geodesic_equations, solve_geodesic
from gr_particle_path_tool.event_horizon import detect_event_horizon

class GRParticlePathTool:
    def __init__(self, root):
        self.root = root
        self.root.title("GR Particle Path Research Tool")
        
        # Use ttk for themed widgets
        self.style = ttk.Style()
        self.style.configure("TFrame", background="lightpink")
        self.style.configure("TLabel", background="lightpink", font=("Arial", 10))
        self.style.configure("TButton", font=("Arial", 10))
        self.style.configure("TRadiobutton", background="lightpink", font=("Arial", 10))

        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky="nsew")

        self.metric_var = tk.StringVar(value="Schwarzschild")
        metrics = [
            "Schwarzschild", "Kerr", "Reissner-Nordström", "Kerr-Newman",
            "FLRW", "de Sitter", "Anti-de Sitter", "Minkowski",
            "Vaidya", "Göttingen", "Bertotti-Robinson"
        ]

        metric_frame = ttk.LabelFrame(main_frame, text="Metric Selection", padding="10")
        metric_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nw")

        for i, metric in enumerate(metrics):
            ttk.Radiobutton(metric_frame, text=metric, variable=self.metric_var, value=metric).grid(row=i, column=0, sticky="w")

        param_frame = ttk.LabelFrame(main_frame, text="Parameters and Initial Conditions", padding="10")
        param_frame.grid(row=0, column=1, padx=10, pady=10, sticky="ne")

        self.mass_entry = self.add_parameter_entry(param_frame, "Mass (M):", "1.0", 0)
        self.a_entry = self.add_parameter_entry(param_frame, "Spin (a):", "0.0", 1)
        self.charge_entry = self.add_parameter_entry(param_frame, "Charge (Q):", "0.0", 2)
        self.r_entry = self.add_parameter_entry(param_frame, "Initial r:", "10.0", 3)
        self.theta_entry = self.add_parameter_entry(param_frame, "Initial θ (radians):", f"{np.pi / 2}", 4)
        self.phi_entry = self.add_parameter_entry(param_frame, "Initial φ (radians):", "0.0", 5)
        self.dr_dtau_entry = self.add_parameter_entry(param_frame, "Initial dr/dτ:", "0.0", 6)
        self.dtheta_dtau_entry = self.add_parameter_entry(param_frame, "Initial dθ/dτ:", "0.0", 7)
        self.dphi_dtau_entry = self.add_parameter_entry(param_frame, "Initial dφ/dτ:", "0.05", 8)

        action_frame = ttk.Frame(main_frame, padding="10")
        action_frame.grid(row=1, column=0, columnspan=2, pady=10)

        ttk.Button(action_frame, text="Run Simulation", command=self.run_simulation).grid(row=0, column=0, padx=10)
        ttk.Button(action_frame, text="Check Stability", command=self.check_stability).grid(row=0, column=1, padx=10)
        ttk.Button(action_frame, text="Parameter Sweep", command=self.parameter_sweep).grid(row=0, column=2, padx=10)
        ttk.Button(action_frame, text="Save Simulation", command=self.save_simulation).grid(row=0, column=3, padx=10)
        ttk.Button(action_frame, text="Load Simulation", command=self.load_simulation).grid(row=0, column=4, padx=10)
        ttk.Button(action_frame, text="Export Data", command=self.export_data).grid(row=0, column=5, padx=10)
        ttk.Button(action_frame, text="Event Horizon Detection", command=self.detect_event_horizon).grid(row=1, column=0, padx=10)
        ttk.Button(action_frame, text="Comparative Analysis", command=self.comparative_analysis).grid(row=1, column=1, padx=10)
        ttk.Button(action_frame, text="Gravitational Wave Signals", command=self.gravitational_wave_signals).grid(row=1, column=2, padx=10)
        ttk.Button(action_frame, text="User-defined Metrics", command=self.user_defined_metrics).grid(row=1, column=3, padx=10)
        ttk.Button(action_frame, text="Advanced Visualization", command=self.advanced_visualization).grid(row=1, column=4, padx=10)
        ttk.Button(action_frame, text="Batch Processing", command=self.batch_processing).grid(row=1, column=5, padx=10)
        ttk.Button(action_frame, text="Connect to Astronomy Databases", command=self.connect_astronomy_db).grid(row=2, column=0, padx=10)
        ttk.Button(action_frame, text="Machine Learning Integration", command=self.ml_integration).grid(row=2, column=1, padx=10)
        ttk.Button(action_frame, text="Detailed Documentation", command=self.show_documentation).grid(row=2, column=2, padx=10)
        ttk.Button(action_frame, text="Customizable Interface", command=self.customize_interface).grid(row=2, column=3, padx=10)
        ttk.Button(action_frame, text="Simulation Playback", command=self.simulation_playback).grid(row=2, column=4, padx=10)

    def add_parameter_entry(self, frame, label_text, default_value, row):
        ttk.Label(frame, text=label_text).grid(column=0, row=row, sticky="w")
        entry = ttk.Entry(frame)
        entry.grid(column=1, row=row, padx=10, pady=5, sticky="w")
        entry.insert(0, default_value)
        return entry

    def get_metric(self):
        metric_name = self.metric_var.get()
        if metric_name == "Schwarzschild":
            return schwarzschild_metric
        elif metric_name == "Kerr":
            return kerr_metric
        elif metric_name == "Reissner-Nordström":
            return reissner_nordstrom_metric
        elif metric_name == "Kerr-Newman":
            return kerr_newman_metric
        elif metric_name == "FLRW":
            return flrw_metric
        elif metric_name == "de Sitter":
            return de_sitter_metric
        elif metric_name == "Anti-de Sitter":
            return anti_de_sitter_metric
        elif metric_name == "Minkowski":
            return minkowski_metric
        elif metric_name == "Vaidya":
            return vaidya_metric
        elif metric_name == "Göttingen":
            return gottingen_metric
        elif metric_name == "Bertotti-Robinson":
            return bertotti_robinson_metric
        else:
            raise ValueError("Unknown metric")

    def run_simulation(self):
        try:
            M = float(self.mass_entry.get())
            r0 = float(self.r_entry.get())
            θ0 = float(self.theta_entry.get())
            dot_r0 = float(self.dr_dtau_entry.get())
            dot_theta0 = float(self.dtheta_dtau_entry.get())
            h = 0.01  # Integration step size

            # Calculate constants L and E
            L_value = np.sqrt(M) * r0 / np.sqrt(r0 - 3 * M)
            E_value = np.sqrt(dot_r0**2 + ((r0**2 + r0**4 * dot_theta0**2 + L_value**2 * (1 / np.sin(θ0))**2) * (1 - 2 * M / r0)) / r0**2)

            # Initial values for [ρ, ρ_dot, θ, θ_dot]
            y0 = np.array([r0, dot_r0, θ0, dot_theta0])
            t_span = (0, 200)
            dt = 0.01
            t_values = np.arange(t_span[0], t_span[1], dt)

            # Solving using RK4
            sol = np.zeros((len(t_values), len(y0)))
            sol[0] = y0

            for i in range(1, len(t_values)):
                sol[i] = self.rk4_step(self.system_of_eqs, sol[i-1], t_values[i-1], dt, L_value)

            # Interpolation functions
            rho_s = lambda tau: np.interp(tau, t_values, sol[:, 0])
            theta_s = lambda tau: np.interp(tau, t_values, sol[:, 2])

            # Equation for φ
            sol_phi = np.zeros(len(t_values))
            sol_phi[0] = 0  # Initial condition for φ

            for i in range(1, len(t_values)):
                sol_phi[i] = self.rk4_step(self.phi_eqs, sol_phi[i-1], t_values[i-1], dt, rho_s, theta_s, L_value)

            # Interpolation function for φ
            phi_s = lambda tau: np.interp(tau, t_values, sol_phi)

            # Plotting the 3D orbit
            rho_values = np.array([rho_s(tau) for tau in t_values])
            theta_values = np.array([theta_s(tau) for tau in t_values])
            phi_values = np.array([phi_s(tau) for tau in t_values])

            x_vals = rho_values * np.sin(theta_values) * np.cos(phi_values)
            y_vals = rho_values * np.sin(theta_values) * np.sin(phi_values)
            z_vals = rho_values * np.cos(theta_values)

            self.plot_geodesic(x_vals, y_vals, z_vals)

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def system_of_eqs(self, t, y, L):
        ρ, ρ_dot, θ, θ_dot = y

        csc_θ = 1 / np.sin(θ)
        cot_θ = np.cos(θ) / np.sin(θ)

        θ_ddot = (L**2 * cot_θ * csc_θ**2 - 2 * ρ**3 * θ_dot * ρ_dot) / ρ**4
        ρ_ddot = (L**2 * csc_θ**2 + ρ**4 * θ_dot**2 - 3 * L**2 * csc_θ**2 / ρ - ρ - 3 * ρ**3 * θ_dot**2) / ρ**3

        return np.array([ρ_dot, ρ_ddot, θ_dot, θ_ddot])

    def phi_eqs(self, t, y, rho_s, theta_s, L):
        φ = y
        ρ_val = rho_s(t)
        θ_val = theta_s(t)
        sin_θ_val = np.sin(θ_val)
        dφ_dt = L / (ρ_val**2 * sin_θ_val**2)
        return dφ_dt

    def rk4_step(self, func, y, t, dt, *args):
        k1 = dt * func(t, y, *args)
        k2 = dt * func(t + dt/2, y + k1/2, *args)
        k3 = dt * func(t + dt/2, y + k2/2, *args)
        k4 = dt * func(t + dt, y + k3, *args)
        return y + (k1 + 2*k2 + 2*k3 + k4) / 6

    def plot_geodesic(self, x_vals, y_vals, z_vals):
        fig = Figure(figsize=(5, 5), dpi=100)
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x_vals, y_vals, z_vals, label="Geodesic Path")
        ax.scatter([0], [0], [0], color='black', s=100)  # Black hole in the center

        max_range = np.array([x_vals.max() - x_vals.min(), y_vals.max() - y_vals.min(), z_vals.max() - z_vals.min()]).max() / 2.0
        mid_x = (x_vals.max() + x_vals.min()) * 0.5
        mid_y = (y_vals.max() + y_vals.min()) * 0.5
        mid_z = (z_vals.max() + z_vals.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('3D Orbit around a Schwarzschild Black Hole')
        ax.legend()

        canvas = FigureCanvasTkAgg(fig, master=self.root)
        canvas.draw()
        canvas.get_tk_widget().grid(row=3, column=0, columnspan=2, pady=10)

    def check_stability(self):
        messagebox.showinfo("Stability Check", "Stability check feature coming soon.")

    def parameter_sweep(self):
        messagebox.showinfo("Parameter Sweep", "Parameter sweep feature coming soon.")

    def save_simulation(self):
        filepath = filedialog.asksaveasfilename(defaultextension=".json")
        if filepath:
            with open(filepath, 'w') as f:
                f.write(self.solution.to_json())

    def load_simulation(self):
        filepath = filedialog.askopenfilename()
        if filepath:
            with open(filepath, 'r') as f:
                data = f.read()
                self.solution = pd.read_json(data)

    def export_data(self):
        filepath = filedialog.asksaveasfilename(defaultextension=".csv")
        if filepath:
            df = pd.DataFrame(self.solution.y.T, columns=['t', 'r', 'theta', 'phi', 'dr/dtau', 'dtheta/dtau', 'dphi/dtau'])
            df.to_csv(filepath, index=False)

    def detect_event_horizon(self):
        metric = self.metric_var.get().lower().replace("-", "_")
        try:
            horizon = detect_event_horizon(metric, float(self.mass_entry.get()), float(self.a_entry.get()), float(self.charge_entry.get()))
            messagebox.showinfo("Event Horizon Detection", f"Event Horizon: {horizon}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def comparative_analysis(self):
        messagebox.showinfo("Comparative Analysis", "Comparative analysis feature coming soon.")

    def gravitational_wave_signals(self):
        messagebox.showinfo("Gravitational Wave Signals", "Gravitational wave signals feature coming soon.")

    def user_defined_metrics(self):
        messagebox.showinfo("User-defined Metrics", "User-defined metrics feature coming soon.")

    def advanced_visualization(self):
        messagebox.showinfo("Advanced Visualization", "Advanced visualization feature coming soon.")

    def batch_processing(self):
        messagebox.showinfo("Batch Processing", "Batch processing feature coming soon.")

    def connect_astronomy_db(self):
        messagebox.showinfo("Connect to Astronomy Databases", "Connect to astronomy databases feature coming soon.")

    def ml_integration(self):
        messagebox.showinfo("Machine Learning Integration", "Machine learning integration feature coming soon.")

    def show_documentation(self):
        messagebox.showinfo("Detailed Documentation", "Detailed documentation feature coming soon.")

    def customize_interface(self):
        messagebox.showinfo("User Customizable Interface", "User customizable interface feature coming soon.")

    def simulation_playback(self):
        messagebox.showinfo("Simulation Playback Control", "Simulation playback control feature coming soon.")

if __name__ == "__main__":
    root = tk.Tk()
    app = GRParticlePathTool(root)
    root.mainloop()
