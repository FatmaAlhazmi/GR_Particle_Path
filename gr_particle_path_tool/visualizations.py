# gr_particle_path_tool/visualizations.py
import matplotlib.pyplot as plt

def plot_geodesic(solution, coords):
    plt.figure(figsize=(10, 6))
    for i, coord in enumerate(coords):
        plt.plot(solution.t, solution.y[i], label=f'{coord}')
    plt.xlabel('Proper time')
    plt.ylabel('Coordinate values')
    plt.legend()
    plt.title('Geodesic Path')
    plt.grid(True)
    plt.show()

