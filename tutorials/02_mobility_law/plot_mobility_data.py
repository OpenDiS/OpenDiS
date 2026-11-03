"""
Script to plot the mobility data used for the tutorial example.
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import json


def load_mobility_data(mobility_file):
    """Load the mobility data.
    
    Parameters
    ----------
    mobility_file : str
        Path to the json file containing the mobility data.
    """

    def nested_numpy_decoder(dct):
        """Convert JSON objects marked as NumPy arrays back into arrays."""
        if "__ndarray__" in dct:
            return np.array(dct["__ndarray__"])
        return dct
    with open(mobility_file, "r", encoding="utf-8") as f:
        data = json.load(f, object_hook=nested_numpy_decoder)

    return data


def plot_mobility_data(data):
    """Plot edge and screw stress-velocity data by temperature.

    Parameters
    ----------
    data : dict
        Nested mapping with ``"edge"`` and ``"screw"`` sections. Temperature
        entries contain two-column stress-velocity arrays.
    """

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    
    # Edge mobility
    data_edge = data["edge"]
    for tempK, shear_vel in ((k, v) for k, v in data_edge.items() if k.endswith("K")):
        ax[0].plot(shear_vel[:,0], shear_vel[:,1], 'o-', label=tempK)
    ax[0].set_xlabel('Shear stress (MPa)')
    ax[0].set_ylabel('Velocity (m/s)')
    ax[0].legend()

    # Screw mobility
    data_screw = data["screw"]
    for tempK, shear_vel in ((k, v) for k, v in data_screw.items() if k.endswith("K")):
        ax[1].plot(shear_vel[:,0], shear_vel[:,1], 'o-', label=tempK)
    ax[1].set_xlabel('Shear stress (MPa)')
    ax[1].set_ylabel('Velocity (m/s)')
    ax[1].legend()

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Load mobility data
    data = load_mobility_data("mobility_data.json")
    # Plot mobility data
    plot_mobility_data(data)
