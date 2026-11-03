"""
Example of a script to fit mobility data to existing mobility laws.
"""

import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple

# Import opendis modules
opendis_paths = ['../../python', '../../lib', '../../core/pydis/python', '../../core/exadis/python/']
[sys.path.append(os.path.abspath(path)) for path in opendis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from framework.mobility_base import MobilityLaw_Base
from pydis import MobilityLaw as MobilityLaw_PyDiS
try:
    import pyexadis
    from pyexadis_base import ExaDisNet, NodeConstraints, SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration, Collision, Topology, Remesh
except ImportError:
    raise ImportError('Cannot import pyexadis')

from plot_mobility_data import load_mobility_data, plot_mobility_data
from test_stress_velocity import stress_velocity

from dataclasses import dataclass
from lm_fitter import Parameters


def init_loop_config(Lbox, burg, plane, radius, N=20):
    """Create a circular glissile dislocation loop in a cubic cell.

    The crystal is rotated so that the Burgers vector is aligned with the
    simulation y-axis and the glide plane normal with the z-axis.

    Parameters
    ----------
    Lbox : float
        Side length of the cubic simulation cell.
    burg : numpy.ndarray
        Burgers vector.
    plane : numpy.ndarray
        Glide plane normal.
    radius : float
        Radius of the circular loop.
    N : int, optional
        Number of nodes used to discretize the loop.

    Returns
    -------
    tuple
        The DisNetManager object and the crystal orientation matrix.
    """
    print("init_circular_loop: radius = %f, N = %d" % (radius, N))
    plane = plane / np.linalg.norm(plane)
    b = burg / np.linalg.norm(burg)
    y = np.cross(plane, b)
    y = y / np.linalg.norm(y)
    linedir = b
    x = np.cross(linedir, plane)
    x = x / np.linalg.norm(x)
    Rorient = np.array([x, linedir, plane])

    burg = np.matmul(Rorient, burg)
    plane = np.matmul(Rorient, plane)
    linedir = np.matmul(Rorient, linedir)
    # linedir = -np.matmul(Rorient, linedir)

    theta = np.arange(N)*2.0*np.pi/N
    nodes = np.vstack([radius*np.cos(theta)+Lbox/2, radius*np.sin(theta)+Lbox/2, np.zeros_like(theta)+Lbox/2]).T
    segs  = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg, nodes[(i+1)%N,:3]-nodes[i,:3])
        countzero = not np.any(pn)
        if countzero:
            pn = plane
        pn = pn / np.linalg.norm(pn)
        segs[i,:] = np.concatenate(([i, (i+1)%N], burg, pn))
    
    cell = pyexadis.Cell(Lbox)
    N = DisNetManager(ExaDisNet(cell, nodes, segs))
    return N, Rorient


@dataclass
class MobilityParams:
    mobility_law: str
    crystal: str
    parameters: Parameters


"""Mobility-law profiles available to the fitting.

Each key is a supported mobility-law name and maps to a ``MobilityParams``
profile containing the ExaDiS law name, its crystal structure, and the initial
mobility parameter values with their bounds. Parameter flags identify screw, 
edge, and temperature-dependent terms so that ``get_mobility_params`` can 
select which parameters should vary for a given dataset.
"""
MOBILITY_PARAMS = {
    "FCC_0": MobilityParams(
        mobility_law="FCC_0",
        crystal="fcc",
        parameters=(
            Parameters()
            # Screw fitting parameters
            .add('Mscrew', value=10000.0, min=0.0, flags=("screw")) # Screw mobility, 1/(Pa.s)
            # Edge fitting parameters
            .add('Medge', value=10000.0, min=0.0, flags=("edge")) # Edge mobility, 1/(Pa.s)
            # Other parameters
            .add('vmax', value=3400.0, min=0.0, max=1e4) # Dislocation velocity capping, m/s
        )
    ),
    "FCC_0_FRIC": MobilityParams(
        mobility_law="FCC_0_FRIC",
        crystal="fcc",
        parameters=(
            Parameters()
            # Screw fitting parameters
            .add('Mscrew', value=10000.0, min=0.0, flags=("screw")) # Screw mobility, 1/(Pa.s)
            .add('Fscrew', value=0.1e9, min=0.0, max=10e9, flags=("screw")) # Screw friction stress, Pa
            # Edge fitting parameters
            .add('Medge', value=10000.0, min=0.0, flags=("edge")) # Edge mobility, 1/(Pa.s)
            .add('Fedge', value=0.0, min=0.0, max=1e9, flags=("edge")) # Edge friction stress, Pa
            # Other parameters
            .add('vmax', value=3400.0, min=0.0, max=1e4) # Dislocation velocity capping, m/s
        )
    ),
    "BCC_0B": MobilityParams(
        mobility_law="BCC_0B",
        crystal="bcc",
        parameters=(
            Parameters()
            # Screw fitting parameters
            .add('Mscrew', value=20.0, min=0.0, flags=("screw")) # Screw mobility, 1/(Pa.s)
            .add('Fscrew', value=0.1e9, min=0.0, max=10e9, flags=("screw")) # Screw friction stress, Pa
            # Edge fitting parameters
            .add('Medge', value=10000.0, min=0.0, flags=("edge")) # Edge mobility, 1/(Pa.s)
            .add('Fedge', value=0.0, min=0.0, max=1e9, flags=("edge")) # Edge friction stress, Pa
            # Other parameters
            .add('Mclimb', value=1e-4, vary=False) # Climb mobility, 1/(Pa.s)
            .add('vmax', value=3400.0, min=0.0, max=1e4) # Dislocation velocity capping, m/s
        )
    ),
    "BCC_0C": MobilityParams(
        mobility_law="BCC_0C",
        crystal="bcc",
        parameters=(
            Parameters()
            # Screw fitting parameters
            .add('Mscrew', value=2600.0, min=1.0, flags=("screw")) # Linear screw mobility, 1/(Pa.s)
            .add('tau0', value=0.5e9, min=0.0, max=10e9, flags=("screw")) # Peierls stress, Pa
            .add('p', value=0.8, min=0.1, max=2.0, flags=("screw")) # Enthalpy exponent, unitless
            .add('q', value=1.2, min=0.1, max=2.0, flags=("screw")) # Enthalpy exponent, unitless
            .add('H0', value=1.0, min=1e-2, max=10.0, flags=("screw")) # Enthalpy pre-factor, eV
            .add('T0', value=1e10, vary=False, flags=("screw")) # Transition temperature, K
            .add('tauv', value=0.0, min=-5e9, max=5e9, flags=("screw")) # Shift stress, Pa
            .add('B1screw', value=0.0, min=-1e10, max=1e10, flags=("screw", "temperature")) # Screw temperature dependant drag
            # Edge fitting parameters
            .add('Medge', value=20000.0, min=1.0, flags=("edge")) # Linear edge mobility, 1/(Pa.s)
            .add('B1edge', value=0.0, flags=("edge", "temperature")) # Edge temperature dependant drag
            .add('Dedge', value=2e-10, flags=("edge")) # Edge non-linear drag
            .add('D1edge', value=1e-13, flags=("edge", "temperature")) # Edge temperature dependant non-linear drag
            .add('v0edge', value=1500.0, min=0.0, flags=("edge")) # Edge non-linear velocity, m/s
            # Other parameters
            .add('Mclimb', value=1e-4, vary=False) # Climb mobility, 1/(Pa.s)
            .add('vmax', value=3000, min=1.0, max=1e4) # Dislocation velocity capping, m/s
        )
    ),
    "BCC_NL": MobilityParams(
        mobility_law="BCC_NL",
        crystal="bcc",
        parameters=(
            Parameters()
            # Screw fitting parameters
            .add('Peierls', value=0.5e9, min=0.0, max=10e9, flags=("screw")) # Peierls stress, Pa
            .add('Bscrew', value=1.0/2600.0, min=1e-6, max=1.0, flags=("screw")) # Linear screw drag, Pa.s
            .add('Rscrew', value=1.0, flags=("screw")) # unitless
            .add('PeierlsTATasym', value=0.0, vary=False, flags=("screw")) # T/AT coefficient, unitless
            # Edge fitting parameters
            .add('B0edge', value=1.0/2600.0, min=1e-6, max=1.0, flags=("edge")) # Linear edge drag, Pa.s
            .add('B1edge', value=7.7e-7, flags=("edge", "temperature")) # Edge temperature dependant drag
            # Other parameters
            .add('Bclimb', value=3e-1, vary=False) # Climb drag, Pa.s
            .add('vmax', value=3000, min=1.0, max=1e4) # Dislocation velocity capping, m/s
        )
    ),
}


def get_mobility_params(mobility_law, dislocation_type, dataset):
    """Prepare the parameters that will be varied during mobility fitting.

    Parameters associated with the other dislocation character are fixed. If
    the data contain only one temperature, temperature-dependent parameters
    are fixed as well.

    Parameters
    ----------
    mobility_law : str
        Name of a mobility law defined in ``MOBILITY_PARAMS``.
    dislocation_type : {"edge", "screw", "all"}
        Dislocation character represented by the data being fitted.
    dataset : sequence of tuple
        ``(theta, tempK, data)`` entries, where ``data`` contains stress and
        velocity columns.

    Returns
    -------
    tuple
        Crystal structure name and an independent ``Parameters`` instance
        configured for the requested fit.
    """

    if mobility_law not in MOBILITY_PARAMS:
        raise KeyError(f"Unknown mobility {mobility_law}")

    profile = MOBILITY_PARAMS[mobility_law]
    crystal = profile.crystal
    params = profile.parameters.copy()

    # Adjust parameters based on dislocation type
    # If we are fitting edge dislocations, then we should fix
    # all parameters that are screw related, and vice-versa
    for name, param in params.parameters.items():
        if (dislocation_type == "edge" and "screw" in param.flags) or \
            (dislocation_type == "screw" and "edge" in param.flags):
            params.fix(name)

    # Adjust parameters based on temperature fitting
    # (if we have multiple temperatures in the dataset)
    temp_fitting = (len(set([T for _, T, _ in dataset])) > 1)
    for name, param in params.parameters.items():
        if "temperature" in param.flags and not temp_fitting:
            params.fix(name)

    return crystal, params


def mobility_residual(params, crystal, burgmag, mobility_law, dataset):
    """Calculate velocity residuals for the non-linear mobility fit.

    Parameters
    ----------
    params : Parameters
        Current mobility-law parameters supplied by the fitter.
    crystal : str
        Crystal structure used by the mobility law.
    burgmag : float
        Burgers vector magnitude, in meters.
    mobility_law : str
        Name of the mobility law to evaluate.
    dataset : sequence of tuple
        ``(theta, tempK, data)`` entries. Each data array contains shear stress
        in MPa and measured velocity in m/s.

    Returns
    -------
    numpy.ndarray
        Concatenated velocity residuals (errors).
    """
    # Evaluate nodal velocities and corresponding residuals
    resids = []
    for theta, tempK, data in dataset:

        shear = np.linspace(0.0, data[:,0].max(), 100)
        vels = stress_velocity(
            crystal, burgmag,
            shear, theta, [tempK],
            mobility_law=mobility_law,
            mobparams=params.valuesdict(),
        )
        vels = np.array(vels).squeeze()

        vel_interp = np.interp(data[:,0], vels[:,0], vels[:,1])
        resid = vel_interp - data[:,1]
        resids.append(resid)

    return np.concatenate(resids)


def plot_mobility_fit(crystal, burgmag, mobility_law, dataset, fitparams):
    """Plot measured mobility data together with the fitted model curves.

    Parameters
    ----------
    crystal : str
        Crystal structure used by the mobility law.
    burgmag : float
        Burgers vector magnitude in meters.
    mobility_law : str
        Name of the fitted mobility law.
    dataset : sequence of tuple
        ``(theta, tempK, data)`` entries to plot.
    fitparams : dict
        Fitted parameters passed to the mobility law.
    """

    fig, ax = plt.subplots()
    for i, (theta, tempK, data) in enumerate(dataset):
        ax.plot(data[:,0], data[:,1], 'o', color=f"C{i}", label=f"{tempK}K")

        shear = np.linspace(0.0, data[:,0].max(), 100)
        vels = stress_velocity(
            crystal, burgmag,
            shear, theta, [tempK],
            mobility_law=mobility_law,
            mobparams=fitparams,
        )
        vels = np.array(vels).squeeze()
        ax.plot(vels[:,0], vels[:,1], color=f"C{i}")

    ax.set_xlabel('Shear stress (MPa)')
    ax.set_ylabel('Velocity (m/s)')
    ax.legend()
    plt.show()


def fit_MD_mobility_example():
    """Fit an ExaDiS mobility law to molecular-dynamics reference data.

    The example loads ``mobility_data.json``, selects the model parameters to
    vary, performs a nonlinear least-squares fit, and plots the resulting
    stress-velocity curves against the reference data.
    """
    pyexadis.initialize()

    """
    Load the mobility data
    """
    data = load_mobility_data("mobility_data.json")
    #plot_mobility_data(data)

    """
    Fit the mobility data
    """
    burgmag = 2.85e-10
    mobility_law = "BCC_0C"
    dislocation_type = "edge"
    theta = 0.0 if dislocation_type == "screw" else 90.0

    # Data to be fitted simultaneously
    # The dataset object must be a list of tuple (theta, tempK, stress_vel) entries,
    # where stress_vel contains stress and velocity columns.
    dataset = [
        (theta,  10.0, data[dislocation_type]["10K"]),
        (theta, 100.0, data[dislocation_type]["100K"]),
        (theta, 200.0, data[dislocation_type]["200K"]),
        (theta, 300.0, data[dislocation_type]["300K"]),
    ]

    # Get the mobility fitting parameters from pre-existing ExaDiS mobility laws
    crystal, params = get_mobility_params(mobility_law, dislocation_type, dataset)

    # One can make manual adjustments to the default mobility 
    # parameters to help the fitter get better fits to the data
    if mobility_law == "BCC_0C" and dislocation_type == "edge":
        # fit adjustements for BCC_0C edge
        params.set_value("vmax", 3000.0)
        params.fix("vmax")

    elif mobility_law == "BCC_0C" and dislocation_type == "screw":
        # fit adjustements for BCC_0C edge
        params.set_value('tau0', value=1e9)
        params.set_value('p', value=0.8)
        params.set_value('q', value=1.2)
        params.set_value('H0', value=1.0)
        params.set_value('tauv', value=-4000.0e6)
        params.set_value('B1screw', value=-3e-7)
        params.fix('tau0')
        params.fix('p')
        params.fix('q')
        params.fix('tau0')
        params.fix('B1screw')

    # Print parameters
    print(params)

    print('Fitting...')
    from lm_fitter import minimize
    fit = minimize(mobility_residual, params, args=(crystal, burgmag, mobility_law, dataset))
    fitparams = fit["params"]
    print('fitsuccess',fit["success"])
    print('RMSE',fit["rmse"])
    print('fitparams =',fit["params"])

    """
    Plot and check results
    """
    plot_mobility_fit(crystal, burgmag, mobility_law, dataset, fitparams)


    if not sys.flags.interactive:
        pyexadis.finalize()


def test_loop():
    """Run a dislocation-loop simulation using example fitted parameters.

    This demonstration initializes a circular BCC loop, configures the force,
    mobility, integration, collision, topology, and remeshing models, and then
    evolves and visualizes the network under an applied stress.
    """
    
    pyexadis.initialize()

    Lbox = 1000.0
    burg = 1.0/np.sqrt(3.0)*np.array([1,1,1])
    plane = 1.0/np.sqrt(2.0)*np.array([-1,1,0])
    N, Rorient = init_loop_config(Lbox, burg, plane, radius=0.1*Lbox, N=20)

    vis = VisualizeNetwork(view_init={'elev': -90, 'azim': 90})
    
    state = {
        "crystal": 'bcc',
        "burgmag": 2.85e-10,
        "mu": 50e9,
        "nu": 0.3,
        "a": 1.0,
        "maxseg": 0.05*Lbox,
        "minseg": 0.01*Lbox,
        "rann": 0.5,
        "maxdt": 1e-12,
    }
    state["Rorient"] = Rorient
    
    calforce  = CalForce(force_mode='LineTension', state=state)

    fitparams = {
        'Mscrew': 2242.2894859343523, 'tau0': 1000000000.0, 'p': 0.8, 'q': 1.2, 'H0': 1.0398232052433396,
        'T0': 10000000000.0, 'tauv': -4613328317.887041, 'B1screw': -3e-07, 'Medge': 205928.62042551287,
        'B1edge': 1.6217730093630888e-07, 'Dedge': 6.464330227652489e-10, 'D1edge': 9.670353453850635e-13,
        'v0edge': 449.5897059965622, 'Mclimb': 0.0001, 'vmax': 3000.0
    }
    mobility  = MobilityLaw(mobility_law='BCC_0C', state=state, **fitparams)

    timeint   = TimeIntegration(integrator='Trapezoid', force=calforce, mobility=mobility, state=state)
    collision = Collision(collision_mode='Retroactive', state=state)
    topology  = Topology(split_mode='TopologySerial', state=state, force=calforce, mobility=mobility)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)
    
    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint, 
                          collision=collision, topology=topology, remesh=remesh, vis=vis,
                          max_step=1000, loading_mode='stress', state=state,
                          applied_stress=np.array([0.0, 0.0, 0.0, -1e9, 0.0, 0.0]),
                          print_freq=10, plot_freq=10, plot_pause_seconds=0.0001,
                          write_freq=100, write_dir='test_loop_results')
    sim.run(N, state)
    
    pyexadis.finalize()


if __name__ == "__main__":
    fit_MD_mobility_example()
    #test_loop()
