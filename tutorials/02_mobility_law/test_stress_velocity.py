"""
Example of a script to test mobility laws.
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
try:
    import pyexadis
    from pyexadis_base import ExaDisNet
    from pyexadis_base import CalForce, MobilityLaw
    from pyexadis_utils import insert_infinite_line
except ImportError:
    raise ImportError('Cannot import pyexadis')

from user_mob_bccglide import My_MobilityLaw

def init_line_config(Lbox, burg, plane, theta, maxseg):
    """init_line_config
    Initialize an infinite straight line along the y-direction
    with glide plane normal oriented along the z-direction.
    Returns the line configuration and appropriate crystal 
    orientation matrix.
    """
    # Crystal orientation
    plane = plane / np.linalg.norm(plane)
    b = burg / np.linalg.norm(burg)
    y = np.cross(plane, b)
    y = y / np.linalg.norm(y)
    linedir = np.cos(theta*np.pi/180.0)*b+np.sin(theta*np.pi/180.0)*y
    x = np.cross(linedir, plane)
    x = x / np.linalg.norm(x)
    Rorient = np.array([x, linedir, plane])
    
    burg = np.matmul(Rorient, burg)
    plane = np.matmul(Rorient, plane)
    linedir = np.matmul(Rorient, linedir)
    # linedir = -np.matmul(Rorient, linedir)
    
    # Create the cell and insert the line
    cell = pyexadis.Cell(Lbox)
    nodes, segs = [], []
    origin = np.array(cell.origin) + 0.5*np.sum(np.array(cell.h), axis=0)
    nodes, segs = insert_infinite_line(cell, nodes, segs, burg, plane, origin, linedir=linedir, maxseg=maxseg)
    
    N = DisNetManager(ExaDisNet(cell, nodes, segs))
    return N, Rorient


def stress_velocity(
    crystal,
    burgmag,
    shear,
    theta,
    tempK=300.0,
    mobility_law='FCC_0',
    mobparams=None
):
    """Compute a straight dislocation's velocity for a range shear stresses.

    The calculation is repeated for every value in either ``theta`` or
    ``tempK`` (one of them must be a list). For each case, an infinite straight
    dislocation is initialized, the applied force and nodal velocities are
    evaluated, and the mean glide speed is recorded for each shear value.

    Parameters
    ----------
    crystal : str
        Crystal structure, typically ``"bcc"`` or ``"fcc"``.
    burgmag : float
        Burgers vector magnitude, in meters.
    shear : numpy.ndarray
        Shear stresses to evaluate, in MPa.
    theta : float or list of float
        Dislocation character angle or angles, in degrees.
    tempK : float or list of float, optional
        Temperature or list of temperatures, in kelvin.
    mobility_law : str or type, optional
        ExaDiS mobility-law name or a ``MobilityLaw_Base`` subclass.
    mobparams : dict, optional
        Parameters passed to the mobility-law.

    Returns
    -------
    list
        One stress-velocity data per angle or temperature. Each data set
        contains ``[shear_stress, average_velocity]`` pairs.
    """
    theta_is_list = isinstance(theta, list)
    temp_is_list = isinstance(tempK, list)
    if theta_is_list == temp_is_list:
        raise ValueError("Exactly one of 'theta' or 'tempK' must be a list")

    base_mobparams = dict(mobparams or {})
    all_results = []

    Lbox = 1000.0
    maxseg = 0.1*Lbox

    state = {
        "crystal": crystal,
        "burgmag": burgmag,
        "mu": 50.0e9,
        "nu": 0.3,
        "a": 1.0,
        "maxseg": maxseg,
        "minseg": 0.2*maxseg,
        "rann": 0.5,
        "rtol": 0.25,
    }

    calforce = CalForce(force_mode='LINE_TENSION_MODEL', state=state)

    if temp_is_list:
        target_name = 'tempK'
        target = tempK
    else:
        target_name = 'theta'
        target = theta

    for val in target:
        theta_val = val if target_name == 'theta' else theta
        tempK_val = val if target_name == 'tempK' else tempK
        current_mobparams = base_mobparams.copy()

        if crystal == 'bcc':
            # 1/2<111>{110} dislocation
            burg = 1.0/np.sqrt(3.0)*np.array([1.,1.,1.])
            plane = np.array([1.,-1.,0.])
        else:
            # 1/2<110>{111} dislocation
            burg = 1.0/np.sqrt(2.0)*np.array([1.,1.,0.])
            plane = np.array([1.,-1.,1.])

        # Dislocation line
        N, Rorient = init_line_config(Lbox, burg, plane, theta_val, maxseg)
        state["Rorient"] = Rorient

        def is_python_mobility(mobility_law) -> bool:
            try:
                return issubclass(mobility_law, MobilityLaw_Base)
            except TypeError:
                return False

        if is_python_mobility(mobility_law):
            # Python-based mobility laws
            mobility = mobility_law(state=state, **current_mobparams)
        else:
            # exadis mobility laws
            if mobility_law in ['BCC_0B', 'FCC_0', 'SimpleGlide']:
                mobility = MobilityLaw(mobility_law=mobility_law, state=state, **current_mobparams)
            elif mobility_law in ['BCC_0C', 'BCC_NL']:
                current_mobparams["tempK"] = tempK_val
                mobility = MobilityLaw(mobility_law=mobility_law, state=state, **current_mobparams)
            else:
                raise ValueError(f"Unsupported mobility {mobility_law}")

        results = []
        for s in range(shear.size):
            # Applied stress
            applied_stress = np.zeros(6)
            applied_stress[4] = np.sin(theta_val*np.pi/180.0)*shear[s]*1e6 # xz
            applied_stress[3] = np.cos(theta_val*np.pi/180.0)*shear[s]*1e6 # yz
            state["applied_stress"] = applied_stress
            
            # Compute velocity
            state = calforce.NodeForce(N, state)
            state = mobility.Mobility(N, state)
            
            # Save average velocity
            vel = state["nodevels"]
            #avg_vel = np.mean(vel[:,0])*state["burgmag"]
            avg_vel = np.mean(np.sqrt(vel[:,0]**2+vel[:,2]**2))*state["burgmag"]
            #avg_vel = np.mean(np.linalg.norm(vel, axis=1))*state["burgmag"]
            results.append([shear[s], avg_vel])

        all_results.append(results)

    return all_results


def test_stress_velocity():
    """
    Test the implementation of mobility laws by probing the stress-velocity
    relation. This is done by computing the velocity of an infinite straight
    dislocation line as a function of shear stress for different character angles.
    """
    pyexadis.initialize()

    burgmag = 2.85e-10

    # Mobility law
    #mobility_law = 'SimpleGlide'
    #mobility_law = 'FCC_0'
    mobility_law = 'BCC_0B'
    #mobility_law = 'BCC_NL'
    #mobility_law = 'BCC_0C'
    #mobility_law = My_MobilityLaw
    
    # Mobility law parameters
    if mobility_law == 'SimpleGlide':
        crystal = ''
        mobparams = {"Medge": 2000.0, "Mscrew": 1000.0}
        linear_mob = True
    elif mobility_law == 'FCC_0':
        crystal = 'fcc'
        mobparams = {"Medge": 60000.0, "Mscrew": 40000.0}
        linear_mob = True
    elif mobility_law == 'BCC_0B':
        crystal = 'bcc'
        mobparams = {"Medge": 2600.0, "Mscrew": 20.0, "Mclimb": 1e-4, "vmax": 3400.0}
        linear_mob = True
    elif mobility_law == 'BCC_NL':
        crystal = 'bcc'
        mobparams = {
            "Peierls": 1.2e9, "Bscrew": 4.6e-4, "Rscrew": 1.0, "PeierlsTATasym": 0.0,
            "B0edge": 0.0, "B1edge": 7.7e-7, "Bclimb": 3.0e-1, "tempK": 300.0, "vmax": -1.0
        }
        linear_mob = False
    elif mobility_law == 'BCC_0C':
        crystal = 'bcc'
        mobparams = {
            "Medge": 2000.0, "Mscrew": 2000.0, "Mclimb": 1e-4, "FscrewTATasym": 0.0,
            "tau0": 0.5e9, "p": 0.75, "q": 1.1, "H0": 1.06, "T0": -1.0, "tauv": 0.0,
            "B1screw": 0.0, "B1edge": 0.0, "Dedge": 0.0, "D1edge": 0.0, "v0edge": -1.0,
            "vmax": -1.0, "tempK": 300.0
        }
        linear_mob = False
    elif mobility_law == My_MobilityLaw:
        crystal = 'bcc'
        mobparams = {"mob": 1000.0, "mob_edge": 1000.0, "mob_screw": 100.0}
        linear_mob = True
    else:
        raise ValueError('invalid mobility law %s' % mobility_law)
        
    theta = [0.0, 30.0, 60.0, 90.0] # character angles in degrees
    shear = np.linspace(0.0, 1000.0, 20+1) # shear stresses MPa

    vel_results = stress_velocity(
        crystal, burgmag,
        shear, theta,
        mobility_law=mobility_law,
        mobparams=mobparams
    )
    
    # Plot results
    add_linear_fit = linear_mob
    
    fig, ax = plt.subplots()
    
    for i in range(len(vel_results)):
        results = np.array(vel_results[i])
        ax.plot(results[:,0], results[:,1], 'o' if add_linear_fit else 'o-', label=r'$\theta$ = %.1f$^o$'%theta[i])
        
        if add_linear_fit:
            pfit = np.polyfit(results[:,0], results[:,1], 1)
            Mfit = pfit[0]/1e6/burgmag
            ax.plot(results[:,0], np.polyval(pfit, results[:,0]), '-', color='C%d'%i, label=r'fit M = %.2f (Pa.s)$^{-1}$'%Mfit)
    
    ax.set_xlabel('Shear stress (MPa)')
    ax.set_ylabel('Velocity (m/s)')
    ax.legend()
    ax.set_title(f'Mobility law: {mobility_law}')
    plt.show()

    if not sys.flags.interactive:
        pyexadis.finalize()


if __name__ == "__main__":
    test_stress_velocity()
