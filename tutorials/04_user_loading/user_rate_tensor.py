import os, sys
import numpy as np

# Import pyexadis
pyexadis_path = '../../core/exadis/python/'
if not pyexadis_path in sys.path: sys.path.append(pyexadis_path)
try:
    import pyexadis
    from pyexadis_base import ExaDisNet, DisNetManager, SimulateNetwork
except ImportError:
    raise ImportError('Cannot import pyexadis')


class My_SimulateNetwork(SimulateNetwork):
    """ User-defined simulation driver to impose a strain-rate tensor loading
    """
    
    def __init__(self, *args, **kwargs) -> None:
        super(My_SimulateNetwork, self).__init__(*args, **kwargs)
        
        self.strain_rate_tensor = kwargs.get("strain_rate_tensor")
        state = kwargs.get("state")
        self.MU, self.NU = state["mu"], state["nu"]
        self.LA = 2*self.MU*self.NU/(1-2*self.NU)
    
    # Override step_update_response() function to apply strain-rate tensor loading
    def step_update_response(self, N: DisNetManager, state: dict):
        """step_update_response: update applied stress and rotation if needed
        """
        if self.loading_mode == 'strain_rate_tensor':
            
            # get values of plastic strain, plastic spin, and density computed internally in exadis
            dEp, dWp, state["density"] = N.get_disnet(ExaDisNet).net.get_plastic_strain()
            dEp = np.array(dEp).ravel()[[0,4,8,5,2,1]] # xx,yy,zz,yz,xz,xy
            dWp = np.array(dWp).ravel()[[5,2,1]] # yz,xz,xy
            state["dEp"] = dEp
            state["dWp"] = dWp
            
            # update strain and stress states based on strain rate tensor
            dE = self.strain_rate_tensor * state["dt"]
            dEe = dE - dEp # elastic strain
            dstress = self.LA*np.sum(dEe[0:3])*np.array([1,1,1,0,0,0]) + 2*self.MU*dEe
            
            # increment stress and strain tensors
            state["applied_stress"] += dstress
            state["Etot"] += dE
            
            def von_mises(T):
                S = np.array(T[[0,5,4,5,1,3,4,3,2]]).reshape(3,3)
                Sdev = S - np.trace(S)/3.0*np.eye(3)
                return np.sqrt(3.0/2.0*np.dot(Sdev.ravel(), Sdev.ravel()))
            
            # store strain and stress values used for the output (e.g. von Mises)
            state["strain"] = von_mises(state["Etot"])
            state["stress"] = von_mises(state["applied_stress"])
            
        else:
            # call base class function
            super().step_update_response(N, state)
            
        return state
