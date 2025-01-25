"""@package docstring
DisNetManager: class for managing multiple implementations of dislocation network

Implements synchronization between different implementations of DisNet
"""
import numpy as np

class DisNetManager:
    """Class for managing multiple implementations of dislocation network

    Implements synchronization between different implementations of DisNet
    """
    def __init__(self, disnet=None):
        self.disnet_dict = {}
        if disnet is not None and not isinstance(disnet, type):
            self.disnet_dict[type(disnet)] = disnet
        else:
            raise ValueError("DisNetManager: user need to provide a disnet object (not class)")
        self._last_active_type = list(self.disnet_dict)[0] if self.disnet_dict else None
        self._active_type = None

    def add_disnet(self, disnet=None, cell=None, cell_list=None):
        """Add DisNet object of disnet_type
        """
        if disnet is not None:
            if cell is not None or cell_list is not None:
                raise ValueError("add_disnet: cell or cell_list should not be provided if disnet is not None")

        self.disnet_dict[type(disnet)] = disnet
        self._last_active_type = type(disnet)

    def synchronize_disnet(self, disnet_src, disnet_des):
        """Synchronize DisNet between disnet_src and disnet_des
        """
        if disnet_src is None:
            raise ValueError("synchronize_disnet: disnet_src is None")
        if disnet_des is None:
            raise ValueError("synchronize_disnet: disnet_des is None")
        if disnet_src == disnet_des:
            raise ValueError("synchronize_disnet: disnet_src and disnet_des are the same")

        if disnet_src in self.disnet_dict:
            G_src = self.disnet_dict[disnet_src]
        else:
            raise ValueError("synchronize_disnet: disnet_src not found")

        if disnet_des in self.disnet_dict:
            G_des = self.disnet_dict[disnet_des]
        else:
            # call constructor to create default DisNet object
            G_des = disnet_des()
            self.add_disnet(G_des)

        G_des.import_data(G_src.export_data())

        self._last_active_type = disnet_des

    def get_disnet(self, disnet_type=None):
        """Get DisNet object of disnet_type
        """
        if not disnet_type is None:
            self._active_type = disnet_type
        else:
            disnet_type = self._last_active_type
        if self._last_active_type is not None and self._active_type is not None and self._last_active_type != self._active_type:
            self.synchronize_disnet(self._last_active_type, self._active_type)
        return self.disnet_dict[disnet_type]
    
    def get_active_type(self):
        """Return the type of DisNet that is active
        """
        return self._active_type

    def add_nodes_links_from_list(self, rn, links, disnet_type=None):
        """Add nodes and links from list
        """
        G = self.get_disnet(disnet_type)
        G.add_nodes_links_from_list(rn, links)
        self._last_active_type = disnet_type

    def export_data(self):
        """Export DisNet data
        """
        G = self.get_disnet()
        return G.export_data()

    def import_data(self, data):
        """import DisNet data
        """
        G = self.get_disnet()
        return G.import_data(data)

    def write_json(self, filename):
        """Write DisNetManager data to JSON file
        """
        data = self.export_data()
        data['version'] = '1.0'
        data['nodes_attr'] = ['domain', 'index', 'x', 'y', 'z', 'constraint']
        data['segs_attr'] = ['node1', 'node2', 'bx', 'by', 'bz', 'nx', 'ny', 'nz']

        import json
        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        with open(filename, 'w') as f:
            json.dump(data, f, cls=NumpyEncoder, indent=4)

    def read_json(self, filename):
        """Read DisNetManager data to JSON file
        """
        import json
        with open(filename, 'r') as f:
            data = json.load(f)
        if data['version'] != '1.0':
            raise ValueError("read_json: version not supported")
        data['cell']['h'] = np.array(data['cell']['h'])
        data['cell']['origin'] = np.array(data['cell']['origin'])
        data['nodes'] = np.array(data['nodes'])
        data['segs']  = np.array(data['segs'])
        self.import_data(data)

    @property
    def G(self):
        """Return graph of DisNet
        """
        G = self.get_disnet()
        return G

    @property
    def cell(self):
        """Return cell of DisNet
        """
        G = self.get_disnet()
        return G.cell
        
    def num_nodes(self):
        """Return total number of nodes in DisNet
        """
        G = self.get_disnet()
        return G.num_nodes()

    def num_segments(self):
        """Return total number of segments in DisNet
        """
        G = self.get_disnet()
        return G.num_segments()

    def is_sane(self, disnet_type=None):
        """Check if DisNet is sane
        """
        G = self.get_disnet(disnet_type)
        return G.is_sane()
