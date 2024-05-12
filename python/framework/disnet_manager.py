"""@package docstring
DisNetManager: class for managing multiple implementations of dislocation network

Implements synchronization between different implementations of DisNet
"""

class DisNetManager:
    """Class for managing multiple implementations of dislocation network

    Implements synchronization between different implementations of DisNet
    """
    def __init__(self, disnet_dict: dict={}, disnet=None):
        if type(disnet_dict) is not dict:
            raise ValueError("DisNetManager: disnet_dict should be a dictionary")
        self.disnet_dict = disnet_dict
        if disnet_dict == {} and disnet is not None:
            self.disnet_dict[type(disnet)] = disnet
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

    @property
    def G(self):
        """Return networkx graph of DisNet
        """
        G = self.get_disnet()
        return G

    @property
    def cell(self):
        """Return cell of DisNet
        """
        G = self.get_disnet()
        return G.cell

    def is_sane(self, disnet_type=None):
        """Check if DisNet is sane
        """
        G = self.get_disnet(disnet_type)
        return G.is_sane()
