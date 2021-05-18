import numpy as np
import matplotlib.pyplot as plt
from ase import Atom, Atoms, cluster
from ase.visualize import view
from ase.md.analysis import DiffusionCoefficient
from ase.io.trajectory import Trajectory
from ase.lattice.cubic import FaceCenteredCubic

from ase.geometry.dimensionality import isolate_components

def get_components(atoms, kcutoff, viewer = False):
    """
    A function to isolate individual components from the input structure and
    view if desired.
    
    Parameters
    Accepts:
    atoms           (Atoms object) ASE atoms object of desired structure
    kcutoff         (int) Cutoff of components to consider. Ideally, this should be the
                            Lennard-Jones cutoff value used.
    viewer          (bool) When True, returns the view of each component.
                            Do not set True unless you only have a few components.
    
    Returns:
    comp           (collections.defaultdict)
    """
    comp = isolate_components(atoms, kcutoff = kcutoff)
    print("counts:", [(k, len(v)) for k, v in sorted(comp.items())])
    if viewer == True:
        for dim, components in comp.items():
            for atoms in components:
                view(atoms, block=True)
    return comp

def get_radial_distance(atoms, bins = 100):
    """
    A histogram wrapperfunction for getting the RDF from the ASE atoms
    object based off of the np.histogram function
    
    Parameters
    Accepts:
    atoms           (Atoms object) ASE atoms object of desired structure
    bins            (int, array) desired bins as described in np.histogram
    
    Returns
    dist_counts     (array) RDF counts
    bin_centers     (array) Bin centers
    """
    dists = atoms.get_all_distances()
    dist_counts, dist_bins = np.histogram(dists, bins, 
                                          (np.min(dists), np.max(dists)))
    w = dist_bins[1]-dist_bins[0]
    bin_centers = dist_bins[:-1]+w/2
    if bins[0] == 0:
        dist_counts[0] = 0
    return dist_counts, bin_centers

def get_integrated_rd(comp, bins, return_fig = False):
    """
    Returns integrated radial distances over all components in a list of components.
    
    Useful when trying to break component sizes up by dimension or size to determine
    if there is a difference in radial distances. Otherwise, the same as get_radial_distance 
    if run on all of the components within a frame
    
    Parameters
    Accepts:
    comp           (collections.defaultdict) A collection of atoms objects isolated from
                          the reaction cell
    bins           (int, array) desired bins as described in np.histogram
    return_fig
    Returns:
    normalized_counts
    
    """
    counter = 0
    for dim, components in comp.items():
        for atoms in components:
            counts, bins_out = get_rdf(atoms, bins)
            if counter == 0:
                int_counts = counts
                counter += 1
            else:    
                int_counts += counts
    
    normalized_counts = int_counts / np.max(int_counts)
    if return_fig == True:
        fig,ax = plt.subplots(figsize = (8,5))
        ax.plot(bins_out, normalized_counts, c = 'b', lw = 2)
        ax.set_xlabel('Distance (Angstroms)',size =15)
        ax.set_ylabel('Normalized counts', size=15)
        ax.set_title('Radial Distribution Function- Clusters', size = 18)
        ax.set_xlim(0,30)
    return normalized_counts, bins_out

def get_cluster_sizes(comp):
    """
    Returns a nd.array of cluster sizes based on the components isolated.
    Parameters
    Accepts:
    comp          (collections.defaultdict) A collection of atoms objects isolated
                         from the reaction cell
    
    Returns:
    cluster_sizes (nd.array) Cluster sizes of input components
    """
    v = [(k, len(v)) for k, v in sorted(comp.items())]
    cluster_sizes = np.zeros([v[0][1],1])
    iter = 0
    for dim, components in comp.items():
        for atoms in components:
            cluster_sizes[iter, 0] = len(atoms.positions)
            iter += 1
    return cluster_sizes