from .calculate_rmsd import *
from rxcclib.utils.cclib.utils import PeriodicTable

def getrmsd(xyz1, xyz2, nohydrogen=False):

    p_atoms = [PeriodicTable.element[x] for x in xyz1.fchk.atomnos]
    p_all = np.array(xyz1.fchk.atomcoords[-1])
    q_atoms = [PeriodicTable.element[x] for x in xyz2.fchk.atomnos]
    q_all = np.array(xyz2.fchk.atomcoords[-1])


    if np.count_nonzero(p_atoms != q_atoms):
        exit("Atoms not in the same order")

    P = p_all
    Q = q_all

    if nohydrogen:
        not_hydrogens = np.where(p_atoms != 'H')
        P = p_all[not_hydrogens]
        Q = q_all[not_hydrogens]

    # elif args.remove_idx:
    #     N, = p_atoms.shape
    #     index = list(range(N))
    #     index = set(index) - set(args.remove_idx)
    #     index = list(index)
    #     P = p_all[index]
    #     Q = q_all[index]

    # elif args.add_idx:
    #     P = p_all[args.add_idx]
    #     Q = q_all[args.add_idx]


    # Calculate 'dumb' RMSD
    normal_rmsd = rmsd(P, Q)

    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    # if args.output:
    #     U = kabsch(P, Q)
    #     p_all -= Pc
    #     p_all = np.dot(p_all, U)
    #     write_coordinates(p_atoms, p_all, title="{} translated".format(args.structure_a))
    #     quit()

    return quaternion_rmsd(P, Q)

