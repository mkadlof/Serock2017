import re


def point_reader(fname):
    """Read the points from PDB file format

        Args:
            fname (string) filename of single chromatin model in pdb file format

        Returns:
            (list) List of three floats tuples representing points in euclidean R^3
    """
    atoms = [i.strip() for i in open(fname) if re.search('^(ATOM|HETATM)', i)]
    points = []
    for i in atoms:
        x = float(i[30:38])
        y = float(i[38:46])
        z = float(i[46:54])
        points.append((x, y, z))
    return points


def save_points_as_xyz(points, filename, fmt='chimera', **kwargs):
    """Saves points as xyz file

        Args:
            points (list): three elements tuples
            filename (str): self-explanatory
            fmt (str): available formats:
                xyz         3 column tab separated
                idxyz       4 column tab separated (first column is an index)
                chimera     kwargs: molecule_name
    """
    prefix = ''
    atoms = ''
    suffix = ''
    n = len(points)
    if fmt == 'xyz':
        for i in range(n):
            x, y, z = points[i]
            atoms += ('{}\t{}\t{}\n'.format(x, y, z))
    elif fmt == 'idxyz':
        for i in range(n):
            x, y, z = points[i]
            atoms += ('{}\t{}\t{}\t{}\n'.format(i + 1, x, y, z))
    elif fmt == 'chimera':
        if kwargs is not None and 'molecule_name' in kwargs:
            molecule_name = kwargs['molecule_name']
        else:
            molecule_name = ''
        prefix = '{}\n{}\n'.format(n, molecule_name)
        for i in range(n):
            x, y, z = points[i]
            atoms += ('C\t{}\t{}\t{}\n'.format(x, y, z))

    with open(filename, 'w') as f:
        f.write(prefix)
        f.write(atoms)
        f.write(suffix)
    print("File {} saved...".format(filename))


def save_points_as_pdb(points, filename, render_connect=True, verbose=True):
    """Save points in PDB file format."""
    atoms = ''
    n = len(points)
    for i in range(n):
        x = points[i][0]
        y = points[i][1]
        try:
            z = points[i][2]
        except IndexError:
            z = 0.0
        atoms += (
            '{0:6}{1:>5}  {2:3}{3:}{4:3} {5:}{6:>4}{7:}   {8:>8.3f}{9:>8.3f}{10:>8.3f}{11:6.2f}{12:6.2f}{13:>12}\n'.
            format(
                "ATOM", i + 1, 'B', ' ', 'BEA', 'A', i + 1, ' ', max(x, -999), max(y, -999), max(z, -999), 0, 0, 'B'))
    connects = ''
    if render_connect:
        if n != 1:
            connects = 'CONECT    1    2\n'
            for i in range(2, n):
                connects += 'CONECT{:>5}{:>5}{:>5}\n'.format(i, i - 1, i + 1)
            connects += 'CONECT{:>5}{:>5}\n'.format(n, n - 1)
    pdb_file_content = atoms + connects
    with open(filename, 'w') as f:
        f.write(pdb_file_content)
    if verbose:
        print("File {} saved...".format(filename))
    return filename


def save_points_as_gro(points, filename, comment="comment"):
    n = len(points)
    x, y, z = zip(*points)
    d = max(x + y + z)
    l = ["{}\n".format(comment), str(n) + '\n']
    for i in range(n):
        x = points[i][0] / 10
        y = points[i][1] / 10
        z = points[i][2] / 10
        w = '{:5}{:5}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}\n'.format(i, "BEA", "B", i + 1, x, y, z)
        l.append(w)
    l.append('{0:5f} {0:5f} {0:5f}'.format(d))
    filename = '{}'.format(filename)
    open(filename, 'w').writelines(l)
    print("File {} saved...".format(filename))
    return filename

