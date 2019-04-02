from pymatgen import Lattice, Structure 
import numpy as np

struct = Structure.from_file('POSCAR.0')
#print(struct.lattice.b)
for item in struct: 
    print(item.coords[0])
