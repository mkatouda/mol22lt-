#!/usr/bin/python

# Copyright (c) 2021 Michio Katouda
# This software is released under the MIT License, see LICENSE.

import sys
import os
import math

def read_mol2(f_mol2):

    # Open mol2 file

    with open(f_mol2) as f:
        lines = f.readlines()

    # Read atom, atomtype, coord, charge

    atomtype = []
    atom = []
    coord = []
    charge = []
    iatom = 0
    flag = False
    for line in lines:
        if line[0:13] == '@<TRIPOS>BOND':
            flag = False
            continue
        elif line[0:13] == '@<TRIPOS>ATOM':
            flag = True
            continue
        elif flag:
            tmp = line.split()
            atomtype.append(tmp[5])
            atom.append(tmp[1])
            coord.append([float(tmp[2]), float(tmp[3]), float(tmp[4])])
            charge.append(float(tmp[8]))
            iatom += 1

    # Read bond index

    bondidx = []
    ibonds = 0
    flag = False
    for line in lines:
        if line[0:13] == '@<TRIPOS>SUBS':
            flag = False
            continue
        elif line[0:13] == '@<TRIPOS>BOND':
            flag = True
            continue
        elif flag:
            tmp = line.split()
            bondidx.append([int(tmp[1])-1, int(tmp[2])-1])
            ibonds += 1

    return atom, atomtype, coord, charge, bondidx


def read_gro(f_gro, f_itp):

    # Open Gromacs itp file

    with open(f_itp) as f:
        lines = f.readlines()

    # Read atom, atomtype, charge

    atomtype = []
    atom = []
    charge = []
    iatom = 0
    flag = False
    for line in lines:
        if line[0] == ';':
            continue
        elif line[0] == '\n':
            flag = False
            continue
        elif line[0:9] == '[ atoms ]':
            flag = True
            continue
        elif flag:
            tmp = line.split()
            atomtype.append(tmp[1])
            atom.append(tmp[4])
            charge.append(float(tmp[6]))
            iatom += 1
    natoms = iatom

    # Read bond index

    bondidx = []
    ibonds = 0
    flag = False
    for line in lines:
        if line[0] == ';':
            continue
        elif line[0] == '\n':
            flag = False
            continue
        elif line[0:9] == '[ bonds ]':
            flag = True
            continue
        elif flag:
            tmp = line.split()
            bondidx.append([int(tmp[0])-1, int(tmp[1])-1])
            ibonds += 1

    # Open Gromacs gro file

    with open(f_gro) as f:
        lines = f.readlines()

    # Read coordinate

    nm2ang = 10.0
    coord = []
    for iatom in range(natoms):
        line = lines[iatom+2].split()
        coord.append([float(line[4])*nm2ang, float(line[5])*nm2ang, float(line[6])*nm2ang])

    return atom, atomtype, coord, charge, bondidx


def write_lt(f_lt, molecule, atom, atomtype, coord, charge, bondidx, itype=0, ff='gaff'):

    natoms = len(atom)
    nbonds = len(bondidx)

    molname = molecule
    cdiff = 0.0
    iatom_bgn = 0
    iatom_end = natoms
    if itype == 1:
        molname = molecule + '-a'
        cdiff = charge[0] / float(natoms - 1)
        iatom_bgn = 1
        iatom_end = natoms
    elif itype == 2:
        molname = molecule + '-b'
        cdiff = charge[natoms - 1] / float(natoms - 1)
        iatom_bgn = 0
        iatom_end = natoms - 1
    elif itype == 3:
        molname = molecule + '-c'
        cdiff = (charge[0] + charge[natoms - 1]) / float(natoms - 2)
        iatom_bgn = 1
        iatom_end = natoms - 1


    # Write Moltemplate LT file

    with open(f_lt, mode='w') as f:

        # Write header

        f.write('import "{}.lt"\n\n'.format(ff))   
        f.write('{} inherits {} {{\n\n'.format(molname, ff.upper()))

        # Write coordinate

        f.write('  write(\'Data Atoms\') {\n')
        for iatom in range(iatom_bgn, iatom_end):
            [cx, cy, cz] = coord[iatom]
            f.write('    $atom:{:6} $mol @atom:{:6} {:>10.6f} {:>12.6f} {:>12.6f} {:>12.6f}\n' \
                    .format(atom[iatom], atomtype[iatom], charge[iatom] + cdiff, cx, cy, cz))
        f.write('  }\n\n')

        # Write bonding information

        f.write('  write(\'Data Bond List\') {\n')
        for ibond in range(nbonds):
            [iatom1, iatom2] = bondidx[ibond]
            if itype == 1:
                if iatom1 == 0 or iatom2 == 0:
                    continue
            elif itype == 2:
                if iatom1 == natoms -1 or iatom2 == natoms -1:
                    continue
            elif itype == 3:
                if iatom1 == 0 or iatom2 ==0:
                    continue
                if iatom1 == natoms -1 or iatom2 == natoms -1:
                    continue
            f.write('    $bond:{:10} $atom:{:6} $atom:{:6}\n' \
                    .format(atom[iatom1]+atom[iatom2], atom[iatom1], atom[iatom2]))
        f.write('  }\n')
        f.write('\n}\n')


def mol22lt_main(f_in, ff):

    molecule = os.path.splitext(os.path.basename(f_in))[0]
    f_in_ext = os.path.splitext(f_in)[1][1:]

    f_mol2 = molecule + '.mol2'
    f_gro = molecule + '.gro'
    f_itp = molecule + '.itp'
    f_lt  = molecule + '.lt'

    if f_in_ext == 'mol2':
        atom, atomtype, coord, charge, bondidx = read_mol2(f_mol2)
    elif f_in_ext == 'gro' or f_in_ext == 'itp':
        atom, atomtype, coord, charge, bondidx = read_gro(f_gro, f_itp)
    else: 
        print('Error: Cannot convert file {}'.format(f_in))
        sys.exit()

    f_lt  = molecule + '.lt'
    itype = 0
    write_lt(f_lt, molecule, atom, atomtype, coord, charge, bondidx, itype, ff)

    f_lt  = molecule + '-a.lt'
    itype = 1
    write_lt(f_lt, molecule, atom, atomtype, coord, charge, bondidx, itype, ff)

    f_lt  = molecule + '-b.lt'
    itype = 2
    write_lt(f_lt, molecule, atom, atomtype, coord, charge, bondidx, itype, ff)

    f_lt  = molecule + '-c.lt'
    itype = 3
    write_lt(f_lt, molecule, atom, atomtype, coord, charge, bondidx, itype, ff)


if __name__ == '__main__':
    f_in = "mol.mol2"
    ff = 'gaff'

    argvs = sys.argv
    argc = len(argvs)
    if(argc > 1):
        f_in = argvs[1]
    if(argc > 2):
        ff = argvs[2]

    mol22lt_main(f_in, ff)
