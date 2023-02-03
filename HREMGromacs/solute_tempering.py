# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2021-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# MIT License                                               #
#############################################################
"""Functions to setup a solute tempering run with plumed
"""

import shutil
import subprocess
import tempfile
from pathlib import Path
import parmed
import copy

import PythonAuxiliaryFunctions.path as _path
import PythonAuxiliaryFunctions.run as _run


def scale_topology_with_plumed(input_top,
                               output_top,
                               scaling_value,
                               plumed='plumed'):
    """Scales the hamiltonian of a given topology with plumed

    In order to make a solute tempering you will need to have a certain number
    of topologies where the solute's hamiltonian is scaled
    This function scales one topology depending on the `scaling_value`
    The topology given as input must have been correctly modified and pre processed
    as explained in the plumed documentation https://www.plumed.org in the solute tempering part

    Parameters
    ------------
    input_top : str or pathlib.Path
        the input topology, the topology must have been preprocessed and modified
        in order for plumed to know which parts of the structure to scale
    output_top : str or pathlib.Path
        the output topology, shall be different from the input one
    scaling_value : float
        how much the solute's hamiltonian shall be scaled (1.0 means it is not scaled
        0.0 there is no potential energy is zero in the caled parts)
    plumed : str or pathlib.Path, optional, default=plumed
        the path to the plumed executable, as a default it will be searched in the
        PATH and in the current directory (in this order)

    Raises
    -----------
    RuntimeError
        if plumed fails to scale the topology
    ValueError
        if `input_top` == `output_top`
    """

    if input_top == output_top:
        raise ValueError(
            'The input and output topologies must be different files')

    plumed = _path.absolute_programpath(plumed)

    command = [str(plumed), 'partial_tempering', f'{scaling_value}']

    with open(input_top, 'r') as input_file_handle:
        with open(output_top, 'w') as output_file_handle:

            r = subprocess.run(command,
                               shell=False,
                               check=False,
                               stdin=input_file_handle,
                               stdout=output_file_handle,
                               stderr=subprocess.PIPE)

    if r.returncode != 0:
        raise RuntimeError(f'Plumed failure\n{r.stderr}')


def edit_preprocessed_top(input_top_file,
                          output_top_file,
                          protein_resSeq_to_scale,
                          ligand_resname,
                          exclude=('HOH', 'WAT', 'SOL')):
    """Edits a preprocessed top file for `scale_topology_with_plumed`

    `scale_topology_with_plumed` needs the preprocessed top file
    to be edited in order to know which residues to scale

    Parameters
    -----------
    input_top_file : str or pathlib.Path
    output_top_file : str or pathlib.Path
    protein_resSeq_to_scale : iterable(int)
        the resSeq (residue numbers) of the residues of the protein to
        scale, if more residues have the same resSeq no check will
        be done and they will all be scaled, often in the gromacs top file
        the resseqs restart from 1 for every molecule so if the protein is not
        the first molecule you have to check carefully
    ligand_resname : str
        as often in the gromacs top file
        the resseqs restart from 1 for every molecule
        for ligands you should give the residue names as written in the top file
    exclude : iterable(str), default=('HOH', 'WAT', 'SOL')
        as often in the gromacs top file
        the resseqs restart from 1 for every molecule
        it is useful to exclude certain molecules like the solvent
    """

    with open(input_top_file, 'r') as f:
        lines = f.readlines()

    #auxiliary variables
    is_atoms = False
    for i, line in enumerate(lines):

        # # pylint: disable=no-else-continue

        #if the line is empty or a comment i can go on with the for loop
        if not line.strip():
            continue
        elif line.strip()[0] == ';':
            continue
        #check if we are in the atoms part
        elif line.strip().replace(' ', '').split(';')[0] == '[atoms]':

            is_atoms = True
            continue
        #end of atoms part
        elif line.strip()[0] == '[':

            is_atoms = False
            continue

        if is_atoms:

            tmp_line = line.strip().split()

            # I am looking for lines like this
            # 1         N1      1    SER      N      1 0.18490000  14.006720
            if len(tmp_line) >= 4:
                if tmp_line[
                        3] not in exclude:  # Is not one of the excluded resnames
                    #heat the needed residues
                    if (int(tmp_line[2].strip())
                            in protein_resSeq_to_scale) or (ligand_resname
                                                            in tmp_line):

                        tmp_line[1] = tmp_line[1] + '_'

                        lines[i] = (' ' * 5).join(tmp_line) + '\n'

    with open(output_top_file, 'w') as f:
        for line in lines:
            f.write(line)


def preprocess_topology(input_top_file,
                        output_top_file,
                        gro_file,
                        mdp_file,
                        gmx_path='gmx'):
    """Uses grompp -pp to precompile the top file (needed for plumed)

    All the given paths must be absolute

    Parameters
    ---------------
    input_top_file : str or Path
    output_top_file : str or Path
    gro_file : str or Path
    mdp_file : str or Path
    gmx_path : str or Path, optional, default=will search gmx in path and working directory
    """

    gmx_path = _path.absolute_programpath(gmx_path)

    commands = [
        f'{gmx_path}', 'grompp', '-f', f'{mdp_file}', '-c', f'{gro_file}',
        '-p', f'{input_top_file}', '-maxwarn', '100', '-pp',
        f'{output_top_file}'
    ]

    _run.subprocess_run(commands,
                        error_string='grompp -pp failed',
                        shell=False)


def geometrical_progression(basis=0.2, denom=7):
    """returns a generator of a geometrical progression

    output(x) = basis^(x/(denom)) for x in [0, +inf)
    This is an useful way to get the values to scale the hamiltonian
    for HREM. A common default is having `basis`=0.2 for a protein-ligand system and 0.1 for
    a ligand-in-vacuum system and if you plan on using N replicas `denom` usually shall be N-1

    Paremeters
    ------------
    basis : float or int, default=0.2
    denom : float or int, default=7
        usually it is the number of HREM replicas minus 1

    Returns
    --------
    generator that returns float
    """

    yield 1.0

    x = 0

    while True:

        x += 1

        yield basis**(x / denom)


def prepare_topologies_for_hrem(top_file,
                                protein_resSeq_to_scale,
                                ligand_resname,
                                gro_file,
                                mdp_file=None,
                                number_of_replicas=8,
                                basis=0.2,
                                gmx_path='gmx',
                                plumed_path='plumed',
                                preprocess_top=True,
                                exclude=('HOH', 'WAT', 'SOL')):
    """High level function to generate the wanted number of scaled topology files for HREM

    This function is a high level interface to many of the functions in this module

    Parameters
    -------------
    top_file : str ot path
        the input topology file (will not be modified)
    protein_resSeq_to_scale : iterable(int)
        the resSeq (residue numbers) of the residues of the protein to
        scale, if more residues have the same resSeq no check will
        be done and they will all be scaled, often in the gromacs top file
        the resseqs restart from 1 for every molecule so if the protein is not
        the first molecule you have to check carefully
    ligand_resname : str
        as often in the gromacs top file
        the resseqs restart from 1 for every molecule
        for ligands you should give the residue names as written in the top file
    gro_file : str or path
        gromacs structure file
    mdp_file : str or path, default=None
        it is always needed except when `preprocess_top` = False
    number_of_replicas : int, default=8
        the number of HREM replicas you want to do
    basis : float, default=0.2
        the basis of the geometrical prograssion
        check `geometrical_progression` in this module
    gmx_path : str or Path, optional, default=gmx
        the path to the gromacs executable, as a default it will be searched in the
        PATH and in the current directory (in this order)
    plumed_path : str or pathlib.Path, optional, default=plumed
        the path to the plumed executable, as a default it will be searched in the
        PATH and in the current directory (in this order)
    preprocess_top : bool, optional, default=True
        if True the .top file will be processed in order to
        remove all the #include statements, this is necessary for plumed
        to work properly, but if the input topology is already preprocessed
        it can be skipped
    exclude : iterable(str), default=('HOH', 'WAT', 'SOL')
        residues not to scale even if they have the right resname
        as often in the gromacs top file
        the resseqs restart from 1 for every molecule
        it is useful to exclude certain molecules like the solvent

    Returns
    -----------
    list(Path)
        a list of the scaled topologies ordered from scaling 1 (not scaled) to scaling `basis`
        they will be enumerated from 0 to `number_of_replicas` - 1
    """

    if preprocess_top:
        tmp_elaborated_top = Path('TMP_elaborated_top.top').resolve()
        preprocess_topology(input_top_file=Path(top_file).resolve(),
                            output_top_file=tmp_elaborated_top,
                            gro_file=Path(gro_file).resolve(),
                            mdp_file=Path(mdp_file).resolve(),
                            gmx_path=gmx_path)

    else:
        shutil.copy(top_file, 'TMP_elaborated_top.top')
        tmp_elaborated_top = Path('TMP_elaborated_top.top').resolve()

    edit_preprocessed_top(input_top_file=tmp_elaborated_top,
                          output_top_file=tmp_elaborated_top,
                          protein_resSeq_to_scale=protein_resSeq_to_scale,
                          ligand_resname=ligand_resname,
                          exclude=exclude)

    output_top_prefix = Path(top_file).name.rsplit('.', 1)[0] + '_scaled_'
    scaling_values_generator = geometrical_progression(
        basis=basis, denom=number_of_replicas - 1)

    output_files = []

    for i in range(number_of_replicas):
        scaling_value = next(scaling_values_generator)

        output_files.append(Path(output_top_prefix + f'{i}.top').resolve())

        scale_topology_with_plumed(input_top=tmp_elaborated_top,
                                   output_top=output_files[-1],
                                   scaling_value=scaling_value,
                                   plumed=plumed_path)

    return output_files


def scale_dihedrals_in_parmed_structure(parmed_structure,
                               scaling_value,
                               atom_list,
                                is_gromacs=False):
    """Scales the dihedrals of the input parmed structure
    It does not change the input structure a copy is done
    
    Parameters
    ------------
    parmed_structure : parmed structure instance
        it will not be modified a copy is done
    scaling_value : float
        new_dihedral = old_dihedral * scaling_value
    atom_list : numpy.array
        the zero indexed atom indexes of the atoms in the
        dihedrals to scale
    is_gromacs : bool, default=False
        because I only know how to do it for gromacs
        if False I will have to transform it in a gromacs
        parmed structure (time consuming)
    
    Returns
    ----------
    the modified parmed structure
    """

    # I only know how to do it for gromacs structures
    if is_gromacs:
        parmed_structure = copy.deepcopy(parmed_structure)
    
    else:
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdirname = Path(tmpdirname)

            parmed_structure.save(str(tmpdirname / "tmp.top"))
            parmed_structure.save(str(tmpdirname / "tmp.pdb"), renumber=False)

            parmed_structure = parmed.load_file(str(tmpdirname / "tmp.top"),
                xyz=str(tmpdirname / "tmp.pdb"))

    for dihedral in parmed_structure.dihedrals:
        if not dihedral.improper:
            if dihedral.atom2.idx in atom_list and dihedral.atom3.idx in atom_list:
                new_type = copy.copy(dihedral.type)

                # Check if iterable
                try:
                    for t in new_type:
                        t.phi_k *= scaling_value

                except TypeError:
                    new_type.phi_k *= scaling_value

                parmed_structure.dihedral_types.append(new_type)
                parmed_structure.dihedral_types.claim()
                dihedral.type = parmed_structure.dihedral_types[-1]

    parmed_structure.dihedrals.changed = True
    parmed_structure.dihedral_types.changed = True

    return parmed_structure


def scale_14_interactions_in_parmed_structure(parmed_structure,
                               scaling_value,
                               atom_list,
                                is_gromacs=False):
    """Scales the 14 interactions of the input parmed structure
    It does not change the input structure a copy is done

    It creates some new atom tipes with _ after the original name
    except if _ is already the last part of the atom name
    
    Parameters
    ------------
    parmed_structure : parmed structure instance
        it will not be modified a copy is done
    scaling_value : float
        new_14_interactions = old_14_interactions * scaling_value
    atom_list : numpy.array
        the zero indexed atom indexes of the atoms in the
        dihedrals to scale
    is_gromacs : bool, default=False
        because I only know how to do it for gromacs
        if False I will have to transform it in a gromacs
        parmed structure (time consuming)
    
    Returns
    ----------
    the modified parmed structure
    """

    # I only know how to do it for gromacs structures
    if is_gromacs:
        parmed_structure = copy.deepcopy(parmed_structure)
    
    else:
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdirname = Path(tmpdirname)

            parmed_structure.save(str(tmpdirname / "tmp.top"))
            parmed_structure.save(str(tmpdirname / "tmp.pdb"), renumber=False)

            parmed_structure = parmed.load_file(str(tmpdirname / "tmp.top"),
                xyz=str(tmpdirname / "tmp.pdb"))

    # Necessary to have [ pairs ] defined explicitly (1-4 interactions)
    parmed_structure.defaults.gen_pairs = "no"


    for atom in parmed_structure.atoms:
        if atom.idx in atom_list:
            atom.atom_type = copy.copy(atom.atom_type)
            atom.epsilon_14 *= scaling_value
            atom.atom_type.epsilon_14 *= scaling_value

            if list(atom.atom_type.name)[-1] != "_":
                atom.atom_type.name += "_"
                
            atom.type = atom.atom_type.name
            atom.name = atom.atom_type.name

    # Sometimes 1-4 interactions are here
    if parmed_structure.adjusts:
        for adjust in parmed_structure.adjusts:
            if adjust.atom1.idx in atom_list and adjust.atom2.idx in atom_list:
                adjust.type = copy.copy(adjust.type)
                # LJ 1-4
                adjust.type.epsilon *= scaling_value
                # Charges 1-4
                adjust.type.chgscale *= scaling_value

                parmed_structure.adjust_types.append(adjust.type)

                parmed_structure.adjust_types.claim()

                parmed_structure.adjusts.claim()

        parmed_structure.adjust_types.changed = True
        parmed_structure.adjusts.changed = True

    parmed_structure.atoms.changed = True

    return parmed_structure


def scale_vdw_and_charges_in_parmed_structure(parmed_structure,
                               scaling_value,
                               atom_list,
                                is_gromacs=False,
                                scale_charges=True,
                                scale_vdw=True):
    """Scales vdw and charges of the input parmed structure
    It does not change the input structure a copy is done
    it does not touch the 14 interactions

    It creates some new atom tipes with _ after the original name
    except if _ is already the last part of the atom name
    
    Parameters
    ------------
    parmed_structure : parmed structure instance
        it will not be modified a copy is done
    scaling_value : float
        new = old * scaling_value
    atom_list : numpy.array
        the zero indexed atom indexes of the atoms in the
        dihedrals to scale
    is_gromacs : bool, default=False
        because I only know how to do it for gromacs
        if False I will have to transform it in a gromacs
        parmed structure (time consuming)
    scale_charges : bool, default=True
    scale_vdw : bool, default=True
    
    Returns
    ----------
    the modified parmed structure
    """

    # I only know how to do it for gromacs structures
    if is_gromacs:
        parmed_structure = copy.deepcopy(parmed_structure)
    
    else:
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdirname = Path(tmpdirname)

            parmed_structure.save(str(tmpdirname / "tmp.top"))
            parmed_structure.save(str(tmpdirname / "tmp.pdb"), renumber=False)
            
            parmed_structure = parmed.load_file(str(tmpdirname / "tmp.top"),
                xyz=str(tmpdirname / "tmp.pdb"))

    for atom in parmed_structure.atoms:
        if atom.idx in atom_list:
            atom.atom_type = copy.copy(atom.atom_type)

            if scale_vdw:
                atom.epsilon *= scaling_value
                atom.atom_type.epsilon *= scaling_value
            
            if scale_charges:
                atom.charge *= scaling_value
                atom.atom_type.charge *= scaling_value
            
            if list(atom.atom_type.name)[-1] != "_":
                atom.atom_type.name += "_"

            atom.type = atom.atom_type.name
            atom.name = atom.atom_type.name
    parmed_structure.atoms.changed = True

    return parmed_structure


def scale_topologies_with_parmed(input_top,
                                output_top,
                                xyz,
                                atom_list,
                                dihedral_scaling=None,
                                interaction14_scaling=None,
                                vdw_scaling=None,
                                q_scaling=None):
    """High level function to scale contributions in a topology file
    It uses parmed, therefore it should work with any supported format
    as input and output but only gromacs topology files have been
    tested
    
    Parameters
    ------------
    input_top : str or Path
        input topology file, any parmed supported format
    output_top : str or Path
        output topology file, any parmed supported format
    xyz : str or Path
        a file contaning the coordinates of the system (like PDB),
        any parmed supported format
    atom_list : numpy.array
        the zero indexed atom indexes of the atoms in the
        dihedrals to scale
    dihedral_scaling : float, default=None
        new = old * scaling_value
    interaction14_scaling : float, default=None
        new = old * scaling_value
    vdw_scaling : float, default=None
        scale van der waals
        new = old * scaling_value
    q_scaling : float, default=None
        scale charges
        new = old * scaling_value
    
    Warning
    --------------
    it should work with any parmed supported format
    as input and output but only gromacs topology files have been
    tested

    If output_top = input_top
    the input one will be overwritten
    """
    parmed_structure = parmed.load_file(str(input_top), xyz=str(xyz))
    
    # I only know how to do it for gromacs structures
    if Path(input_top).suffix != ".top":
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdirname = Path(tmpdirname)

            parmed_structure.save(str(tmpdirname / "tmp.top"))
            parmed_structure.save(str(tmpdirname / "tmp.pdb"), renumber=False)

            parmed_structure = parmed.load_file(str(tmpdirname / "tmp.top"),
                xyz=str(tmpdirname / "tmp.pdb"))

    if dihedral_scaling is not None:
        parmed_structure = scale_dihedrals_in_parmed_structure(parmed_structure,
                                scaling_value=dihedral_scaling,
                                atom_list=atom_list,
                                is_gromacs=True)

    if interaction14_scaling is not None:
        parmed_structure = scale_14_interactions_in_parmed_structure(parmed_structure,
                                scaling_value=interaction14_scaling,
                                atom_list=atom_list,
                                is_gromacs=True)

    if vdw_scaling is not None:
        parmed_structure = scale_vdw_and_charges_in_parmed_structure(parmed_structure,
                                scaling_value=vdw_scaling,
                                atom_list=atom_list,
                                is_gromacs=True,
                                scale_charges=False,
                                scale_vdw=True)

    if q_scaling is not None:
        parmed_structure = scale_vdw_and_charges_in_parmed_structure(parmed_structure,
                                scaling_value=vdw_scaling,
                                atom_list=atom_list,
                                is_gromacs=True,
                                scale_charges=True,
                                scale_vdw=False)
    
    parmed_structure.save(str(output_top), overwrite=True)