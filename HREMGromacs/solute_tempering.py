# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2021-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# MIT License                                               #
#############################################################
"""Functions to setup a solute tempering run with plumed
"""

import subprocess
from pathlib import Path

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
    mdp_file : str or path, default=None
        it is always needed except when `preprocess_top`=False
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
        they will be enumerated from 0 to `number_of_replicas`-1
    """

    if preprocess_top:
        tmp_elaborated_top = Path('TMP_elaborated_top.top').resolve()
        preprocess_topology(input_top_file=Path(top_file).resolve(),
                            output_top_file=tmp_elaborated_top,
                            gro_file=Path(gro_file).resolve(),
                            mdp_file=Path(mdp_file).resolve(),
                            gmx_path=gmx_path)

    else:
        tmp_elaborated_top = gro_file

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
