# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2021-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# MIT License                                               #
#############################################################
"""Module containing the functions to setup the multidir directories for the HREM
"""

import shutil
from pathlib import Path


def make_hrem_multidir_directory(number_of_replicas=8,
                                 dir_name='BATTERY',
                                 exist_ok=False):
    """Makes the multidir directory and sub-directories for HREM

    Parameters
    ----------
    number_of_replicas : int, default=8
    dir_name : str, default=BATTERY
    exist_ok : bool, default=False
        if True won't raise an exception if the direcotry already exists

    Returns
    --------
    pathlib.Path
        the path object of the directory

    Notes
    ----------
    the sirectory will be filled with sub-directories called
    scaled0, ..., scaled `number_of_replicas`-1
    """

    dir_path = Path(dir_name).resolve()
    dir_path.mkdir(parents=True, exist_ok=exist_ok)

    for i in range(number_of_replicas):
        (dir_path / f'scaled{i}').mkdir(parents=True, exist_ok=True)

    return dir_path


def copy_file_in_hrem_multidir(file_name='empty_plumed.dat',
                               dir_name='BATTERY'):
    """copies a file inside all the scaled* direcotries inside a multidir directory

    if the file doesn't exist it will be created on the fly (empty), very useful to copy the
    needed plumed input file for HREM (often empty)

    Parameters
    -----------
    file_name : str or path, default=empty_plumed.dat
    dir_name : str or path, default=BATTERY
    """

    dir_name = Path(dir_name)
    sub_dirs = dir_name.glob('scaled*')

    file_name = Path(file_name)
    file_name.touch()

    for sub_dir in sub_dirs:

        shutil.copy(str(file_name), str(sub_dir))


def make_multiple_hrem_batteries(number_of_batteries,
                                 replicas_per_battery=8,
                                 plumed_file='empty_plumed.dat'):
    """high level function that makes multiple hrem-multidir directories

    often it is not possible to run the number of nanoseconds wanted in
    the maximum wall time of an HPC cluster (often 24h) so you might want
    to run multiple separated runs in parallel
    The directories created will be called BATTERY0, BATTERY1, ...
    andcontaind scaled0, scaled1, ...

    Parameters
    ------------
    number_of_batteries : int
        the number of independent hrem you want to do
    replicas_per_battery : int, default=8
        the number of replicas you will run, and therefore the number of scaled*
        sub-directories
    plumed_file : str or path, default=empty_plumed.dat
        the plumed input file, if it doesn't exist it is created empty

    Returns
    ---------
    list(pathlib.Path)
        the path objects of the various created directories (not the sub-direcotries)

    Notes
    --------
    on a gpu, with a 15000 atoms system a very conservative extimate is
    that you can do 15 nanoseconds in 24h
    on only cpu on the other end you may go with 5 ns on 64 cores with broadwell cpus
    but in any case this are arbitrary and very conservative numbers
    """
    output = []

    for i in range(number_of_batteries):
        output.append(
            make_hrem_multidir_directory(
                number_of_replicas=replicas_per_battery,
                dir_name=f'BATTERY{i}',
                exist_ok=True))

        copy_file_in_hrem_multidir(file_name=plumed_file, dir_name=output[-1])

    return output


def make_multidir_mdrun_string_for_hrem(multidir_directories,
                                        gromacs_path='gmx_mpi',
                                        plumed_file='empty_plumed.dat',
                                        deffnm='HREM',
                                        tpr_file='HREM.tpr',
                                        replex=100):
    """Helper function to make the mdrun string for an HREM (with multidir)

    Parameters
    -----------
    multidir_directories : list(str) or list(path)
        the list of the multidir directories. must be ordered
        the first one is the reference state the last one the most
        scaled replica
    gromacs_path : str, default=gmx_mpi
    plumed_file : str, default=empty_plumed.dat
    deffnm : str, default=HREM
    tpr_file : str, default=HREM.tpr
    replex : int, default=100
        after how many steps a swap shall be attempted

    Returns
    ---------
    str
        it will look something like
        gmx_mpi mdrun -v -plumed empty_plumed.dat -replex 100 -hrex -dlb no
        -multidir BATTERY/scaled0 BATTERY/scaled1 -s HREM.tpr -deffnm HREM
    """

    multidir_str = [str(i) for i in multidir_directories]
    multidir_str = ' '.join(multidir_str)

    multidir_str = (
        f'{gromacs_path} mdrun -v -plumed {plumed_file}'
        f' -replex {replex} -hrex -dlb no -s {tpr_file} -deffnm {deffnm}'
        f' -multidir {multidir_str}')

    return multidir_str
