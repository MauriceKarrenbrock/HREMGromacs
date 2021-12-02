# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2021-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# MIT License                                               #
#############################################################
"""Some classes to easily make mdp files for HREM

Of course this mdp files are nothing more than a suggestion
"""

import FSDAMGromacs.mdp_files as _mdp_files


class HREMSuperclass(_mdp_files.MdpFile):
    """Makes the mdp file for a protein-ligand system

    It is a superclass of `FSDAMGromacs.mdp_files.MdpFile`

    Attributes
    -----------
    _template : list
        Private, list of strings of the template
        that will be written on the mdp file
        ("\n" new line caracters may be missing)

    Parameters
    ----------
    mdp_file : str
        the name of the output mdp file (or the path)
    timestep_ps : float, default=0.002
        the timestep of the MD run (ps)
    number_of_steps : int, default=16000000
        number of MD steps
    temperature : float, default=298.15
        temperature in Kelvin (K)
    COM_pull_goups : list of strings
        the list of gromacs groups (System,Protein,<ligand_resname>,...)
        that will have an harmonic COM-COM (center of mass) constrain
        Usually for unbound state no pulling is needed (default)
        for bound state it is ["Protein", "<ligand residue name>", "<dummy heavy atom>"]
    harmonic_kappa : list
        [ ["group_1", "group_2", harmonic_kappa_value], ... ] (str, str, float)
        it is a nested list containing the couple-couple harmonic kappa value
        for the umbrella COM-COM pulling, a good numbe may be 120
        if you don't want to groups to pull each other set kappa to 0
    pbc_atoms : iterable of int, optional
        an iterable containing the number of each nearest atom to the
        geometric center of the `COM_pull_goups` molecule
        must be as long as `COM_pull_goups`, for small molecules
        you can set it to zero (gromacs will guess it) for big ones (proteins)
        it must be given, if you keep it None it will gess it on any molecule
    constraints : str, optional, default=h-bonds
        which bounds to constrain, the syntax is tha same as the gromacs one
        none, h-bonds, all-bonds

    Methods
    ---------
    execute()
        the only pubblic method, writes the mdp file
        with the right template

    Notes
    -----------
    `COM_pull_goups` in the bound state may need a dummy heavy atom
    (no LJ nor Q but infinite mass) because
    COM COM pulling in gromacs crashes because of poor PBC implemetation if the protein
    crosses the box
    """
    def __init__(self,
                 mdp_file,
                 timestep_ps=0.002,
                 number_of_steps=16000000,
                 temperature=298.15,
                 COM_pull_goups=None,
                 harmonic_kappa=None,
                 pbc_atoms=None,
                 constraints=None):

        if constraints is None:
            constraints = 'h-bonds'

        super().__init__(mdp_file,
                         alchemical_molecule=None,
                         timestep_ps=timestep_ps,
                         number_of_steps=number_of_steps,
                         temperature=temperature,
                         lambda_steps=None,
                         COM_pull_goups=COM_pull_goups,
                         harmonic_kappa=harmonic_kappa,
                         pbc_atoms=pbc_atoms,
                         constraints=constraints)

    def _hook(self):
        self._template += self._create_COMCOM_pulling_strings()


class ProteinLigandHREM(HREMSuperclass):
    """Subclasses `HREMSuperclass`
    """
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self._template = [
            '; VARIOUS PREPROCESSING OPTIONS',
            '; Preprocessor information: use cpp syntax.',
            '; e.g.: -I/home/joe/doe -I/home/mary/roe',
            'include                  =',
            '; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)',
            'define                   =', '', '; RUN CONTROL PARAMETERS',
            'integrator               = md', '; Start time and timestep in ps',
            f'dt                       = {self.timestep_ps}',
            f'nsteps                   = {self.number_of_steps}',
            '; For exact run continuation or redoing part of a run',
            'init-step                = 0',
            '; Part index is updated automatically on checkpointing (keeps files separate)',
            'simulation-part          = 1',
            '; mode for center of mass motion removal',
            'comm-mode                = Linear',
            '; number of steps for center of mass motion removal',
            'nstcomm                  = 100',
            '; group(s) for center of mass motion removal',
            'comm-grps                =', '',
            '; TEST PARTICLE INSERTION OPTIONS',
            'rtpi                     = 0.05', '', '; OUTPUT CONTROL OPTIONS',
            '; Output frequency for coords (x), velocities (v) and forces (f)',
            'nstxout                  = 5000',
            'nstvout                  = 5000',
            'nstfout                  = 100000',
            '; Output frequency for energies to log file and energy file',
            'nstlog                   = 5000',
            'nstcalcenergy            = 100',
            'nstenergy                = 5000',
            '; Output frequency and precision for .xtc file',
            'nstxtcout                = 10000',
            'xtc-precision            = 1000',
            '; This selects the subset of atoms for the .xtc file. You can',
            '; select multiple groups. By default all atoms will be written.',
            'xtc-grps                 =', '; Selection of energy groups',
            'energygrps               = System', '',
            '; NEIGHBORSEARCHING PARAMETERS',
            '; cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)',
            '; nblist update frequency', 'cutoff-scheme            = Verlet',
            'nstlist                  = 20',
            'verlet-buffer-tolerance  = 0.0001',
            '; ns algorithm (simple or grid)',
            'ns_type                  = grid',
            '; Periodic boundary conditions: xyz, no, xy',
            'pbc                      = xyz', 'periodic-molecules       = no',
            '; Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,',
            '; a value of -1 means: use rlist', '; nblist cut-off',
            'rlist                    = 1.0',
            '; long-range cut-off for switched potentials',
            'rlistlong                = -1', '',
            '; OPTIONS FOR ELECTROSTATICS AND VDW',
            '; Method for doing electrostatics',
            'coulombtype              = PME', 'rcoulomb-switch          = 0',
            'rcoulomb                 = 1.0',
            '; Relative dielectric constant for the medium and the reaction field',
            'epsilon-r                = 1', 'epsilon-rf               = 0',
            '; Method for doing Van der Waals',
            'vdw-type                 = Cut-off', '; cut-off lengths',
            'rvdw-switch              = 0', 'rvdw                     = 1.0',
            '; Apply long range dispersion corrections for Energy and Pressure',
            'DispCorr                 = EnerPres',
            '; Extension of the potential lookup tables beyond the cut-off',
            'table-extension          = 1',
            '; Separate tables between energy group pairs',
            'energygrp-table          =',
            '; Spacing for the PME/PPPM FFT grid',
            'fourierspacing           = 0.12',
            '; FFT grid size, when a value is 0 fourierspacing will be used',
            'fourier-nx               = 0', 'fourier-ny               = 0',
            'fourier-nz               = 0', '; EWALD/PME/PPPM parameters',
            'pme-order                = 4', 'ewald-rtol               = 1e-06',
            'ewald-geometry           = 3d', 'epsilon-surface          =',
            'optimize-fft             = no', '',
            '; IMPLICIT SOLVENT ALGORITHM', 'implicit-solvent         = No',
            '', '; OPTIONS FOR WEAK COUPLING ALGORITHMS',
            '; Temperature coupling', 'tcoupl                   = v-rescale',
            'nsttcouple               = -1', 'nh-chain-length          = 1',
            '; Groups to couple separately',
            'tc-grps                  = System',
            '; Time constant (ps) and reference temperature (K)',
            'tau-t                    = 0.1',
            f'ref-t                    = {self.temperature}',
            '; pressure coupling', 'pcoupl                   = c-rescale',
            'pcoupltype               = Isotropic',
            'nstpcouple               = -1',
            '; Time constant (ps), compressibility (1/bar) and reference P (bar)',
            'tau-p                    = 1.0',
            'compressibility          = 4.6e-5',
            'ref-p                    = 1',
            '; Scaling of reference coordinates, No, All or COM',
            'refcoord-scaling         = COM', '',
            '; GENERATE VELOCITIES FOR STARTUP RUN',
            'gen-vel                  = no',
            f'gen-temp                 = {self.temperature}',
            'gen-seed                 = 173529', '', '; OPTIONS FOR BONDS',
            f'constraints              = {self.constraints}',
            '; Type of constraint algorithm',
            'constraint-algorithm     = Lincs',
            '; Do not constrain the start configuration',
            'continuation             = yes',
            '; Use successive overrelaxation to reduce the number of shake iterations',
            'Shake-SOR                = no', '; Relative tolerance of shake',
            'shake-tol                = 0.00001',
            '; Highest order in the expansion of the constraint coupling matrix',
            'lincs-order              = 5',
            '; Number of iterations in the final step of LINCS. 1 is fine for',
            '; normal simulations, but use 2 to conserve energy in NVE runs.',
            '; For energy minimization with constraints it should be 4 to 8.',
            'lincs-iter               = 2',
            '; Lincs will write a warning to the stderr if in one step a bond',
            '; rotates over more degrees than',
            'lincs-warnangle          = 30',
            '; Convert harmonic bonds to morse potentials',
            'morse                    = no'
        ]


class OnlyLigandHREM(HREMSuperclass):
    """Subclasses `HREMSuperclass`

    Ligand in vacuum
    """
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self._template = [
            '; VARIOUS PREPROCESSING OPTIONS',
            '; Preprocessor information: use cpp syntax.',
            '; e.g.: -I/home/joe/doe -I/home/mary/roe',
            'include                  =',
            '; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)',
            'define                   =', '', '; RUN CONTROL PARAMETERS',
            'integrator               = md', '; Start time and timestep in ps',
            f'dt                       = {self.timestep_ps}',
            f'nsteps                   = {self.number_of_steps}',
            '; For exact run continuation or redoing part of a run',
            'init-step                = 0',
            '; Part index is updated automatically on checkpointing (keeps files separate)',
            'simulation-part          = 1',
            '; mode for center of mass motion removal',
            'comm-mode                = Linear',
            '; number of steps for center of mass motion removal',
            'nstcomm                  = 100',
            '; group(s) for center of mass motion removal',
            'comm-grps                =', '',
            '; TEST PARTICLE INSERTION OPTIONS',
            'rtpi                     = 0.05', '', '; OUTPUT CONTROL OPTIONS',
            '; Output frequency for coords (x), velocities (v) and forces (f)',
            'nstxout                  = 5000',
            'nstvout                  = 5000',
            'nstfout                  = 100000',
            '; Output frequency for energies to log file and energy file',
            'nstlog                   = 5000',
            'nstcalcenergy            = 100',
            'nstenergy                = 5000',
            '; Output frequency and precision for .xtc file',
            'nstxtcout                = 10000',
            'xtc-precision            = 1000',
            '; This selects the subset of atoms for the .xtc file. You can',
            '; select multiple groups. By default all atoms will be written.',
            'xtc-grps                 =', '; Selection of energy groups',
            'energygrps               = System', '',
            '; NEIGHBORSEARCHING PARAMETERS',
            '; cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)',
            '; nblist update frequency', 'cutoff-scheme            = Verlet',
            'nstlist                  = 100',
            'verlet-buffer-tolerance  = 0.0001',
            '; ns algorithm (simple or grid)',
            'ns_type                  = grid',
            '; Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,',
            '; a value of -1 means: use rlist', '; nblist cut-off',
            'rlist                    = 4.5',
            '; long-range cut-off for switched potentials',
            'rlistlong                = -1', '',
            '; OPTIONS FOR ELECTROSTATICS AND VDW',
            '; Method for doing electrostatics',
            'coulombtype              = cut-off',
            'rcoulomb-switch          = 0', 'rcoulomb                 = 4.5',
            '; Relative dielectric constant for the medium and the reaction field',
            'epsilon-r                = 1', 'epsilon-rf               = 0',
            '; Method for doing Van der Waals',
            'vdw-type                 = Cut-off', '; cut-off lengths',
            'rvdw-switch              = 0', 'rvdw                     = 4.5',
            '; Apply long range dispersion corrections for Energy and Pressure',
            'DispCorr                 = EnerPres',
            '; Extension of the potential lookup tables beyond the cut-off',
            'table-extension          = 1',
            '; Separate tables between energy group pairs',
            'energygrp-table          =', '; IMPLICIT SOLVENT ALGORITHM',
            'implicit-solvent         = No', '',
            '; OPTIONS FOR WEAK COUPLING ALGORITHMS', '; Temperature coupling',
            'tcoupl                   = v-rescale',
            'nsttcouple               = -1', 'nh-chain-length          = 1',
            '; Groups to couple separately',
            'tc-grps                  = System',
            '; Time constant (ps) and reference temperature (K)',
            'tau-t                    = 0.002',
            f'ref-t                    = {self.temperature}',
            '; pressure coupling', 'pcoupl                   = no',
            '; Scaling of reference coordinates, No, All or COM',
            'refcoord-scaling         = COM', '',
            '; GENERATE VELOCITIES FOR STARTUP RUN',
            'gen-vel                  = yes',
            f'gen-temp                 = {self.temperature}',
            'gen-seed                 = 173529', '', '; OPTIONS FOR BONDS',
            f'constraints              = {self.constraints}',
            '; Type of constraint algorithm',
            'constraint-algorithm     = Lincs',
            '; Do not constrain the start configuration',
            'continuation             = yes',
            '; Use successive overrelaxation to reduce the number of shake iterations',
            'Shake-SOR                = no', '; Relative tolerance of shake',
            'shake-tol                = 0.00001',
            '; Highest order in the expansion of the constraint coupling matrix',
            'lincs-order              = 5',
            '; Number of iterations in the final step of LINCS. 1 is fine for',
            '; normal simulations, but use 2 to conserve energy in NVE runs.',
            '; For energy minimization with constraints it should be 4 to 8.',
            'lincs-iter               = 2',
            '; Lincs will write a warning to the stderr if in one step a bond',
            '; rotates over more degrees than',
            'lincs-warnangle          = 50',
            '; Convert harmonic bonds to morse potentials',
            'morse                    = no'
        ]
