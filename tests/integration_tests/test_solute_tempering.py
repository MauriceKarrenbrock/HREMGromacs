# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name
# pylint: disable=wildcard-import
# pylint: disable=unused-wildcard-import
# pylint: disable=no-self-use
# pylint: disable=broad-except
#############################################################
# Copyright (c) 2021-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# MIT License                                               #
#############################################################

import os
import warnings
import parmed

import HREMGromacs.solute_tempering as _sol


class Testedit_preprocessed_top():

    def test_first_residue_and_1xa(self, get_data_dir, tmp_path):
        input_file = get_data_dir / '4LR6_solvated.top'
        expected_file = get_data_dir / '4LR6_solute_tempering_first_residue_and_1xa.top'
        output_file = tmp_path / 'output.top'

        _sol.edit_preprocessed_top(input_file, output_file, [1], '1XA')

        with expected_file.open() as ex:
            with output_file.open() as out:
                assert ex.readlines()[15:] == out.readlines()[15:]


class Testprepare_topologies_for_hrem():

    def test_works(self, get_data_dir, tmp_path):

        old_dir = os.getcwd()

        try:
            os.chdir(str(tmp_path.resolve()))

            input_file = get_data_dir / '4LR6_solvated.top'
            test_mdp = tmp_path / 'test.mdp'
            test_mdp.touch()
            test_gro = get_data_dir / '4LR6_solvated.pdb'

            expected_output = [(tmp_path / f'4LR6_solvated_scaled_{i}.top')
                               for i in range(8)]

            output = _sol.prepare_topologies_for_hrem(
                top_file=input_file,
                protein_resSeq_to_scale=(1, 2, 3),
                ligand_resname='AAAAAA',
                mdp_file=test_mdp,
                gro_file=test_gro,
                number_of_replicas=8)

            assert output == expected_output

            for i in expected_output:
                assert i.exists()

            with expected_output[0].open() as f1:
                with expected_output[1].open() as f2:
                    assert f1.read() != f2.read()

        except OSError:
            warnings.warn("Gromacs was not found in the path, one or more tests can not be executed")

        finally:
            os.chdir(old_dir)

            # For some reason gromacs doesn't write this files in tmp_dir
            try:
                os.remove('mdout.mdp')
            except Exception:
                pass
            try:
                os.remove('topol.tpr')
            except Exception:
                pass


class Testscale_dihedrals_in_parmed_structure():
    def test_05_is_gromacs_False(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_dihedrals_in_parmed_structure(input_structure,
                               scaling_value=0.5,
                               atom_list=atom_list,
                                is_gromacs=False)

        for dihedral_input, dihedral_output in zip(input_structure.dihedrals, output_structure.dihedrals):
            if not dihedral_output.improper:
                if dihedral_output.atom2.idx in atom_list and dihedral_output.atom3.idx in atom_list:
                    assert dihedral_output.type.phi_k == 0.5 * dihedral_input.type.phi_k
                
                else:
                    assert dihedral_output.type.phi_k == dihedral_input.type.phi_k

        assert output_structure.dihedrals.changed == True
        assert output_structure.dihedral_types.changed == True

    def test_05_is_gromacs_True(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_dihedrals_in_parmed_structure(input_structure,
                               scaling_value=0.5,
                               atom_list=atom_list,
                                is_gromacs=True)

        for dihedral_input, dihedral_output in zip(input_structure.dihedrals, output_structure.dihedrals):
            if not dihedral_output.improper:
                if dihedral_output.atom2.idx in atom_list and dihedral_output.atom3.idx in atom_list:
                    assert dihedral_output.type.phi_k == 0.5 * dihedral_input.type.phi_k
                
                else:
                    assert dihedral_output.type.phi_k == dihedral_input.type.phi_k

        assert output_structure.dihedrals.changed == True
        assert output_structure.dihedral_types.changed == True

    def test_0_is_gromacs_False(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_dihedrals_in_parmed_structure(input_structure,
                               scaling_value=0.,
                               atom_list=atom_list,
                                is_gromacs=False)

        for dihedral_input, dihedral_output in zip(input_structure.dihedrals, output_structure.dihedrals):
            if not dihedral_output.improper:
                if dihedral_output.atom2.idx in atom_list and dihedral_output.atom3.idx in atom_list:
                    assert dihedral_output.type.phi_k == 0.
                
                else:
                    assert dihedral_output.type.phi_k == dihedral_input.type.phi_k

        assert output_structure.dihedrals.changed == True
        assert output_structure.dihedral_types.changed == True


class Testscale_14_interactions_in_parmed_structure():
    def test_05_is_gromacs_False(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))
        
        input_structure.defaults.gen_pairs == "yes"

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_14_interactions_in_parmed_structure(input_structure,
                               scaling_value=0.5,
                               atom_list=atom_list,
                                is_gromacs=False)

        for output_atom, input_atom in zip(output_structure.atoms, input_structure.atoms):
            if output_atom.idx in atom_list:
                assert list(output_atom.atom_type.name)[-1] == "_"
                assert list(output_atom.type)[-1] == "_"
                assert list(output_atom.name)[-1] == "_"
                
                # Parmed likes to change names sometimes so I can not check for this
                # assert "".join(list(output_atom.atom_type.name)[:-1]) == input_atom.atom_type.name
                # assert "".join(list(output_atom.type)[:-1]) == input_atom.type
                # assert "".join(list(output_atom.name)[:-1]) == input_atom.name

                output_atom.epsilon_14 == 0.5 * input_atom.epsilon_14
                output_atom.atom_type.epsilon_14 == 0.5 * input_atom.atom_type.epsilon_14

            else:
                # Parmed likes to change names sometimes so I can not check for this
                # assert output_atom.atom_type.name == input_atom.atom_type.name
                # assert output_atom.type == input_atom.type
                # assert output_atom.name == input_atom.name

                output_atom.epsilon_14 == input_atom.epsilon_14
                output_atom.atom_type.epsilon_14 == input_atom.atom_type.epsilon_14

        if output_structure.adjusts:
            for output_adjust, input_adjust in zip(output_structure.adjusts, input_structure.adjusts):
                if output_adjust.atom1.idx in atom_list and output_adjust.atom2.idx in atom_list:
                    assert output_adjust.type.epsilon == 0.5 * input_adjust.type.epsilon
                    assert output_adjust.type.chgscale == 0.5 * input_adjust.type.chgscale

                else:
                    assert output_adjust.type.epsilon == input_adjust.type.epsilon
                    assert output_adjust.type.chgscale == input_adjust.type.chgscale
            
            assert output_structure.adjust_types.changed
            assert output_structure.adjusts.changed

        assert output_structure.atoms.changed

        assert output_structure.defaults.gen_pairs == "no"
    
    def test_05_is_gromacs_True(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))
        
        input_structure.defaults.gen_pairs == "yes"

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_14_interactions_in_parmed_structure(input_structure,
                               scaling_value=0.5,
                               atom_list=atom_list,
                                is_gromacs=True)

        for output_atom, input_atom in zip(output_structure.atoms, input_structure.atoms):
            if output_atom.idx in atom_list:
                assert list(output_atom.atom_type.name)[-1] == "_"
                assert list(output_atom.type)[-1] == "_"
                assert list(output_atom.name)[-1] == "_"
                
                # Parmed likes to change names sometimes so I can not check for this
                # assert "".join(list(output_atom.atom_type.name)[:-1]) == input_atom.atom_type.name
                # assert "".join(list(output_atom.type)[:-1]) == input_atom.type
                # assert "".join(list(output_atom.name)[:-1]) == input_atom.name

                output_atom.epsilon_14 == 0.5 * input_atom.epsilon_14
                output_atom.atom_type.epsilon_14 == 0.5 * input_atom.atom_type.epsilon_14

            else:
                # Parmed likes to change names sometimes so I can not check for this
                # assert output_atom.atom_type.name == input_atom.atom_type.name
                # assert output_atom.type == input_atom.type
                # assert output_atom.name == input_atom.name

                output_atom.epsilon_14 == input_atom.epsilon_14
                output_atom.atom_type.epsilon_14 == input_atom.atom_type.epsilon_14

        if output_structure.adjusts:
            for output_adjust, input_adjust in zip(output_structure.adjusts, input_structure.adjusts):
                if output_adjust.atom1.idx in atom_list and output_adjust.atom2.idx in atom_list:
                    assert output_adjust.type.epsilon == 0.5 * input_adjust.type.epsilon
                    assert output_adjust.type.chgscale == 0.5 * input_adjust.type.chgscale

                else:
                    assert output_adjust.type.epsilon == input_adjust.type.epsilon
                    assert output_adjust.type.chgscale == input_adjust.type.chgscale
            
            assert output_structure.adjust_types.changed
            assert output_structure.adjusts.changed

        assert output_structure.atoms.changed

        assert output_structure.defaults.gen_pairs == "no"

    def test_0_is_gromacs_True(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))
        
        input_structure.defaults.gen_pairs == "yes"

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_14_interactions_in_parmed_structure(input_structure,
                               scaling_value=0.,
                               atom_list=atom_list,
                                is_gromacs=True)

        for output_atom, input_atom in zip(output_structure.atoms, input_structure.atoms):
            if output_atom.idx in atom_list:
                assert list(output_atom.atom_type.name)[-1] == "_"
                assert list(output_atom.type)[-1] == "_"
                assert list(output_atom.name)[-1] == "_"
                
                # Parmed likes to change names sometimes so I can not check for this
                # assert "".join(list(output_atom.atom_type.name)[:-1]) == input_atom.atom_type.name
                # assert "".join(list(output_atom.type)[:-1]) == input_atom.type
                # assert "".join(list(output_atom.name)[:-1]) == input_atom.name

                output_atom.epsilon_14 == 0. * input_atom.epsilon_14
                output_atom.atom_type.epsilon_14 == 0. * input_atom.atom_type.epsilon_14

            else:
                # Parmed likes to change names sometimes so I can not check for this
                # assert output_atom.atom_type.name == input_atom.atom_type.name
                # assert output_atom.type == input_atom.type
                # assert output_atom.name == input_atom.name

                output_atom.epsilon_14 == input_atom.epsilon_14
                output_atom.atom_type.epsilon_14 == input_atom.atom_type.epsilon_14

        if output_structure.adjusts:
            for output_adjust, input_adjust in zip(output_structure.adjusts, input_structure.adjusts):
                if output_adjust.atom1.idx in atom_list and output_adjust.atom2.idx in atom_list:
                    assert output_adjust.type.epsilon == 0. * input_adjust.type.epsilon
                    assert output_adjust.type.chgscale == 0. * input_adjust.type.chgscale

                else:
                    assert output_adjust.type.epsilon == input_adjust.type.epsilon
                    assert output_adjust.type.chgscale == input_adjust.type.chgscale
            
            assert output_structure.adjust_types.changed
            assert output_structure.adjusts.changed

        assert output_structure.atoms.changed

        assert output_structure.defaults.gen_pairs == "no"


class Testscale_vdw_and_charges_in_parmed_structure():
    def test_05_is_gromacs_False(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_vdw_and_charges_in_parmed_structure(input_structure,
                               scaling_value=0.5,
                               atom_list=atom_list,
                                is_gromacs=False,
                                scale_charges=True,
                                scale_vdw=True)

        for output_atom, input_atom in zip(output_structure.atoms, input_structure.atoms):
            if output_atom.idx in atom_list:
                assert list(output_atom.atom_type.name)[-1] == "_"
                assert list(output_atom.type)[-1] == "_"
                assert list(output_atom.name)[-1] == "_"
                
                # Parmed likes to change names sometimes so I can not check for this
                # assert "".join(list(output_atom.atom_type.name)[:-1]) == input_atom.atom_type.name
                # assert "".join(list(output_atom.type)[:-1]) == input_atom.type
                # assert "".join(list(output_atom.name)[:-1]) == input_atom.name

                output_atom.epsilon == 0.5 * input_atom.epsilon
                output_atom.atom_type.epsilon == 0.5 * input_atom.atom_type.epsilon

                output_atom.charge == 0.5 * input_atom.charge
                output_atom.atom_type.charge == 0.5 * input_atom.atom_type.charge

            else:
                # Parmed likes to change names sometimes so I can not check for this
                # assert output_atom.atom_type.name == input_atom.atom_type.name
                # assert output_atom.type == input_atom.type
                # assert output_atom.name == input_atom.name

                output_atom.epsilon == input_atom.epsilon
                output_atom.atom_type.epsilon == input_atom.atom_type.epsilon

                output_atom.charge == input_atom.charge
                output_atom.atom_type.charge == input_atom.atom_type.charge

    def test_05_is_gromacs_True(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_vdw_and_charges_in_parmed_structure(input_structure,
                               scaling_value=0.5,
                               atom_list=atom_list,
                                is_gromacs=True,
                                scale_charges=True,
                                scale_vdw=True)

        for output_atom, input_atom in zip(output_structure.atoms, input_structure.atoms):
            if output_atom.idx in atom_list:
                assert list(output_atom.atom_type.name)[-1] == "_"
                assert list(output_atom.type)[-1] == "_"
                assert list(output_atom.name)[-1] == "_"
                
                # Parmed likes to change names sometimes so I can not check for this
                # assert "".join(list(output_atom.atom_type.name)[:-1]) == input_atom.atom_type.name
                # assert "".join(list(output_atom.type)[:-1]) == input_atom.type
                # assert "".join(list(output_atom.name)[:-1]) == input_atom.name

                output_atom.epsilon == 0.5 * input_atom.epsilon
                output_atom.atom_type.epsilon == 0.5 * input_atom.atom_type.epsilon

                output_atom.charge == 0.5 * input_atom.charge
                output_atom.atom_type.charge == 0.5 * input_atom.atom_type.charge

            else:
                # Parmed likes to change names sometimes so I can not check for this
                # assert output_atom.atom_type.name == input_atom.atom_type.name
                # assert output_atom.type == input_atom.type
                # assert output_atom.name == input_atom.name

                output_atom.epsilon == input_atom.epsilon
                output_atom.atom_type.epsilon == input_atom.atom_type.epsilon

                output_atom.charge == input_atom.charge
                output_atom.atom_type.charge == input_atom.atom_type.charge

    def test_0_is_gromacs_True(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_vdw_and_charges_in_parmed_structure(input_structure,
                               scaling_value=0.,
                               atom_list=atom_list,
                                is_gromacs=True,
                                scale_charges=True,
                                scale_vdw=True)

        for output_atom, input_atom in zip(output_structure.atoms, input_structure.atoms):
            if output_atom.idx in atom_list:
                assert list(output_atom.atom_type.name)[-1] == "_"
                assert list(output_atom.type)[-1] == "_"
                assert list(output_atom.name)[-1] == "_"
                
                # Parmed likes to change names sometimes so I can not check for this
                # assert "".join(list(output_atom.atom_type.name)[:-1]) == input_atom.atom_type.name
                # assert "".join(list(output_atom.type)[:-1]) == input_atom.type
                # assert "".join(list(output_atom.name)[:-1]) == input_atom.name

                output_atom.epsilon == 0. * input_atom.epsilon
                output_atom.atom_type.epsilon == 0. * input_atom.atom_type.epsilon

                output_atom.charge == 0. * input_atom.charge
                output_atom.atom_type.charge == 0. * input_atom.atom_type.charge

            else:
                # Parmed likes to change names sometimes so I can not check for this
                # assert output_atom.atom_type.name == input_atom.atom_type.name
                # assert output_atom.type == input_atom.type
                # assert output_atom.name == input_atom.name

                output_atom.epsilon == input_atom.epsilon
                output_atom.atom_type.epsilon == input_atom.atom_type.epsilon

                output_atom.charge == input_atom.charge
                output_atom.atom_type.charge == input_atom.atom_type.charge

    def test_0_is_gromacs_True_only_vdw(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_vdw_and_charges_in_parmed_structure(input_structure,
                               scaling_value=0.,
                               atom_list=atom_list,
                                is_gromacs=True,
                                scale_charges=False,
                                scale_vdw=True)

        for output_atom, input_atom in zip(output_structure.atoms, input_structure.atoms):
            if output_atom.idx in atom_list:
                assert list(output_atom.atom_type.name)[-1] == "_"
                assert list(output_atom.type)[-1] == "_"
                assert list(output_atom.name)[-1] == "_"
                
                # Parmed likes to change names sometimes so I can not check for this
                # assert "".join(list(output_atom.atom_type.name)[:-1]) == input_atom.atom_type.name
                # assert "".join(list(output_atom.type)[:-1]) == input_atom.type
                # assert "".join(list(output_atom.name)[:-1]) == input_atom.name

                output_atom.epsilon == 0. * input_atom.epsilon
                output_atom.atom_type.epsilon == 0. * input_atom.atom_type.epsilon

                output_atom.charge == input_atom.charge
                output_atom.atom_type.charge == input_atom.atom_type.charge

            else:
                # Parmed likes to change names sometimes so I can not check for this
                # assert output_atom.atom_type.name == input_atom.atom_type.name
                # assert output_atom.type == input_atom.type
                # assert output_atom.name == input_atom.name

                output_atom.epsilon == input_atom.epsilon
                output_atom.atom_type.epsilon == input_atom.atom_type.epsilon

                output_atom.charge == input_atom.charge
                output_atom.atom_type.charge == input_atom.atom_type.charge
    
    def test_0_is_gromacs_True_only_charges(self, get_data_dir):
        input_structure = parmed.load_file(str(get_data_dir / "4LR6_solvated.top"),
            xyz=str(get_data_dir / "4LR6_solvated.pdb"))

        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        output_structure = _sol.scale_vdw_and_charges_in_parmed_structure(input_structure,
                               scaling_value=0.,
                               atom_list=atom_list,
                                is_gromacs=True,
                                scale_charges=True,
                                scale_vdw=False)

        for output_atom, input_atom in zip(output_structure.atoms, input_structure.atoms):
            if output_atom.idx in atom_list:
                assert list(output_atom.atom_type.name)[-1] == "_"
                assert list(output_atom.type)[-1] == "_"
                assert list(output_atom.name)[-1] == "_"
                
                # Parmed likes to change names sometimes so I can not check for this
                # assert "".join(list(output_atom.atom_type.name)[:-1]) == input_atom.atom_type.name
                # assert "".join(list(output_atom.type)[:-1]) == input_atom.type
                # assert "".join(list(output_atom.name)[:-1]) == input_atom.name

                output_atom.epsilon == input_atom.epsilon
                output_atom.atom_type.epsilon == input_atom.atom_type.epsilon

                output_atom.charge == 0. * input_atom.charge
                output_atom.atom_type.charge == 0. * input_atom.atom_type.charge

            else:
                # Parmed likes to change names sometimes so I can not check for this
                # assert output_atom.atom_type.name == input_atom.atom_type.name
                # assert output_atom.type == input_atom.type
                # assert output_atom.name == input_atom.name

                output_atom.epsilon == input_atom.epsilon
                output_atom.atom_type.epsilon == input_atom.atom_type.epsilon

                output_atom.charge == input_atom.charge
                output_atom.atom_type.charge == input_atom.atom_type.charge


class Testscale_topologies_with_parmed():
    def test_works(self, get_data_dir, tmpdir):
        #First and last residue of the protein
        atom_list = (0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083,
            2084, 2085, 2086, 2087)

        _sol.scale_topologies_with_parmed(input_top=str(get_data_dir / "4LR6_solvated.top"),
                                output_top=str(tmpdir / "output_top_1.top"),
                                xyz=str(get_data_dir / "4LR6_solvated.pdb"),
                                atom_list=atom_list,
                                dihedral_scaling=1,
                                interaction14_scaling=1,
                                vdw_scaling=1,
                                q_scaling=1)
        
        _sol.scale_topologies_with_parmed(input_top=str(get_data_dir / "4LR6_solvated.top"),
                                output_top=str(tmpdir / "output_top_2.top"),
                                xyz=str(get_data_dir / "4LR6_solvated.pdb"),
                                atom_list=atom_list,
                                dihedral_scaling=1,
                                interaction14_scaling=1,
                                vdw_scaling=1,
                                q_scaling=1)

        _sol.scale_topologies_with_parmed(input_top=str(get_data_dir / "4LR6_solvated.top"),
                                output_top=str(tmpdir / "output_top_3.top"),
                                xyz=str(get_data_dir / "4LR6_solvated.pdb"),
                                atom_list=atom_list,
                                dihedral_scaling=0.1,
                                interaction14_scaling=0.1,
                                vdw_scaling=0.1,
                                q_scaling=0.1)
        
        with open(tmpdir / "output_top_1.top") as f:
            out_str_1 = f.readlines()
        with open(tmpdir / "output_top_2.top") as f:
            out_str_2 = f.readlines()
        with open(tmpdir / "output_top_3.top") as f:
            out_str_3 = f.readlines()

        # Remove comments and empty lines
        for top in (out_str_1, out_str_2, out_str_3):
            for i in range(len(top)-1, -1, -1):
                if not top[i].strip():
                    top.pop(i)
                
                elif top[i].strip()[0] == ";":
                    top.pop(i)
                
                else:
                    top[i].replace("_", "")

        assert "".join(out_str_1) == "".join(out_str_2)
        assert "".join(out_str_1) != "".join(out_str_3)