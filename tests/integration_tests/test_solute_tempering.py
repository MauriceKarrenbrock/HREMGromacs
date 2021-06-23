# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name
# pylint: disable=wildcard-import
# pylint: disable=unused-wildcard-import
# pylint: disable=no-self-use
#############################################################
# Copyright (c) 2021-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# MIT License                                               #
#############################################################

import os

import HREMGromacs.solute_tempering as _sol


class Testedit_preprocessed_top():
    def test_second_residue_and_1xa(self, get_data_dir, tmp_path):
        input_file = get_data_dir / '4LR6_solvated.top'
        expected_file = get_data_dir / '4LR6_solute_tempering_first_residue_and_1xa.top'
        output_file = tmp_path / 'output.top'

        _sol.edit_preprocessed_top(input_file, output_file, [1])

        with expected_file.open() as ex:
            with output_file.open() as out:
                assert ex.readlines() == out.readlines()


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

            output = _sol.prepare_topologies_for_hrem(top_file=input_file,
                                                      resSeq_to_scale=(1, 2,
                                                                       3),
                                                      mdp_file=test_mdp,
                                                      gro_file=test_gro,
                                                      number_of_replicas=8)

            assert output == expected_output

            for i in expected_output:
                assert i.exists()

            with expected_output[0].open() as f1:
                with expected_output[1].open() as f2:
                    assert f1.read() != f2.read()

        finally:
            os.chdir(old_dir)

            # For some reason gromacs doesn't write this files in tmp_dir
            os.remove('mdout.mdp')
            os.remove('topol.tpr')
