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

from pathlib import Path
from unittest.mock import MagicMock

import pytest

import HREMGromacs.solute_tempering as _sol


class Testscale_topology_with_plumed():
    def test_ValueError(self):

        with pytest.raises(ValueError):
            _sol.scale_topology_with_plumed('input', 'input', 1.0)

    def test_RuntimeError(self, mocker, tmp_path):

        m_abspath = mocker.patch(
            'PythonAuxiliaryFunctions.path.absolute_programpath',
            return_value='plumed')
        run_output = MagicMock(name='mock_output')
        run_output.returncode.return_value = 1
        m_run = mocker.patch('subprocess.run', return_value=run_output)

        input_top = tmp_path / 'input.top'
        input_top.write_text('a')
        output_top = tmp_path / 'output.top'
        output_top.write_text('a')

        with pytest.raises(RuntimeError):

            _sol.scale_topology_with_plumed(input_top, output_top, 1.0)

        m_abspath.assert_called_once_with('plumed')
        m_run.assert_called_once()


class Testpreprocess_topology():
    def test_works(self, mocker):
        m_run = mocker.patch('PythonAuxiliaryFunctions.run.subprocess_run')
        m_abspath = mocker.patch(
            'PythonAuxiliaryFunctions.path.absolute_programpath',
            return_value='gmx')

        _sol.preprocess_topology(input_top_file='i_top',
                                 output_top_file='o_top',
                                 gro_file='gro',
                                 mdp_file='mdp',
                                 gmx_path='gmx')

        commands = [
            'gmx', 'grompp', '-f', 'mdp', '-c', 'gro', '-p', 'i_top',
            '-maxwarn', '100', '-pp', 'o_top'
        ]

        m_abspath.assert_called_once_with('gmx')
        m_run.assert_called_once_with(commands,
                                      error_string='grompp -pp failed',
                                      shell=False)


class Testgeometrical_progression():
    def test_works(self):
        generator = _sol.geometrical_progression(basis=0.2, denom=7)

        expected_output = (1., 0.7945974047018523, 0.6313850355589192,
                           0.5016969106227039, 0.3986470631277377,
                           0.3167639217533158, 0.25169979012836535, 0.2)

        for expected in expected_output:
            assert next(generator) == expected


class Testprepare_topologies_for_hrem():
    def test_works(self, mocker):
        scaling_values = (1., 0.5, 0.2)

        def _scaling_values_generator():
            for i in scaling_values:
                yield i

        scaling_values_generator = _scaling_values_generator()

        m_preprocess_topology = mocker.patch(
            'HREMGromacs.solute_tempering.preprocess_topology')
        m_edit_preprocessed_top = mocker.patch(
            'HREMGromacs.solute_tempering.edit_preprocessed_top')
        m_geometrical_progression = mocker.patch(
            'HREMGromacs.solute_tempering.geometrical_progression',
            return_value=scaling_values_generator)
        m_scale_topology_with_plumed = mocker.patch(
            'HREMGromacs.solute_tempering.scale_topology_with_plumed')

        input_top = 'test.top'
        expected_output = [
            Path(i).resolve() for i in
            ['test_scaled_0.top', 'test_scaled_1.top', 'test_scaled_2.top']
        ]

        output = _sol.prepare_topologies_for_hrem(top_file=input_top,
                                                  resSeq_to_scale=(1, 2, 3),
                                                  mdp_file='test.mdp',
                                                  gro_file='test.gro',
                                                  number_of_replicas=3)

        assert output == expected_output

        m_preprocess_topology.assert_called_once_with(
            input_top_file=Path(input_top).resolve(),
            output_top_file=Path('TMP_elaborated_top.top').resolve(),
            gro_file=Path('test.gro').resolve(),
            mdp_file=Path('test.mdp').resolve(),
            gmx_path='gmx')

        m_edit_preprocessed_top.assert_called_once_with(
            input_top_file=Path('TMP_elaborated_top.top').resolve(),
            output_top_file=Path('TMP_elaborated_top.top').resolve(),
            resSeq_to_scale=(1, 2, 3))

        m_geometrical_progression.assert_called_once_with(basis=0.2, denom=2)

        for n, i in enumerate(expected_output):
            m_scale_topology_with_plumed.assert_any_call(
                input_top=Path('TMP_elaborated_top.top').resolve(),
                output_top=i,
                scaling_value=scaling_values[n],
                plumed='plumed')
