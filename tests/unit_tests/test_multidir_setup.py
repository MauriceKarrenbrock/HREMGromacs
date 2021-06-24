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

import HREMGromacs.multidir_setup as _multi


class Testmake_multidir_mdrun_string_for_hrem():
    def test_works(self):
        input_list = [
            'BATTERY0/scaled0', 'BATTERY0/scaled1', 'BATTERY0/scaled2',
            'BATTERY0/scaled3'
        ]
        output = _multi.make_multidir_mdrun_string_for_hrem(
            multidir_directories=input_list)

        expected_output = (
            'gmx_mpi mdrun -v -plumed empty_plumed.dat -replex 100 -hrex -dlb no'
            ' -s HREM.tpr -deffnm HREM'
            ' -multidir BATTERY0/scaled0 BATTERY0/scaled1 BATTERY0/scaled2 BATTERY0/scaled3'
        )

        assert output == expected_output
