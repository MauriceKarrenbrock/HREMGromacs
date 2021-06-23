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

import HREMGromacs.multidir_setup as _multi


class Testmake_hrem_multidir_directory():
    def test_works(self, tmp_path):
        output = _multi.make_hrem_multidir_directory(number_of_replicas=8,
                                                     dir_name=tmp_path /
                                                     'BATTERY',
                                                     exist_ok=True)

        assert output == tmp_path / 'BATTERY'
        for i in range(8):
            assert ((tmp_path / 'BATTERY') / f'scaled{i}').exists()


class Testcopy_file_in_hrem_multidir():
    def test_works(self, tmp_path):
        old_dir = os.getcwd()

        try:
            os.chdir(str(tmp_path.resolve()))

            dir_ = _multi.make_hrem_multidir_directory(number_of_replicas=8,
                                                       dir_name=tmp_path /
                                                       'BATTERY',
                                                       exist_ok=True)

            dirs = dir_.glob('scaled*')

            _multi.copy_file_in_hrem_multidir(dir_name=dir_)

            for dd in dirs:
                assert (dd / 'empty_plumed.dat').exists()

        finally:
            os.chdir(old_dir)


class Testmake_multiple_hrem_batteries():
    def test_works(self, tmp_path):
        old_dir = os.getcwd()

        try:
            os.chdir(str(tmp_path.resolve()))

            output_dirs = _multi.make_multiple_hrem_batteries(
                number_of_batteries=12, replicas_per_battery=4)

            expected_output = [tmp_path / f'BATTERY{i}' for i in range(12)]

            assert output_dirs == expected_output

            for dir_ in output_dirs:

                sub_dirs = list(dir_.glob('scaled*'))
                sub_dirs.sort()

                expected_subdirs = [dir_ / f'scaled{i}' for i in range(4)]

                assert sub_dirs == expected_subdirs

                for sub_dir in sub_dirs:
                    assert (sub_dir / 'empty_plumed.dat').exists()

        finally:
            os.chdir(old_dir)
