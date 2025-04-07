import os
from pathlib import Path

import numpy as np
import pytest

from tests.config_for_tests import DEFAULT_CONFIG_PATH, SIMULATION_SCRIPT
from quvac.simulation import quvac_simulation
from quvac.utils import read_yaml, write_yaml


@pytest.fixture(scope="session")
def get_tmp_path(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("simulation")
    return tmp_path


def save_ini(path, ini_data):
    path = os.path.join(str(path), "test_simulation")
    Path(path).mkdir(parents=True, exist_ok=True)

    ini_file = os.path.join(path, "ini.yml")
    write_yaml(ini_file, ini_data)
    return ini_file


def test_simulation(tmp_path):
    ini_data = read_yaml(DEFAULT_CONFIG_PATH)
    ini_file = save_ini(tmp_path, ini_data)
    quvac_simulation(ini_file)

    ini_data["mode"] = "postprocess"
    ini_file = save_ini(tmp_path, ini_data)
    quvac_simulation(ini_file)


def test_simulation_test(tmp_path):
    ini_data = read_yaml(DEFAULT_CONFIG_PATH)
    ini_data["performance"]["test_run"] = True
    ini_file = save_ini(tmp_path, ini_data)
    quvac_simulation(ini_file)


def test_channels(tmp_path):
    ini_data = read_yaml(DEFAULT_CONFIG_PATH)
    ini_data["integrator"]["type"] = "vacuum_emission_channels"
    ini_file = save_ini(tmp_path, ini_data)
    quvac_simulation(ini_file)


def test_precision(tmp_path):
    ini_data = read_yaml(DEFAULT_CONFIG_PATH)
    ini_data["performance"]["precision"] = "float32"
    ini_data["postprocess"]["perp_polarization_type"] = None
    ini_file = save_ini(tmp_path, ini_data)
    quvac_simulation(ini_file)


