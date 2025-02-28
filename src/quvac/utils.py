"""
Useful generic utilities
"""

import pkgutil
import importlib
import inspect
import os
import platform
import resource
import zipfile
from contextlib import contextmanager
from pathlib import Path
import shutil

import numpy as np
import pyfftw
import yaml


def read_yaml(yaml_file):
    with open(yaml_file, "r") as stream:
        try:
            data = yaml.safe_load(stream)
            return data
        except yaml.YAMLError as exc:
            print(exc)
            return exc


def write_yaml(yaml_file, data):
    with open(yaml_file, "w") as outfile:
        yaml.dump(data, outfile, default_flow_style=False)


def format_time(seconds):
    days, seconds = divmod(seconds, 86400)
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    out_str = [
        f"{days:.0f} days" * bool(days),
        f"{hours:.0f} h" * bool(hours),
        f"{minutes:.0f} min" * bool(minutes),
        f"{seconds:.2f} s",
    ]
    return " ".join(out_str)


def format_memory(mem):
    """
    mem: float
        Memory in KB (kilobyte)
    """
    units = "KB MB GB TB".split()
    idx = 0
    while mem > 1024:
        mem /= 1024
        idx += 1
    return f"{mem:.2f} {units[idx]}"


def save_wisdom(ini_file, wisdom_file=None, add_host_name=False):
    if wisdom_file is None:
        wisdom_path = os.path.dirname(ini_file)
        if not os.path.exists(wisdom_path):
            Path(wisdom_path).mkdir(parents=True, exist_ok=True)
        wisdom_name = "fftw-wisdom"
        if add_host_name:
            wisdom_name += platform.node()
        wisdom_file = os.path.join(wisdom_path, wisdom_name)
    else:
        wisdom_path = os.path.dirname(wisdom_file)
        if not os.path.exists(wisdom_path):
            Path(wisdom_path).mkdir(parents=True, exist_ok=True)
    wisdom = pyfftw.export_wisdom()
    with open(wisdom_file, "wb") as f:
        f.write(b"\n".join(wisdom))


def load_wisdom(wisdom_file):
    with open(wisdom_file, "rb") as f:
        wisdom = f.read()
    return tuple(wisdom.split(b"\n"))


def get_maxrss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


def zip_directory_shutil(directory_path, output_path):
    shutil.make_archive(output_path, 'zip', directory_path)


def find_classes_in_package(package_name):
    """Find all class names in a given package."""
    classes = []
    
    # Import the package
    package = importlib.import_module(package_name)
    
    # Recursively find all modules in the package
    for _, module_name, _ in pkgutil.walk_packages(package.__path__, package.__name__ + "."):
        try:
            module = importlib.import_module(module_name)
            
            # Inspect module members and find classes
            for name, obj in inspect.getmembers(module, inspect.isclass):
                # Ensure class belongs to the module (not an imported one)
                if obj.__module__ == module_name:
                    classes.append(f"{module_name}.{name}")
        
        except Exception as e:
            print(f"Skipping {module_name} due to error: {e}")

    return classes
