import sys
import os

import gwaslab.hm_harmonize_sumstats
from .harmonize import fast_checkref


def find_imports_to_reference(reference_module, verbose=False):
    """
    Find all the modules that import the reference_module and return a list of them.
    """

    out = [reference_module]

    # Define the directory of the gwaslab package
    gwaslab_directory = os.path.dirname(gwaslab.__file__)

    # Walk through the directory and find all python files that are not __init__.py
    # and check if they import 
    for root, dirs, files in os.walk(gwaslab_directory):
        for filename in files:
            if filename.endswith('.py') and filename != '__init__.py':
                filepath = os.path.join(root, filename)
                with open(filepath, 'r') as file:
                    for line in file:
                        if f'from {reference_module}' in line:
                            if verbose: print(f'Found reference to {reference_module} in {filepath}: {line.strip()}')
                            out_el = filepath.replace(gwaslab_directory, '').replace(os.sep, '.').strip('.').strip('.py')
                            out_el = f'gwaslab.{out_el}'
                            out.append(out_el)
                            break

    return list(set(out))


def apply_monkeypatch(verbose=False):
    """
    Overwrite the original functions with the new ones.
    """

    # We need to delete some modules from sys.modules to make sure that the new function is used
    modules_to_delete = find_imports_to_reference('gwaslab.g_Sumstats', verbose=verbose) # gwaslab.g_Sumstats is the module that originally imports gwaslab.hm_harmonize_sumstats.checkref
    modules_to_delete.append('gwaslab') # very important!
    modules_to_delete = [m for m in modules_to_delete if m in sys.modules.keys()]

    for m in modules_to_delete:
        if verbose: print(f"Deleting module {m}")
        del sys.modules[m]

    # Overwrite original functions with the new ones
    gwaslab.hm_harmonize_sumstats.checkref = fast_checkref