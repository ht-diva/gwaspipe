import sys
import os

#import gwaslab.hm_harmonize_sumstats
#import gwaslab.io_to_formats
import gwaslab

from .harmonize import fast_checkref
from .format import tofmt

PATCHES = {
    gwaslab.hm_harmonize_sumstats.checkref: fast_checkref,
    gwaslab.io_to_formats.tofmt: tofmt
}

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


def clean_imports(verbose=False):
    """
    Deletes some modules from sys.modules to make sure that the new functions are used.
    """
    # modules_to_delete = find_imports_to_reference('gwaslab.g_Sumstats', verbose=verbose) # gwaslab.g_Sumstats is the module that originally imports gwaslab.hm_harmonize_sumstats.checkref
    # modules_to_delete.append('gwaslab') # very important!
    # modules_to_delete = set([m for m in modules_to_delete if m in sys.modules.keys()])

    # for m in modules_to_delete:
    #     if verbose: print(f"Deleting module {m}")
    #     del sys.modules[m]

    to_keep = [
        'gwaslab.g_Log', 'gwaslab.g_vchange_status', 'gwaslab.bd_config', 'gwaslab.g_version', 'gwaslab.bd_download',
        'gwaslab.bd_common_data', 'gwaslab.qc_check_datatype', 'gwaslab.util_in_fill_data', 'gwaslab.qc_fix_sumstats', 'gwaslab.hm_harmonize_sumstats',
        'gwaslab.io_preformat_input', 'gwaslab.bd_get_hapmap3', 'gwaslab.io_to_formats'
    ]
    if 'gwaslab' in to_keep: to_keep.remove('gwaslab') # make sure we also delete 'gwaslab' from sys.modules

    for k in sys.modules.copy():
        if k.startswith('gwaslab') and k not in to_keep:
            if verbose: print(f"Deleting module {k}")
            del sys.modules[k]


def apply_monkeypatch(verbose=False):
    """
    Overwrite the original functions with the new ones.
    """

    # Clean the imports
    clean_imports(verbose=verbose)

    # Apply the patches
    gwaslab.hm_harmonize_sumstats.checkref = fast_checkref
    gwaslab.io_to_formats.tofmt = tofmt