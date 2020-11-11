import os
import shlex
from subprocess import check_call
from pathlib import Path

from random import sample


def shell_command(command):
    check_call(shlex.split(command_string))

def dock(receptor_path, molec_path, config, subsample_perc=0.01):
    # prepare paths
    

    # random subsampling

    # dock every molecule