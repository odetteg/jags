import sys
from pathlib import Path
from constants.common import *

rule get_map_cmds:
    output:
        str(cmds_dir) + "/map.txt"
    run:
        bwa_cmds(ref=ref_genome)