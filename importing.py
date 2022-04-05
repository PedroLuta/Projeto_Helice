import numpy as nup
import math
import matplotlib.pyplot as plt
import random
random.seed()
from numba import njit, jit
import time
import os
import subprocess as sp
from math import inf as infinite
#import concurrent.futures
#from smt.sampling_methods import LHS
import json
import sys
from math import inf as infinite

root_dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root_dir_path + r'\Libraries\Airfoils_library')
sys.path.append(root_dir_path + r'\Libraries\Optimization')
sys.path.append(root_dir_path + r'\Libraries\Prop_design')
sys.path.append(root_dir_path + r'\Libraries\Simulation')
sys.path.append(root_dir_path + r'\Libraries\Auxiliary')

import Auxiliary

from numpy.core.fromnumeric import size
from numpy.linalg.linalg import LinAlgError

pi = nup.pi
euler = nup.e



