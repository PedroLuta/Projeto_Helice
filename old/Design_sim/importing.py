import numpy as nup
import math
import matplotlib.pyplot as plt
import random
from numba import njit, jit
import time
import os
import subprocess as sp
from math import inf as infinite
#from smt.sampling_methods import LHS
import json

from numpy.core.fromnumeric import size
from numpy.linalg.linalg import LinAlgError

pi = nup.pi
euler = nup.e