#! /usr/bin/env python3

import numpy as np
import sys


if __name__=="__main__":
    eps_r = sys.argv[1]
    Z0 = sys.argv[2]

    eps_eff = get_effective_eps(eps_r)
    