#!/usr/bin/env python3
import pstats
from pstats import SortKey
p = pstats.Stats('profile.out')
p.strip_dirs().sort_stats('cumulative').print_stats(10)
