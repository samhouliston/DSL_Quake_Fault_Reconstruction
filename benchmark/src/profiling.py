import numpy as np
import cProfile

from main import main

profiler = cProfile.Profile()
profiler.enable()
main()
profiler.disable()
profiler.dump_stats("main_og_boxlib_parallel.stats")
