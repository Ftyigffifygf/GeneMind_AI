# genemind_ai/app/modules/performance_utils.py

from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor

# Caching decorator for functions that return the same output for the same input
def cached(maxsize=128):
    def decorator(func):
        return lru_cache(maxsize=maxsize)(func)
    return decorator

# Thread pool for parallel execution of CPU-bound tasks
executor = ThreadPoolExecutor()

def run_in_executor(func, *args, **kwargs):
    """Runs a function in a thread pool executor."""
    return executor.submit(func, *args, **kwargs)


