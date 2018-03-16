"""Creates a decorator which can be used to timeout any
long-running function.

Adapted from the "Function Timeout" in the Python Decorator
Library([1]_).

.. [1] https://wiki.python.org/moin/PythonDecoratorLibrary#Function_Timeout

"""
from functools import wraps
import signal


class TimeoutError(Exception):
    pass

def timeout(seconds=300, error_message='Function call timed out.'):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        @wraps(func)
        def wrapper(*args, **kwargs):
            # Set the signal handler and default alarm time of 300 s
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                # Disable the alarm
                signal.alarm(0)
            return result

        return wrapper

    return decorator

