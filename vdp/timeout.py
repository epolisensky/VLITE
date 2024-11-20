"""Adapted from "Function Timeout" in the Python
Decorator Library 
(https://wiki.python.org/moin/PythonDecoratorLibrary#Function_Timeout).

"""
from functools import wraps
import signal


class TimeoutError(Exception):
    """Exception raised when function runs too long."""
    pass

#updated to 60 min
def timeout(seconds=36000, error_message='Function call timed out.'):
    """Creates a decorator which can be used to timeout any
    long-running function.

    Parameters
    ----------
    seconds : float, optional
        Number of seconds a function is permitted to run
        before a TimeoutError exception is raised. Default
        is 3600 seconds (60 minutes).
    error_message : str, optional
        Error message to be printed when the TimeoutError
        exception is raised.

    Returns
    -------
    decorator
        Decorator function which can wrap around any other
        function to terminate if it runs for too long.
    """
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

