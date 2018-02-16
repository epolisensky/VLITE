class Error(Exception):
    """Base class for all exceptions."""
    pass


class ConfigError(Error):
    """Raised when there is an error in the configuration file."""
    pass
