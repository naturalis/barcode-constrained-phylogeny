import logging


def get_formatted_logger(name, config):
    """
    Instantiates a python `logging` object for the specified name and verbosity level.
    The function further configures this object so that the logging message it emits
    are formatted to include a time stamp.
    :param name: name for the logger
    :param config: bactria.Config object
    :return:
    """
    # Create a logger
    requested_logger = logging.getLogger(name)

    # Set the log level
    requested_logger.setLevel(config.get('log_level'))

    # Create a console handler
    handler = logging.StreamHandler()

    # Define the format for the log messages,
    # including the date and time stamp
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    # Set the date format
    date_format = '%Y-%m-%d %H:%M:%S'

    # Create a formatter using the specified format and date format
    formatter = logging.Formatter(log_format, datefmt=date_format)

    # Set the formatter for the handler
    handler.setFormatter(formatter)

    # Add the handler to the logger
    requested_logger.addHandler(handler)

    return requested_logger

