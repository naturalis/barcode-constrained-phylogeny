import logging


def get_formatted_logger(name, level):
    # Create a logger
    logger = logging.getLogger(name)

    # Set the log level
    logger.setLevel(level)

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
    logger.addHandler(handler)

    return logger
