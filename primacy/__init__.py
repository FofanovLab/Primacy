import click

def error(msg, logger=False):
    """Prints an error message to stderr and logs."""
    click.secho(msg, fg='red', err=True)
    if logger:
        logger.error(msg)


def warn(msg, logger=False):
    '''Prints a warning message to stderr.'''
    click.secho(msg, fg='yellow')
    if logger:
        logger.warning(msg)
    


def info(msg, logger=False):
    click.secho(msg, fg='green')
    if logger:
        logger.info(msg)