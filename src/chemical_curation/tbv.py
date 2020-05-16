import click
from chemical_curation import curate

@click.group(name='tbv')
def cli():
    """
    CLI to curate chemical data using the Trust But Verify process.
    """
    pass


@click.command(name='curate')
@click.option('-o', '--output-dir', 'output_dir',
              help='The directory in which the output files should be placed')
@click.option('-t', '--target', 'targets', default=["ic50", "ki", "kd"],
              help='Can be invoked multiple times to specify multiple targets',
              show_default=True)
@click.option('-r', '--review-threshold', 'review_threshold',
              help='', default=10, show_default=True)
@click.option('-v', '--verbosity', 'verbosity',
              help='Level of debugging to show', default='INFO',
              type=click.Choice(['DEBUG','INFO','WARNING','ERROR','CRITICAL']))
@click.argument('filenames', nargs=-1) # nargs=-1 allows any number of filenames to be passed
def cli_curate(filenames, output_dir, targets, review_threshold, verbosity):
    """
    Run the entire curation process end-to-end on every file in FILENAMES.
    """
    
    curate.main(filenames, output_dir, targets, review_threshold, verbosity)

    
cli.add_command(cli_curate)
