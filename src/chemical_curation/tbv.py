import click
from chemical_curation import curate

@click.group(name='tbv')
def cli():
    """
    CLI to curate chemical data.
    """
    pass


@click.command(name='curate')
@click.option('-o', '--output_dir', 'output_dir',
              help='The directory in which the output files should be placed')
@click.option('-t', '--target', 'targets', default=["ic50", "ki", "kd"],
              help='Can be invoked multiple times to specify multiple targets',
              show_default=True)
@click.argument('filenames', nargs=-1) # nargs=-1 allows an unspecified number of filenames to be passed
def cli_curate(filenames, output_dir, targets):
    """
    Run the entire curation process end-to-end on every file in FILENAMES.
    """
    curate.main(filenames, output_dir, targets)

    
cli.add_command(cli_curate)
