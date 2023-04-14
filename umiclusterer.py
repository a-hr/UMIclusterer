import click

from src import UMIclusterer


@click.command()
@click.argument(
    "bam",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--target_regions",
    "-t",
    type=click.Path(exists=True),
    help="Path to the file containing the target regions.",
    required=True,
)
@click.option(
    "--saf",
    "-s",
    is_flag=True,
    help="Target regions are in SAF format.",
    required=False,
)
@click.option(
    "--debug",
    "-d",
    is_flag=True,
    help="Enables debug mode.",
    required=False,
)
def main(bam, target_regions, saf, debug):
    UMIclusterer(bam, target_regions, saf, debug)


if __name__ == "__main__":
    main()
