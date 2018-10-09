#!/usr/bin/env python
"""
HISTORY
2018-10-09 modified by T.Enoto
2018-10-01 created by T.Enoto 
"""

import click 
import giantradiopulse.radio 

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
	if ctx.invoked_subcommand is None:
		print(ctx.get_help())

@cli.command(help="Convert GTI text file to fitsfile.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outfitsfile", type=click.Path(), default=None)
def convert_gti_txt2fits(file_path,outfitsfile):
	gtifile = giantradiopulse.radio.open_gti(file_path)
	gtifile.writeAsFitsFormat(outfitsfile=outfitsfile)

@cli.command(help="Convert GiantRadioPulse list text file to fitsfile.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outfitsfile", type=click.Path(), default=None)
def convert_grp_txt2fits(file_path,outfitsfile):
	grpfile = giantradiopulse.radio.open_gitantradiopulse(file_path)
	grpfile.writeAsFitsFormat(outfitsfile=outfitsfile)

def main():
	cli()

if __name__ == "__main__":
	main()
