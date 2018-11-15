#!/usr/bin/env python
"""
HISTORY
2018-10-24 created by T.Enoto 
"""

import click 
import giantradiopulse.process

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
	if ctx.invoked_subcommand is None:
		print(ctx.get_help())

@cli.command(help="Prepare X-ray and radio fitsfiles for barycentric and add phase.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outdir", type=click.Path(), default='out/crabgrp')
def prepare_datafiles(file_path,outdir):
	process_manager = giantradiopulse.process.ProcessManager(
		file_path,outdir=outdir)
	process_manager.prepare_datafiles()

@cli.command(help="Run correlation study.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outdir", type=click.Path(), default='out/crabgrp')
def run_correlation_study(file_path,outdir):
	process_manager = giantradiopulse.process.ProcessManager(
		file_path,outdir=outdir)
	process_manager.run_correlation_study()

@cli.command(help="Add observation.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outdir", type=click.Path(), default='out/crabgrp')
def add_observations(file_path,outdir):
	process_manager = giantradiopulse.process.ProcessManager(
		file_path,outdir=outdir)
	process_manager.add_observations()

@cli.command(help="Generate light curves.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outdir", type=click.Path(), default='out/crabgrp')
def generate_lightcurves(file_path,outdir):
	process_manager = giantradiopulse.process.ProcessManager(
		file_path,outdir=outdir)
	process_manager.generate_lightcurves()	

@cli.command(help="Generate background spectrum.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outdir", type=click.Path(), default='out/crabgrp')
def generate_bgdspec(file_path,outdir):
	process_manager = giantradiopulse.process.ProcessManager(
		file_path,outdir=outdir)
	process_manager.generate_bgdspec()	

@cli.command(help="Generate source spectrum.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outdir", type=click.Path(), default='out/crabgrp')
def generate_srcspec(file_path,outdir):
	process_manager = giantradiopulse.process.ProcessManager(
		file_path,outdir=outdir)
	process_manager.generate_srcspec()	

@cli.command(help="Add GRP flag.")
@click.argument("file_path", type=click.Path(exists=True))
@click.argument("indir", type=click.Path(exists=True))
@click.option("--outdir", type=click.Path(), default='out/crabgrp')
def add_grp_flag_to_xrayevents(file_path,indir,outdir):
	process_manager = giantradiopulse.process.ProcessManager(
		file_path,outdir=outdir)
	process_manager.add_grp_flag_to_xrayevents(indir)	

def main():
	cli()

if __name__ == "__main__":
	main()
