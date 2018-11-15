#!/usr/bin/env python
"""
HISTORY
2018-11-15 created by T.Enoto 
"""

import click 
import giantradiopulse.xrayevents

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
	if ctx.invoked_subcommand is None:
		print(ctx.get_help())

@cli.command(help="Analysis for combined X-ray event fitsfiles.")
@click.argument("xrayevt_list_file", type=click.Path(exists=True))
@click.option("--outdir",default='out/merge')
def run(xrayevt_list_file,outdir):
	xrayevent = giantradiopulse.xrayevents.XrayEventList(xrayevt_list_file,outdir=outdir)
	xrayevent.run()

def main():
	cli()

if __name__ == "__main__":
	main()
