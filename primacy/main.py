import click
import logging
import json
import os
from sequence import Sequence
from collection import get_primer_collection


@click.group()
def main():
    pass

@click.group()
def gui():
    pass

@gui.command()
@click.argument('config', type=click.Path(exists=True))
@click.argument('outpath', type=click.Path(exists=False), default=os.getcwd())
@click.option('--threads', '-t', default=2, type=click.IntRange(min=1))
def primer_collection(config, threads, outpath=os.getcwd()):
    with open(config, 'r') as json_in:
        config_obj = json.loads(json_in.read())
    background_paths = \
        config_obj['primer_collection']['params']['background_seq']
    tm_params = config_obj['primer_collection']['params']['pcr_salts']
    seq_obj_dict = {seq_id: Sequence(
        seq_id=seq_id,
        primer_len_min=seq_info['primer_len_range']['min'],
        primer_len_max=seq_info['primer_len_range']['max'],
        **seq_info)
        for seq_id, seq_info in config_obj['sequences'].items()}
    get_primer_collection(
        seq_obj_dict, tm_params, background_paths, outpath, threads)


if __name__ == '__main__':
    gui()