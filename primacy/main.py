import click
import logging
import json
import os
from sequence import Sequence
from collection import get_primer_collection
from primer_score import get_primer_scores_mp
from utils import read_config, get_json_obj, write_json_obj


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
    config_obj = read_config(config)
    background_paths = \
        config_obj['primer_collection']['params']['background_seq']
    tm_params = config_obj['primer_collection']['params']['pcr_salts']
    seq_obj_dict = {seq_id: Sequence(
        seq_id=seq_id,
        primer_len_min=seq_info['primer_len_range']['min'],
        primer_len_max=seq_info['primer_len_range']['max'],
        **seq_info)
        for seq_id, seq_info in config_obj['sequences'].items()}
    # Run primer collection and return dictionary of output files
    output = get_primer_collection(
        seq_obj_dict, tm_params,
        background_paths, os.path.abspath(outpath), threads)
    # modify config json file with path to primer_collection for each sequence
    for seq_id, outpath in output.items():
        config_obj['sequences'][seq_id]['outfile'] = outpath
    write_json_obj(config_obj, config)


@gui.command()
@click.argument('config', type=click.Path(exists=True))
@click.option('--threads', '-t', default=2, type=click.IntRange(min=1))
def primer_score(config, threads):
    config_obj = read_config(config)
    primer_obj_paths = [
        seq_obj['outfile'] for 
            seq_obj in config_obj['sequences'].values()]
    get_primer_scores_mp(
        primer_obj_paths, config_obj['primer_scores']['params'], threads)


if __name__ == '__main__':
    gui()