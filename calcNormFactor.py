#!/usr/bin/env python
'''
Usage: calcNormFactor.py [options] -t <tgt_genome> -s <spk_genome> [-q <mq>] <sample_tbl>

Options:
    -t <tgt_genome>     target genome, e.g. ce10
    -s <spk_genome>     spike-in genome, e.g. cb3
    -q <mq>             mapq threshold [default: 10]
    <sample_tbl>        sample table, contains "sample", "factor", "input", "strain"
    --debug             print debug info
'''

from __future__ import print_function
import sys
import signal
import logging
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import os.path
import itertools
import pandas as pd
import numpy as np

def extract_number(filename):
    with open(filename) as fh:
        for line in fh:
            return line.rstrip()

def main(args):
    tgt_genome = args['t']
    spk_genome = args['s']
    mq = args['q']
    sample_tbl = pd.read_table(args['sample_tbl'], sep='\t')

    chip_samples = sample_tbl[sample_tbl.factor!='input']['sample']
    chip_tgt_total = []
    chip_tgt_blacklist = []
    chip_tgt_peak = []
    chip_spk_total = []
    chip_spk_blacklist = []
    chip_spk_peak = []
    input_tgt_total = []
    input_tgt_blacklist = []
    input_spk_total = []
    input_spk_blacklist = []
    chip_spk_peak_height = []

    for s in chip_samples:
        s_input = sample_tbl[sample_tbl['sample']==s]['input'].item()
        chip_FNs = ['{}/{}.{}.mq{}.{}.readCount'.format(s,s,genome,mq,typ)
                for genome,typ in itertools.product((tgt_genome,spk_genome), ('total','blacklist','peak'))]
        input_FNs = ['{}/{}.{}.mq{}.{}.readCount'.format(s_input,s_input,genome,mq,typ)
                for genome,typ in itertools.product((tgt_genome,spk_genome), ('total','blacklist'))]
        chip_spk_pkht_fn = '{}/{}.{}.mq{}.avgPeakHeight'.format(s,s,spk_genome,mq)
        for fn in chip_FNs + input_FNs + [chip_spk_pkht_fn]:
            assert os.path.exists(fn), '{} does not exist'.format(fn)
        # Contents of 'numbers' are:
        #  chip_tgt_total,  chip_tgt_blacklist, chip_tgt_peak,
        #  chip_spk_total,  chip_spk_blacklist, chip_spk_peak,
        #  input_tgt_total, input_tgt_blacklist,
        #  input_spk_total, input_spk_blacklist,
        #  chip_spk_peak_height 
        numbers = [extract_number(fn) for fn in chip_FNs + input_FNs + [chip_spk_pkht_fn]]
        chip_tgt_total.append(int(numbers[0]))
        chip_tgt_blacklist.append(int(numbers[1]))
        chip_tgt_peak.append(int(numbers[2]))
        chip_spk_total.append(int(numbers[3]))
        chip_spk_blacklist.append(int(numbers[4]))
        chip_spk_peak.append(int(numbers[5]))
        input_tgt_total.append(int(numbers[6]))
        input_tgt_blacklist.append(int(numbers[7]))
        input_spk_total.append(int(numbers[8]))
        input_spk_blacklist.append(int(numbers[9]))
        chip_spk_peak_height.append(float(numbers[10]))

    dat = pd.DataFrame({
                'sample':chip_samples,
                'chip_tgt_total':chip_tgt_total,
                'chip_tgt_blacklist':chip_tgt_blacklist,
                'chip_tgt_peak':chip_tgt_peak,
                'chip_spk_total':chip_spk_total,
                'chip_spk_blacklist':chip_spk_blacklist,
                'chip_spk_peak':chip_spk_peak,
                'input_tgt_total':input_tgt_total,
                'input_tgt_blacklist':input_tgt_blacklist,
                'input_spk_total':input_spk_total,
                'input_spk_blacklist':input_spk_blacklist,
                'chip_spk_peak_height':chip_spk_peak_height})

    dat['norm_factor'] = 1e6 / (dat.chip_spk_total*dat.input_tgt_total/dat.input_spk_total)
    dat['norm_factor_rmblack'] = 1e6 / ((dat.chip_spk_total-dat.chip_spk_blacklist)*
            (dat.input_tgt_total-dat.input_tgt_blacklist)/(dat.input_spk_total-dat.input_spk_blacklist))
    dat['norm_factor_peak'] = 1e6 / (dat.chip_spk_peak*(dat.input_tgt_total-dat.input_tgt_blacklist)/(dat.input_spk_total-dat.input_spk_blacklist))
    dat['norm_factor_bg'] = 1e6 / ((dat.chip_spk_total-dat.chip_spk_blacklist-dat.chip_spk_peak)
            *(dat.input_tgt_total-dat.input_tgt_blacklist)/(dat.input_spk_total-dat.input_spk_blacklist))
    dat['norm_factor_snr'] = 1e6 * (1+1/(dat.chip_spk_peak_height-1)) / (dat.chip_spk_total*dat.input_tgt_total/dat.input_spk_total)
    dat['norm_factor_rmblack_snr'] = 1e6 * (1+1/(dat.chip_spk_peak_height-1)) / (
            (dat.chip_spk_total-dat.chip_spk_blacklist)*(dat.input_tgt_total-dat.input_tgt_blacklist)/(dat.input_spk_total-dat.input_spk_blacklist))
    dat['norm_factor_peak_snr'] = 1e6 * (1+1/(dat.chip_spk_peak_height-1)) / (
            dat.chip_spk_peak*(dat.input_tgt_total-dat.input_tgt_blacklist)/(dat.input_spk_total-dat.input_spk_blacklist))
    dat['norm_factor_bg_snr'] = 1e6 * (1+1/(dat.chip_spk_peak_height-1)) / (
            (dat.chip_spk_total-dat.chip_spk_blacklist-dat.chip_spk_peak)*(dat.input_tgt_total-dat.input_tgt_blacklist)/(dat.input_spk_total-dat.input_spk_blacklist))
    dat[['sample',
        'chip_tgt_total','chip_tgt_blacklist','chip_tgt_peak',
        'chip_spk_total','chip_spk_blacklist','chip_spk_peak',
        'input_tgt_total','input_tgt_blacklist',
        'input_spk_total','input_spk_blacklist',
        'chip_spk_peak_height',
        'norm_factor','norm_factor_rmblack','norm_factor_peak','norm_factor_bg',
        'norm_factor_snr','norm_factor_rmblack_snr','norm_factor_peak_snr','norm_factor_bg_snr']].to_csv(sys.stdout, sep='\t',header=True,index=False)


if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)
    args = {k.lstrip('-<').rstrip('>'):args[k] for k in args}
    try:
        if args.get('debug'):
            logLevel = logging.DEBUG
        else:
            logLevel = logging.WARN
        logging.basicConfig(
                level=logLevel,
                format='%(asctime)s; %(levelname)s; %(funcName)s; %(message)s',
                datefmt='%y-%m-%d %H:%M:%S')
        if args.get('prof'):
            import cProfile
            cProfile.run('main(args)')
        else:
            main(args)
    except KeyboardInterrupt:
        logging.warning('Interrupted')
        sys.exit(1)
