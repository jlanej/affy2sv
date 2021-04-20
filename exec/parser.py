#!/usr/bin/python2.7
# encoding: utf-8


import argparse
import os
import pandas
import sys


# Headers
headers = {
        'sds': 'Chromosome Display StartIndex MarkerCount MinSignal MaxSignal MedianCnState HomFrequency HetFrequency Mosaicism LOH MedianSignal',
        'cnds': 'ProbeSetName Chromosome Position Log2Ratio WeightedLog2Ratio SmoothSignal',
        'apds': 'ProbeSetName Chromosome Position AllelePeaks0 AllelePeaks1',
        'mabsds': 'Index SCAR',
        'cn_ds': 'SegmentID Chromosome StartPosition StopPosition MarkerCount MeanMarkerDistance State Confidence',
        'lohds': 'SegmentID Chromosome StartPosition StopPosition MarkerCount MeanMarkerDistance LOH Confidence',
        'cnnlohds': 'SegmentID Chromosome StartPosition StopPosition MarkerCount MeanMarkerDistance CNNeutralLOH Confidence',
        'genotype': 'Index Call Confidence ForcedCall ASignal BSignal SignalStrength Contrast',
}

sel_full = ('cnds', 'genotype', 'apds')
sel_short = ('cnds', 'genotype')

def set_sel(short=False):
    global sel
    if short:
        sel = sel_short
    else:
        sel = sel_full

def mach(line):
    for key in headers:
        if line == headers[key] and key in sel:
            return key
    return None





def main():

    # Step 1: Parsing arguments
    # --------------------------------------------------------------------------
    # Here we got the program's arguments and its values.

    parser = argparse.ArgumentParser(description='Parser and filter for cychp for CytoScan HD and CytoSCan 750K to be sued with affy2sv.')
    parser.add_argument('-i', '--input_file', dest='input_file', help='Input CYCHP file.')
    parser.add_argument('-o', '--output_path', dest='output_path', help='Output directory.')
    parser.add_argument('-s','--short', dest='short', help='Does not take care of APDS.',
        action='store_true', default=False)
    #parser.add_argument('-c','--clean', dest='clean', help='Remove non ".f." files.',
    #    action='store_true', default=False)
    
    args = parser.parse_args()
	
    if not os.path.isfile(args.input_file):
        print '[ERROR]: Given "input file" argument is not a file.'
        sys.exit(1)

    if not os.path.isdir(args.output_path):
        print '[ERROR]: Provided output path is not a directory.'
        sys.exit(1)

    set_sel(args.short)

    # Step 2
    # --------------------------------------------------------------------------
    # The idea is to copy to different files the tables into the 
    # cydh.cycho.txt files. So we I read all the lines and we detect when we 
    # enter into a tables looking for the '#' character. When we enter into a 
    # table we create a new file and we copy the incoming lines until we 
    # find a new '#' character.

    with open(args.input_file, 'r') as fi:
        out = os.path.join(args.output_path, os.path.basename(fi.name))
        write = False
        nt, p, = 1, '#'
        for line in fi:
            if p.startswith('#') and not line.startswith('#'):
                write = True
                tn = mach(line.strip().replace('\t', ' '))
                if tn is None:
                    fo = None
                else:
                    fo = open(out + '.' + tn + '.txt', 'w')
            if not p.startswith('#') and line.startswith('#'):
                write = False
                fo.close() if fo is not None else None
            if write and fo is not None:
                fo.write(line)
            p = line

    # Step 3
    # --------------------------------------------------------------------------
    # Load the files into pandas and filter using '.ix'. Then store the result.

    if not args.short:
        apds = pandas.read_csv(out + '.apds.txt', sep='\t', low_memory=False)
    cnds = pandas.read_csv(out + '.cnds.txt', sep='\t', low_memory=False)
    gnds = pandas.read_csv(out + '.genotype.txt', sep='\t', low_memory=False)


    names = cnds.ix[ gnds.Index ]['ProbeSetName']
    names = names.reset_index()['ProbeSetName']
    gnds['ProbeSetName'] = names

    if not args.short:
        apds = apds.set_index('ProbeSetName')
    gnds = gnds.set_index('ProbeSetName')
    cnds = cnds.set_index('ProbeSetName')
    
    if not args.short:
        apds = apds.groupby(apds.index).first()
    gnds = gnds.groupby(gnds.index).first()
    cnds = cnds.groupby(cnds.index).first()

    if not args.short:
        im = list(set.intersection(set(gnds.index), set(cnds.index), set(apds.index)))
    else:
        im = list(set.intersection(set(gnds.index), set(cnds.index)))
    s = pandas.Series(im)

    if not args.short:
        apds_f = apds.ix[s]
    cnds_f = cnds.ix[s]
    gnds_f = gnds.ix[s]
    
    if not args.short:
        apds_f['ProbeSetName'] = pandas.Series(im, index=apds_f.index)
    gnds_f['ProbeSetName'] = pandas.Series(im, index=gnds_f.index)
    cnds_f['ProbeSetName'] = pandas.Series(im, index=cnds_f.index)

    if not args.short:
        apds_f.to_csv(str.replace(out, ".txt", "") + '.apds.f.txt', sep='\t', index=False)
    gnds_f.to_csv(str.replace(out, ".txt", "") + '.genotype.f.txt', sep='\t', index=False)
    cnds_f.to_csv(str.replace(out, ".txt", "") + '.cnds.f.txt', sep='\t', index=False)
    
    #if args.clean:
    # if not args.short:
        # os.remove(out + '.apds.txt')
    # os.remove(out + '.genotype.txt')
    # os.remove(out + '.cnds.txt')


if __name__ == '__main__':
    main()