#!/usr/bin/env python

from multiprocessing import Pool
import pandas as pd
import argparse
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(description='mappabilityCal')
    parser.add_argument('-m', '--map', help="each position's mappability")
    parser.add_argument(
        '-b', '--bed', help="windows distribution among each chromosome")
    args = parser.parse_args()
    return args


def mappability(chromosome):
    mapp = {}
    with open(args.map, 'r') as f:
        for i in f:
            line = i.strip().split()
            if line[0] == chromosome:
                position = int(line[1])
                value = line[3]
                mapp[position] = value
    bed = []
    with open(args.bed, 'r') as f:
        for i in f:
            line = i.strip().split()
            if line[0] == chromosome:
                start = int(line[1])
                end = int(line[2])
                bed.append(line)
    result = []
    for i in bed:
        value = 0
        start = int(i[1])
        end = int(i[2])
        ori = start
        while start < end:
            if start in mapp:
                value += float(mapp[start])
                start += 1
            else:
                start += 1
        bin = [chromosome, ori, end, value/(ori - end)]
        result.append(bin)
    return result


if __name__ == '__main__':
    args = get_args()
    chrs = []
    for i in range(1, 23):
        chromosome = 'chr' + str(i)
        chrs.append(chromosome)
    chrs.append('chrX')
    chrs.append('chrY')
    pool = Pool(24)
    resultList = {}
    for chromosome in chrs:
        result = pool.apply_async(mappability, args=(chromosome,))
        resultList[chromosome]=result
    pool.close()
    pool.join()
    with open('mappability.txt', 'w') as f:
        for chromosome in chrs:
            results=resultList[chromosome].get()
            for bin in results:
                f.write('\t'.join(map(str, bin))+'\n')
