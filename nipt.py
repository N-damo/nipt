import argparse
import subprocess
from scipy import stats
import pandas as pd
from argparse import ArgumentParser


def main():
    parser = ArgumentParser(description='gc biased correction')
    parser.add_argument('-bin', '--bin_size', help='bin_size', default=50000)
    parser.add_argument('-b', '--bam', help='bam file')
    parser.add_argument('-c', '--chr_size', help='chromosome size,each row is chr\\tsize',default='hg38.chroms.size')
    parser.add_argument('-m', '--map', help="each position's mappability")
    parser.add_argument('-g', '--gc_correct_r',
                        help="provide loess.R script", default='loess.R')
    parser.add_argument('-w', '--wilcox_test_r',
                        help="provide wilcox.R script", default='wilcox.R')
    args = parser.parse_args()
    return args


class correct(object):

    def __init__(self, bam, bin_size):
        self.bam = bam
        self.bin_size = bin_size
        self.filename = self.bam.split('/')[-1].replace('.bam', '')
        self.write_sh()
        self.test()
        self.write_r()

    def write_sh(self):
        cmd1 = 'bedtools makewindows -g {} -w {} > hg38_{}.bed ;'.format(
            args.chr_size, self.bin_size, self.bin_size)
        cmd2 = 'bedtools nuc -fi /share/data1/genome/hs38DH.fa  -bed hg38_{}.bed | cut -f 1-3,5 >  {}.gc.bed ;'.format(
            self.bin_size, self.bin_size)
        cmd3 = 'bedtools coverage  -a hg38_{}.bed -b {} > {}.counts ;'.format(
            self.bin_size, self.bam, self.filename)
        cmd4 = 'python mappability.py -m {} -b hg38_{}.bed ;'.format(
            args.map, self.bin_size)
        cmd5 = "paste <(grep -v '#' {}.gc.bed) <(cut -f4 mappability.txt) <(cut -f4 {}.counts) |sed '1i chr\\tstart\\tend\\tGC\\tMP\\tRC'| awk '$6 != 0 {}' > niptest ;".format(
            self.bin_size, self.filename, '{print}')
        cmd6 = "R --no-save <{}>f".format(args.gc_correct_r)
        with open('correct.sh', 'w') as f:
            f.write('\n'.join([cmd1, cmd2, cmd3, cmd4, cmd5, cmd6])+'\n')
            f.write('wait;')
        cmd = 'bash correct.sh'
        subprocess.call(cmd, shell=True)
        return None

    def test(self):
        df = pd.read_table('nipt.txt', header=0)
        df = df[df['GC'] != 0]
        df2 = df[~df['chr'].isin(['chr13', 'chr18', 'chr21'])]
        chr_dict = {}
        all = df2['RC'].sum()
        with open('basic.txt', 'w') as f:
            f.write('chr\treads\tpercent\n')
            for group in df.groupby('chr'):
                chr = group[0]
                rc_sum = group[1]['RC'].sum()
                chr_dict[chr] = rc_sum
                p = rc_sum/float(all)
                print chr, rc_sum, p
                f.write('\t'.join([chr, str(rc_sum), str(p)])+'\n')
        for trisomy in ['chr13', 'chr18', 'chr21']:
            chr = df[df['chr'] == trisomy]
            chr.reset_index(inplace=True)
            chr = chr.drop('index', axis=1)
            artificial = chr.copy()
            for i in range(len(artificial)):
                gc = artificial.loc[i, 'GC']
                _ = list(df2[df2['GC'] == gc]['RC'])
                if len(_) >= 3:
                    artificial.loc[i, 'RC'] = sum(
                        _[:3])/3.0  # stats.mode(_)[0][0]
                else:
                    artificial.loc[i, 'RC'] = ''
            _ = artificial[artificial['RC'] == ''].index
            artificial = artificial.drop(_)
            chr = chr.drop(_)
            artificial['chr'] = 'artificial' + trisomy
            artificial.to_csv('artificial' + trisomy,
                              index=False, header=True, sep=',')
            chr.to_csv(trisomy, index=False, header=True, sep=',')
        return None

    def write_r(self):
        cmd = "R --no-save <{}>f".format(args.wilcox_test_r)
        subprocess.call(cmd, shell=True)
        return None


if __name__ == '__main__':
    args = main()
    bin_size = args.bin_size
    correct(args.bam, bin_size)
