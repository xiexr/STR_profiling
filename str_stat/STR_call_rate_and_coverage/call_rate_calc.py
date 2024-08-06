from pathlib import Path
import re
import pickle
import sys


def process_files(pathway,file_pattern, output_pickle):
    totalSTRposdt = {}
    pathmode = re.compile(r'\d+')
    p = Path(pathway)

    for files in p.glob(file_pattern):
        with open(files, 'r') as r:
            for line in r:
                if line.startswith('\n'):
                    continue
                elif line.startswith(('hipSTR', 'lobSTR','gatk')):
                    pos = int(pathmode.search(line.strip().split(',')[0]).group()[2:])
                    chrom = int(line.strip().split(',')[0][3:5])
                    totalSTRposdt.setdefault(chrom, {}).setdefault(pos, '')

    with open(output_pickle, 'wb') as w:
        pickle.dump(totalSTRposdt, w, protocol=pickle.HIGHEST_PROTOCOL)


def calculate_variant_frequencies(input_pickle, input_file):
    variantdt = {}
    Lastdt = {}

    with open(input_pickle, 'rb') as r:
        STRposdt = pickle.load(r)

    with open(input_file, 'r') as r:
        chrom = int(re.search(r'\d+', input_file).group())

        for line in r:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    fline = line.strip().split('\t')[9:]
                    for index, each in enumerate(fline):
                        variantdt.setdefault(index, each)
            else:
                fline = line.strip().split('\t')
                name = int(line.strip().split('\t')[1])
                variant = fline[9:]

                if name in STRposdt.get(chrom, {}):
                    for index, each in enumerate(variant):
                        if '.' in each.split(':')[0]:
                            continue
                        else:
                            Lastdt.setdefault(variantdt.get(index), 0)
                            Lastdt[variantdt.get(index)] += 1

    with open(f'callrateof{chrom}.txt', 'w') as w:
        for k, v in Lastdt.items():
            w.write(f'{k}\t{v}\n')


if __name__ == "__main__":
    process_files('/hipSTR','hipchrom*pos.csv', 'totalhipSTRpos.pkl')
    process_files('/lobSTR','lobchrom*pos.csv', 'totallobSTRpos.pkl')
    process_files('/gatk', 'gatkchrom*pos.csv', 'totallobgatkpos.pkl')

    _, input_file = sys.argv
    calculate_variant_frequencies('totalhipSTRpos.pkl', input_file_chrom)
    calculate_variant_frequencies('totallobSTRpos.pkl', input_file_chrom)
    calculate_variant_frequencies('totalgatkSTRpos.pkl', input_file_chrom)

