from pathlib import Path
import re

p = Path()
pattern = re.compile(r'\d+')


def process_and_write_entries(tempdt, w):
    minnum = 999999999
    if 'onlyTRUE' in tempdt:
        w.write(tempdt['onlyTRUE'])
    elif 'secondtrue' in tempdt:
        for k, v in tempdt.items():
            if k == 'secondtrue':
                for each in v:
                    if each[0] < minnum:
                        minnum = each[0]
                        mink = k
            else:
                if k[0] < minnum:
                    minnum = k[0]
                    mink = k
        if mink == 'secondtrue':
            for index, each in enumerate(tempdt.get(mink)):
                if each[0] == minnum:
                    w.write(tempdt['secondtrue'][index][1])
        else:
            w.write(tempdt.get(mink))
    else:
        templs = sorted(tempdt.keys(), key=lambda x: x[0])
        if len(templs) > 1:
            if templs[1][0] - templs[0][0] < 100000:
                templs1 = sorted(tempdt.keys(), key=lambda x: x[1])
                w.write(tempdt[templs1[0]])
            else:
                if templs:
                    w.write(tempdt[templs[0]])
        else:
            if templs:
                w.write(tempdt[templs[0]])


def process_files(file_pattern):
    for files in p.glob(file_pattern):
        with open(files, 'r') as r:
            with open(f'{files.stem}_output_new.txt', 'w') as w:
                print(files.name)
                tempdt = {}
                current_id = None
                for line in r:
                    fline = line.strip().split('\t')
                    query_id = fline[0]
                    before_chrom = fline[0].split('-')[0]
                    after_chrom = str(int(fline[1][3:]))

                    if current_id is None:
                        current_id = query_id

                    if query_id != current_id:
                        process_and_write_entries(tempdt, w)
                        tempdt = {}
                        current_id = query_id

                    if before_chrom == after_chrom:
                        nippos = int(pattern.search(fline[0].split('-')[-1]).group(0)[2:])
                        zspos = int(fline[8])
                        abspos = abs(nippos - zspos)
                        if float(fline[2]) == 100 and int(fline[3]) == 100:
                            tempdt.setdefault('onlyTRUE', line)
                        elif float(fline[2]) == 100 and int(fline[7]) == 100:
                            tempdt.setdefault('secondtrue', []).append((abspos, line))
                        else:
                            if int(fline[7]) == 100:
                                tempdt.setdefault((float(abspos), float(fline[-2])), line)
                process_and_write_entries(tempdt, w)
                print(f'{files.name} ending')



process_files('*_MH63_100bp.txt')
process_files('*_ZS97_100bp.txt')


from pathlib import Path

def filter_positions(totalline, threshold=10000):
    ls = []
    dt = {}
    final_result = []

    for line in totalline:
        fline = line.strip().split('\t')
        ls.append((fline[0], int(fline[8])))

    for i in range(1, len(ls)):
        current_line = ls[i][1]
        previous_line = ls[i - 1][1]
        if abs(previous_line - current_line) > threshold:
            dt.setdefault(ls[i-1][0], '')

    for line in totalline:
        fline = line.strip().split('\t')
        if fline[0] not in dt:
            final_result.append(line)


    dt.clear()
    ls.clear()

    for line in final_result:
        fline = line.strip().split('\t')
        ls.append((fline[0], int(fline[8])))

    final_result_second = []
    for i in range(1, len(ls)):
        current_line = ls[i][1]
        previous_line = ls[i - 1][1]
        if abs(previous_line - current_line) > threshold:
            dt.setdefault(ls[i-1][0], '')


    for line in final_result:
        fline = line.strip().split('\t')
        if fline[0] not in dt:
            final_result_second.append(line)

    return final_result_second


p = Path()

for files in p.glob('*_100bp_output_new.txt'):
    print(files.name)
    final_output_file = f'{files.stem}_filter.txt'

    with open(files, 'r') as r:
        totalline = r.readlines()

    filtered_lines = filter_positions(totalline)

    with open(final_output_file, 'w') as w:
        for line in filtered_lines:
            w.write(line)

    print(f'{files.name} filtering completed')