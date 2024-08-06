import sys


_, file_name = sys.argv

def count_repeats(string, key):
    count = 0
    while key * (count + 1) in string:
        count += 1
    return count

results = {}


with open(file_name, 'r') as file:
    for line in file:
        if line.startswith('>'):
            parts = line.strip().split('-')
            title = line.strip('>')
            period = parts[-1]
        else:
            repeats = count_repeats(line.strip(), period)
            results.setdefault(parts[0], []).append(repeats)


output_file = f'{file_name[:-4]}rptime.txt'
with open(output_file, 'w') as output:
    output.write('STRname\trefrptime\tthreekrptime\n')
    for key, values in results.items():
        if values[0] != values[1]:
            print(key)
        output.write(f'{key}\t{values[0]}\t{values[1]}\n')
