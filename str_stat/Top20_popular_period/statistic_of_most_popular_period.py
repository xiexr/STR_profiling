from collections import Counter
tempstatls=[]
with open('Last_result.txt','r') as r:
    for line in r:
        if line.startswith('\n'):
            continue
        elif line.startswith('period'):
            fline=line.strip().split(',')[1]
            tempstatls.append(fline)
tempstatdt=Counter(tempstatls)
with open('mostpopularperiod.txt','w') as w:
    for k,v in tempstatdt.items():
        w.write(f'{k}\t{v}\n')