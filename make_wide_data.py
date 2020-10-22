f = open('diacl_qualitative_coding.tsv','r')
text = f.read()
f.close()


text = [l.split('\t') for l in text.split('\n')]


states = {}
for l in text:
    if l[0] not in states.keys():
        states[l[0]]={}
    states[l[0]][l[1]]=l[2]


chars = sorted(set([l for k in states.keys() for l in states[k].keys()]))


wide = [['Language']+chars]
for l in states.keys():
    line = [l.replace('(','').replace(')','')]
    for c in chars:
        if c in states[l].keys():
            line.append(states[l][c])
        else:
            line.append('NA')
    wide.append(line)


f = open('diacl_qualitative_wide.tsv','w')
for l in wide:
    print('\t'.join(l),file=f)


f.close()