with open("UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem.txt") as f:
    data = f.read()
lines = [s.strip().split() for s in data.splitlines()]
#print(lines[0])
genotype = []
chromosome = []
position1 = []
position2 = []
mutation_type = []
length = []
allele = []
mutation = []
left_reference = []
right_reference = []
left_consensus = []
right_consensus = []
trinucleotide = []
for line in lines:
    genotype.append(line[4])
    chromosome.append('chr' + str(line[6]))
    position1.append(line[8])
    position2.append(line[9])
    mutation_type.append(line[10])
    length.append(line[11])
    allele.append(line[12])
    mutation.append(line[13])
    if len(line) > 22:
        left_reference.append(line[23])
        right_reference.append(line[24])
        left_consensus.append(line[25])
        right_consensus.append(line[26])
        if line[27].isnumeric():
            trinucleotide.append(line[24])
        else:
            trinucleotide.append(line[27])
    else:
        left_reference.append('NA')
        right_reference.append('NA')
        left_consensus.append('NA')
        right_consensus.append('NA')
        trinucleotide.append('NA')
f1 = open("UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem_insertions.bed", 'w+')
f1a = open("WT_DIPs_sorted_anz2_tandem_insertions.bed", 'w+')
f1b = open("rad16_DIPs_sorted_anz2_tandem_insertions.bed", 'w+')
f1c = open("rad26_DIPs_sorted_anz2_tandem_insertions.bed", 'w+')
f1d = open("rad30_DIPs_sorted_anz2_tandem_insertions.bed", 'w+')
f2 = open("UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem_deletions.bed", 'w+')
f2a = open("WT_DIPs_sorted_anz2_tandem_deletions.bed", 'w+')
f2b = open("rad16_DIPs_sorted_anz2_tandem_deletions.bed", 'w+')
f2c = open("rad26_DIPs_sorted_anz2_tandem_deletions.bed", 'w+')
f2d = open("rad30_DIPs_sorted_anz2_tandem_deletions.bed", 'w+')

def WriteChr(chr_name):
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(1, len(lines)):
        if mutation_type[x] == 'DIP':
            if chromosome[x] == chr_name:
                chromosome1_dict[int(position1[x])] = x
    for key1 in chromosome1_dict.keys():
        chromosome1.append(int(key1))
    chromosome1.sort()
    for y in range(0, len(chromosome1)):
        
        file_data.append(chromosome[chromosome1_dict[chromosome1[y]]] + '\t' + position1[chromosome1_dict[chromosome1[y]]] + '\t' + position2[chromosome1_dict[chromosome1[y]]] + '\t' + trinucleotide[chromosome1_dict[chromosome1[y]]] + '\t' + allele[chromosome1_dict[chromosome1[y]]] + '\t' + mutation[chromosome1_dict[chromosome1[y]]] + '\t' + genotype[chromosome1_dict[chromosome1[y]]])
    for fd in range(0, len(file_data)):
        if '-' in allele[chromosome1_dict[chromosome1[fd]]]:
            f1.write(file_data[fd])
            f1.write('\n')
            if genotype[chromosome1_dict[chromosome1[fd]]] == 'WT':
                f1a.write(file_data[fd])
                f1a.write('\n')
            if  genotype[chromosome1_dict[chromosome1[fd]]] == 'rad16D/rad16D':
                f1b.write(file_data[fd])
                f1b.write('\n')
            if  genotype[chromosome1_dict[chromosome1[fd]]] == 'rad26D/rad26D':
                f1c.write(file_data[fd])
                f1c.write('\n')
            if  genotype[chromosome1_dict[chromosome1[fd]]] == 'rad30D/rad30D':
                f1d.write(file_data[fd])
                f1d.write('\n')

        if '-' in mutation[chromosome1_dict[chromosome1[fd]]]:
            f2.write(file_data[fd])
            if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                f2.write('\t+')
            else:
                f2.write('\t-')
            f2.write('\n')
            if genotype[chromosome1_dict[chromosome1[fd]]] == 'WT':
                f2a.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f2a.write('\t+')
                else:
                    f2a.write('\t-')
                f2a.write('\n')
            if  genotype[chromosome1_dict[chromosome1[fd]]] == 'rad16D/rad16D':
                f2b.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f2b.write('\t+')
                else:
                    f2b.write('\t-')
                f2b.write('\n')
            if  genotype[chromosome1_dict[chromosome1[fd]]] == 'rad26D/rad26D':
                f2c.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f2c.write('\t+')
                else:
                    f2c.write('\t-')
                f2c.write('\n')
            if  genotype[chromosome1_dict[chromosome1[fd]]] == 'rad30D/rad30D':
                f2d.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f2d.write('\t+')
                else:
                    f2d.write('\t-')
                f2d.write('\n')


WriteChr('chrI')
WriteChr('chrII')
WriteChr('chrIII')
WriteChr('chrIV')
WriteChr('chrV')
WriteChr('chrVI')
WriteChr('chrVII')
WriteChr('chrVIII')
WriteChr('chrXI')
WriteChr('chrX')
WriteChr('chrXI')
WriteChr('chrXII')
WriteChr('chrXIII')
WriteChr('chrXIV')
WriteChr('chrXV')
WriteChr('chrXVI')

countsi = {}
countsd = {}


for x in range(0, len(lines)):
    if mutation_type[x] == 'DIP':
        if '-' in allele[x]:
            if length[x] not in countsi.keys():
                countsi[length[x]] = 1
            else:
                countsi[length[x]] = countsi[length[x]] + 1
        if '-' in mutation[x]:
            if length[x] not in countsd.keys():
                countsd[length[x]] = 1
            else:
                countsd[length[x]] = countsd[length[x]] + 1

f3 = open('indel_counts', 'w+')
f3.write('Insertions')
f3.write('\n')
for key in countsi.keys():
    f3.write(str(key) + ": " + str(countsi[key]))
    f3.write('\n')
f3.write('Deletions')
f3.write('\n')
for key in countsd.keys():
    f3.write(str(key) + ": " + str(countsd[key]))
    f3.write('\n')
f3.close()

f1.close()
f1a.close()
f1b.close()
f1c.close()
f1d.close()

f2a = open("WT_DIPs_sorted_anz2_tandem_deletions.bed")
f2b = open("rad16_DIPs_sorted_anz2_tandem_deletions.bed")
f2c = open("rad26_DIPs_sorted_anz2_tandem_deletions.bed")
f2d = open("rad30_DIPs_sorted_anz2_tandem_deletions.bed")

def CountStrands(f):
    strand = []
    data = f.read()
    lines = [s.strip().split() for s in data.splitlines()]
    for line in lines:
        strand.append(line[7])
    return strand

strands = CountStrands(f2a)
strands16 = CountStrands(f2b)
strands26 = CountStrands(f2c)
strands30 = CountStrands(f2d)



f4 = open('stands.txt', 'w+')
counts = {}
counts['WT+'] = 0
counts['WT-'] = 0
counts['rad16+'] = 0
counts['rad16-'] = 0
counts['rad26+'] = 0
counts['rad26-'] = 0
counts['rad30+'] = 0
counts['rad30-'] = 0
for a in range(0, len(strands)):
    if strands[a] == '+':
        counts['WT+'] = counts['WT+'] + 1
    else:
        counts['WT-'] = counts['WT-'] + 1
for b in range(0, len(strands16)):
    if strands16[b] == '+':
        counts['rad16+'] = counts['rad16+'] + 1
    else:
        counts['rad16-'] = counts['rad16-'] + 1
for c in range(0, len(strands26)):
    if strands26[c] == '+':
        counts['rad26+'] = counts['rad26+'] + 1
    else:
        counts['rad26-'] = counts['rad26-'] + 1
for d in range(0, len(strands30)):
    if strands30[d]== '+':
        counts['rad30+'] = counts['rad30+'] + 1
    else:
        counts['rad30-'] = counts['rad30-'] + 1

for key in counts.keys():
    f4.write(key + ': ' + str(counts[key]))
    f4.write('\n')
f4.close()

f2.close()
f2a.close()
f2b.close()
f2c.close()
f2d.close()
numbers = []
strain = []
with open('SraRunTable.txt') as r:
    rundata = r.read()
Sralines = [s.strip().split(',') for s in rundata.splitlines()]
for line in Sralines:
    numbers.append(line[0])
    print(numbers[len(numbers) - 1])
    strain.append(line[26])
    print(strain[len(strain) - 1])
print(strain[0])

fo = open('accession_numbers.txt', 'w+')
fo2 = open('strains.txt', 'w+')
for x in range(1, len(numbers)):
    if 'WT' not in strain[x]:
        fo.write(numbers[x])
        fo.write('\n')
for y in range(1, len(strain)):
    if 'WT' not in strain[y]:
        fo2.write(strain[y])
        fo2.write('\n')

fo.close()
fo2.close()
