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
f1 = open("UV_yeast_Merged_SNP_DIPs_sorted_tandem_insertions.bed", 'w+')
f1a = open("WT_DIPs_sorted_tandem_insertions.bed", 'w+')
f1b = open("rad16_DIPs_sorted_tandem_insertions.bed", 'w+')
f1c = open("rad26_DIPs_sorted_tandem_insertions.bed", 'w+')
f1d = open("rad30_DIPs_sorted_tandem_insertions.bed", 'w+')
f2 = open("UV_yeast_Merged_SNP_DIPs_sorted_tandem_deletions.bed", 'w+')
f2a = open("WT_DIPs_sorted_tandem_deletions.bed", 'w+')
f2b = open("rad16_DIPs_sorted_tandem_deletions.bed", 'w+')
f2c = open("rad26_DIPs_sorted_tandem_deletions.bed", 'w+')
f2d = open("rad30_DIPs_sorted_tandem_deletions.bed", 'w+')

def WriteChr(chr_name, f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
             trinucleotide, allele, mutation, genotype, mutation_type):
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(1, len(mutation_type)):
        if mutation_type[x] == 'DIP' or mutation_type[x] == 'Complex DIP':
            if chromosome[x] == chr_name:
                chromosome1_dict[int(position1[x])] = x
    for key1 in chromosome1_dict.keys():
        chromosome1.append(int(key1))
    chromosome1.sort()
    for y in range(0, len(chromosome1)):
        file_data.append(chromosome[chromosome1_dict[chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[chromosome1[y]]]) - 1) + '\t' + position2[chromosome1_dict[chromosome1[y]]] + '\t' + trinucleotide[chromosome1_dict[chromosome1[y]]] + '\t' + allele[chromosome1_dict[chromosome1[y]]] + '\t' + mutation[chromosome1_dict[chromosome1[y]]])
    
    for fd in range(0, len(file_data)):
        if '-' not in mutation[chromosome1_dict[chromosome1[fd]]]:
            f1.write(file_data[fd])
            f1.write('\n')
            if genotype[chromosome1_dict[chromosome1[fd]]] == 'WT':
                f1a.write(file_data[fd])
                f1a.write('\n')
            if 'rad16D' in  genotype[chromosome1_dict[chromosome1[fd]]]:
                f1b.write(file_data[fd])
                f1b.write('\n')
            if 'rad26D' in genotype[chromosome1_dict[chromosome1[fd]]]:
                f1c.write(file_data[fd])
                f1c.write('\n')
            if 'rad30D' in  genotype[chromosome1_dict[chromosome1[fd]]]:
                f1d.write(file_data[fd])
                f1d.write('\n')

        if '-' not in allele[chromosome1_dict[chromosome1[fd]]]:
            f2.write(file_data[fd])
            if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                f2.write('\t+')
            elif allele[chromosome1_dict[chromosome1[fd]]] == 'A' or allele[chromosome1_dict[chromosome1[fd]]] == 'G':
                f2.write('\t-')
            else:
                f2.write('\tNA')
            f2.write('\n')
            if genotype[chromosome1_dict[chromosome1[fd]]] == 'WT':
                f2a.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f2a.write('\t+')
                elif allele[chromosome1_dict[chromosome1[fd]]] == 'A' or allele[chromosome1_dict[chromosome1[fd]]] == 'G':
                    f2a.write('\t-')
                else:
                    f2a.write('\tNA')
                f2a.write('\n')
            if 'rad16D' in  genotype[chromosome1_dict[chromosome1[fd]]]:
                f2b.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f2b.write('\t+')
                elif allele[chromosome1_dict[chromosome1[fd]]] == 'A' or allele[chromosome1_dict[chromosome1[fd]]] == 'G':
                    f2b.write('\t-')
                else:
                    f2b.write('\tNA')
                f2b.write('\n')
            if 'rad26D' in genotype[chromosome1_dict[chromosome1[fd]]]:
                f2c.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f2c.write('\t+')
                elif allele[chromosome1_dict[chromosome1[fd]]] == 'A' or allele[chromosome1_dict[chromosome1[fd]]] == 'G':
                    f2c.write('\t-')
                else:
                    f2c.write('\tNA')
                f2c.write('\n')
            if 'rad30D' in  genotype[chromosome1_dict[chromosome1[fd]]]:
                f2d.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f2d.write('\t+')
                elif allele[chromosome1_dict[chromosome1[fd]]] == 'A' or allele[chromosome1_dict[chromosome1[fd]]] == 'G':
                    f2d.write('\t-')
                else:
                    f2d.write('\tNA')
                f2d.write('\n')


WriteChr('chrI', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrII', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrIII', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrIV', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrIX', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrV', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrVI', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrVII', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrVIII', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrX', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXI', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXII', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXIII', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXIV', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXV', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXVI', f1, f1a, f1b, f1c, f1d, f2, f2a, f2b, f2c, f2d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)

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

f2a = open("WT_DIPs_sorted_tandem_deletions.bed")
f2b = open("rad16_DIPs_sorted_tandem_deletions.bed")
f2c = open("rad26_DIPs_sorted_tandem_deletions.bed")
f2d = open("rad30_DIPs_sorted_tandem_deletions.bed")
f2e = open('WT_singlebase_deletions.bed', 'w+')
f2f = open('rad16_singlebase_deletions.bed', 'w+')
f2g = open('rad26_singlebase_deletions.bed', 'w+')
f2h = open('rad30_singlebase_deletions.bed', 'w+')
f2i = open('WT_TriNuc_singlebase_deletions.bed', 'w+')
f2j = open('rad16_TriNuc_singlebase_deletions.bed', 'w+')
f2k = open('rad26_TriNuc_singlebase_deletions.bed', 'w+')
f2l = open('rad30_TriNuc_singlebase_deletions.bed', 'w+')

def GetReverseComplement(str):
    dict = {}
    dict['A'] = 'T'
    dict['C'] = 'G'
    dict['G'] = 'C'
    dict['T'] = 'A'
    str2 = ''
    for i in range(0, len(str)):
        str2 = str2 + dict[str[i]]
    return str2



def FilterBaseNumber(file1, file2):
    data = file1.read()
    lines = [s.strip().split() for s in data.splitlines()]
    for line in lines:
        if line[5] != 'NA':
            if line[4] == 'T' or line[4] == 'A':
                file2.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + 'NTN' + '\t' + '-' + '\t' + line[6])
                file2.write('\n')
            if line[4] == 'C' or line[4] =='G':
                file2.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + 'NCN' + '\t' + '-' + '\t' + line[6])
                file2.write('\n')
FilterBaseNumber(f2a, f2e)
FilterBaseNumber(f2b, f2f)
FilterBaseNumber(f2c, f2g)
FilterBaseNumber(f2d, f2h)

f2a = open("WT_DIPs_sorted_tandem_deletions.bed")
f2b = open("rad16_DIPs_sorted_tandem_deletions.bed")
f2c = open("rad26_DIPs_sorted_tandem_deletions.bed")
f2d = open("rad30_DIPs_sorted_tandem_deletions.bed")

def MakeTriNucFile(file1, file2):
    data = file1.read()
    lines = [s.strip().split() for s in data.splitlines()]
    for line in lines:
        if line[6] != 'NA' and line[3] != 'NA':
            if line[4] == 'A' or line[4] == 'G':
                file2.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + GetReverseComplement(line[3]) + '\t' + '-' + '\t' + line[6])
                file2.write('\n')
            if line[4] == 'C' or line[4] =='T':
                file2.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + '-' + '\t' + line[6])
                file2.write('\n')
MakeTriNucFile(f2a, f2i)
MakeTriNucFile(f2b, f2j)
MakeTriNucFile(f2c, f2k)
MakeTriNucFile(f2d, f2l)

f2.close()
f2a.close()
f2b.close()
f2c.close()
f2d.close()
f2e.close()
f2f.close()
f2g.close()
f2h.close()
f2i.close()
f2j.close()
f2k.close()
f2l.close()

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

with open("All_UniqueMutations_04_08_2020_tandem.txt") as f3:
    data3 = f3.read()
lines = [s.strip().split() for s in data3.splitlines()]
genotype2 = []
photoreactivation = []
batch = []
chromosome2 = []
position1_2 = []
position2_2 = []
mutation_type2 = []
allele2 = []
mutation2 = []
trinucleotide2 = []
for line in lines:
    if line[4] == '0':
        if line[5] == '1':
            genotype2.append(line[1])
            photoreactivation.append(line[4])
            batch.append(line[5])
            chromosome2.append('chr' + line[7])
            position1_2.append(line[8])
            position2_2.append(line[9])
            mutation_type2.append(line[10])
            allele2.append(line[12])
            mutation2.append(line[14])
def GetSequence(file):
    f1 = open(file)
    sequence1 = f1.read()
    sequence1 = str(sequence1.strip().upper())
    sequence1 = ''.join(sequence1.split('\n'))
    return sequence1

for pos in range(0, len(position2_2)):
    if chromosome2[pos] == 'chrI':
        sequence = GetSequence('chr1.fa')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrII':
        sequence = GetSequence('chr2.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrIII':
        sequence = GetSequence('chr3.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrIV':
        sequence = GetSequence('chr4.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrV':
        sequence = GetSequence('chr5.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrVI':
        sequence = GetSequence('chr6.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrVII':
        sequence = GetSequence('chr7.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrVIII':
        sequence = GetSequence('chr8.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrIX':
        sequence = GetSequence('chr9.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrX':
        sequence = GetSequence('chr10.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrXI':
        sequence = GetSequence('chr11.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrXII':
        sequence = GetSequence('chr12.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrXIII':
        sequence = GetSequence('chr13.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrXIV':
        sequence = GetSequence('chr14.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrXV':
        sequence = GetSequence('chr15.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    elif chromosome2[pos] == 'chrXVI':
        sequence = GetSequence('chr16.txt')
        trinucleotide2.append(sequence[int(position1_2[pos]) - 2: int(position2_2[pos]) + 1])
    else:
        trinucleotide2.append('NA')

f3a = open('UniqueMutations_04_08_2020_tandem_insertions.txt', 'w+')
f3b = open('WT_UniqueMutations_04_08_2020_tandem_insertions.txt', 'w+')
f3c = open('rad16_UniqueMutations_04_08_2020_tandem_insertions.txt', 'w+')
f3d = open('rad26_UniqueMutations_04_08_2020_tandem_insertions.txt', 'w+')
f3e = open('rad30_UniqueMutations_04_08_2020_tandem_insertions.txt', 'w+')
f3f = open('UniqueMutations_04_08_2020_tandem_deletions.txt', 'w+')
f3g = open('WT_UniqueMutations_04_08_2020_tandem_deletions.txt', 'w+')
f3h = open('rad16_UniqueMutations_04_08_2020_tandem_deletions.txt', 'w+')
f3i = open('rad26_UniqueMutations_04_08_2020_tandem_deletions.txt', 'w+')
f3j = open('rad30_UniqueMutations_tandem_deletions.txt', 'w+')

WriteChr('chrI', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrII', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrIII', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrIV', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrIX', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrV', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrVI', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrVII', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrVIII', f3a, f3b, f3c,  f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrX', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrXI', f3a, f3b, f3c,  f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrXII', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrXIII', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrXIV', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrXV', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)
WriteChr('chrXVI', f3a, f3b, f3c, f3d, f3e, f3f, f3g, f3h, f3i, f3j, chromosome2, position1_2, position2_2,
             trinucleotide2, allele2, mutation2, genotype2, mutation_type2)

f3a.close()
f3b.close()
f3c.close()
f3d.close()
f3e.close()
f3f.close()
f3g.close()
f3h.close()
f3i.close()
f3j.close()

f3j = open('rad30_UniqueMutations_tandem_deletions.txt')
f3k = open('rad30_UniqueMutations_singlebase_deletions.bed', 'w+')
FilterBaseNumber(f3j, f3k)
f3j.close()
f3k.close()
f3j = open('rad30_UniqueMutations_tandem_deletions.txt')
f3l = open('rad30_UniqueMutations_singlebase_TriNuc_deletions.bed', 'w+')

MakeTriNucFile(f3j, f3l)
f3j.close()
f3l.close()


def ContextCounts(file1, file2):
    context = []
    TS = []
    NTS = []
    intergenic = []
    data = file1.read()
    lines = [s.strip().split() for s in data.splitlines()]
    for l in range (1, len(lines)):
        context.append(lines[l][0])
        TS.append(int(lines[l][3]))
        NTS.append(int(lines[l][4]))
        intergenic.append(int(lines[l][5]))
    dipyrimidine_count = 0
    other_count = 0
    genic_count = 0
    intergenic_count = 0
    for a in range(1, len(TS)):
        genic_count = genic_count + TS[a] + NTS[a]
        intergenic_count = intergenic_count + intergenic[a]
        if 'TT' in context[a] or 'CC' in context[a]:
            dipyrimidine_count = dipyrimidine_count + TS[a] + NTS[a] + intergenic[a]
        else:
            other_count = other_count + TS[a] + NTS[a] + intergenic[a]
    file2.write('Genic Count: ' + str(genic_count))
    file2.write('\n')
    file2.write('Intergenic Count: ' + str(intergenic_count))
    file2.write('\n')
    file2.write('Ratio: ' + str(genic_count/intergenic_count))
    file2.write('\n')
    file2.write('Dipyrimidine Count: ' + str(dipyrimidine_count))
    file2.write('\n')
    file2.write('Non-dipyrimidine Count: ' + str(other_count))
    file2.write('\n')
    file2.write('Ratio: ' + str(dipyrimidine_count/other_count))

f4 = open('WT_deletion_stranddata.txt')
f4b = open('WT_deletion_strandresults.txt', 'w+')
f5 = open('rad16_deletion_stranddata.txt')
f5b = open('rad16_deletion_strandresults.txt', 'w+')
f6 = open('rad26_deletion_stranddata.txt')
f6b = open('rad26_deletion_strandresults.txt', 'w+')
f7 = open('rad30_deletion_stranddata.txt')
f7b = open('rad30_deletion_strandresults.txt', 'w+')
f8 = open('rad30_deletion_stranddata2.txt')
f8b = open('rad30_deletion_strandresults2.txt', 'w+')

ContextCounts(f4, f4b)
ContextCounts(f5, f5b)
ContextCounts(f6, f6b)
ContextCounts(f7, f7b)
ContextCounts(f8, f8b)

f4.close()
f4b.close()
f5.close()
f5b.close()
f6.close()
f6b.close()
f7.close()
f7b.close()
f8.close()
f8b.close()



