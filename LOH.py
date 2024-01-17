import random

def GetSequence(file):
    f1 = open(file)
    sequence1 = f1.read()
    sequence1 = str(sequence1.strip().upper())
    sequence1 = ''.join(sequence1.split('\n'))
    return sequence1

def GetReverseComplement(str):
    dict = {}
    dict['A'] = 'T'
    dict['C'] = 'G'
    dict['G'] = 'C'
    dict['T'] = 'A'
    str2 = ''
    if str == '-':
        return '-'
    else:
        for i in range(0, len(str)):
            str2 = str2 + dict[str[i]]
        return str2[::-1]


with open("UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem.txt") as f:
    data = f.read()
lines = [s.strip().split() for s in data.splitlines()]
genotype = []
chromosome = []
position1 = []
position2 = []
mutation_type = []
length = []
allele = []
mutation = []
frequency = []
left_reference = []
right_reference = []
left_consensus = []
right_consensus = []
trinucleotide = []
isolate_list = {}
isolates = []
lines.pop(0)
for line in lines:
    if line[0] not in isolate_list.values():
        isolate_list[int(len(isolates))] = line[0]
    isolates.append(line[0])
    genotype.append(line[4])
    chromosome.append('chr' + str(line[6]))
    position1.append(line[8])
    position2.append(line[9])
    mutation_type.append(line[10])
    length.append(line[11])
    allele.append(line[12])
    mutation.append(line[13])
    if len(line) >= 22: 
        if '-' not in  line[18] and line[18].isalpha() == False:
            if float(line[18]) >= 90:
                frequency.append('Homozygous')
            else:
                frequency.append('Heterozygous')
        elif float(line[19]) >= 90:
            frequency.append('Homozygous')
        else:
            frequency.append('Heterozygous')
        left_reference.append(line[23])
        right_reference.append(line[24])
        left_consensus.append(line[25])
        right_consensus.append(line[26])
        if line[27].isalpha() == False:
            trinucleotide.append(line[24])
        else:
            trinucleotide.append(line[27])
    else:
        frequency.append('NA')
        left_reference.append('NA')
        right_reference.append('NA')
        left_consensus.append('NA')
        right_consensus.append('NA')
        trinucleotide.append('NA')

def WriteChr(chr_name, f1 , f2, f3, f4, f5, chromosome, position1, position2,
             trinucleotide, allele, mutation, genotype, mutation_type):
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(0, len(mutation_type)):
        if mutation_type[x] == 'SNP' :
            if chromosome[x] == chr_name:
                chromosome1_dict[int(position1[x])] = x
    for key1 in chromosome1_dict.keys():
        chromosome1.append(int(key1))
    chromosome1.sort()
    for y in range(0, len(chromosome1)):
        file_data.append(chromosome[chromosome1_dict[chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[chromosome1[y]]]) - 1) + '\t' + position2[chromosome1_dict[chromosome1[y]]] + '\t' + trinucleotide[chromosome1_dict[chromosome1[y]]] + '\t' + allele[chromosome1_dict[chromosome1[y]]] + '\t' + mutation[chromosome1_dict[chromosome1[y]]])
    
    for fd in range(0, len(file_data)):
        if frequency[chromosome1_dict[chromosome1[fd]]] != 'NA' and frequency[chromosome1_dict[chromosome1[fd]]] == 'Homozygous':
            f1.write(file_data[fd])
            if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                f1.write('\t+')
            elif allele[chromosome1_dict[chromosome1[fd]]] == 'A' or allele[chromosome1_dict[chromosome1[fd]]] == 'G':
                f1.write('\t-')
            else:
                f1.write('\tNA')
            f1.write('\t' + genotype[chromosome1_dict[chromosome1[fd]]])
            f1.write('\n')
        if genotype[chromosome1_dict[chromosome1[fd]]] == 'WT':
            if frequency[chromosome1_dict[chromosome1[fd]]] != 'NA' and frequency[chromosome1_dict[chromosome1[fd]]] == 'Homozygous':
                f2.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f2.write('\t+')
                elif allele[chromosome1_dict[chromosome1[fd]]] == 'A' or allele[chromosome1_dict[chromosome1[fd]]] == 'G':
                    f2.write('\t-')
                else:
                    f2.write('\tNA')
                f2.write('\t' + genotype[chromosome1_dict[chromosome1[fd]]])
                f2.write('\n')
        if 'rad16D' in genotype[chromosome1_dict[chromosome1[fd]]]:
            if frequency[chromosome1_dict[chromosome1[fd]]] != 'NA' and frequency[chromosome1_dict[chromosome1[fd]]] == 'Homozygous':
                f3.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f3.write('\t+')
                elif allele[chromosome1_dict[chromosome1[fd]]] == 'A' or allele[chromosome1_dict[chromosome1[fd]]] == 'G':
                    f3.write('\t-')
                else:
                    f3.write('\tNA')
                f3.write('\t' + genotype[chromosome1_dict[chromosome1[fd]]])
                f3.write('\n')
        if 'rad26D' in genotype[chromosome1_dict[chromosome1[fd]]]:
            if frequency[chromosome1_dict[chromosome1[fd]]] != 'NA' and frequency[chromosome1_dict[chromosome1[fd]]] == 'Homozygous':
                f4.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f4.write('\t+')
                elif allele[chromosome1_dict[chromosome1[fd]]] == 'A' or allele[chromosome1_dict[chromosome1[fd]]] == 'G':
                    f4.write('\t-')
                else:
                    f4.write('\tNA')
                f4.write('\t' + genotype[chromosome1_dict[chromosome1[fd]]])
                f4.write('\n')
        if 'rad30D' in genotype[chromosome1_dict[chromosome1[fd]]]:
            if frequency[chromosome1_dict[chromosome1[fd]]] != 'NA' and frequency[chromosome1_dict[chromosome1[fd]]] == 'Homozygous':
                f5.write(file_data[fd])
                if allele[chromosome1_dict[chromosome1[fd]]] == 'C' or allele[chromosome1_dict[chromosome1[fd]]] == 'T':
                    f5.write('\t+')
                elif allele[chromosome1_dict[chromosome1[fd]]] == 'A' or allele[chromosome1_dict[chromosome1[fd]]] == 'G':
                    f5.write('\t-')
                else:
                    f5.write('\tNA')
                f5.write('\t' + genotype[chromosome1_dict[chromosome1[fd]]])
                f5.write('\n')
        
        


f6 = open("LOH_substitutions.bed")
f7 = open('LOH_substitution_counts.txt', 'w+')
def LOHCounts(f1, f2):
    WTcounts = 0
    rad16counts = 0
    rad26counts = 0
    rad30counts = 0
    data1 = f1.read()
    lines1 = [s.strip().split() for s in data1.splitlines()]
    for line in lines1:
        if line[7] == 'WT':
            WTcounts = WTcounts + 1
        if 'rad16D' in line[7]:
            rad16counts = rad16counts + 1
        if 'rad26D' in line[7]:
            rad26counts = rad26counts + 1
        if 'rad30D' in line[7]:
            rad30counts = rad30counts + 1
    f2.write('WT counts = ' + str(WTcounts))
    f2.write('\n')
    f2.write('WT LOH frequency = '+ str(WTcounts/ 12199))
    f2.write('\n')
    f2.write('rad16 counts= ' + str(rad16counts))
    f2.write('\n')
    f2.write('rad16 LOH frequency= ' + str(rad16counts/23600))
    f2.write('\n')
    f2.write('rad26 counts= ' + str(rad26counts))
    f2.write('\n')
    f2.write('rad26 LOH frequency= ' + str(rad26counts/9948))
    f2.write('\n')
    f2.write('rad30 counts= ' + str(rad30counts))
    f2.write('rad30 LOH frequency= ' + str(rad30counts/16684))
    return [WTcounts, rad16counts, rad26counts, rad30counts]

LOH_counts = LOHCounts(f6, f7)
WTcounts = LOH_counts[0]
rad16counts = LOH_counts[1]
rad26counts = LOH_counts[2]
rad30counts = LOH_counts[3]

with open("WT_ogg1_UVA_All_mutations_comb_unique.txt") as f2:
    data2 = f2.read()
lines2 = [s.strip().split() for s in data2.splitlines()]
isolate_list_2 = {}
isolates2 = []
UV = []
genotype2 = []
chromosome2 = []
mutation_type_2 = []
position1_2 = []
position2_2 = []
allele2 = []
mutation2 = []
frequency2 = []
trinucleotide2 = []
lines2.pop(0)
for line in lines2:
    if line[0] not in isolate_list_2.values():
        isolate_list_2[int(len(isolates2))] = line[0]
    isolates2.append(line[0])
    genotype2.append(line[1])
    UV.append(line[2])
    chromosome2.append(str(line[5]))
    position1_2.append(int(line[6]) - 1)
    position2_2.append(int(line[6]))
    mutation_type_2.append(line[11])
    allele2.append(line[13])
    mutation2.append(line[14])
    frequency2.append(line[15])

for pos2 in range(0, len(position2_2)):
        if chromosome2[pos2] == 'chrI':
            sequence = GetSequence('chr1.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrII':
            sequence = GetSequence('chr2.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrIII':
            sequence = GetSequence('chr3.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrIV':
            sequence = GetSequence('chr4.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrV':
            sequence = GetSequence('chr5.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrVI':
            sequence = GetSequence('chr6.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrVII':
            sequence = GetSequence('chr7.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrVIII':
            sequence = GetSequence('chr8.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrIX':
            sequence = GetSequence('chr9.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrX':
            sequence = GetSequence('chr10.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrXI':
            sequence = GetSequence('chr11.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrXII':
            sequence = GetSequence('chr12.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrXIII':
            sequence = GetSequence('chr13.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrXIV':
            sequence = GetSequence('chr14.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrXV':
            sequence = GetSequence('chr15.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        elif chromosome2[pos2] == 'chrXVI':
            sequence = GetSequence('chr16.txt')
            trinucleotide2.append(sequence[int(position1_2[pos2]) - 1: int(position2_2[pos2]) + 1])
        else:
            trinucleotide2.append('NA')

   
    
WT_counts = 0
WT_LOH = 0
ogg1_counts = 0
ogg1_LOH = 0
UV_WT_counts = 0
UV_WT_LOH = 0
UV_ogg1_counts = 0
UV_ogg1_LOH = 0

f8 = open('ogg1_LOH_substitutions.bed', 'w+')
f9 = open('ogg1_LOH_indels.bed', 'w+')
f10 = open('ogg1_LOH_singlebase_deletions.bed', 'w+')

for i in range (0, len(genotype2)):
    if UV[i] == 'NO':
        if genotype2[i] == 'WT':
            WT_counts = WT_counts + 1
            if frequency2[i] == 'Homozygous':
                WT_LOH = WT_LOH + 1
        elif genotype2[i] == 'OGG':
            ogg1_counts = ogg1_counts + 1
            if frequency2[i] == 'Homozygous':
                ogg1_LOH = ogg1_LOH + 1
    elif UV[i] == 'UV':
        if genotype2[i] == 'WT':
            UV_WT_counts = UV_WT_counts + 1
            if frequency2[i] == 'Homozygous':
                UV_WT_LOH = UV_WT_LOH + 1
        elif genotype2[i] == 'OGG':
            UV_ogg1_counts = UV_ogg1_counts + 1
            if frequency2[i] == 'Homozygous':
                UV_ogg1_LOH = UV_ogg1_LOH + 1
                if mutation_type_2[i] == 'SNV':
                    f8.write(chromosome2[i] + '\t' + str(position1_2[i]) + '\t' + str(position2_2[i]) +
                         '\t' + allele2[i] + '\t' + trinucleotide2[i] + '\t' + mutation2[i])
                    f8.write('\n')
                if mutation_type_2 == 'Insertion' or mutation_type_2 == 'Deletion':
                    f9.write(chromosome2[i]  + '\t' + str(position2_2[i]) +
                         '\t' + allele2[i] + '\t' + mutation2[i])
                    f9.write('\n')
                if mutation_type_2 == 'Deletion':
                    if len(str(allele2[i])) == 1:
                        f10.write(chromosome2[i]  + '\t' + str(position1_2[i])  + '\t' + str(position2_2[i]) +
                         '\t' + allele2[i] + '\t' + mutation2[i])
                        f10.write('\n')


WT_frequency = WT_LOH/ WT_counts
print(WT_frequency)
ogg1_frequency = ogg1_LOH/ ogg1_counts
print(ogg1_frequency)
    
UV_WT_frequency = UV_WT_LOH/ UV_WT_counts
print(UV_WT_frequency)
UV_ogg1_frequency = UV_ogg1_LOH/UV_ogg1_counts
print(UV_ogg1_frequency)

print(len(chromosome))
print(len(position2))
print(len(frequency))
def GetLOHNeighbors(isolate_list, isolates, frequency, chromosome, position):
    LOH_neighbors = {}
    other_neighbors = {}
    prev_key = 0
    del isolate_list[0]
    isolate_list[int(len(position))] = 'end'
    for key in isolate_list.keys():
        #if prev_key == key - 1:
            #if frequency[key - 1] == 'Homozygous':
                #LOH_neighbors[key - 1] = 0
            #if frequency[key - 1] == 'Heterozygous':
                #other_neighbors[key - 1] = 0
        #else:
            if frequency[prev_key] == 'Homozygous':
                if chromosome[prev_key] == chromosome[prev_key + 1] and frequency[prev_key + 1] == 'Homozygous':
                    LOH_neighbors[prev_key] = 1
                else:
                    LOH_neighbors[prev_key] = 0
            else:
                if chromosome[prev_key] == chromosome[prev_key + 1] and frequency[prev_key + 1] == 'Homozygous':
                    other_neighbors[prev_key] = 1
                else:
                    other_neighbors[prev_key] = 0
    
            a = prev_key + 1
            b = key - 2
            if b >= a:
                for w in range(a,b + 1):
                    if frequency[w] == 'Homozygous':
                        if chromosome[w] == chromosome[w + 1] and chromosome[w] == chromosome[w - 1] and frequency[w - 1] == 'Homozygous' and frequency[w + 1] == 'Homozygous':
                            LOH_neighbors[w] = 2

                        elif frequency[w - 1] == 'Homozygous' and chromosome[w-1] == chromosome[w]:
                            LOH_neighbors[w] = 1
                        elif frequency[w + 1] == 'Homozygous' and chromosome[w+1] == chromosome[w]:
                            LOH_neighbors[w] = 1
                
                        else:
                            LOH_neighbors[w] = 0
                    else:
                        if chromosome[w] == chromosome[w + 1] and chromosome[w] == chromosome[w - 1] and frequency[w - 1] == 'Homozygous' and frequency[w + 1] == 'Homozygous':
                            other_neighbors[w] = 2

                        elif frequency[w - 1] == 'Homozygous' and chromosome[w-1] == chromosome[w]:
                            other_neighbors[w] = 1
                        elif frequency[w + 1] == 'Homozygous' and chromosome[w+1] == chromosome[w]:
                            other_neighbors[w] = 1
                
                        else:
                            other_neighbors[w] = 0
            if frequency[b + 1]  == 'Homozygous':
                if frequency[b ] == 'Homozygous' and chromosome[b + 1] == chromosome[b]:
                    LOH_neighbors[b + 1] =1

                else:
                    LOH_neighbors[b + 1] = 0
            else:
                if frequency[b] == 'Homozygous' and chromosome[b + 1] == chromosome[b]:
                    other_neighbors[b + 1] =1

                else:
                    other_neighbors[b + 1] = 0
            prev_key = key

    return [LOH_neighbors, other_neighbors]


neighbors = GetLOHNeighbors(isolate_list_2, isolates2, frequency2, chromosome2, position2_2)
LOH_neighbors = neighbors[0]
other_neighbors = neighbors[1]


neighbors2 = GetLOHNeighbors(isolate_list, isolates, frequency, chromosome, position2)
LOH_neighbors2 = neighbors2[0]
other_neighbors2 = neighbors2[1]


sum = 0
sum2 = 0
sum3 = 0
sum4 = 0
for key in LOH_neighbors.keys():
    if UV[key] == 'NO':
        if genotype2[key] == 'WT':
            sum = sum + LOH_neighbors[key]
        elif genotype2[key] == 'OGG':
            sum2 = sum2 + LOH_neighbors[key]
    elif UV[key] == 'UV':
        if genotype2[key] == 'WT':
            sum3 = sum3 + LOH_neighbors[key]
        elif genotype2[key] == 'OGG':
            sum4 = sum4 + LOH_neighbors[key]

print(sum/WT_LOH)
print(sum2/ogg1_LOH)
print(sum3/UV_WT_LOH)
print(sum4/UV_ogg1_LOH)

sum5 = 0
sum6 = 0
sum7 = 0
sum8 = 0
      
for key in other_neighbors.keys():
    if UV[key] == 'NO':
        if genotype2[key] == 'WT':
            sum5 = sum5 + other_neighbors[key]
        elif genotype2[key] == 'OGG':
            sum6 = sum6 + other_neighbors[key]
    elif UV[key] == 'UV':
        if genotype2[key] == 'WT':
            sum7 = sum7 + other_neighbors[key]
        elif genotype2[key] == 'OGG':
            sum8 = sum8 + other_neighbors[key]
    
        
print(sum5/WT_counts)
print(sum6/ogg1_counts)
print(sum7/UV_WT_counts)
print(sum8/UV_ogg1_counts)


sum9 = 0
sum10 = 0
sum11 = 0
sum12 = 0
sum13 = 0
sum14 = 0
sum15 = 0
sum16 = 0

for key in LOH_neighbors2.keys():
    if genotype[key] == 'WT':
        sum9 = sum9 + LOH_neighbors2[key]
    elif 'rad16D' in  genotype[key]:
        sum10 = sum10 + LOH_neighbors2[key]
    elif 'rad26D' in genotype[key]:
        sum11 = sum11 + LOH_neighbors2[key]
    elif 'rad30D' in genotype[key]:
        sum12 = sum12 + LOH_neighbors2[key]

print(sum9/WTcounts)
print(sum10/rad16counts)
print(sum11/rad26counts)
print(sum12/rad30counts)

for key in other_neighbors2.keys():
    if genotype[key] == 'WT':
        sum13 = sum13 + other_neighbors2[key]
    elif 'rad16D' in genotype[key]:
        sum14 = sum14 + other_neighbors2[key]
    elif 'rad26D' in genotype[key]:
        sum15 = sum15 + other_neighbors2[key]
    elif 'rad30D' in genotype[key]:
        sum16 = sum16 + other_neighbors2[key]
    
        
print(sum13/(12199 - WTcounts))
print(sum14/(23600 - rad16counts))
print(sum15/(9948 - rad26counts))
print(sum16/(16684 - rad30counts))

def WriteData(chr_name, f1, f2, genotype, chromosome, position1, position2,
             trinucleotide, allele, mutation, mutation_type, frequency):
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(0, len(mutation_type)):
        if chromosome[x] == chr_name and trinucleotide[x] != 'NA' and UV[x] == 'UV' and frequency[x] == 'Homozygous':
                chromosome1_dict[int(position1[x])] = x
    for key1 in chromosome1_dict.keys():
        chromosome1.append(int(key1))
    chromosome1.sort()
    for y in range(0, len(chromosome1)):
        if mutation_type[chromosome1_dict[chromosome1[y]]] == 'SNV' or mutation_type[chromosome1_dict[chromosome1[y]]] == 'Deletion':
            if allele[chromosome1_dict[chromosome1[y]]] == 'A' or allele[chromosome1_dict[chromosome1[y]]] == 'G':
                file_data.append(chromosome[chromosome1_dict[chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[chromosome1[y]]])) + '\t' + str(int(position2[chromosome1_dict[chromosome1[y]]]))+ '\t' + GetReverseComplement(trinucleotide[chromosome1_dict[chromosome1[y]]])  + '\t' + GetReverseComplement(mutation[chromosome1_dict[chromosome1[y]]]) + '\t-')
            elif allele[chromosome1_dict[chromosome1[y]]] == 'C' or allele[chromosome1_dict[chromosome1[y]]] == 'T':
                file_data.append(chromosome[chromosome1_dict[chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[chromosome1[y]]]) ) + '\t' + str(int(position2[chromosome1_dict[chromosome1[y]]])) + '\t' + trinucleotide[chromosome1_dict[chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[chromosome1[y]]] + '\t+')
    for fd in range(0, len(file_data)):
        if genotype[chromosome1_dict[chromosome1[fd]]] == 'OGG':
            if mutation_type[chromosome1_dict[chromosome1[fd]]] == 'SNV':
                f1.write(file_data[fd])
                f1.write('\n')
            if mutation_type[chromosome1_dict[chromosome1[fd]]] == 'Deletion' and len(allele[chromosome1_dict[chromosome1[fd]]]) == 1:
                f2.write(file_data[fd])
                f2.write('\n')



f11 = open('ogg1_LOH_TriNuc_substitutions.bed', 'w+')
f12 = open('ogg1_LOH_TriNuc_deletions.bed', 'w+')


f1 = open("LOH_substitutions.bed", 'w+')
f2 = open("LOH_WT_substitutions.bed", 'w+')
f3 = open("LOH_rad16_substitutions.bed", 'w+')
f4 = open("LOH_rad26_substitutions.bed", 'w+')
f5 = open("LOH_rad30_substitutions.bed", 'w+')

WriteChr('chrI', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrII', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrIII', f1, f2, f3, f4, f5,chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrIV', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrIX', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrV', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrVI', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrVII', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrVIII', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrX', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXI', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXII', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXIII', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXIV', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXV', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)
WriteChr('chrXVI', f1, f2, f3, f4, f5, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type)


WriteData('chrI', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrII', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrIII', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrIV', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrIX', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrV', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrVI', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrVII', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrVIII', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrX', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrXI', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrXII', f11, f12, genotype2,  chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrXIII', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrXIV', f11, f12, genotype2,  chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrXV', f11, f12, genotype2,  chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)
WriteData('chrXVI', f11, f12, genotype2, chromosome2, position1_2, position2_2, trinucleotide2, allele2, mutation2, mutation_type_2, frequency2)

f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()
f7.close()
f8.close()
f9.close()
f10.close()
f11.close()
f12.close()

f6.close()
f7.close()