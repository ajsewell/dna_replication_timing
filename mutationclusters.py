def GetReverseComplement(str):
    if '-' in str:
        return str
    dict = {}
    dict['A'] = 'T'
    dict['C'] = 'G'
    dict['G'] = 'C'
    dict['T'] = 'A'
    dict['-'] = '-'
    str2 = ''
    for i in range(0, len(str)):
        str2 = str2 + dict[str[i]]
    return str2[::-1]

def GetSequence(file):
    f1 = open(file)
    sequence1 = f1.read()
    sequence1 = str(sequence1.strip().upper())
    sequence1 = ''.join(sequence1.split('\n'))
    return sequence1




with open("UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem.txt") as f:
    data = f.read()
lines = [s.strip().split() for s in data.splitlines()]
genotypes = {'WT': 0, 'rad16':0, 'rad26':0, 'rad30': 0}
WT_chromosomes = {}
rad16_chromosomes = {}
rad26_chromosomes = {}
rad30_chromosomes = {}
WT_indels = 0
rad16_indels = 0
rad26_indels = 0
rad30_indels = 0
genotype = []
background = []
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
strand = []
isolate_list = {}
lines.pop(0)
for line in lines:
   
    if 'Mito' not in line[6] and 'ySR128' not in line[3]:
    
        genotype.append(line[4])
    
        mutation_type.append(line[10])
        if 'WT' in line[4]:
            genotypes['WT'] = genotypes['WT'] + 1
        if 'rad16' in line[4]:
            genotypes['rad16'] = genotypes['rad16'] + 1
        if 'rad26' in line[4]:
            genotypes['rad26'] = genotypes['rad26'] + 1
        if 'rad30' in line[4]:
            genotypes['rad30'] = genotypes['rad30'] + 1
        chromosome.append('chr' + str(line[6]))
        if 'WT' in line[4]:
            if str(line[6]) in WT_chromosomes.keys():
                WT_chromosomes[str(line[6])] = WT_chromosomes[str(line[6])] + 1

            else:
                WT_chromosomes[str(line[6])] = 1
            if line[10] == "DIP" or line[10] == "Complex DIP":
                WT_indels = WT_indels + 1
        if 'rad16' in line[4]:
            if str(line[6]) in rad16_chromosomes.keys():
                rad16_chromosomes[str(line[6])] = rad16_chromosomes[str(line[6])] + 1
            else:
                rad16_chromosomes[str(line[6])] = 1
            if line[10] == "DIP" or line[10] == "Complex DIP":
                rad16_indels = rad16_indels + 1
        if 'rad26' in line[4]:
            if str(line[6]) in rad26_chromosomes.keys():
                rad26_chromosomes[str(line[6])] = rad26_chromosomes[str(line[6])] + 1
            else:
                rad26_chromosomes[str(line[6])] = 1
            if line[10] == "DIP" or line[10] == "Complex DIP":
                rad26_indels = rad26_indels + 1
        if 'rad30' in line[4]:
            if str(line[6]) in rad30_chromosomes.keys():
                rad30_chromosomes[str(line[6])] = rad30_chromosomes[str(line[6])] + 1
            else:
                rad30_chromosomes[str(line[6])] = 1
            if line[10] == "DIP" or line[10] == "Complex DIP":
                rad30_indels = rad30_indels + 1
        position1.append(line[8])
        position2.append(line[9])
        length.append(line[11])
        allele.append(line[12])
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
                if line[24][1] == "G" or line[24][1] == "A":
                    trinucleotide.append(GetReverseComplement(line[24]))
                    strand.append('-')
                    if '-' not in line[13]:
                        mutation.append(GetReverseComplement(line[13]))
                    else:
                        mutation.append(line[13])
                else:
                    trinucleotide.append(line[24])
                    mutation.append(line[13])
                    strand.append('+')
            else:
                if line[27][1] == "A" or line[27][1] == "G":
                    trinucleotide.append(GetReverseComplement(line[27]))
                    strand.append('-')
                    if '-' not in line[13]:
                        mutation.append(GetReverseComplement(line[13]))
                    else:
                        mutation.append(line[13])
                else:
                    trinucleotide.append(line[27])
                    mutation.append(line[13])
                    strand.append('+')

    
        else:
            frequency.append('NA')
            left_reference.append('NA')
            right_reference.append('NA')
            left_consensus.append('NA')
            right_consensus.append('NA')
            trinucleotide.append('NA')
        if (line[0], line[4]) not in isolate_list.keys():
            isolate_list[(line[0],  line[4])] = [(line[9], 'chr' + str(line[6]), trinucleotide[len(trinucleotide) - 1], mutation[len(mutation) - 1], strand[len(strand) - 1], mutation_type[len(mutation_type) - 1], line[0], line[10])]
        else:
            isolate_list[(line[0],  line[4])].append((line[9], 'chr'+ str(line[6]),  trinucleotide[len(trinucleotide) - 1], mutation[len(mutation) - 1], strand[len(strand) - 1], mutation_type[len(mutation_type) - 1], line[0], line[10]))

print(genotypes)
print(WT_indels/genotypes['WT'])
print(rad16_indels/genotypes['rad16'])
print(rad26_indels/genotypes['rad26'])
print(rad30_indels/genotypes['rad30'])


cluster_file = open('clusters.txt', 'w+')


def FindClusters( isolate_list, strain):
    lengths = {}
    clusters = {}
    homopolymers = {}
    cluster_numbers = {}
    for key in isolate_list.keys():
        number = 0
        if strain in key[1]:
            if isolate_list[key][0][1] == isolate_list[key] [1][1] and (int(isolate_list[key][1][0]) <= int((isolate_list[key][ 0][0])) + 10):
                clusters[isolate_list[key][0]] = [isolate_list[key][0]]
                clusters[isolate_list[key][1]] = [isolate_list[key][1]]
            if isolate_list[key][len (isolate_list[key]) - 2][1] == isolate_list[key] [len (isolate_list[key]) - 1][1] and (int(isolate_list[key][len (isolate_list[key]) - 1][0]) <= int((isolate_list[key][ len (isolate_list[key]) - 2][0])) + 10):
                if isolate_list[key][len(isolate_list[key]) - 1] not in clusters.keys():
                    clusters[isolate_list[key][len(isolate_list[key]) - 1]] = [isolate_list[key][len(isolate_list[key]) - 1]]
                if isolate_list[key][len(isolate_list[key]) - 2] not in clusters.keys():
                    clusters[isolate_list[key][len(isolate_list[key]) - 2]] = [isolate_list[key][len(isolate_list[key]) - 2]]
            for j in range (0,len (isolate_list[key]) - 1):
                #flag = True
                #n = 1
                #while flag == True:
                #if j <= len(isolate_list[key]) - 1:
                    if (isolate_list[key][j][1] == isolate_list[key][j - 1][1] and (int(isolate_list[key][j][0]) <= int((isolate_list[key][j - 1][0])) + 10)) or (isolate_list[key][j + 1][1] == isolate_list[key][j][1] and (int(isolate_list[key][j + 1][0]) <= int((isolate_list[key][j][0])) + 10)): #or ((isolate_list[key][j + n][6], int(isolate_list[key][j + n][0]), isolate_list[key][j + n][1]) in homopolymers.keys() and (int(isolate_list[key][j + n][0]) <= int((isolate_list[key][j + n - 1][0])) + 10 + homopolymers[key])) : 
                            if isolate_list[key][j] not in clusters.keys():# and isolate_list[key][j][7] == "DIP":
                                #if j == 0 or isolate_list[key][j - n] not in clusters.keys():
                                clusters[isolate_list[key][j]] = [isolate_list[key][j]]
                                    #number = number + 1
            
                                    #if isolate_list[key][j + n][7] == "DIP":
                                    #clusters[isolate_list[key][j]].append(isolate_list[key][j + n])
                                    
                                      
                                    #else:
                                    #n = n + 1
                                #else:
                                    #flag = False
                            #else:#elif isolate_list[key][j + n] == "DIP":
                               # clusters[isolate_list[key][j]].append(isolate_list[key][j + n])
                                    #number = number + 1
                                #n = n + 1

                           
                        #else:
                            #flag = False
                    #else:
                        #flag = False
        #cluster_numbers[key] = number
    for key in clusters.keys():
        if len(clusters[key]) in lengths.keys():
            lengths[len(clusters[key])] = lengths[len(clusters[key])] + 1
        else:
            lengths[len(clusters[key])] = 1
        if len(clusters[key]) > 5:
            print(str(clusters[key][0]) + ": " + str(len(clusters[key])))
    for key in lengths.keys():
        print(str(key) + ": " + str(lengths[key]))
    for key in cluster_numbers.keys():
        print(str(key) + ": " + str(cluster_numbers[key]/ len(isolate_list[key])))
    new_clusters = {}
    for key in clusters.keys():
        for i in range(0, len(clusters[key])):
            if (clusters[key][i][6], clusters[key][i][1], clusters[key][i][0]) not in new_clusters.keys():
                new_clusters[(clusters[key][i][6], clusters[key][i][1], clusters[key][i][0])] = [(clusters[key][i][6], clusters[key][i][1], clusters[key][i][0])]



    return clusters




def WriteChr(chr_name, clusters, file):
    
    counts = 0
    
    chromosome1 =[]

    for key in clusters.keys():
        if key[1] == chr_name:
            for a in range(0, len(clusters[key])):
                chromosome1.append(((int(clusters[key][a][0])), clusters[key][a][2], clusters[key][a][3], clusters[key][a][4], clusters[key][a][6]))
                counts = counts + 1

    chromosome1.sort()
    for y in range(0, len(chromosome1)):
        if len(chromosome1[y][1]) == 3 and '-' not in chromosome1[y][1]:
            file.write(chr_name  + '\t'   + str(chromosome1[y][0] - 1) + '\t' + str(chromosome1[y][0]) + '\t' + chromosome1[y][1] + '\t'+ chromosome1[y][2] + '\t'+ chromosome1[y][3])
            file.write('\n')  
    return counts

clusters = FindClusters( isolate_list, 'WT')
clusters16 = FindClusters(isolate_list, 'rad16')
clusters26 = FindClusters(isolate_list, 'rad26')
clusters30 = FindClusters(isolate_list, 'rad30')
for key in clusters:
    for a in  range(0, len(clusters[key])):
        cluster_file.write(str(clusters[key][a][0]) + str(clusters[key][a][1]))
    cluster_file.write('\n')

file = open('clusters.bed', 'w+')
file16 = open('clusters16.bed', 'w+')
file30 = open('clusters30.bed', 'w+')
file26 = open('clusters26.bed', 'w+')

WT_totals = []
rad16_totals = []
rad26_totals = []
rad30_totals = []

for key in WT_chromosomes.keys():
    WT_totals.append(WT_chromosomes[key])
for key in rad16_chromosomes.keys():
    rad16_totals.append(rad16_chromosomes[key])
for key in rad26_chromosomes.keys():
    rad26_totals.append(rad26_chromosomes[key])
for key in rad30_chromosomes.keys():
    rad30_totals.append(rad30_chromosomes[key])



count_file = open('cluster_counts', 'w+')
count_file16 = open('cluster_counts16', 'w+')
count_file26 = open('cluster_counts26', 'w+')
count_file30 = open('cluster_counts30', 'w+')
def GetCounts(clusters, file, count_file, totals):
    counts = []
    counts.append(WriteChr('chrI', clusters, file))
    counts.append(WriteChr('chrII', clusters, file))
    counts.append(WriteChr('chrIII', clusters, file))
    counts.append(WriteChr('chrIV', clusters, file))
    counts.append(WriteChr('chrIX', clusters, file))
    counts.append(WriteChr('chrV', clusters, file))
    counts.append(WriteChr('chrVI', clusters, file))
    counts.append(WriteChr('chrVII', clusters, file))
    counts.append(WriteChr('chrVIII', clusters, file))
    counts.append(WriteChr('chrX', clusters, file))
    counts.append(WriteChr('chrXI', clusters, file))
    counts.append(WriteChr('chrXII', clusters, file))
    counts.append(WriteChr('chrXIII', clusters, file))
    counts.append(WriteChr('chrXIV', clusters, file))
    counts.append(WriteChr('chrXV', clusters, file))
    counts.append(WriteChr('chrXVI', clusters, file))
    file.close()
    for i in range(0, len(counts)):
        count_file.write('Chromosome ' + str(i) + ':' + str(counts[i]/totals[i]))
        count_file.write('\n')
    count_file.close()

GetCounts(clusters, file, count_file, WT_totals)
GetCounts(clusters16, file16, count_file16, rad16_totals)
GetCounts(clusters26, file26, count_file26, rad26_totals)
GetCounts(clusters30, file30, count_file30, rad30_totals)

cluster_file.close()

indel_file = open('indel_list', 'w+')

with open("UVB_New_Combined_All_v2_format_filter.txt") as f2:
    data2 = f2.read()
lines2 = [s.strip().split() for s in data2.splitlines()]
lines2.pop(0)
genotypes2 = {'WT': 0, 'Rad16':0, 'Rad26':0, 'Rad30': 0}
WT_indels_2 = 0
rad16_indels_2 = 0
rad26_indels_2 = 0
rad30_indels_2 = 0
isolates2 = []
genotype2 = []
chromosome2 = []
position1_2 = []
position2_2 = []
mutation_type_2 = []
length2 = []
allele2 = []
mutation2 = []
frequency2 = []
trinucleotide2 = []
strand2 = []
isolate_list2 = {}
WT_chromosomes2 = {}
rad16_chromosomes2 = {}
rad26_chromosomes2 = {}
rad30_chromosomes2 = {}
for line in lines2:
    if line[13] != 'Replacement' and line[13] != 'MNV' and line[14] != 'Replacement' and line[14] != 'MNV' and 'Mito' not in line[5] and 'Mito' not in line[4]:
        isolates2.append(line[0])
        genotype2.append(line[1])
        if line[2] != 'RP':
            mutation_type_2.append(line[13])
        else:
            mutation_type_2.append(line[14])
        if 'WT' in line[1]:
            genotypes2['WT'] = genotypes2['WT'] + 1
        if 'Rad16' in line[1]:
            genotypes2['Rad16'] = genotypes2['Rad16'] + 1
        if 'Rad26' in line[1]:
            genotypes2['Rad26'] = genotypes2['Rad26'] + 1
        if 'Rad30' in line[1]:
            genotypes2['Rad30'] = genotypes2['Rad30'] + 1
        if line[2] != 'RP':
            chromosome2.append(line[4])
        else:
            chromosome2.append(line[5])
        if 'WT' in line[1]:
            if chromosome2[len(chromosome2) - 1] in WT_chromosomes2.keys():
                WT_chromosomes2[chromosome2[len(chromosome2) - 1]] = WT_chromosomes2[chromosome2[len(chromosome2) - 1]] + 1

            else:
                WT_chromosomes2[chromosome2[len(chromosome2) - 1]] = 1
            if line[13] == "Deletion" or line[13] == "Insertion" or line[14] == "Deletion" or line[14] == "Insertion":
                WT_indels_2 = WT_indels_2 + 1
        if 'Rad16' in line[1]:
            if chromosome2[len(chromosome2) - 1] in rad16_chromosomes2.keys():
                rad16_chromosomes2[chromosome2[len(chromosome2) - 1]] = rad16_chromosomes2[chromosome2[len(chromosome2) - 1]] + 1
            else:
                rad16_chromosomes2[chromosome2[len(chromosome2) - 1]] = 1
            if line[13] == "Deletion" or line[13] == "Insertion" or line[14] == "Deletion" or line[14] == "Insertion":
                rad16_indels_2 = rad16_indels_2 + 1
                indel_file.write(line[13])
                indel_file.write('\n')
        if 'Rad26' in line[1]:
            if chromosome2[len(chromosome2) - 1] in rad26_chromosomes2.keys():
                rad26_chromosomes2[chromosome2[len(chromosome2) - 1]] = rad26_chromosomes2[chromosome2[len(chromosome2) - 1]] + 1
            else:
                rad26_chromosomes2[chromosome2[len(chromosome2) - 1]] = 1
            if line[13] == "Deletion" or line[13] == "Insertion" or line[14] == "Deletion" or line[14] == "Insertion":
                rad26_indels_2 = rad26_indels_2 + 1
        if 'Rad30' in line[1]:
            if chromosome2[len(chromosome2) - 1] in rad30_chromosomes2.keys():
                rad30_chromosomes2[chromosome2[len(chromosome2) - 1]] = rad30_chromosomes2[chromosome2[len(chromosome2) - 1]] + 1
            else:
                rad30_chromosomes2[chromosome2[len(chromosome2) - 1]] = 1
            if line[13] == "Deletion" or line[13] == "Insertion" or line[14] == "Deletion" or line[14] == "Insertion":
                rad30_indels_2 = rad30_indels_2 + 1
        if line[2] != 'RP':
            position1_2.append(int(line[6]) - 1)
            position2_2.append(int(line[6]))
            allele2.append(line[14])
            if line[14] == 'C' or line[14] == 'T':
                mutation2.append(line[15])
                strand2.append('+')
            else:
                mutation2.append(GetReverseComplement(line[15]))
                strand2.append('-')
        else:
            position1_2.append(int(line[7]) - 1)
            position2_2.append(int(line[7]))
            allele2.append(line[15])
            if line[15] == 'C' or line[15] == 'T':
                mutation2.append(line[16])
                strand2.append('+')
            else:
                mutation2.append(GetReverseComplement(line[16]))
                strand2.append('-')


for pos in range(0, len(position2_2)):
    if mutation_type_2[pos] == 'SNV' or mutation_type_2[pos] == 'Deletion' and len(allele2[pos]) == 1:
        if chromosome2[pos] == 'chrI':
            sequence = GetSequence('chr1.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))

        elif chromosome2[pos] == 'chrII':
            sequence = GetSequence('chr2.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrIII':
            sequence = GetSequence('chr3.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrIV':
            sequence = GetSequence('chr4.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrV':
            sequence = GetSequence('chr5.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrVI':
            sequence = GetSequence('chr6.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrVII':
            sequence = GetSequence('chr7.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrVIII':
            sequence = GetSequence('chr8.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrIX':
            sequence = GetSequence('chr9.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrX':
            sequence = GetSequence('chr10.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrXI':
            sequence = GetSequence('chr11.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrXII':
            sequence = GetSequence('chr12.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrXIII':
            sequence = GetSequence('chr13.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrXIV':
            sequence = GetSequence('chr14.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrXV':
            sequence = GetSequence('chr15.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        elif chromosome2[pos] == 'chrXVI':
            sequence = GetSequence('chr16.txt')
            if strand2[pos] == '+':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
        else:
            trinucleotide2.append('NA')
    else:
        trinucleotide2.append('NA')
        


for i in range (0, len(isolates2)):
    if (isolates2[i], genotype2[i]) not in isolate_list2.keys():
        isolate_list2[(isolates2[i], genotype2[i])] = [(position2_2[i], chromosome2[i], trinucleotide2[i], mutation2[i], strand2[i], mutation_type_2[i], isolates2[i])]
    else:
        isolate_list2[(isolates2[i], genotype2[i])].append((position2_2[i], chromosome2[i], trinucleotide2[i], mutation2[i], strand2[i], isolates2[i] ,mutation_type_2[i], isolates2[i]))

print(genotypes2)

cluster_file2 = open('clusters_UVB_singlebase.txt', 'w+')
clusters2 = FindClusters( isolate_list2, 'WT')
clusters16_2 = FindClusters(isolate_list2, 'Rad16' )
clusters26_2 = FindClusters(isolate_list2, 'Rad26' )
clusters30_2 = FindClusters(isolate_list2, 'Rad30' )
for key in clusters2:
    for a in  range(0, len(clusters2[key])):
        cluster_file2.write(str(clusters2[key][a][0]) + str(clusters2[key][a][1]))
    cluster_file2.write('\n')

indel_file.close()
file2 = open('clusters_UVB.bed', 'w+')
file16_2 = open('clusters16_UVB.bed', 'w+')
file30_2 = open('clusters30_UVB.bed', 'w+')
file26_2 = open('clusters26_UVB.bed', 'w+')

WT_totals2 = []
rad16_totals2 = []
rad26_totals2 = []
rad30_totals2 = []

count_file2 = open('cluster_counts_UVB', 'w+')
count_file16_2 = open('cluster_counts16_UVB', 'w+')
count_file26_2 = open('cluster_counts26_UVB', 'w+')
count_file30_2 = open('cluster_counts30_UVB', 'w+')

for key in WT_chromosomes2.keys():
    WT_totals2.append(WT_chromosomes2[key])
for key in rad16_chromosomes2.keys():
    rad16_totals2.append(rad16_chromosomes2[key])
for key in rad26_chromosomes2.keys():
    rad26_totals2.append(rad26_chromosomes2[key])
for key in rad30_chromosomes2.keys():
    rad30_totals2.append(rad30_chromosomes2[key])


GetCounts(clusters2, file2, count_file2, WT_totals2)
GetCounts(clusters16_2, file16_2, count_file16_2, rad16_totals2)
GetCounts(clusters26_2, file26_2, count_file26_2, rad26_totals2)
GetCounts(clusters30_2, file30_2, count_file30_2, rad30_totals2)

cluster_file2.close()
file2.close()
file16_2.close()
file30_2.close()
file26_2.close()


