def GetReverseComplement(str):
    if str == 'NA':
        return 
    dict = {}
    dict['A'] = 'T'
    dict['C'] = 'G'
    dict['G'] = 'C'
    dict['T'] = 'A'
    str2 = ''
    for i in range(0, len(str)):
        str2 = str2 + dict[str[i]]
    return str2[::-1]

with open("WT_UV_bothruns_Muts_allSNVs_sorted.bed") as f:
    data = f.read()
with open("rad16_UV_bothruns_Muts_allSNVs_sorted.bed") as f16:
    data16 = f16.read()
with open("rad30_UV_bothruns_Muts_allSNVs_sorted.bed") as f30:
    data30 = f30.read()
with open("WT_UVA_Muts_allSNVs_sorted.bed") as fA:
    dataUVA = fA.read()
with open("ogg1_UVA_Muts_allSNVs_sorted.bed") as fOg:
    dataOg = fOg.read()
with open("rad26_UV_bothruns_Muts_allSNVs_sorted.bed") as f26:
    data26 = f26.read()
with open("WT_UVB_Muts_allSNVs_sorted.bed") as fB:
    dataUVB = fB.read()
with open("Rad16_UVB_Muts_allSNVs_sorted.bed") as fB16:
    dataB16 = fB16.read()
with open('Rad26_UVB_Muts_allSNVs_sorted.bed') as fB26:
    dataB26 = fB26.read()
with open('Rad30_UVB_Muts_allSNVs_sorted.bed') as fB30:
    dataB30 = fB30.read()
with open("WT_UV_run1_Muts_allSNVs_sorted.bed") as fr1:
    datar1 = fr1.read()
with open("WT_UV_run2_Muts_allSNVs_sorted.bed") as fr2:
    datar2 = fr2.read()

#extract columns from bed files to make arrays
# with matching indices   
def getChromosome(data):       
    lines = [s.strip().split() for s in data.splitlines()]
    chr_counts = {}
    chromosome = []
    for i in range (0, len(lines)):
        chromosome.append(lines[i][0])
        if lines[i][0] not in chr_counts.keys():
            chr_counts[lines[i][0]] = 1
        else:
            chr_counts[lines[i][0]] = chr_counts[lines[i][0]] + 1
    for key in chr_counts.keys():
        print(key + ": " + str(chr_counts[key]))
    return chromosome
chromosome = getChromosome(data)
chromosome16 = getChromosome(data16)
chromosome30 = getChromosome(data30)
chromosomeA = getChromosome(dataUVA)
chromosomeOg = getChromosome(dataOg)
chromosome26 = getChromosome(data26)
chromosomeB = getChromosome(dataUVB)
chromosomeB16 = getChromosome(dataB16)
chromosomeB26 = getChromosome(dataB26)
chromosomeB30 = getChromosome(dataB30)
chromosomer1 = getChromosome(datar1)
chromosomer2 = getChromosome(datar2)

def getPosition1(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    position1 = []
    for j in range (0, len(lines)):
        position1.append(lines[j][1])
        position1[j] = float(position1[j])
    return position1
position1 = getPosition1(data)
position1_16 = getPosition1(data16)
position1_30 = getPosition1(data30)
position1_A = getPosition1(dataUVA)
position1_Og = getPosition1(dataOg)
position1_26 = getPosition1(data26)
position1_B = getPosition1(dataUVB)
position1_B16 = getPosition1(dataB16)
position1_B26 = getPosition1(dataB26)
position1_B30 = getPosition1(dataB30)
position1_r1 = getPosition1(datar1)
position1_r2 = getPosition1(datar2)


def getPosition2(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    position2 = []
    for k in range (0, len(lines) ):
        position2.append(lines[k][2])
        position2[k] = float(position2[k])
    return position2
position2 = getPosition2(data)
position2_16 = getPosition2(data16)
position2_30 = getPosition2(data30)
position2_A = getPosition2(dataUVA)
position2_Og = getPosition2(dataOg)
position2_26 = getPosition2(data26)
position2_B = getPosition2(dataUVB)
position2_B16 = getPosition2(dataB16)
position2_B26 = getPosition2(dataB26)
position2_B30 = getPosition2(dataB30)
position2_r1 = getPosition2(datar1)
position2_r2 = getPosition2(datar2)


def getTrinucleotide(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    trinucleotide = []
    for m in range (0, len(lines)):
        trinucleotide.append(lines[m][3].replace('-', ''))
    return trinucleotide
trinucleotide = getTrinucleotide(data)
trinucleotide16 = getTrinucleotide(data16)
trinucleotide30 = getTrinucleotide(data30)
trinucleotideA = getTrinucleotide(dataUVA)
trinucleotideOg = getTrinucleotide(dataOg)
trinucleotide26 = getTrinucleotide(data26)
trinucleotideB = getTrinucleotide(dataUVB)
trinucleotideB16 = getTrinucleotide(dataB16)
trinucleotideB26 = getTrinucleotide(dataB26)
trinucleotideB30 = getTrinucleotide(dataB30)
trinucleotider1 = getTrinucleotide(datar1)
trinucleotider2 = getTrinucleotide(datar2)

with open("Raghuraman_replication_timing_saccer3.bed") as t:
    data2 = t.read()
lines2 = [u.strip().split() for u in data2.splitlines()]
chr = []
for v in range (0, len(lines2)):
    chr.append(lines2[v][0])
       
    
start = []
for w in range (0, len(lines2)):
    start.append(lines2[w][2])

#scores need to be sorted in order to 
#separate into early, late, and middle replication times
#but unsorted to assign times to mutations
#since indices must match
scores = []
sorted_scores = []
for x in range (0, len(lines2)):
    scores.append(lines2[x][3])
    sorted_scores.append(lines2[x][3])
    
    

for sc in range(0, len(scores)):
    scores[sc] = float(scores[sc])

for ss in range(0, len(sorted_scores)):
    sorted_scores[ss] = float(sorted_scores[ss])

sorted_scores.sort()

#Define 1/3 of genome as early replicating,
#1/3 as middle replicating, and 1/3 as late 
#Separate further into 3 bins each containing 1/9
#of the genome
early = {}
early1 = {}
early2 = {}
early3 = {}
middle = {}
middle1 = {}
middle2 = {}
middle3 = {}
late = {}
late1 =  {}
late2 = {}
late3 = {}

fe = open('early_regions.bed', 'w+')
fm = open('middle_regions.bed', 'w+')
fl = open('late_regions.bed', 'w+')



for y in range(0, len(scores)):
    if y <= len(scores)/3:
        early[y] = sorted_scores[y]
        if y <= len(scores)/9:
            early1[y] = sorted_scores[y]
        elif y <= (2* len(scores))/ 9:
            early2[y] = sorted_scores[y]
        elif y <= len(scores)/3:
            early3[y] = sorted_scores[y]
early_1 = 0
early1_1 = early1[int(len(scores)/9)]
early2_1 = early2[int((2* len(scores))/ 9)]
early_2 = early[len(early) - 1]
for z in range(0, len(scores)):
    if z >=  (2 * len(scores))/3 and z <= len(scores):
        late[z] = sorted_scores[z]
        if z >= (8 * len(scores))/9:
            late3[z] = sorted_scores[z]
        elif z >= (7 * len(scores))/9:
             late2[z] = sorted_scores[z]
        elif z >=  (2 * len(scores))/3:
            late1[z] = sorted_scores[z] 
       
late_1 = late[int(2 * len(scores)/3 + 1)]
late2_1 = late2[int((7 * len(scores))/9)+ 1]
late3_1 = late3[int((8 * len(scores))/9) + 1]
late_2 = late[int(len(scores)- 1)]
for q in range(0, len(scores) ):
    if q > len(scores)/3 and q < 2 * len(scores)/3:
        middle[q] = sorted_scores[q]
        if q > len(scores)/3 and q <= (4 * len(scores))/9:
            middle1[q] = sorted_scores[q]
        elif q > len(scores)/3 and q <= (5 * len(scores))/9:
            middle2[q] = sorted_scores[q]
        elif q > len(scores)/3 and q < (2 * len(scores))/3:
            middle3[q] = sorted_scores[q]
middle_1 = early_2 
middle1_1 = middle1[int((4 * len(scores))/9)]
middle2_1 = middle2[int((5 * len(scores))/9)]
middle_2 = late_1


for s in range(0, len(scores)):
    if scores[s] <= early_2:
        fe.write(chr[s] + '\t' + str(int(start[s]) - 1) + '\t' + str(start[s]) + '\t' + str(scores[s]))
        fe.write('\n')
    if scores[s] <= middle_2:
        fm.write(chr[s] + '\t' + str(int(start[s]) - 1) + '\t' + str(start[s]) + '\t' + str(scores[s]))
        fm.write('\n')
    else:
        fl.write(chr[s] + '\t' + str(int(start[s]) - 1) + '\t' + str(start[s]) + '\t' + str(scores[s]))
        fl.write('\n')
fe.close()
fm.close()
fl.close()

#percentage of early, middle, and late replicating regions in chromosome
def ChromosomePercentages(chr_name):
    ch_total = 0
    ch_early = 0
    ch_middle = 0
    ch_late = 0
    for s6 in range(0, len(scores)):
        if chr[s6] == chr_name:
            if scores[s6] <= middle_1:
                ch_total = ch_total + 1
                ch_early = ch_early + 1
            if scores[s6] <= late_1:
                ch_total = ch_total + 1
                ch_middle = ch_middle + 1
            else:
                ch_total = ch_total + 1
                ch_late = ch_late + 1
    early_percentage = (ch_early / ch_total) * 100
    print(early_percentage)
    middle_percentage = (ch_middle / ch_total) * 100
    print(middle_percentage)
    late_percentage = (ch_late / ch_total) * 100
    print(late_percentage)
    return ([early_percentage, middle_percentage, late_percentage])
percentages = ChromosomePercentages('chrIX')
print(percentages)


f1 = open("replication_time_definition", 'w+')
f1.write('Replication Time Definitions')
f1.write('\n')
f1.write('Early=' + str(early_1) +'-' + str(early_2))
f1.write('\n')
f1.write('Bin1=' + str(early_1) +'-' + str(early1_1))
f1.write('\n')
f1.write('Bin2=' + str(early1_1) +'-' + str(early2_1))
f1.write('\n')
f1.write('Bin3=' + str(early2_1) +'-' + str(early_2))
f1.write('\n')
f1.write('Middle=' + str(middle_1) +'-' + str(middle_2))
f1.write('\n')
f1.write('Bin4=' + str(middle_1) +'-' + str(middle1_1))
f1.write('\n')
f1.write('Bin5=' + str(middle1_1) +'-' + str(middle2_1))
f1.write('\n')
f1.write('Bin6=' + str(middle2_1) +'-' + str(middle_2))
f1.write('\n')
f1.write('Late=' + str(late_1) +'-' + str(late_2))
f1.write('\n')
f1.write('Bin7=' + str(late_1) +'-' + str(late2_1))
f1.write('\n')
f1.write('Bin8=' + str(late2_1) +'-' + str(late3_1))
f1.write('\n')
f1.write('Bin9=' + str(late3_1) +'-' + str(late_2))

f1.close()

#Determines replication time at which a mutation occurs
def calculateTime(time1, time2, pos, pos1, pos2):
    fraction = (pos -pos1)/(pos2 - pos1)
    time = time1 + (time2- time1) * fraction
    return float(time)

for s1 in range(0, len(start)):
    start[s1] = float(start[s1])

#Goes through position and start arrays in parallel
#to assign replication times for each mutation based on scores
def AssignTimes(chromosome, position2, start, scores):
    times = []
    po = 0
    st = 0

    while po < len(position2)  and st <= len(start) - 2:
        
        if position2[po] >= start[st] and position2[po] <= start[st + 1] and chromosome[po] == chr[st]:
            time1 = (scores[st])
            time2 = (scores[st + 1])
            pos = (position2[po])
            pos1 = (start[st])
            pos2 = (start[st + 1])
            time = calculateTime(time1, time2, pos, pos1, pos2)
            times.append(time)
            po = po + 1

        elif chromosome[po] < chr[st]:
            po = po + 1
            times.append('NA')
        elif chromosome[po] > chr[st]:
            st = st + 1
    
        elif position2[po] >= start[st] :
            st = st + 1
        
        elif position2[po] <= start[st]:
            po = po + 1
            times.append('NA')
    return times

times = AssignTimes(chromosome, position2, start, scores)
times16 = AssignTimes(chromosome16, position2_16, start, scores)
times30 = AssignTimes(chromosome30 , position2_30, start, scores)
times_A = AssignTimes(chromosomeA, position2_A, start, scores)
times_Og = AssignTimes(chromosomeOg, position2_Og, start, scores)
times26 = AssignTimes(chromosome26, position2_26, start, scores)
times_B = AssignTimes(chromosomeB, position2_B, start, scores)
timesB16 = AssignTimes(chromosomeB16, position2_B16, start, scores)
timesB26 = AssignTimes(chromosomeB26, position2_B26, start, scores)
timesB30 = AssignTimes(chromosomeB30, position2_B30, start, scores)
timesr1 = AssignTimes(chromosomer1, position2_r1, start, scores)
timesr2 = AssignTimes(chromosomer2, position2_r2, start, scores)

#Numbers and percentages of mutations in early,middle,
#and late replication times. Adds to dictionaries and separates
#into three files
def CountMutations(filename, times, chromosome, chr_name):
    early_mut = {}
    early1_mut = {}
    early2_mut = {}
    early3_mut = {}
    middle_mut = {}
    middle1_mut = {}
    middle2_mut = {}
    middle3_mut = {}
    late_mut = {}
    late1_mut = {}
    late2_mut = {}
    late3_mut = {}
    chr_early = {}
    chr_middle = {}
    chr_late = {}

    for tm in range(0, len(times) ):
        if times[tm] != 'NA' and times[tm ] >= 0 and times[tm] < early_2:
            early_mut[tm] = times[tm]
            if times[tm ] >= 0 and times[tm] < early1_1:
                early1_mut[tm] = times[tm]
            if times[tm] >= early1_1 and times[tm] < early2_1:
                early2_mut[tm] = times[tm]
            if times[tm] >= early2_1 and times[tm] <early_2:
                early3_mut[tm] = times[tm]
            if chromosome[tm] == chr_name:
                chr_early[tm] = times[tm]
        if times[tm] != 'NA' and times[tm ] >= middle_1 and times[tm] < middle_2:
            middle_mut[tm] = times[tm]
            if times[tm] >= middle_1 and times[tm] < middle1_1:
                middle1_mut[tm] = times[tm]
            if times[tm] >= middle1_1 and times[tm] < middle2_1:
                middle2_mut[tm] = times[tm]
            if times[tm] >= middle2_1 and times[tm] < middle_2:
                middle3_mut[tm] = times[tm]
            if chromosome[tm] == chr_name:
                chr_middle[tm] = times[tm]
        if times[tm] != 'NA' and times[tm ] >= late_1 :
            late_mut[tm] = times[tm]
            if times[tm] >= late_1 and times[tm] < late2_1:
                late1_mut[tm] = times[tm]
            if times[tm] >= late2_1 and times[tm] < late3_1:
                late2_mut[tm] = times[tm]
            if times[tm] >= late3_1:
                late3_mut[tm] = times[tm]
            if chromosome[tm] == chr_name:
                chr_late[tm] = times[tm]
    total = len(early_mut)  + len(middle_mut) + len(late_mut)
    total_chr = len(chr_early)  + len(chr_middle) + len(chr_late)
    if total != 0:
        percent_early = (len(early_mut)/ total) * 100
        percent_middle = (len(middle_mut)/ total) * 100
        percent_late = (len(late_mut)/ total) * 100
    else:
        percent_early = 0
        percent_middle = 0
        percent_late = 0
    if total_chr != 0:
        percent_early_chr = (len(chr_early)/ total_chr) * 100
        percent_middle_chr = (len(chr_middle)/ total_chr) * 100
        percent_late_chr = (len(chr_late)/ total_chr) * 100
    else:
        percent_early_chr = 'NA'
        percent_middle_chr = 'NA'
        percent_late_chr = 'NA'
    f1 = open(filename, 'w+')
    f1.write('Early Mutation Number: '+ str(len(early_mut)) + '  Percentage: ' + str(percent_early) + '%')
    f1.write('\n')
    f1.write('Middle Mutation Number: '+ str(len(middle_mut)) + '  Percentage: ' + str(percent_middle) + '%')
    f1.write('\n')
    f1.write('Late Mutation Number: ' + str(len(late_mut)) + '  Percentage: ' + str(percent_late) + '%') 
    f1.write('\n')
    f1.write('Early Mutation Number for '+ chr_name + ': ' + str(len(chr_early)) + '  Percentage: ' + str(percent_early_chr) + '%')
    f1.write('\n')
    f1.write('Middle Mutation Number for '+ chr_name + ': ' + str(len(chr_middle)) + '  Percentage: ' + str(percent_middle_chr) + '%')
    f1.write('\n')
    f1.write('Late Mutation Number for '+ chr_name + ': ' + str(len(chr_late)) + '  Percentage: ' + str(percent_late_chr) + '%')


    f1.close()
    return [early_mut, middle_mut, late_mut, early1_mut, early2_mut, early3_mut, middle1_mut, middle2_mut, middle3_mut, late1_mut, late2_mut, late3_mut]

muts = CountMutations("WT_mutations", times, chromosome, 'chrXVI')
early_mut = muts[0]
middle_mut = muts[1]
late_mut = muts[2]
muts_16= CountMutations("rad16_mutations", times16, chromosome16, 'chrXVI')
early_mut_16 = muts_16[0]
middle_mut_16 = muts_16[1]
late_mut_16 = muts_16[2]
muts_30 = CountMutations("rad30_mutations", times30, chromosome30, 'chrI')
early_mut_30 = muts_30[0]
middle_mut_30 = muts_30[1]
late_mut_30 = muts_30[2]
muts_A = CountMutations("UVA_mutations", times_A, chromosomeA, 'chrI')
early_mut_A = muts_A[0]
middle_mut_A = muts_A[1]
late_mut_A = muts_A[2]
muts_Og = CountMutations("ogg1_mutations", times_Og, chromosomeOg, 'chrI')
early_mut_Og = muts_Og[0]
middle_mut_Og = muts_Og[1]
late_mut_Og = muts_Og[2]
muts_26= CountMutations("rad26_mutations", times26, chromosome26, 'chrI')
early_mut_26 = muts_26[0]
middle_mut_26 = muts_26[1]
late_mut_26 = muts_26[2]
muts_B = CountMutations("UVB_mutations", times_B, chromosomeB, 'chrXVI')
early_mut_B = muts_B[0]
middle_mut_B = muts_B[1]
late_mut_B = muts_B[2]
muts_B16 = CountMutations("UVB_rad16_mutations", timesB16, chromosomeB16, 'chrXVI')
early_mut_B16 = muts_B16[0]
middle_mut_B16 = muts_B16[1]
late_mut_B16 = muts_B16[2]
muts_B26 = CountMutations("UVB_rad26_mutations", timesB26, chromosomeB26, 'chrI')
early_mut_B26 = muts_B26[0]
middle_mut_B26 = muts_B26[1]
late_mut_B26 = muts_B26[2]
muts_B30 = CountMutations("UVB_rad30_mutations", timesB30, chromosomeB30, 'chrI')
early_mut_B30 = muts_B30[0]
middle_mut_B30 = muts_B30[1]
late_mut_B30 = muts_B30[2]
muts_r1 = CountMutations("WT_r1_mutations", timesr1, chromosomer1, 'chrVI')
early_mut_r1 = muts_r1[0]
middle_mut_r1 = muts_r1[1]
late_mut_r1 = muts_r1[2]
muts_r2 = CountMutations("WT_r2_mutations", timesr2, chromosomer2, 'chrVI')
early_mut_r2 = muts_r2[0]
middle_mut_r2 = muts_r2[1]
late_mut_r2 = muts_r2[2]

#Calculates average times for each of nine bins for use in linear model
def AverageTimes(muts, file):
    early1_mut = muts[3]
    early2_mut = muts[4]
    early3_mut = muts[5]
    middle1_mut = muts[6]
    middle2_mut = muts[7]
    middle3_mut = muts[8]
    late1_mut = muts[9]
    late2_mut = muts[10]
    late3_mut = muts[11]
    sum1 = 0
    for mut1 in early1_mut.keys():
        sum1 = sum1 + early1_mut[mut1]
    average1 = sum1/len(early1_mut)
    sum2 = 0
    for mut2 in early2_mut.keys():
        sum2 = sum2 + early2_mut[mut2]
    average2 = sum2/len(early2_mut)
    sum3 = 0
    for mut3 in early3_mut.keys():
        sum3 = sum3 + early3_mut[mut3]
    average3 = sum3/len(early3_mut)
    sum4 = 0
    for mut4 in middle1_mut.keys():
        sum4 = sum4 + middle1_mut[mut4]
    average4 = sum4/len(middle1_mut)
    sum5 = 0
    for mut5 in middle2_mut.keys():
        sum5 = sum5 + middle2_mut[mut5]
    average5 = sum5/len(middle2_mut)
    sum6 = 0
    for mut6 in middle3_mut.keys():
        sum6 = sum6 + middle3_mut[mut6]
    average6 = sum6/len(middle3_mut)
    sum7 = 0
    for mut7 in late1_mut.keys():
        sum7 = sum7 + late1_mut[mut7]
    average7 = sum7/len(late1_mut)
    sum8 = 0
    for mut8 in late2_mut.keys():
        sum8 = sum8 + late2_mut[mut8]
    average8 = sum8/len(late2_mut)
    sum9 = 0
    for mut9 in late3_mut.keys():
        sum9 = sum9 + late3_mut[mut9]
    average9 = sum9/len(late3_mut)
    f = open(file, 'w+')
    f.write('Bin1: ' + str(average1))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(early1_mut)))
    f.write('\n')
    f.write('Bin2: ' + str(average2))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(early2_mut)))
    f.write('\n')
    f.write('Bin3: ' + str(average3))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(early3_mut)))
    f.write('\n')
    f.write('Bin4: ' + str(average4))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(middle1_mut)))
    f.write('\n')
    f.write('Bin5: ' + str(average5))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(middle2_mut)))
    f.write('\n')
    f.write('Bin6: ' + str(average6))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(middle3_mut)))
    f.write('\n')
    f.write('Bin7: ' + str(average7))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(late1_mut)))
    f.write('\n')
    f.write('Bin8: ' + str(average8))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(late2_mut)))
    f.write('\n')
    f.write('Bin9: ' + str(average9))
    f.write('\n')
    f.write('Mutation Number: ' + str(len(late3_mut)))
    f.close()
AverageTimes(muts, 'average_times.txt')
AverageTimes(muts_16, 'average_times_rad16.txt')
AverageTimes(muts_26, 'average_times_rad26.txt')
AverageTimes(muts_30, 'average_times_rad30.txt')
AverageTimes(muts_B, 'average_times_UVB.txt')
AverageTimes(muts_B16, 'average_times_rad16_UVB.txt')
AverageTimes(muts_B26, 'average_times_rad26_UVB.txt')
AverageTimes(muts_B30, 'average_times_rad30_UVB.txt')
AverageTimes(muts_A, 'average_times_UVA.txt')
AverageTimes(muts_Og, 'average_times_ogg1.txt')

#Uses times array to add a column for replication times to bed files
def AddTimeColumn(filename, file1, file2, file3, data, times, early_mut, middle_mut, late_mut):
    fi = open(filename, 'w+')
    fa = open(file1 , 'w+')
    fb = open(file2, 'w+')
    fc = open(file3, 'w+')   
    e = 0   
    for d in data.splitlines():
        d.strip()
        if e < len(times):
     
            d = d + '\t' + str(times[e])
        
        
            if e in early_mut.keys():
                fa.write(d)
                fa.write('\n')
            if e in middle_mut.keys():
                fb.write(d)
                fb.write('\n')
            if e in late_mut.keys():
                fc.write(d)
                fc.write('\n')
            e = e + 1
        #else:
            #d = d + ' ' + 'NA'
        
        fi.write(d)
        fi.write('\n')
    
   
    fi.close()
    fa.close()
    fb.close()
    fc.close()
AddTimeColumn("WT_UV_bothruns_allSNVs_sorted.bed", "early_mutations.bed", "middle_mutations.bed", "late_mutations.bed",data, times, early_mut, middle_mut, late_mut )
AddTimeColumn("rad16_UV_bothruns_allSNVs_sorted.bed", "early_mutations_rad16.bed", "middle_mutations_rad16.bed", "late_mutations_rad16.bed",data16, times16, early_mut_16, middle_mut_16, late_mut_16 )
AddTimeColumn("rad30_UV_bothruns_allSNVs_sorted.bed", "early_mutations._rad30bed", "middle_mutations_rad30.bed", "late_mutations_rad30.bed",data30, times30, early_mut_30, middle_mut_30, late_mut_30 )
AddTimeColumn("WT_UVA_allSNVs_sorted.bed", "early_mutations_UVA.bed", "middle_mutations_UVA.bed", "late_mutations_UVA.bed",dataUVA, times_A, early_mut_A, middle_mut_A, late_mut_A)
AddTimeColumn("ogg1_UVA_allSNVs_sorted.bed", "early_mutations_ogg1.bed", "middle_mutations_ogg1.bed", "late_mutations_ogg1.bed",dataOg, times_Og, early_mut_Og, middle_mut_Og, late_mut_Og)
AddTimeColumn("rad26_UV_bothruns_allSNVs_sorted.bed", "early_mutations_rad26.bed", "middle_mutations_rad26.bed", "late_mutations_rad26.bed",data26, times26, early_mut_26, middle_mut_26, late_mut_26 )
AddTimeColumn("WT_UVB_allSNVs_sorted.bed", "early_mutations_UVB.bed", "middle_mutations_UVB.bed", "late_mutations_UVB.bed",dataUVB, times_B, early_mut_B, middle_mut_B, late_mut_B)
AddTimeColumn("Rad16_UVB_allSNVs_sorted.bed", "early_mutations_rad16UVB.bed", "middle_mutations_rad16UVB.bed", "late_mutations_rad16UVB.bed",dataB16, timesB16, early_mut_B16, middle_mut_B16, late_mut_B16)
AddTimeColumn("Rad26_UVB_allSNVs_sorted.bed", "early_mutations_rad26UVB.bed", "middle_mutations_rad26UVB.bed", "late_mutations_rad26UVB.bed",dataB26, timesB26, early_mut_B26, middle_mut_B26, late_mut_B26)
AddTimeColumn("Rad30_UVB_allSNVs_sorted.bed", "early_mutations_rad30UVB.bed", "middle_mutations_rad30UVB.bed", "late_mutations_rad30UVB.bed",dataB30, timesB30, early_mut_B30, middle_mut_B30, late_mut_B30)
AddTimeColumn("WT_UV_run1_allSNVs_sorted.bed", "early_mutations_r1.bed", "middle_mutations_r1.bed", "late_mutations_r1.bed",datar1, timesr1, early_mut_r1, middle_mut_r1, late_mut_r1 )
AddTimeColumn("WT_UV_run2_allSNVs_sorted.bed", "early_mutations_r2.bed", "middle_mutations_r2.bed", "late_mutations_r2.bed",datar2, timesr2, early_mut_r2, middle_mut_r2, late_mut_r2 )

#Calculates number of mutations for each codon
def MutationRate(file1, trinucleotides, muts, file2):
    mutationRates = {}
    base_percentages = {}
    mutationRates_early = {}
    mutationRates_middle = {}
    mutationRates_late = {}
    base_percentages_early = {}
    base_percentages_middle = {}
    base_percentages_late = {}
    fm = open(file1)
    datam = fm.read()
    fm2 = open(file2, 'w+')
    for line in datam.splitlines():
        line = datam.strip().split()
    for d in range(0, len(line)):
        mutationRates[line[d]] = 0
        base_percentages[line[d]] = 0
        mutationRates_early[line[d]] = 0
        mutationRates_middle[line[d]] = 0
        mutationRates_late[line[d]] = 0
        base_percentages_early[line[d]] = 0
        base_percentages_middle[line[d]] = 0
        base_percentages_late[line[d]] = 0
    for key in mutationRates.keys():
        for t in range (0, len(trinucleotides)):
            if trinucleotides[t][1] == 'C' or trinucleotides[t][1] == 'T':
                if trinucleotides[t] == key:
                    mutationRates[key] = mutationRates[key] + 1
                    base_percentages[key] = (mutationRates[key] / len(trinucleotides))
                if trinucleotides[t] == key and t in muts[0].keys():
                    mutationRates_early[key] = mutationRates_early[key] + 1
                    base_percentages_early[key] = (mutationRates_early[key] / len(muts[0].keys()))
                if trinucleotides[t] == key and t in muts[1].keys():
                    mutationRates_middle[key] = mutationRates_middle[key] + 1
                    base_percentages_middle[key] = (mutationRates_middle[key] / len(muts[1].keys()))
                if trinucleotides[t] == key and t in muts[2].keys():
                    mutationRates_late[key] = mutationRates_late[key] + 1
                    base_percentages_late[key] = (mutationRates_late[key] / len(muts[2].keys()))
            if trinucleotides[t][1] == 'A' or trinucleotides[t][1] == 'G':
                if GetReverseComplement(trinucleotides[t]) == key:
                    mutationRates[key] = mutationRates[key] + 1
                    base_percentages[key] = (mutationRates[key] / len(trinucleotides))
                    if t in muts[0].keys():
                        mutationRates_early[key] = mutationRates_early[key] + 1
                        base_percentages_early[key] = (mutationRates_early[key] / len(muts[0].keys()))
                    if t in muts[1].keys():
                        mutationRates_middle[key] = mutationRates_middle[key] + 1
                        base_percentages_middle[key] = (mutationRates_middle[key] / len(muts[1].keys()))
                    if t in muts[2].keys():
                        mutationRates_late[key] = mutationRates_late[key] + 1
                        base_percentages_late[key] = (mutationRates_late[key] / len(muts[2].keys()))
             
    for key in mutationRates.keys():
        fm2.write(key + ": Number = " + str(mutationRates[key])+ ' Percentage =' + str(base_percentages[key] * 100) + "%")
        fm2.write('\n')
    fm2.write('Early Replication: ')
    fm2.write('\n')
    for key in mutationRates_early.keys():
        fm2.write(key + ": Number = " + str(mutationRates_early[key])+ ' Percentage =' + str(base_percentages_early[key] * 100) + "%")
        fm2.write('\n')
    fm2.write('Middle Replication: ')
    fm2.write('\n')
    for key in mutationRates_middle.keys():
        fm2.write(key + ": Number = " + str(mutationRates_middle[key])+ ' Percentage =' + str(base_percentages_middle[key] * 100) + "%")
        fm2.write('\n')
    fm2.write('Late Replication: ')
    fm2.write('\n')
    for key in mutationRates_late.keys():
        fm2.write(key + ": Number = " + str(mutationRates_late[key])+ ' Percentage =' + str(base_percentages_late[key] * 100) + "%")
        fm2.write('\n')
    return(mutationRates)

mutationRates = MutationRate('mutationrates.txt', trinucleotide, muts, 'mutationrate_results')
mutationRates16 = MutationRate('mutationrates.txt', trinucleotide16, muts_16, 'mutationrate16_results')
mutationRatesr1 = MutationRate('mutationrates.txt', trinucleotider1, muts_r1, 'mutationrate_r1_results')
mutationRatesr2 = MutationRate('mutationrates.txt', trinucleotider2, muts_r2, 'mutationrate_r2_results')
mutationRates30 = MutationRate('mutationrates.txt', trinucleotide30, muts_30, 'mutationrate30_results')
mutationRatesUVA = MutationRate('mutationrates.txt', trinucleotideA,muts_A, 'mutationrateUVA_results')  
mutationRatesOg = MutationRate('mutationrates.txt', trinucleotideOg, muts_Og, 'mutationrateOgg1_results')
mutationRates26 = MutationRate('mutationrates.txt', trinucleotide26,muts_26, 'mutationrate26_results') 
mutationRatesUVB = MutationRate('mutationrates.txt', trinucleotideB, muts_B, 'mutationrateUVB_results')
mutationRatesB16 = MutationRate('mutationrates.txt', trinucleotideB16, muts_B16, 'mutationrateUVB16_results')
mutationRatesB26 = MutationRate('mutationrates.txt', trinucleotideB26, muts_B26, 'mutationrateUVB26_results')
mutationRatesB30 = MutationRate('mutationrates.txt', trinucleotideB30, muts_B30, 'mutationrateUVB30_results')


rep_times = []
rep_times.append(early_1)
rep_times.append(early_2)
rep_times.append(middle_1)
rep_times.append(middle_2)
rep_times.append(late_1)
rep_times.append(late_2)


def ChromosomeAverageTimes(chromosome, file):
    value1 = float(0)
    value2 = float(0)
    value3 = float(0)
    value4 = float(0)
    value5 = float(0)
    value6 = float(0)
    value7 = float(0)
    value8 = float(0)
    value9 = float(0)
    value10 = float(0)
    value11 = float(0)
    value12 = float(0)
    value13 = float(0)
    value14 = float(0)
    value15 = float(0)
    value16 = float(0)
    valuea = 0
    valueb = 0
    valuec = 0
    valued = 0
    valuee = 0
    valuef = 0
    valueg = 0
    valueh = 0
    valuei = 0
    valuej = 0
    valuek = 0
    valuel = 0

    valuem = 0
    valuen = 0
    valueo = 0
    valuep = 0
    for t in range(0, len(scores)):
        if chromosome[t] == 'chrI':
            value1 = value1 + float(scores[t])
            valuea = valuea + 1
        if chromosome[t] == 'chrII':
            value2 = value2 + float(scores[t])
            valueb = valueb + 1
        if chromosome[t] == 'chrIII':
            value3 = value3 + float(scores[t])
            valuec = valuec + 1
        if chromosome[t] == 'chrIV':
            value4 = value4 + float(scores[t])
            valued = valued + 1
        if chromosome[t] == 'chrV':
            value5 = value5 + float(scores[t])
            valuee = valuee + 1
        if chromosome[t] == 'chrVI' :
            value6 = value6 + float(scores[t])
            valuef = valuef + 1
        if chromosome[t] == 'chrVII':
            value7 = value7 + float(scores[t])
            valueg = valueg + 1
        if chromosome[t] == 'chrVIII':
            value8 = value8 + float(scores[t])
            valueh = valueh + 1
        if chromosome[t] == 'chrIX':
            value9 = value9 + float(scores[t])
            valuei = valuei + 1
        if chromosome[t] == 'chrX':
            value10 = value10 + float(scores[t])
            valuej = valuej + 1
        if chromosome[t] == 'chrXI':
            value11 = value11 + float(scores[t])
            valuek = valuek + 1
        if chromosome[t] == 'chrXII':
            value12 = value12 + float(scores[t])
            valuel = valuel + 1
        if chromosome[t] == 'chrXIII':
            value13 = value13 + float(scores[t])
            valuem = valuem + 1
        if chromosome[t] == 'chrXIV':
            value14 = value14 + float(scores[t])
            valuen = valuen + 1
        if chromosome[t] == 'chrXV':
            value15 = value15 + float(scores[t])
            valueo = valueo + 1
        if chromosome[t] == 'chrXVI':
            value16 = value16 + float(scores[t])
            valuep = valuep + 1
 
    ft = open(file, 'w+')
    ft.write('Chromosome 1: ' + str(value1/valuea))
    ft.write('\n')
    ft.write('Chromosome 2: ' + str(value2/valueb))
    ft.write('\n')
    ft.write('Chromosome 3: ' + str(value3/valuec))
    ft.write('\n')
    ft.write('Chromosome 4: ' + str(value4/valued))
    ft.write('\n')
    ft.write('Chromosome 5: ' + str(value5/valuee))
    ft.write('\n')
    ft.write('Chromosome 6: ' + str(value6/valuef))
    ft.write('\n')
    ft.write('Chromosome 7: ' + str(value7/valueg))
    ft.write('\n')
    ft.write('Chromosome 8: ' + str(value8/valueh))
    ft.write('\n')
    ft.write('Chromosome 9: ' + str(value9/valuei))
    ft.write('\n')
    ft.write('Chromosome 10: ' + str(value10/valuej))
    ft.write('\n')
    ft.write('Chromosome 11: ' + str(value11/valuek))
    ft.write('\n')
    ft.write('Chromosome 12: ' + str(value12/valuel))
    ft.write('\n')
    ft.write('Chromosome 13: ' + str(value13/valuem))
    ft.write('\n')
    ft.write('Chromosome 14: ' + str(value14/valuen))
    ft.write('\n')
    ft.write('Chromosome 15: ' + str(value15/valueo))
    ft.write('\n')
    ft.write('Chromosome 16: ' + str(value16/valuep))
    ft.write('\n')

ChromosomeAverageTimes(chr, 'chromosome_average_times.txt')
        

early_TrinucCounts = {}
middle_TrinucCounts = {}
late_TrinucCounts = {}

dicts = [early_TrinucCounts, middle_TrinucCounts, late_TrinucCounts]

#Counts numbers of trinucleotides in early, middle, and late replicating regions
#based on input from yeast genome assembly R64
def TrinucleotideContext(file, name, start, scores, early, middle, dicts, num):
    early_TrinucCounts = dicts[0]
    middle_TrinucCounts = dicts[1]
    late_TrinucCounts = dicts[2]
    
    
    chr_scores = []
    ft = open(file)
    sequence = ft.read()
    sequence = str(sequence.strip().upper())
    sequence = ''.join(sequence.split('\n'))
    for a in range(0, len(scores)):
        if chr[a] == name:
            chr_scores.append(scores[a])
    if chr_scores[0] <= early:
        for x in range(0, 250):
            if num == 3:
                if sequence[x:x+num][1] == 'A' or sequence[x:x+num][1] == 'G':
                    seq1 = GetReverseComplement(sequence[x:x+num])
                else:
                    seq1 = sequence[x:x+num]
            else:
                seq1 = sequence[x:x+num]
                
            if seq1 not in early_TrinucCounts.keys():
                early_TrinucCounts[seq1] = 1
        
            else:
                early_TrinucCounts[seq1] = early_TrinucCounts[seq1] + 1
    elif chr_scores[0] <= middle:
        for y in range(0, 250):
            if num == 3:
                if sequence[y:y+num][1] == 'A' or sequence[y:y+num][1] == 'G':
                    seq2 = GetReverseComplement(sequence[y:y+num])
                else:
                    seq2 = sequence[y:y+num]
            else:
                seq2 = sequence[y:y+num]
            if seq2 not in middle_TrinucCounts.keys():
                middle_TrinucCounts[seq2] = 1
            else:
                middle_TrinucCounts[seq2] = middle_TrinucCounts[seq2] + 1
    else:
        for z in range(0, 250):
            print(z)
            if num == 3:
                if sequence[z:z+num][1] == 'A' or sequence[z:z+num][1] == 'G':
                    seq3 = GetReverseComplement(sequence[z:z+num])
                else:
                    seq3 = sequence[z:z+num]
            else:
                seq3 = sequence[z:z+num]

            if seq3 not in late_TrinucCounts.keys():
                late_TrinucCounts[seq3] = 1
            else:
                late_TrinucCounts[seq3] = late_TrinucCounts[seq3] + 1
            
    for i in range(1, len(scores)):
        if chr[i] == name:
            if scores[i] <= early:
                for j in range(int(start[i - 1]), int(start[i ])):
                    print(j)
                    if num == 3:
                        if sequence[j:j+num][1] == 'A' or sequence[j:j+num][1] == 'G':
                            seq4 = GetReverseComplement(sequence[j:j+num])
                        else:
                            seq4 = sequence[j:j+num]
                    else:
                        seq4 = sequence[j:j+num]
                    if seq4 not in early_TrinucCounts.keys():
                        early_TrinucCounts[seq4] = 1
                    else:
                        early_TrinucCounts[seq4] = early_TrinucCounts[seq4] + 1

    
            elif scores[i] > early and scores[i] <= middle:
                for L in range(int(start[i - 1]), int(start[i ])):
                    print(L)
                    if num == 3:
                        if sequence[L:L+num][1] == 'A' or sequence[L:L+num][1] == 'G':
                            seq5 = GetReverseComplement(sequence[L:L+num])
                        else:
                            seq5 = sequence[L:L+num]
                    else:
                        seq5 = sequence[L:L+num]
                    if seq5 not in middle_TrinucCounts.keys():
                        middle_TrinucCounts[seq5] = 1
                    else:
                        middle_TrinucCounts[seq5] = middle_TrinucCounts[seq5] + 1
    
            elif scores[i] > middle:
                for n in range(int(start[i -1]), int(start[i])):
                    print(n)
                    if num == 3:
                        if sequence[n:n+num][1] == 'A' or sequence[n:n+num][1] == 'G':
                            seq6 = GetReverseComplement(sequence[n:n+num])
                        else:
                            seq6 = sequence[n:n+num]
                    else:
                        seq6 = sequence[n:n+num]
                    if seq6 not in late_TrinucCounts.keys():
                        late_TrinucCounts[seq6] = 1
                    else:
                        late_TrinucCounts[seq6] = late_TrinucCounts[seq6] + 1
    ft.close()
    return [early_TrinucCounts, middle_TrinucCounts, late_TrinucCounts]
dict1 = TrinucleotideContext('chr1.fa', 'chrI', start, scores, early_2, middle_2, dicts, 3)
dict2 = TrinucleotideContext('chr2.txt', 'chrII', start, scores, early_2, middle_2, dict1, 3)
dict3 = TrinucleotideContext('chr3.txt', 'chrIII', start, scores, early_2, middle_2, dict2, 3)
dict4 = TrinucleotideContext('chr4.txt', 'chrIV', start, scores, early_2, middle_2, dict3, 3)
dict5 = TrinucleotideContext('chr5.txt', 'chrV', start, scores, early_2, middle_2, dict4, 3)
dict6 = TrinucleotideContext('chr6.txt', 'chrVI', start, scores, early_2, middle_2, dict5, 3)
dict7 = TrinucleotideContext('chr7.txt', 'chrVII', start, scores, early_2, middle_2, dict6, 3)
dict8 = TrinucleotideContext('chr8.txt', 'chrVIII', start, scores, early_2, middle_2, dict7, 3)
dict9 = TrinucleotideContext('chr9.txt', 'chrIX', start, scores, early_2, middle_2, dict8, 3)
dict10 = TrinucleotideContext('chr10.txt', 'chrX', start, scores, early_2, middle_2, dict9, 3)
dict11 = TrinucleotideContext('chr11.txt', 'chrXI', start, scores, early_2, middle_2, dict10, 3)
dict12 = TrinucleotideContext('chr12.txt', 'chrXII', start, scores, early_2, middle_2, dict11, 3)
dict13 = TrinucleotideContext('chr13.txt', 'chrXIII', start, scores, early_2, middle_2, dict12, 3)
dict14 = TrinucleotideContext('chr14.txt', 'chrXIV', start, scores, early_2, middle_2, dict13, 3)
dict15 = TrinucleotideContext('chr15.txt', 'chrXV', start, scores, early_2, middle_2, dict14, 3)
counts = TrinucleotideContext('chr16.txt', 'chrXVI', start, scores, early_2, middle_2, dict15, 3)

early_Dinuc = {}
middle_Dinuc = {}
late_Dinuc = {}

dicts2 = [early_Dinuc, middle_Dinuc, late_Dinuc]
dict1b = TrinucleotideContext('chr1.txt', 'chrI', start, scores, early_2, middle_2, dicts2, 2)
dict2b = TrinucleotideContext('chr2.txt', 'chrII', start, scores, early_2, middle_2, dict1b, 2)
dict3b = TrinucleotideContext('chr3.txt', 'chrIII', start, scores, early_2, middle_2, dict2b, 2)
dict4b = TrinucleotideContext('chr4.txt', 'chrIV', start, scores, early_2, middle_2, dict3b, 2)
dict5b = TrinucleotideContext('chr5.txt', 'chrV', start, scores, early_2, middle_2, dict4b, 2)
dict6b = TrinucleotideContext('chr6.txt', 'chrVI', start, scores, early_2, middle_2, dict5b, 2)
dict7b = TrinucleotideContext('chr7.txt', 'chrVII', start, scores, early_2, middle_2, dict6b, 2)
dict8b = TrinucleotideContext('chr8.txt', 'chrVIII', start, scores, early_2, middle_2, dict7b, 2)
dict9b = TrinucleotideContext('chr9.txt', 'chrIX', start, scores, early_2, middle_2, dict8b, 2)
dict10b = TrinucleotideContext('chr10.txt', 'chrX', start, scores, early_2, middle_2, dict9b, 2)
dict11b = TrinucleotideContext('chr11.txt', 'chrXI', start, scores, early_2, middle_2, dict10b, 2)
dict12b = TrinucleotideContext('chr12.txt', 'chrXII', start, scores, early_2, middle_2, dict11b, 2)
dict13b = TrinucleotideContext('chr13.txt', 'chrXIII', start, scores, early_2, middle_2, dict12b, 2)
dict14b = TrinucleotideContext('chr14.txt', 'chrXIV', start, scores, early_2, middle_2, dict13b, 2)
dict15b = TrinucleotideContext('chr15.txt', 'chrXV', start, scores, early_2, middle_2, dict14b, 2)
counts2 = TrinucleotideContext('chr16.txt', 'chrXVI', start, scores, early_2, middle_2, dict15b, 2)



early_Trinuc = {}
middle_Trinuc = {}
late_Trinuc = {}

chr_dicts = [early_Trinuc, middle_Trinuc, late_Trinuc]

chr1_counts =   TrinucleotideContext('chr1.txt', 'chrI', start, scores, early_2, middle_2, chr_dicts, 3)  
                                                   
def FindExpected(mutationRates, count, file):
    early_totals = count[0] 
    middle_totals = count[1] 
    late_totals = count[2] 
    expected_early = 0
    expected_middle = 0
    expected_late = 0
    f = open(file, 'w+')
    for key in mutationRates.keys():
        
    #mutated over total for each codon
        if key in late_totals.keys() and key in early_totals.keys():
            mutationRates[key] = mutationRates[key]/(count[0][key] + count[1][key] + count[2][key])
        elif key not in late_totals.keys():
            mutationRates[key] = mutationRates[key]/(count[0][key] + count[1][key]) 
        else:
            mutationRates[key] = mutationRates[key]/(count[2][key] + count[1][key])
        if key in early_totals.keys():    
            expected_early = expected_early + (early_totals[key] * mutationRates[key])
        if key in middle_totals.keys():
            expected_middle = expected_middle + (middle_totals[key] * mutationRates[key])
        if key in late_totals.keys():
            expected_late = expected_late + (late_totals[key] * mutationRates[key])
    total = expected_early + expected_middle + expected_late
    f.write('Early')
    f.write('\n')
    for key in early_totals.keys():
       
        f.write(key + ': ' + str((early_totals[key] )))
        f.write('\n')
       
    f.write('Middle')
    f.write('\n')
    for key in middle_totals.keys():
        f.write(key + ': ' + str((middle_totals[key])))
        f.write('\n')
       
    f.write('Late')
    f.write('\n')
    for key in late_totals.keys():
        
        f.write(key + ': ' + str((late_totals[key])))
        f.write('\n')
        
    print(total) 
    f.close()     

    percentage_early = expected_early/total
    percentage_middle = expected_middle/total
    percentage_late = expected_late/total
    print(percentage_early)
    print(percentage_middle)
    print(percentage_late)
    print(expected_early)
    print(expected_middle)
    print(expected_late)
    return[percentage_early, percentage_middle, percentage_late]

FindExpected(mutationRates, counts, 'Rates.txt')
FindExpected(mutationRates16, counts, 'Rates16.txt')
FindExpected(mutationRatesr1, counts, 'Rates_r1.txt')
FindExpected(mutationRatesr2, counts, 'Rates_r2.txt')
FindExpected(mutationRates30, counts, 'Rates30.txt')
FindExpected(mutationRates26, counts, 'Rates26.txt')
FindExpected(mutationRatesOg, counts, 'RatesOg.txt')
FindExpected(mutationRatesUVA, counts, 'RatesA.txt')
FindExpected(mutationRatesUVB, counts, 'Rates_UVB.txt')
FindExpected(mutationRatesB16, counts, 'Rates_UVB16.txt')
FindExpected(mutationRatesB26, counts, 'Rates_UVB26.txt')
FindExpected(mutationRatesB30,  counts, 'Rates_UVB30.txt')

expected_percentages = FindExpected(mutationRates, chr1_counts, 'Rates_chr1.txt')
expected_percentages16 = FindExpected(mutationRates16, chr1_counts, 'Rates16_chr1.txt')
expected_percentagesUVB = FindExpected(mutationRatesUVB, chr1_counts, 'Rates_UVB_chr1.txt')
expected_percentages26 = FindExpected(mutationRates26, chr1_counts, 'Rates26_chr1.txt')
expected_percentagesB26 = FindExpected(mutationRatesB26, chr1_counts, 'Rates26_UVB_chr1.txt')
expected_percentagesB16 = FindExpected(mutationRatesB16, chr1_counts, 'Rates_UVB16_chr1.txt')

def MutationNumbers(file,file2,chromosome, expected_percentages, chr_name, muts):
    f1 = open(file)
    sequence1 = f1.read()
    sequence1 = str(sequence1.strip().upper())
    sequence1 = ''.join(sequence1.split('\n'))
    count1 = 0
    for c in range(0, len(chromosome)):
        if chromosome[c] == chr_name:
            count1 = count1 + 1
    density1 = count1/len(sequence1)

    expected_early_1 = count1 * expected_percentages[0]
    expected_middle_1 = count1 * expected_percentages[1]
    expected_late_1 = count1 * expected_percentages[2]
    f2 = open(file2, 'w+')
    f2.write(chr_name)
    f2.write('\n')
    f2.write('Early: ' + str(len(muts[0])))
    f2.write('\n')
    f2.write('Middle: ' + str(len(muts[1])))
    f2.write('\n')
    f2.write('Late: ' + str(len(muts[2])))
    f2.write('\n')
    f2.write('Total: ' + str(len(muts[0]) + len(muts[1]) + len(muts[2])))
    f2.write('\n')
    f2.write('Expected Early: ' + str(int(expected_early_1)))
    f2.write('\n')
    f2.write('Expected Middle: ' + str(int(expected_middle_1)))
    f2.write('\n')
    f2.write('Expected Late Percentage: ' + str(expected_percentages[2]))
    f2.write('\n')
    f2.write('Expected Late: ' + str(int(expected_late_1)))
    f2.write('\n')
    f2.write('Mutation Density: ' + str(density1))
    f1.close()
    f2.close()
    
    
MutationNumbers('chr1.txt', 'chr1_value.txt', chromosome, expected_percentages,'chrI', muts)
MutationNumbers('chr1.txt', 'chr1_value_rad16.txt', chromosome16, expected_percentages16,'chrI', muts_16)
MutationNumbers('chr1.txt', 'chr1_value_UVB', chromosomeB, expected_percentagesUVB,'chrI', muts_B)
MutationNumbers('chr1.txt', 'chr1_value_UVBrad16', chromosomeB16, expected_percentagesB16,'chrI', muts_B16)
MutationNumbers('chr1.txt', 'chr1_value_rad26', chromosome26, expected_percentages26,'chrI', muts_26)
MutationNumbers('chr1.txt', 'chr1_value_UVBrad26', chromosomeB26, expected_percentagesB26,'chrI', muts_B26)