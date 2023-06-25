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
with open("WT_UV_run1_Muts_allSNVs_sorted.bed") as fr1:
    datar1 = fr1.read()
with open("WT_UV_run2_Muts_allSNVs_sorted.bed") as fr2:
    datar2 = fr2.read()
with open("WT_UV_run1_Muts_allDNVs_sorted.bed") as fd1:
    datad1 = fd1.read()
with open("WT_UV_run2_Muts_allDNVs_sorted.bed") as fd2:
    datad2 = fd2.read()
with open("UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem_insertions.bed") as fin:
    datain = fin.read()
with open('WT_DIPs_sorted_anz2_tandem_insertions.bed') as fina:
    dataina = fina.read()
with open('rad16_DIPs_sorted_anz2_tandem_insertions.bed') as finb:
    datainb = finb.read()
with open("UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem_deletions.bed") as fdel:
    datadel = fdel.read()
with open('WT_DIPs_sorted_anz2_tandem_deletions.bed') as fdela:
    datadela = fdela.read()
with open('rad16_DIPs_sorted_anz2_tandem_deletions.bed') as fdelb:
    datadelb = fdelb.read()

#extract columns from bed files to make arrays
# with matching indices   
def getChromosome(data):       
    lines = [s.strip().split() for s in data.splitlines()]
    chromosome = []
    for i in range (0, len(lines)):
        chromosome.append(lines[i][0])
    return chromosome
chromosome = getChromosome(data)
chromosome16 = getChromosome(data16)
chromosome30 = getChromosome(data30)
chromosomeA = getChromosome(dataUVA)
chromosomeOg = getChromosome(dataOg)
chromosome26 = getChromosome(data26)
chromosomeB = getChromosome(dataUVB)
chromosomeB16 = getChromosome(dataB16)
chromosomer1 = getChromosome(datar1)
chromosomer2 = getChromosome(datar2)
chromosomed1 = getChromosome(datad1)
chromosomed2 = getChromosome(datad2)
chromosomein = getChromosome(datain)
chromosomeina = getChromosome(dataina)
chromosomeinb = getChromosome(datainb)
chromosomedel = getChromosome(datadel)
chromosomedela = getChromosome(datadela)
chromosomedelb = getChromosome(datadelb)
    
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
position1_r1 = getPosition1(datar1)
position1_r2 = getPosition1(datar2)
position1_d1 = getPosition1(datad1)
position1_d2 = getPosition1(datad2)
position1_in = getPosition1(datain)
position1_ina = getChromosome(dataina)
position1_inb = getChromosome(datainb)
position1_del = getPosition1(datadel)
position1_dela = getPosition1(datadela)
position1_delb = getPosition1(datadelb)

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
position2_r1 = getPosition2(datar1)
position2_r2 = getPosition2(datar2)
position2_d1 = getPosition2(datad1)
position2_d2 = getPosition2(datad2)
position2_in = getPosition2(datain)
position2_ina = getPosition2(dataina)
position2_inb = getPosition2(datainb)
position2_del = getPosition2(datadel)
position2_dela = getPosition2(datadela)
position2_delb = getPosition2(datadelb)

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
trinucleotider1 = getTrinucleotide(datar1)
trinucleotider2 = getTrinucleotide(datar2)
trinucleotidein = getTrinucleotide(datain)
trinucleotideina = getTrinucleotide(dataina)
trinucleotideinb = getTrinucleotide(datainb)
trinucleotidedel = getTrinucleotide(datadel)
trinucleotidedela = getTrinucleotide(datadela)
trinucleotidedelb = getTrinucleotide(datadelb)

tridel = []
tridela = []
tridelb = []

for tri in range(0,len(trinucleotidein)):
    trinucleotidein[tri]= str(trinucleotidein[tri][0]) + str(trinucleotidein[tri][len(trinucleotidein[tri]) - 1])
for tria in range(0,len(trinucleotideina)):
    trinucleotideina[tria] = str(trinucleotideina[tria][0]) + str(trinucleotideina[tria][len(trinucleotideina[tria]) - 1])
for trib in range(0,len(trinucleotidedelb)):
    trinucleotideinb[trib] = str(trinucleotideinb[trib][0]) + str(trinucleotideinb[trib][len(trinucleotideinb[trib])- 1])

def getMononucleotide(data):
    lines = [s.strip().split() for s in data.splitlines()]  
    mononucleotide = []
    for n in range (0, len(lines)):
        mononucleotide.append(lines[n][4])
    return mononucleotide
mononucleotide = getMononucleotide(data)
mononucleotide16 = getMononucleotide(data16)

#sign = []
#for p in range (0, len(lines) ):
    #sign.append(lines[p][5])

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
print(early1_1)
early2_1 = early2[int((2* len(scores))/ 9)]
print(early2_1)
early_2 = early[len(early) - 1]
print(early_2)
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
print(late_1)
late2_1 = late2[int((7 * len(scores))/9)+ 1]
print(late2_1)
late3_1 = late3[int((8 * len(scores))/9) + 1]
print(late3_1)
late_2 = late[int(len(scores)- 1)]
print(late_2)
for q in range(0, len(scores) ):
    if q > len(scores)/3 and q < 2 * len(scores)/3:
        middle[q] = sorted_scores[q]
        if q > len(scores)/3 and q <= (4 * len(scores))/9:
            middle1[q] = sorted_scores[q]
        elif q > len(scores)/3 and q <= (5 * len(scores))/9:
            middle2[q] = sorted_scores[q]
        elif q > len(scores)/3 and q < (2 * len(scores))/3:
            middle3[q] = sorted_scores[q]
middle_1 = early_2 #middle[len(scores)/4 + 1]
middle1_1 = middle1[int((4 * len(scores))/9)]
print(middle1_1)
middle2_1 = middle2[int((5 * len(scores))/9)]
print(middle2_1)
middle_2 = late_1


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
percentages = ChromosomePercentages('chrXVI')




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
        #if position2[po] >= start[st] and position2[po] < start2[st] and chromosome[po] == chr[st]:
            #time = scores[st] #calculateTime(time1, time2, pos, pos1, pos2)
            #times.append(time)
            #po = po + 1
        #elif position2[po] > start2[st] and position2[po] <= start[st + 1] and chromosome[po] == chr[st]:
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
timesr1 = AssignTimes(chromosomer1, position2_r1, start, scores)
timesr2 = AssignTimes(chromosomer2, position2_r2, start, scores)
timesd1 = AssignTimes(chromosomed1, position2_d1, start, scores)
timesd2 = AssignTimes(chromosomed2, position2_d2, start, scores)
timesin = AssignTimes(chromosomein, position2_in, start, scores)
timesina = AssignTimes(chromosomeina, position2_ina, start, scores)
timesinb = AssignTimes(chromosomeinb, position2_inb, start, scores)
timesdel = AssignTimes(chromosomedel, position2_del, start, scores)
timesdela = AssignTimes(chromosomedela, position2_dela, start, scores)
timesdelb = AssignTimes(chromosomedelb, position2_delb, start, scores)

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
    percent_early = (len(early_mut)/ total) * 100
    percent_early_chr = (len(chr_early)/ total_chr) * 100
    percent_middle = (len(middle_mut)/ total) * 100
    percent_middle_chr = (len(chr_middle)/ total_chr) * 100
    percent_late = (len(late_mut)/ total) * 100
    percent_late_chr = (len(chr_late)/ total_chr) * 100
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

muts = CountMutations("WT_mutations", times, chromosome, 'chrVI')
muts_16= CountMutations("rad16_mutations", times16, chromosome16, 'chrVI')
muts_30 = CountMutations("rad30_mutations", times30, chromosome30, 'chrVI')
muts_A = CountMutations("UVA_mutations", times_A, chromosomeA, 'chrVI')
muts_Og = CountMutations("ogg1_mutations", times_Og, chromosomeOg, 'chrVI')
muts_26= CountMutations("rad26_mutations", times26, chromosome26, 'chrVI')
muts_B = CountMutations("UVB_mutations", times_B, chromosomeB, 'chrVI')
muts_B16 = CountMutations("UVB_rad16_mutations", timesB16, chromosomeB16, 'chrVI')
muts_r1 = CountMutations("WT_r1_mutations", timesr1, chromosomer1, 'chrVI')
muts_r2 = CountMutations("WT_r2_mutations", timesr2, chromosomer2, 'chrVI')
muts_d1 = CountMutations("WT_d1_mutations", timesd1, chromosomed1, 'chrVI')
early_mut_d1 = muts_d1[0]
middle_mut_d1 = muts_d1[1]
late_mut_d1 = muts_d1[2]
muts_d2 = CountMutations("WT_d2_mutations", timesd2, chromosomed2, 'chrVI')
early_mut_d2 = muts_d2[0]
middle_mut_d2 = muts_d2[1]
late_mut_d2 = muts_d2[2]
muts_in = CountMutations("insertion_mutations", timesin, chromosomein, 'chrVI')
early_mut_in = muts_in[0]
middle_mut_in = muts_in[1]
late_mut_in = muts_in[2]
muts_ina = CountMutations("WT_insertion_mutations", timesina, chromosomeina, 'chrVI')
early_mut_ina = muts_ina[0]
middle_mut_ina = muts_ina[1]
late_mut_ina = muts_ina[2]
muts_inb = CountMutations("rad16_insertion_mutations", timesinb, chromosomeinb, 'chrVI')
early_mut_inb = muts_inb[0]
middle_mut_inb = muts_inb[1]
late_mut_inb = muts_inb[2]
muts_del = CountMutations("deletion_mutations", timesdel, chromosomedel, 'chrVI')
early_mut_del = muts_del[0]
middle_mut_del = muts_del[1]
late_mut_del = muts_del[2]
muts_dela = CountMutations("WT_deletion_mutations", timesdela, chromosomedela, 'chrVI')
early_mut_dela = muts_dela[0]
middle_mut_dela = muts_dela[1]
late_mut_dela = muts_dela[2]
muts_delb = CountMutations("rad16_deletion_mutations", timesdelb, chromosomedelb, 'chrVI')
early_mut_delb = muts_delb[0]
middle_mut_delb = muts_delb[1]
late_mut_delb = muts_delb[2]

def CountMutations_ch(times, chromosome, name, rep_times):
    early_2 = rep_times[1]
    middle_1 = rep_times[2]
    middle_2 = rep_times[3]
    late_1 = rep_times[4]
    late_2 = rep_times[5]
    chr_early = {}
    chr_middle = {}
    chr_late = {}

    for tm in range(0, len(times) ):
        if times[tm] != 'NA' and times[tm ] >= 0 and times[tm] < early_2:
            if chromosome[tm] == name:
                chr_early[tm] = times[tm]
        if times[tm] != 'NA' and times[tm ] >= middle_1 and times[tm] < middle_2:
            if chromosome[tm] == name:
                chr_middle[tm] = times[tm]
        if times[tm] != 'NA' and times[tm ] >= late_1:
            if chromosome[tm] == name:
                chr_late[tm] = times[tm]
    return ([chr_early, chr_middle, chr_late])
#muts2 = CountMutations_ch(times, chromosome, 'chrXVI', rep_times )
#muts2_16 = CountMutations_ch(times16, chromosome16, 'chrXVI', rep_times)
#muts2_30 = CountMutations_ch(times30, chromosome30, 'chrXVI', bins)
#muts2_A = CountMutations_ch(times_A, chromosomeA, 'chrXVI', bins)
#muts2_og = CountMutations_ch(times_Og, chromosomeOg, 'chrXVI', bins)
#muts2_26 = CountMutations_ch(times26, chromosome26, 'chrXVI', bins)
#muts2_B = CountMutations_ch(times_B, chromosomeB, 'chrXVI', rep_times)
#muts2_B16 = CountMutations_ch(timesB16, chromosomeB16, 'chrXVI', rep_times)


def ClassifyMutations(muts, trinucleotide, mononucleotide, file):
    early_muts = muts[0]
    middle_muts = muts[1]
    late_muts = muts[2]
    
    earlycount1 = 0
    middlecount1 = 0
    latecount1 = 0
    earlycount2 = 0
    middlecount2 = 0
    latecount2 = 0
    earlycount3 = 0
    middlecount3 = 0
    latecount3 = 0
    earlycount4 = 0
    middlecount4 = 0
    latecount4 = 0
    earlycount5 = 0
    middlecount5 = 0
    latecount5 = 0
    earlycount6 = 0
    middlecount6 = 0
    latecount6 = 0
    earlycount7 = 0
    middlecount7 = 0
    latecount7 = 0
    earlycount8 = 0
    middlecount8 = 0
    latecount8 = 0
    earlycount9 = 0
    middlecount9 = 0
    latecount9 = 0
    earlycount10 = 0 
    middlecount10 = 0
    latecount10 = 0
    earlycount11 = 0
    middlecount11 = 0
    latecount11 = 0
    earlycount12 = 0
    middlecount12 = 0
    latecount12 = 0
    for t in range(0, len(trinucleotide) ):
        if trinucleotide[t][1] == 'A'and mononucleotide[t] == 'C':
                if t in early_muts.keys():
                    earlycount1 = earlycount1 + 1
                if t in middle_muts.keys():
                    middlecount1 = middlecount1 + 1
                if t in late_muts.keys():
                    latecount1 = latecount1 + 1
        if trinucleotide[t][1] == 'A'and mononucleotide[t] == 'G':
            if t in early_muts.keys():
                    earlycount2 = earlycount2 + 1
            if t in middle_muts.keys():
                middlecount2 = middlecount2 + 1
            if t in late_muts.keys():
                latecount2 = latecount2 + 1
        if trinucleotide[t][1] == 'A'and mononucleotide[t] == 'T':
            if t in early_muts.keys():
                    earlycount3 = earlycount3 + 1
            if t in middle_muts.keys():
                middlecount3 = middlecount3 + 1
            if t in late_muts.keys():
                latecount3 = latecount3 + 1
        if trinucleotide[t][1] == 'C'and mononucleotide[t] == 'A':
            if t in early_muts.keys():
                    earlycount4 = earlycount4 + 1
            if t in middle_muts.keys():
                middlecount4 = middlecount4 + 1
            if t in late_muts.keys():
                latecount4 = latecount4 + 1
        if trinucleotide[t][1] == 'C'and mononucleotide[t] == 'G':
            if t in early_muts.keys():
                    earlycount5 = earlycount5 + 1
            if t in middle_muts.keys():
                middlecount5 = middlecount5 + 1
            if t in late_muts.keys():
                latecount5 = latecount5 + 1
        if trinucleotide[t][1] == 'C'and mononucleotide[t] == 'T':
            if t in early_muts.keys():
                earlycount6 = earlycount6 + 1
            if t in middle_muts.keys():
                middlecount6 = middlecount6 + 1
            if t in late_muts.keys():
                latecount6 = latecount6 + 1
        if trinucleotide[t][1] == 'G'and mononucleotide[t] == 'A':
            if t in early_muts.keys():
                earlycount7 = earlycount7 + 1
            if t in middle_muts.keys():
                middlecount7 = middlecount7 + 1
            if t in late_muts.keys():
                latecount7 = latecount7 + 1
        if trinucleotide[t][1] == 'G'and mononucleotide[t] == 'C':
            if t in early_muts.keys():
                earlycount8 = earlycount8 + 1
            if t in middle_muts.keys():
                middlecount8 = middlecount8 + 1
            if t in late_muts.keys():
                latecount8 = latecount8 + 1
        if trinucleotide[t][1] == 'G'and mononucleotide[t] == 'T':
            if t in early_muts.keys():
                earlycount9 = earlycount9 + 1
            if t in middle_muts.keys():
                middlecount9 = middlecount9 + 1
            if t in late_muts.keys():
                latecount9 = latecount9 + 1
        if trinucleotide[t][1] == 'T'and mononucleotide[t] == 'A':
            if t in early_muts.keys():
                earlycount10 = earlycount10 + 1
            if t in middle_muts.keys():
                middlecount10 = middlecount10 + 1
            if t in late_muts.keys():
                latecount10 = latecount10 + 1
        if trinucleotide[t][1] == 'T'and mononucleotide[t] == 'C':
            if t in early_muts.keys():
                earlycount11 = earlycount11 + 1
            if t in middle_muts.keys():
                middlecount11 = middlecount11 + 1
            if t in late_muts.keys():
                latecount11 = latecount11 + 1
        if trinucleotide[t][1] == 'T'and mononucleotide[t] == 'G':
            if t in early_muts.keys():
                earlycount12 = earlycount12 + 1
            if t in middle_muts.keys():
                middlecount12 = middlecount12 + 1
            if t in late_muts.keys():
                latecount12 = latecount12 + 1
    f = open(file, 'w+')
    f.write('A>C')
    f.write('\n')
    f.write('Early: '+ str(earlycount1))
    if earlycount1 + middlecount1 + latecount1 > 0:   
            f.write(' Percentage: ' + str(earlycount1/(earlycount1 + middlecount1 + latecount1)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount1))
    if earlycount1 + middlecount1 + latecount1 > 0:
        f.write(' Percentage: ' + str(middlecount1/(earlycount1 + middlecount1 + latecount1)))
    f.write('\n')
    f.write('Late: '+ str(latecount1)) 
    if earlycount1 + middlecount1 + latecount1 > 0:
        f.write(' Percentage: ' + str(latecount1/(earlycount1 + middlecount1 + latecount1)))
    f.write('\n')
    f.write('A>G')
    f.write('\n')
    f.write('Early: '+ str(earlycount2))
    if earlycount2 + middlecount2 + latecount2 >0:
        f.write(' Percentage: ' + str(earlycount2/(earlycount2 + middlecount2 + latecount2)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount2)) 
    if earlycount2 + middlecount2 + latecount2 >0:
        f.write(' Percentage: ' + str(middlecount2/(earlycount2 + middlecount2 + latecount2)))
    f.write('\n')
    f.write('Late: '+ str(latecount2)) 
    if earlycount2 + middlecount2 + latecount2 >0:
        f.write(' Percentage: ' + str(latecount2/(earlycount2 + middlecount2 + latecount2)))
    f.write('\n')
    f.write('A>T')
    f.write('\n')
    f.write('Early: '+ str(earlycount3))
    if earlycount3 + middlecount3 + latecount3 > 0: 
        f.write(' Percentage: ' + str(earlycount3/(earlycount3 + middlecount3 + latecount3)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount3)) 
    if earlycount3 + middlecount3 + latecount3 > 0:
        f.write(' Percentage: ' + str(middlecount3/(earlycount3 + middlecount3 + latecount3)))
    f.write('\n')
    f.write('Late: '+ str(latecount3)) 
    if earlycount3 + middlecount3 + latecount3 > 0:
        f.write(' Percentage: ' + str(latecount3/(earlycount3 + middlecount3 + latecount3)))
    f.write('\n')
    f.write('C>A')
    f.write('\n')
    f.write('Early: '+ str(earlycount4)) 
    if earlycount4 + middlecount4 + latecount4 > 0:
        f.write(' Percentage: ' + str(earlycount4/(earlycount4 + middlecount4 + latecount4)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount4)) 
    if earlycount4 + middlecount4 + latecount4 > 0:
        f.write(' Percentage: ' + str(middlecount4/(earlycount4 + middlecount4 + latecount4)))
    f.write('\n')
    f.write('Late: '+ str(latecount4)) 
    if earlycount4 + middlecount4 + latecount4 > 0:
        f.write(' Percentage: ' + str(latecount4/(earlycount4 + middlecount4 + latecount4)))
    f.write('\n')
    f.write('C>G')
    f.write('\n')
    f.write('Early: '+ str(earlycount5)) 
    if earlycount5 + middlecount5 + latecount5 > 0:
        f.write(' Percentage: ' + str(earlycount5/(earlycount5 + middlecount5 + latecount5)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount5)) 
    if earlycount5 + middlecount5 + latecount5 > 0:
        f.write(' Percentage: ' + str(middlecount5/(earlycount5 + middlecount5 + latecount5)))
    f.write('\n')
    f.write('Late: '+ str(latecount5)) 
    if earlycount5 + middlecount5 + latecount5 > 0:
        f.write(' Percentage: ' + str(latecount5/(earlycount5 + middlecount5 + latecount5)))
    f.write('\n')
    f.write('C>T')
    f.write('\n')
    f.write('Early: '+ str(earlycount6)) 
    if earlycount6 + middlecount6 + latecount6 > 0:
        f.write(' Percentage: ' + str(earlycount6/(earlycount6 + middlecount6 + latecount6)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount6)) 
    if earlycount6 + middlecount6 + latecount6 > 0:
        f.write(' Percentage: ' + str(middlecount6/(earlycount6 + middlecount6 + latecount6)))
    f.write('\n')
    f.write('Late: '+ str(latecount6)) 
    if earlycount6 + middlecount6 + latecount6 > 0:
        f.write(' Percentage: ' + str(latecount6/(earlycount6 + middlecount6 + latecount6)))
    f.write('\n')
    f.write('G>A')
    f.write('\n')
    f.write('Early: '+ str(earlycount7)) 
    if earlycount7 + middlecount7 + latecount7 > 0:
        f.write(' Percentage: ' + str(earlycount7/(earlycount7 + middlecount7 + latecount7)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount7)) 
    if earlycount7 + middlecount7 + latecount7 > 0:
        f.write(' Percentage: ' + str(middlecount7/(earlycount7 + middlecount7 + latecount7)))
    f.write('\n')
    f.write('Late: '+ str(latecount7)) 
    if earlycount7 + middlecount7 + latecount7 > 0:
        f.write(' Percentage: ' + str(latecount7/(earlycount7 + middlecount7 + latecount7)))
    f.write('\n')
    f.write('G>C')
    f.write('\n')
    f.write('Early: '+ str(earlycount8)) 
    if earlycount8 + middlecount8 + latecount8 > 0:
        f.write(' Percentage: ' + str(earlycount8/(earlycount8 + middlecount8 + latecount8)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount8)) 
    if earlycount8 + middlecount8 + latecount8 > 0:
        f.write(' Percentage: ' + str(middlecount8/(earlycount8 + middlecount8 + latecount8)))
    f.write('\n')
    f.write('Late: '+ str(latecount8)) 
    if earlycount8 + middlecount8 + latecount8 > 0:
        f.write(' Percentage: ' + str(latecount8/(earlycount8 + middlecount8 + latecount8)))
    f.write('\n')
    f.write('G>T')
    f.write('\n')
    f.write('Early: '+ str(earlycount9)) 
    if earlycount9 + middlecount9 + latecount9 > 0:
        f.write(' Percentage: ' + str(earlycount9/(earlycount9 + middlecount9 + latecount9)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount9)) 
    if earlycount9 + middlecount9 + latecount9 > 0:
        f.write(' Percentage: ' + str(middlecount9/(earlycount9 + middlecount9 + latecount9)))
    f.write('\n')
    f.write('Late: '+ str(latecount9)) 
    if earlycount9 + middlecount9 + latecount9 > 0:
        f.write(' Percentage: ' + str(latecount9/(earlycount9 + middlecount9 + latecount9)))
    f.write('\n')
    f.write('T>A')
    f.write('\n')
    f.write('Early: '+ str(earlycount10)) 
    if earlycount10 + middlecount10 + latecount10 > 0:
        f.write(' Percentage: ' + str(earlycount10/(earlycount10 + middlecount10 + latecount10)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount10)) 
    if earlycount10 + middlecount10 + latecount10 > 0:
        f.write(' Percentage: ' + str(middlecount10/(earlycount10 + middlecount10 + latecount10)))
    f.write('\n')
    f.write('Late: '+ str(latecount10)) 
    if earlycount10 + middlecount10 + latecount10 > 0:
        f.write(' Percentage: ' + str(latecount10/(earlycount10 + middlecount10 + latecount10)))
    f.write('\n')
    f.write('T>C')
    f.write('\n')
    f.write('Early: '+ str(earlycount11)) 
    if earlycount11 + middlecount11 + latecount11 > 0:
        f.write(' Percentage: ' + str(earlycount11/(earlycount11 + middlecount11 + latecount11)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount11)) 
    if earlycount11 + middlecount11 + latecount11 > 0:
        f.write(' Percentage: ' + str(middlecount11/(earlycount11 + middlecount11 + latecount11)))
    f.write('\n')
    f.write('Late: '+ str(latecount11)) 
    if earlycount11 + middlecount11 + latecount11 > 0:
        f.write(' Percentage: ' + str(latecount11/(earlycount11 + middlecount11 + latecount11)))
    f.write('\n')
    f.write('T>G')
    f.write('\n')
    f.write('Early: '+ str(earlycount12)) 
    if earlycount12 + middlecount12 + latecount12 > 0:
        f.write(' Percentage: ' + str(earlycount12/(earlycount12 + middlecount12 + latecount12)))
    f.write('\n')
    f.write('Middle: '+ str(middlecount12)) 
    if earlycount12 + middlecount12 + latecount12 > 0:
        f.write(' Percentage: ' + str(middlecount12/(earlycount12 + middlecount12 + latecount12)))
    f.write('\n')
    f.write('Late: '+ str(latecount12)) 
    if earlycount12 + middlecount12 + latecount12 > 0:
        f.write(' Percentage: ' + str(latecount12/(earlycount12 + middlecount12 + latecount12)))
    f.write('\n')
    


ClassifyMutations(muts, trinucleotide, mononucleotide, 'base_changes')
ClassifyMutations(muts_16, trinucleotide16, mononucleotide16, 'base_changes_rad16')

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
AverageTimes(muts_r1, 'average_times_r1.txt')
AverageTimes(muts_r2, 'average_times_r2.txt')
AverageTimes(muts_d1, 'average_times_d1.txt')
AverageTimes(muts_d2, 'average_times_d2.txt')

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
#AddTimeColumn("WT_UV_bothruns_allSNVs_sorted.bed", "early_mutations.bed", "middle_mutations.bed", "late_mutations.bed",data, times, early_mut, middle_mut, late_mut )
#AddTimeColumn("rad16_UV_bothruns_allSNVs_sorted.bed", "early_mutations_rad16.bed", "middle_mutations_rad16.bed", "late_mutations_rad16.bed",data16, times16, early_mut_16, middle_mut_16, late_mut_16 )
#AddTimeColumn("rad30_UV_bothruns_allSNVs_sorted.bed", "early_mutations._rad30bed", "middle_mutations_rad30.bed", "late_mutations_rad30.bed",data30, times30, early_mut_30, middle_mut_30, late_mut_30 )
#AddTimeColumn("WT_UVA_allSNVs_sorted.bed", "early_mutations_UVA.bed", "middle_mutations_UVA.bed", "late_mutations_UVA.bed",dataUVA, times_A, early_mut_A, middle_mut_A, late_mut_A)
#AddTimeColumn("ogg1_UVA_allSNVs_sorted.bed", "early_mutations_ogg1.bed", "middle_mutations_ogg1.bed", "late_mutations_ogg1.bed",dataOg, times_Og, early_mut_Og, middle_mut_Og, late_mut_Og)
#AddTimeColumn("rad26_UV_bothruns_allSNVs_sorted.bed", "early_mutations_rad26.bed", "middle_mutations_rad26.bed", "late_mutations_rad26.bed",data26, times26, early_mut_26, middle_mut_26, late_mut_26 )
#AddTimeColumn("WT_UVB_allSNVs_sorted.bed", "early_mutations_UVB.bed", "middle_mutations_UVB.bed", "late_mutations_UVB.bed",dataUVB, times_B, early_mut_B, middle_mut_B, late_mut_B)
#AddTimeColumn("Rad16_UVB_allSNVs_sorted.bed", "early_mutations_rad16UVB.bed", "middle_mutations_rad16UVB.bed", "late_mutations_rad16UVB.bed",dataB16, timesB16, early_mut_B16, middle_mut_B16, late_mut_B16)
#AddTimeColumn("WT_UV_run1_allSNVs_sorted.bed", "early_mutations_r1.bed", "middle_mutations_r1.bed", "late_mutations_r1.bed",datar1, timesr1, early_mut_r1, middle_mut_r1, late_mut_r1 )
#AddTimeColumn("WT_UV_run2_allSNVs_sorted.bed", "early_mutations_r2.bed", "middle_mutations_r2.bed", "late_mutations_r2.bed",datar2, timesr2, early_mut_r2, middle_mut_r2, late_mut_r2 )
#AddTimeColumn("WT_UV_run1_allDNVs_sorted.bed", "early_mutations_d1.bed", "middle_mutations_d1.bed", "late_mutations_d1.bed",datad1, timesd1, early_mut_d1, middle_mut_d1, late_mut_d1 )
#AddTimeColumn("WT_UV_run2_allDNVs_sorted.bed", "early_mutations_d2.bed", "middle_mutations_d2.bed", "late_mutations_d2.bed",datad2, timesd2, early_mut_d2, middle_mut_d2, late_mut_d2 )
AddTimeColumn("UV_insertions_sorted.bed", "early_mutations_insertions.bed", "middle_mutations_insertions.bed", "late_mutations_insertions.bed",datain, timesin, early_mut_in, middle_mut_in, late_mut_in )
AddTimeColumn("WT_insertions_sorted.bed", "WT_early_mutations_insertions.bed", "WT_middle_mutations_insertions.bed", "WT_late_mutations_insertions.bed",dataina, timesina, early_mut_ina, middle_mut_ina, late_mut_ina )
AddTimeColumn("rad16_insertions_sorted.bed", "rad16_early_mutations_insertions.bed", "rad16_middle_mutations_insertions.bed", "rad_16_late_mutations_insertions.bed",datainb, timesinb, early_mut_inb, middle_mut_inb, late_mut_inb )
#AddTimeColumn("UV__deletions_sorted.bed", "early_mutations_deletions.bed", "middle_mutations_deletions.bed", "late_mutations_deletions.bed",datadel, timesdel, early_mut_del, middle_mut_del, late_mut_del )
#AddTimeColumn("WT__deletions_sorted.bed", "WT_early_mutations_deletions.bed", "WT_middle_mutations_deletions.bed", "WT_late_mutations_deletions.bed",datadela, timesdela, early_mut_dela, middle_mut_dela, late_mut_dela )
#AddTimeColumn("rad16__deletions_sorted.bed", "rad16_early_mutations_deletions.bed", "rad16_middle_mutations_deletions.bed", "rad16_late_mutations_deletions.bed",datadelb, timesdelb, early_mut_delb, middle_mut_delb, late_mut_delb )

#Calculates number of mutations for each codon
def MutationRate(file1, trinucleotides, muts, file2):
    mutationRates = {}
    base_percentages = {}
    mutationRates_late = {}
    base_percentages_late = {}
    fm = open(file1)
    datam = fm.read()
    fm2 = open(file2, 'w+')
    for line in datam.splitlines():
        line = datam.strip().split()
    for d in range(0, len(line)):
        mutationRates[line[d]] = 0
        base_percentages[line[d]] = 0
        mutationRates_late[line[d]] = 0
        base_percentages_late[line[d]] = 0
    for key in mutationRates.keys():
        for t in range (0, len(trinucleotides)):
            if trinucleotides[t] == key:
                mutationRates[key] = mutationRates[key] + 1
                base_percentages[key] = (mutationRates[key] / len(trinucleotides))
            if trinucleotides[t] == key and t in muts[2].keys():
                mutationRates_late[key] = mutationRates_late[key] + 1
                base_percentages_late[key] = (mutationRates_late[key] / len(trinucleotides))
             
    for key in mutationRates.keys():
        fm2.write(key + ": Number = " + str(mutationRates[key])+ ' Percentage =' + str(base_percentages[key] * 100) + "%")
        fm2.write('\n')
    fm2.write('Late Replication: ')
    fm2.write('\n')

    for key in mutationRates_late.keys():
        fm2.write(key + ": Number = " + str(mutationRates_late[key])+ ' Percentage =' + str(base_percentages_late[key] * 100) + "%")
        fm2.write('\n')
    return(mutationRates)

#mutationRates = MutationRate('mutationrates.txt', trinucleotide, muts, 'mutationrate_results')
#mutationRates16 = MutationRate('mutationrates.txt', trinucleotide16, muts_16, 'mutationrate16_results')
#mutationRatesr1 = MutationRate('mutationrates.txt', trinucleotider1, muts_r1, 'mutationrate_r1_results')
#mutationRatesr2 = MutationRate('mutationrates.txt', trinucleotider2, muts_r2, 'mutationrate_r2_results')
#mutationRates30 = MutationRate('mutationrates.txt', trinucleotide30, 'mutationrate30_results')
#mutationRatesUVA = MutationRate('mutationrates.txt', trinucleotideA, 'mutationrateUVA_results')  
#mutationRatesOg = MutationRate('mutationrates.txt', trinucleotideOg, 'mutationrateOgg1_results')
#mutationRates26 = MutationRate('mutationrates.txt', trinucleotide26,muts_26, 'mutationrate26_results') 
#mutationRatesUVB = MutationRate('mutationrates.txt', trinucleotideB, 'mutationrateUVB_results')
#mutationRatesB16 = MutationRate('mutationrates.txt', trinucleotideB16, 'mutationrateUVB16_results')
mutationRatesin = MutationRate('di_mutationrates.txt', trinucleotidein,muts_in, 'mutationrateinsertion_results') 
mutationRatesina = MutationRate('di_mutationrates.txt', trinucleotideina ,muts_ina, 'WT_mutationrateinsertion_results') 
mutationRatesinb = MutationRate('di_mutationrates.txt', trinucleotideinb ,muts_inb, 'rad16_mutationrateinsertion_results') 
#mutationRatesdel = MutationRate('mutationrates.txt', trinucleotidedel, muts_del, 'mutationratedeletion_results') 
#mutationRatesdela = MutationRate('mutationrates.txt', trinucleotidedela, muts_dela, 'WT_mutationratedeletion_results') 
#mutationRatesdelb = MutationRate('mutationrates.txt', trinucleotidedelb, muts_delb, 'rad16_mutationratedeletion_results') 

def MutationRate_chr(file1, trinucleotides, chr_name , file2):
    mutationRates = {}
    fm = open(file1)
    datam = fm.read()
    fm2 = open(file2, 'w+')
    for line in datam.splitlines():
        line = datam.strip().split()
    for d in range(0, len(line)):
        mutationRates[line[d]] = 0
    for key in mutationRates.keys():
        for t in range (0, len(trinucleotides)):
            #print(key)
            if trinucleotides[t] == key and chromosome[t] == chr_name:
                mutationRates[key] = mutationRates[key] + 1 
    for key in mutationRates.keys():
        fm2.write(key + ": Number = " + str(mutationRates[key]))
        fm2.write('\n')
    return(mutationRates)
 

rep_times = []
rep_times.append(early_1)
rep_times.append(early_2)
rep_times.append(middle_1)
rep_times.append(middle_2)
rep_times.append(late_1)
rep_times.append(late_2)

def MakeChromosomeFile(file1, filea, fileb, filec, data, chromosome, number, times, muts):
    early = muts[0]
    middle = muts[1]
    late = muts[2]
    f2 = open(file1, 'w+')
    fa = open(filea, 'w+')
    fb = open(fileb, 'w+')
    fc = open(filec, 'w+')
    e = 0
    for d in data.splitlines():
        if e < len(times):
            if chromosome[e] == number:
                d = d + '\t' + str(times[e])
                f2.write(d)
                f2.write('\n')
                if e in early.keys():
                    fa.write(d)
                    fa.write('\n')
                if e in middle.keys():
                    fb.write(d)
                    fb.write('\n')
                if e in late.keys():
                    fc.write(d)
                    fc.write('\n')

        e = e + 1
    f1.close()
    f2.close()
    fa.close()
    fb.close()
    fc.close()


#MakeChromosomeFile('chr16_WT_muts','chr16_WT_early', 'chr16_WT_middle', 'chr16_WT_late', data, chromosome, 'chrXVI', times, muts2)
#MakeChromosomeFile('chr16_muts_rad16','chr16_early_rad16', 'chr16_middle_rad16', 'chr16_late_rad16', data16, chromosome16, 'chrXVI', times16, muts2_16)
#MakeChromosomeFile('chr16_muts_rad30','chr16_early_rad30', 'chr16_middle_rad30', 'chr16_late_rad30', data30, chromosome30, 'chrXVI', times30, muts2_30)
#MakeChromosomeFile('chr16_muts_UVA','chr16_early_UVA', 'chr16_middle_UVA', 'chr16_late_UVA', dataUVA, chromosomeA, 'chrXVI', times_A, muts2_A)
#MakeChromosomeFile('chr16_muts_ogg1','chr16_early_ogg1', 'chr16_middle_ogg1', 'chr16_late_ogg1', dataOg, chromosomeOg, 'chrXVI', times_Og, muts2_og)
#MakeChromosomeFile('chr16_muts_rad26','chr16_early_rad26', 'chr16_middle_rad26', 'chr16_late_rad26', data26, chromosome26, 'chrXVI', times26, muts2_26)
#MakeChromosomeFile('chr16_muts_UVB','chr16_early_UVB', 'chr16_middle_UVB', 'chr16_late_UVB', dataUVB, chromosomeB, 'chrXVI', times_B, muts2_B)
#MakeChromosomeFile('chr16_muts_UVB_rad16','chr16_early_UVB_rad16', 'chr16_middle_UVB_rad6', 'chr16_late_UVB_rad16', dataB16, chromosomeB16, 'chrXVI', timesB16, muts2_B16)

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
        for x in range(0, 249):
            if sequence[x:x+num] not in early_TrinucCounts.keys():
                early_TrinucCounts[sequence[x:x+num]] = 1
        
            else:
                early_TrinucCounts[sequence[x:x+num]] = early_TrinucCounts[sequence[x:x+3]] + 1
    elif chr_scores[0] <= middle:
        for y in range(0, 249):
            if sequence[y:y+num] not in middle_TrinucCounts.keys():
                middle_TrinucCounts[sequence[y:y+num]] = 1
            else:
                middle_TrinucCounts[sequence[y:y+num]] = middle_TrinucCounts[sequence[y:y+num]] + 1
    else:
        for z in range(0, 249):
            if sequence[z:z+num] not in late_TrinucCounts.keys():
                late_TrinucCounts[sequence[z:z+num]] = 1
            else:
                late_TrinucCounts[sequence[z:z+num]] = late_TrinucCounts[sequence[z:z+num]] + 1
    for i in range(1, len(scores)):
        if chr[i] == name:
            if scores[i] <= early:
                for j in range(int(start[i - 1]), int(start[i ])):
                    if sequence[j:j+num] not in early_TrinucCounts.keys():
                        early_TrinucCounts[sequence[j:j+num]] = 1
                    else:
                        early_TrinucCounts[sequence[j:j+num]] = early_TrinucCounts[sequence[j:j+num]] + 1

    for k in range(1, len(scores)):
        if chr[k] == name:
            if scores[k] > early and scores[k] <= middle:
                for L in range(int(start[k - 1]), int(start[k ])):
                    if sequence[L:L+num] not in middle_TrinucCounts.keys():
                        middle_TrinucCounts[sequence[L:L+num]] = 1
                    else:
                        middle_TrinucCounts[sequence[L:L+num]] = middle_TrinucCounts[sequence[L:L+num]] + 1
    for m in range(1, len(scores)):
        if chr[m] == name:
            if scores[m] > middle:
                for n in range(int(start[m -1]), int(start[m])):
                    if sequence[n:n+num] not in late_TrinucCounts.keys():
                        late_TrinucCounts[sequence[n:n+num]] = 1
                    else:
                        late_TrinucCounts[sequence[n:n+num]] = late_TrinucCounts[sequence[n:n+num]] + 1
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
#counts = TrinucleotideContext('chr16.txt', 'chrXVI', start, scores, early_2, middle_2, dict15, 3)

early_Dinuc = {}
middle_Dinuc = {}
late_Dinuc = {}

dicts2 = [early_Dinuc, middle_Dinuc, late_Dinuc]
dict1b = TrinucleotideContext('chr1.fa', 'chrI', start, scores, early_2, middle_2, dicts2, 2)
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

#chr16_counts =   TrinucleotideContext('chr16.txt', 'chrXVI', start, scores, early_2, middle_2, chr_dicts, 3)  
                                                   
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
        f.write(key + ': ' + str((early_totals[key] * mutationRates[key])))
        f.write('\n')
        f.write('Percentage: ' + str(((early_totals[key] * mutationRates[key])/(early_totals[key] + middle_totals[key] + late_totals[key]))*100))
        f.write('\n')
    f.write('Middle')
    f.write('\n')
    for key in middle_totals.keys():
        f.write(key + ': ' + str((middle_totals[key] * mutationRates[key])))
        f.write('\n')
        f.write('Percentage: ' + str(((middle_totals[key] * mutationRates[key])/(early_totals[key] + middle_totals[key] + late_totals[key]))*100))
        f.write('\n')
    f.write('Late')
    f.write('\n')
    for key in late_totals.keys():
        f.write(key + ': ' + str((late_totals[key] * mutationRates[key])))
        f.write('\n')
        f.write('Percentage: ' + str(((late_totals[key] * mutationRates[key])/(early_totals[key] + middle_totals[key] + late_totals[key]))*100))
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

#FindExpected(mutationRates, counts, 'Rates.txt')
#FindExpected(mutationRates16, counts, 'Rates16.txt')
#FindExpected(mutationRatesr1, counts, 'Rates_r1.txt')
#FindExpected(mutationRatesr2, counts, 'Rates_r2.txt')
#FindExpected(mutationRates30, counts)
#FindExpected(mutationRates26, counts, 'Rates26.txt')
#FindExpected(mutationRatesOg, counts)
#FindExpected(mutationRatesUVA, counts)
#FindExpected(mutationRatesUVB, counts)
#FindExpected(mutationRatesB16, counts)

#expected_percentages = FindExpected(mutationRates, chr16_counts)
#expected_percentages16 = FindExpected(mutationRates16, chr16_counts)
#expected_percentagesUVB = FindExpected(mutationRatesUVB, chr16_counts)
#expected_percentagesB16 = FindExpected(mutationRatesB16, chr16_counts)
FindExpected(mutationRatesin, counts2, 'Rates_insertions.txt')
FindExpected(mutationRatesina, counts2, 'WT_Rates_insertions.txt')
FindExpected(mutationRatesinb, counts2, 'rad16_Rates_insertions.txt')
#FindExpected(mutationRatesdel, counts, 'Rates_deletions.txt')
#FindExpected(mutationRatesdela, counts, 'WT_Rates_deletions.txt')
#FindExpected(mutationRatesdelb, counts, 'rad16_Rates_deletions.txt')

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
    f2.write('Expected Late: ' + str(int(expected_late_1)))
    f2.write('\n')
    f2.write('Mutation Density: ' + str(density1))
    f1.close()
    f2.close()
    
    
#MutationNumbers('chr16.txt', 'chr16_value.txt', chromosome, expected_percentages,'chrXVI', muts2)
#MutationNumbers('chr16.txt', 'chr16_value_rad16.txt', chromosome16, expected_percentages16,'chrXVI', muts2_16)
#MutationNumbers('chr16.txt', 'chr16_value_UVB', chromosomeB, expected_percentagesUVB,'chrXVI', muts2_B)
#MutationNumbers('chr16.txt', 'chr16_value_UVBrad16', chromosomeB16, expected_percentagesB16,'chrXVI', muts2_B16)