#!/usr/bin/python
import pandas as pd

def aligned(line, col):
    mirna = line.split('\t')[col]
    if line.split('_')[0] == mirna:
        return 'yes'
    return 'no'


def same(mirna, fastq_name):
    if fastq_name.split('_')[0] == mirna:
        return 'yes'
    return 'no'


def _sam(fn, tool, writer, data, dataseq, lendata, quals):
    check, save = {}, {}
    with open(fn) as in_handle:
        for line in in_handle:
            cols = line.split("\t")
            if line.startswith("@"):
                continue
            if cols[2].find("hsa") == -1:
                continue
            if cols[0] not in check:
                check[cols[0]] = 100
            score = [len(flag) for flag in cols if flag.startswith('MD')]
            if score:
                score = score[0]
            else:
                score = 0
            if (cols[2].find("hsa") > -1) and score < check[cols[0]]:
                check[cols[0]] = score
                slot = cols[2].split("-")
                add = cols[0].find("add:null")
                mut = cols[0].find("mut:null")
                save[cols[0]] = "%s\t%s\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                cols[0], cols[2], cols[0].split("_")[0], add, mut,
                lendata[cols[0]], quals[cols[0]], aligned(line, 2), tool)

    for k in save:
        writer.write(save[k])
    # if not aligned, print them here
    for k in data.keys():
        if k not in check:
            add = k.find("add:null")
            mut = k.find("mut:null")
            writer.write("%s\t%s\tno\tNA\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            k, data[k], add, mut, lendata[k], quals[k], 'NA', tool))


def read_sim(file_name):
    sim = open(file_name, 'r')
    data = {}
    dataseq = {}
    lendata = {}
    quals = {}
    name = ""
    qual_info = False
    for line in sim:
        line = line.strip()
        if line.startswith("@"):
            name = line.replace("@", "")
            slot = name.split("_")[0].split("-")
            data[name] = "-".join(slot[0:3])
        elif line.startswith("+"):
            qual_info = True
        elif not qual_info:
            lendata[name] = len(line)
            # print name
            dataseq[line] = name
        else:
            quals[name] = line.rstrip()
            qual_info = False
    sim.close()
    return data, dataseq, lendata, quals


def write_dict(name, dict_to_write):
    writer = open(name, 'w')
    for key, value in dict_to_write.items():
        writer.write(str(key) + "\t" + str(value) + "\n")
    writer.close()


def write_dicts(base_name, data, dataseq, lendata, quals):
    write_dict(base_name + "_data", data)
    write_dict(base_name + "_dataseq", dataseq)
    write_dict(base_name + "_lendata", lendata)
    write_dict(base_name + "_quals", quals)


def microrazers(input_name, results, tool, data, dataseq, lendata, quals):
    check = {}
    mir = open(input_name, 'r')
    for line in mir:
        cols = line.split("\t")
        if (cols[4].find("hsa") >= 0) and cols[0] not in check:
            slot = cols[2].split("-")
            add = cols[0].find("add:null")
            mut = cols[0].find("mut:null")
            check[cols[0]] = 1
            results.write("%s\t%s\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            cols[0], cols[4], cols[0].split("_")[0], add, mut,
            lendata[cols[0]], quals[cols[0]], aligned(line, 4), tool))
    mir.close()

    # if not aligned, print them here
    for k in data:
        if k not in check:
            add = k.find("add:null")
            mut = k.find("mut:null")
            results.write("%s\t%s\tno\tNA\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            k, data[k], add, mut, lendata[k], quals[k], 'NA', tool))


def real_value(dict, key):
    if key in dict:
        return dict[key]
    shortest = ""
    for k in dict:
        if k.startswith(key) and (shortest == "" or len(key) < len(shortest)):
            shortest = k
    return dict[k]


def mirexpress(input_name, results, tool, data, dataseq, lendata, quals):
    check = {}
    mir = open(input_name, 'r')
    ref = ""
    for line in mir:
        if line.startswith("hsa"):
            ref = line.strip()
            next(mir)
            next(mir)
            # mir.next()
            continue
        cols = line.split()
        if (line.find("hsa") < 0) and (line.find("1") >= 0) and cols[
            0] not in check:
            name = real_value(dataseq, cols[0])
            slot = name.split("-")
            add = name.find("add:null")
            mut = name.find("mut:null")
            check[cols[0]] = 1
            results.write("%s\t%s\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            name, ref, name.split("_")[0], add, mut, lendata[name],
            quals[name], same(ref, name), tool))
    mir.close()

    # if not aligned, print them here
    for k in dataseq.keys():
        if k not in check:
            name = dataseq[k]
            add = name.find("add:null")
            mut = name.find("mut:null")
            results.write("%s\t%s\tno\tNA\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            name, data[name], add, mut, lendata[name], quals[name], 'NA',
            tool))




def quagmir(simFile,resFile, file):

    df = pd.read_csv(simFile, header=None)
    readsDF = pd.DataFrame({  'NAME'     : df[0][df.index % 4 == 0 ].tolist(),
                       'SEQUENCE' : df[0][df.index % 4 == 1 ].tolist(),
                       'QUALS'    : df[0][df.index % 4 == 3 ].tolist()})

    readsDF.NAME = readsDF.NAME.apply(lambda x: x.replace('@',''))

    resultsDF = pd.DataFrame.from_csv(resFile, sep='\t').reset_index()
    resultsDF = resultsDF[['MIRNA','SEQUENCE', 'READS', 'LEN_READ']]

    allDF = pd.merge(readsDF, resultsDF, how='outer', on=['SEQUENCE'])

    allDF['TOOL'] = 'QuagmiR'
    allDF['NAME_MIR'] = allDF.NAME.apply(lambda x: x.split("_")[1])
    allDF['ADD'] = allDF.NAME.apply(lambda x: x.find('add:null'))
    allDF['MUT'] = allDF.NAME.apply(lambda x: x.find('mut:null'))

    aligned = allDF[pd.notnull(allDF.MIRNA)]

    notAligned = allDF[pd.notnull(allDF.MIRNA)==False].copy()
    notAligned['LEN_READ'] = notAligned.SEQUENCE.apply(lambda x: len(x))
    notAligned['MIRNA'] = notAligned.NAME.apply(lambda x: x.split('_')[1])
    notAligned['ALIGNED']  = 'NA'
    notAligned['GOOD_ALIGNED']  = 'no'

    goodAligned = aligned[aligned.apply(lambda x: x.MIRNA.startswith(x.NAME_MIR)>-1, axis=1)].copy()
    goodAligned['ALIGNED'] = 'yes'
    goodAligned['GOOD_ALIGNED'] = 'yes'
    notAligned['GOOD_ALIGNED']  = 'no'
    badAligned = aligned[aligned.apply(lambda x: x["MIRNA"].startswith(x.NAME_MIR)==-1, axis=1)]
    badAligned['GOOD_ALIGNED']  = 'no'


    df = goodAligned.append(badAligned)
    df = df.append(notAligned)

    columns = ['NAME', 'NAME_MIR', 'GOOD_ALIGNED', 'MIRNA', 'ADD', 'MUT', 'LEN_READ', 'QUALS', 'ALIGNED', 'TOOL']
    with open(file, 'a') as f:
        df.sort_index()[columns].to_csv(f, sep='\t', header=False, index=False)



results = open('results.txt', 'w')

s3_1 = read_sim('simulations/3p_1.fastq')
s3_2 = read_sim('simulations/3p_2.fastq')
s5_1 = read_sim('simulations/5p_1.fastq')
s5_2 = read_sim('simulations/5p_2.fastq')

star_3p_2 = read_sim("simulations/star_3p_2.fastq")
star_5p_2 = read_sim("simulations/star_5p_2.fastq")

# write_dicts("aaa", *star_5p_2)
# write_dicts("s3_1", *s3_1)

microrazers('outputs/test_3p-ran5_1v.razers', results, 'microrazers_3p_1', *s3_1)
microrazers('outputs/test_3p-ran5_2v.razers', results, 'microrazers_3p_2', *s3_2)
microrazers('outputs/test_5p-ran5_1v.razers', results, 'microrazers_5p_1', *s5_1)
microrazers('outputs/test_5p-ran5_2v.razers', results, 'microrazers_5p_2', *s5_2)

# write_dicts("s3_1", *s3_1)

mirexpress('outputs/_2_test_3p-ran5_1v.alignment_result', results, 'mirexpress_3p_1', *s3_1)
mirexpress('outputs/_1_test_3p-ran5_2v.alignment_result', results, 'mirexpress_3p_2', *s3_2)
mirexpress('outputs/_1_test_5p-ran5_1v.alignment_result', results, 'mirexpress_5p_1', *s5_1)
mirexpress('outputs/_1_test_5p-ran5_2v.alignment_result', results, 'mirexpress_5p_2', *s5_2)

_sam('outputs/_1_test_3p-ran5_1v.sam', 'razer3_3p_ran5_1v', results, *s3_1)
_sam('outputs/_1_test_3p-ran5_2v.sam', 'razer3_3p_ran5_2v', results, *s3_2)
_sam('outputs/_2_test_5p-ran5_1v.sam', 'razer3_5p_ran5_1v', results, *s5_1)
_sam('outputs/_1_test_5p-ran5_2v.sam', 'razer3_5p_ran5_2v', results, *s5_2)

_sam('outputs/test_3p-ran5_1v.fastq.Aligned.out.sorted.sam', 'star_3p_1v', results, *s3_1)
_sam('outputs/star_3p_2.sam', 'star_3p_2v', results, *star_3p_2)
_sam('outputs/test_5p-ran5_1v.fastq.Aligned.out.sorted.sam', 'star_5p_1v', results, *s5_1)
_sam('outputs/star_5p_2.sam', 'star_5p_2v', results, *star_5p_2)

results.close()

quagmir('simulations/3p_1.fastq','outputs/_2_test_3p-ran5_1v.fastq.isomir.sequence_info.tsv', 'results.txt')
quagmir('simulations/3p_2.fastq','outputs/_2_test_3p-ran5_2v.fastq.isomir.sequence_info.tsv', 'results.txt')
