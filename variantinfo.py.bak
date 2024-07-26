import myvariant
import pymysql
import pymysql.cursors
import re


def rev_comp(seq, cap=True, rev=True):
    if rev:
        revseq = seq[::-1]
    else:
        revseq = seq
    if cap:
        revseq = revseq.upper()
    revcomp = ''
    for item in revseq:
        if item == "T":
            revcomp += "A"
        elif item == "C":
            revcomp += "G"
        elif item == "A":
            revcomp += "T"
        elif item == "G":
            revcomp += "C"
        else:
            revcomp += item
    return revcomp


def exon_mapping(cdstart, txstart, cdend, starts, ends, strand):
    coding_starts = []
    coding_ends = []
    gap = cdstart - txstart
    i = 0
    temp_end = gap
    while i<len(starts) and i<len(ends):
        temp_start = temp_end
        temp_end = temp_start + abs(int(ends[i]) - int(starts[i]))
        coding_starts += [temp_start +1]
        i += 1
    coding_starts[0] = 0
    return coding_starts


def get_gen(chgvs, exon_coord, ends, chromo, strand):
    neg_intron = 0
    pos_intron = 0
    if "-" in chgvs:
        nums = re.findall("\d+", chgvs)
        neg_intron = int(re.findall("\d+",chgvs.split("-")[1])[0])
        tail = chgvs[chgvs.index(nums[len(nums)-1])+len(nums[len(nums)-1]):]
        chgvs = chgvs.split("-")[0]
    elif "+" in chgvs:
        nums = re.findall("\d+", chgvs)
        pos_intron = int(re.findall("\d+",chgvs.split("+")[1])[0])
        tail = chgvs[chgvs.index(nums[len(nums)-1])+len(nums[len(nums)-1]):]
        chgvs = chgvs.split("+")[0]
    else:
        chgvs = chgvs.split("+")[0]
        nums = re.findall("\d+", chgvs)
        tail = chgvs[chgvs.index(nums[len(nums)-1])+len(nums[len(nums)-1]):]
    gen = []
    for item in nums:
        i = 0
        while i<len(exon_coord)-1:
            if int(item) >= int(exon_coord[i]) and int(item) < int(exon_coord[i+1]):
                diff = int(exon_coord[i+1]) - int(item)
                gen += [int(ends[i]) + diff]
            i += 1
    if strand == "-":
        gen[0] = gen[0] + neg_intron - pos_intron
    if strand == "+":
        gen[0] = gen[0] - neg_intron + pos_intron
    hgvs = chromo+":g."+str(gen[0])+rev_comp(tail, rev=False)
    return hgvs


def get_info(gene_symbol, hg, database='hg19'):
    ##db_index = input("""
    ##Select database:
    ##0: hg38
    ##1: hg19
    ##""")
    ##
    ##database = ['hg38','hg19'][db_index]
    conn = pymysql.connect(host='genome-mysql.cse.ucsc.edu',
                                 user='genomep',
                                 password='password',
                                 db=database)

    cur = conn.cursor()

    #gene_symbol = raw_input("Please enter the gene symbol: ")
    gene_symbol = gene_symbol.upper()
    statement1 = "select name from refGene where name2 = '%s'" % gene_symbol
    cur.execute(statement1)
    temp = cur.fetchall()
    fuzzy = False
    if len(temp[0]) == 0:
        fuzzy = True
        temp = []
        i = 1
        while i<len(gene_symbol):
            temp1 = gene_symbol[0:len(gene_symbol)-i] + "_" + gene_symbol[len(gene_symbol)-i+1:]
            test = "select name2 from refGene where name2 like '%s'" % temp1
            cur.execute(test)
            holder = cur.fetchall()
            for item in holder:
                if len(item) > 1:
                    for subitem in item:
                        print subitem
                        temp += [subitem[0]]
                else:
                    temp += [item[0]]
            i += 1
    if len(temp[0]) == 0:
        i = 3
        while i < len(gene_symbol):
            temp1 = gene_symbol[0:i] + "%"
            test = "select name2 from refGene where name2 like '%s'" % temp1
            cur.execute(test)
            holder = cur.fetchall()
            temp += [holder]
            i += 1

    genes = []
    for item in temp:
        if item not in genes:
            genes += [item]

    if fuzzy:
        print "Did you mean?:"
        for index, item in enumerate(genes):
            print index,":", item
        print "none : none"
        gene_index = raw_input("Enter your choice: ")
        try:
            gene_index = int(gene_index)
            statement1 = "select name2 from refGene where name2 = '%s'" % genes[gene_index]
            cur.execute(statement1)
            temp = cur.fetchall()
        except:
            temp = []
            i = 0
            while i < len(gene_symbol)+1:
                temp1 = gene_symbol[0:i+2] + "%"
                test = "select name2 from refGene where name2 like '%s'" % temp1
                cur.execute(test)
                holder = cur.fetchall()
                temp += [holder]
                i += 1
            genes = []
            for item in temp[0]:
                if item[0] not in genes:
                    genes += [item[0]]
            print "Did you mean?:"
            for index, item in enumerate(genes):
                print index,":", item
            print "none : none"
            gene_index = raw_input("Enter your choice: ")
            try:
                gene_index = int(gene_index)
                statement1 = "select name from refGene where name2 = '%s'" % genes[gene_index]
                cur.execute(statement1)
                temp = cur.fetchall()
            except:
                print "Damn"
        
    transcripts = []
    for index, item in enumerate(temp):
        transcripts += [item[0]]
    if len(transcripts) > 1:
        print index,": ",item[0]
        trans_choice = input("Please select a transcript from the list above: ")
    else:
        trans_choice = 0
    transcript = transcripts[trans_choice]
    statement2 = "select * from refGene where name = '%s'" % transcript
    cur.execute(statement2)
    info = cur.fetchall()
    info = info[0]


    txstart=info[4]
    txend=info[5]
    cdstart=info[6]
    cdend=info[7]
    chrom=info[2]
    exonstarts = info[9].split(",")
    exonends = info[10].split(",")
    while '' in exonstarts:
        exonstarts.remove('')
    while '' in exonends:
        exonends.remove('')
    strand = info[3]
    if strand == "-":
        temp = exonstarts
        temp2 = exonends
        exonstarts = temp2[::-1]
        exonends = temp[::-1]
        
        temp = txstart
        temp2 = txend
        txstart = temp2
        txend = temp

        temp = cdstart
        temp2 = cdend
        cdstart = temp2
        cdend = temp


    mv = myvariant.MyVariantInfo()
    #hg=raw_input("Enter variant: ")
    #recreate coding mapping
    es = exon_mapping(cdstart, txstart, cdend, exonstarts, exonends, strand)
    genomic = get_gen(hg, es, exonends, chrom, strand)

    query = "'%s' AND %s" % (hg, gene_symbol)
    cv = mv.query(query)

    hits = []
    for item in cv['hits']:
        if item['_id'] == genomic:
            hits += [item]
    return hits, genomic

#convert the variant to genomic and match the "_id" from myvariant hits
##or just do not specify field and search based only on "_id"


##aa
##genename
##fathmm-mkl
    ##coding_pred
##alt
##clinvar
    ##clinsig
##interpro_domain
##mutationtaster
    ##pred
##rsid
##phylo
##phastcons
##ref
##ancestral_allele
##gerp++
##ensembl
##lrt
        ##pred
##chrom
##dann
##hg38
##genocanyon
##hg18
##hg19
##siphy_29way
