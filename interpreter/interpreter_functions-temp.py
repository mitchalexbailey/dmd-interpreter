import os,sys
import subprocess
import re


##Functions for input interpretation
def exoninput(mut): #entr): #entr variable for use with GUI input
	"""
	Returns a boolean. 
		ex_input: boolean of whether or not input is in exon notation.

	Determines whether or not the input corresponds to an exon input (i.e. exon numbers instead of nucleotide positions).
	""" 
	ex_input = False
	if "X" in mut.upper():
		ex_input = True
	return ex_input





def parse(inp, mutype, exon_positions, ex_input):
	"""
	Returns two lists of ints:
		nums: list of numbers in input corresponding to cDNA positions.
		introns: list of two lists:
			intron_positive: list of [[cdna position, intron position downstream from cdna]]
			intron_negative: list of [[cdna position, intron position upstream from cdna]]

	Takes the user input and interprets the numbers.
	Distinguishes between exon positions and intron positions.
	Returns the positions to num_list global variable (nums2 local) and the intron positions to intron_list global (intron local).
	"""
	i = 0
	temp = []
	intron_positive = []
	intron_negative = []
	nums = []
	introns = []
	strnums = []
	num = re.findall("[-+]?\d+[\.]?\d*", inp)
	if len(num) == 0:
		return nums, introns
	if mutype == "Invalid mutation":
		num = re.findall("\d+",inp)
		for item in num:
			nums += [int(item)]
		return nums, introns
	if not ex_input:
		while i<len(num):
			if '+' not in num[i] and '-' not in num[i]:
				if num[i] not in strnums:
					nums += [int(num[i])]
					strnums += [num[i]]
			test_neg = num[i],'-?'
			test_pos = num[i],'+?'
			test_neg = ''.join(test_neg)
			test_pos = ''.join(test_pos)
			if '+' in num[i]:
				intron_positive += [int(num[i-1]),int(num[i])]
			if '-' in num[i]:
				intron_negative += [int(num[i-1]),int(num[i])]
			if test_pos in inp:
				intron_positive += [int(num[i]),'?']
			if test_neg in inp:
				intron_negative += [int(num[i]),'?']
			i+=1
	if ex_input:
		num = re.findall("\d+", inp)
		for char in num:
			temp += [int(char)]
		if max(temp) < 80 and min(temp) > 0:
			ind = min(temp) - 1
			start = exon_positions[ind][1]
			if len(temp)>1:
				ind = max(temp) - 1
				end = exon_positions[ind][2]
			else:
				end = exon_positions[ind][2]
			nums += [min(start,end)]
			nums += [max(start,end)]
		else:
			pass
	introns = [intron_positive, intron_negative]
	return nums, introns




def type_inp(inp, ex_input):
	"""
	Returns two strings.
		muttype: string of mutation type
		search_term: string used to build HGVS notation.

	Looks for key character combinations in keyboard input to interpret what the mutation type is.
	If no keys are found, invalid mutation type is returned, and all other processes are aborted (user receives new prompt after being told the input was invalid).
	"""
	type_specified = False
	inp = inp.upper()
	if "DEL" in inp and "INS" in inp:
		muttype = "Deletion/Insertion"
		search_term = "delins"
		type_specified = True
		return muttype, search_term
	if "DEL" in inp:
		muttype = "Deletion"
		search_term = "del"
		type_specified = True
		return muttype, search_term
	if "DUP" in inp:
		muttype = "Duplication"
		search_term = 'dup'
		return muttype, search_term
	if "INS" in inp:
		muttype = "Insertion"
		search_term = "ins"
		type_specified = True
		return muttype, search_term
	if ">" in inp or "POINT" in inp or "TRANS" in inp:
		muttype = "Point mutation"
		search_term = ">"
		type_specified = True
		return muttype, search_term
	if "TO" in inp and not ex_input:
		muttype = "Point mutation"
		search_term = ">"
		type_specified = True
		return muttype, search_term
	else:
		muttype = "Invalid mutation"
		search_term = ""
		type_specified = False
		return muttype, search_term



def intron_count(ilst, nlst):
	"""
	Returns one boolean:
		result: true/false for intron only variant
	
	Iterates through ilist (list of characters presumed to be intron positions by parse()) to indicate the number of intron positions versus the number of exon-only positions.
	If there are only intron positions present, then intron_count() returns true, indicating the user input is a variant only occuring in the intron.
	"""
	result = False
	intron_count = 0
	for item in ilst:
		for subitem in item:
			if subitem != "?":
				intron_count += 1
	if len(ilst) == 0 and len(nlst) == 0:
		return result
	elif intron_count > 0 and len(nlst) != 0:
		if len(nlst) == 1:
			result = True
		elif len(nlst) != 1 and intron_count > 1:
			if (max(nlst)-min(nlst)) == 1:
				result = True
			else:
				pass
		else:
			pass
	elif len(nlst) == 0 and (len(ilst[0]) != 0 or len(ilst[1]) != 0):
		result = True
	else:
		pass
	return result




def HGVS(positions, mtype, search_term, leng, ex_input, intron_list, pm_change, pm_pre, genomic=False):
	"""
	Returns a single string:
		hgvs: string of standard notation for mutation.

	Converts the user input into proper HGVS.
	This will also work if incorrect HGVS has been entered.
	"""
	hgvs = 'c.'
	if genomic:
		hgvs = 'ChrX:'
		pm_change = complement(pm_change)
		pm_pre = complement(pm_pre)
	if mtype == 'Deletion/Insertion':
		try:
			leng = leng[0]-leng[1]
			if leng < 0:
				leng = leng*(-1)
		except:
			leng = leng
	if len(positions)!=0:
		if leng > 30 and ex_input:
			notation = str(min(positions)), "-?_" ,str(max(positions)),"+?", search_term
		if not genomic:
			if not ex_input and len(positions) > 1:
				notation = str(min(positions)),"_",str(max(positions)), search_term, pm_change
			if not ex_input and (len(intron_list[0])+len(intron_list[1]))>2:
				if len(intron_list[0]) ==0 and len(intron_list[1])==0:
					notation = str(min(positions)),"_",str(max(positions)), pm_pre, search_term, pm_change
				if len(intron_list[0]) != 0 or len(intron_list[1]) != 0:
					if min(positions) in intron_list[0]:
						if max(positions) in intron_list[0] and len(intron_list[0])>1:
							notation = str(min(positions)),"+",str(intron_list[0][1]),"_",str(max(positions)),"+",str(intron_list[0][3]), pm_pre, search_term, pm_change
						if max(positions) in intron_list[1]:
							notation = str(min(positions)),"+",intron_list[0][1],"_",str(max(positions)),"-",str(intron_list[0][3]), pm_pre, search_term, pm_change
						if max(positions) not in intron_list[0] and max(positions) not in intron_list[1]:
							notation = str(min(positions)),"+",intron_list[0][1],"_",str(max(positions)), pm_pre, search_term, pm_change
					if min(positions) in intron_list[1]:
						if max(positions) in intron_list[0]:
							notation = str(min(positions)),"-",str(intron_list[0][1]),"_",str(max(positions)),"+",str(intron_list[1][1]), pm_pre, search_term, pm_change
						if max(positions) in intron_list[1] and len(intron_list[1])>1:
							notation = str(min(positions)),"-",intron_list[0][1],"_",str(max(positions)),"-",str(intron_list[0][3]), pm_pre, search_term, pm_change
						if max(positions) not in intron_list[0] and max(positions) not in intron_list[1]:
							notation = str(min(positions)),"-",intron_list[0][1],"_",str(max(positions)), pm_pre, search_term, pm_change
					if max(positions) in intron_list[0] and min(positions) not in intron_list[0] and min(positions) not in intron_list[1]:
						if max(positions) in intron_list[0]:
							notation = str(min(positions)),"_",str(max(positions)),"+",str(intron_list[0][3]), pm_pre, search_term, pm_change
						if max(positions) in intron_list[1]:
							notation = str(min(positions)),"_",str(max(positions)),"-",str(intron_list[0][3]), pm_pre, search_term, pm_change
			if mtype == "Deletion/Insertion":
				notation = str(min(positions)),"_",str(max(positions)),"delins", pm_change
			if len(positions) == 1 and not genomic:
				if (len(intron_list[0])+len(intron_list[1])) == 1:
					if len(intron_list[0]) == 2:
						notation = str(positions[0]),"+",str(intron_list[0][1]), search_term, pm_change
					if (len(intron_list[1]))==2:
						notation = str(positions[0]),str(intron_list[1][1]), search_term, pm_change
				if (len(intron_list[0])+len(intron_list[1])) == 0:
					notation = str(positions[0]), search_term, pm_change
			if len(positions) == 1 and genomic:
				notation = str(positions[0]), search_term, pm_change
			if mtype == "Point mutation":
				if len(intron_list[0]) ==0 and len(intron_list[1]) == 0:
					notation = str(positions[0]),pm_pre,search_term, pm_change
				elif len(intron_list[0]) != 0 and not genomic:
					notation = str(positions[0]),"+",str(intron_list[0][1]),pm_pre,search_term, pm_change
				elif len(intron_list[1]) != 0 and not genomic:
					notation = str(positions[0]),str(intron_list[1][1]),pm_pre,search_term, pm_change
				else:
					print "uncaught issue in HGVS()"
		else:
			if len(positions) == 1:
				notation = str(positions[0]),pm_pre,search_term, pm_change
			elif len(positions) > 1:
				notation = str(min(positions)), "_", str(max(positions)),pm_pre,search_term,pm_change
			else:
				print "uncaught issue in genomic HGVS()"
	notation = ''.join(notation)
	hgvs += notation
	return hgvs






def get_length(mutation, mutation_type, positions, mutation_inp, xon, dna=''):
	"""
	Returns one int and two strings:
		len_mut: integer, length of mutated sequence
		point_change: string, the mutant nucleotide(s)
		pre_change: string, the wild type nucleotide(s)

	Determines length of mutated DNA, returns: length of mutation (len_mut), the mutated nucleotides (point_change), the wildtype nucleotides at the mutated loci (pre_change)
	If numbers are available from input (ranges), they are subtracted to find the exact length.
	If there is no range, but rather the inserted sequence (ex. 456insTTT), then the number of characters following the key word is used to determine length.
	If the mutation type is a Deletion/Insertion, a list will be returned for len_mut (item 0) [deleted, inserted, total] as [int, int, str]
	"""
	len_mut = 0
	point_change = ''
	pre_change = ''
	if "Deletion/Insertion at Nucleotide position" in mutation_inp or "Point mutation at Nucleotide position" in mutation_inp or "Insertion at Nucleotide position" in mutation_inp:
		return 0,'', ''
	if mutation_type == "Deletion/Insertion":
		if "del" in mutation and "ins" in mutation:
			len_dele = (max(positions)+1) - min(positions)
		i = 0
		dele = ''
		ins = ''
		while i<len(mutation):
			if mutation[i] == ("l" or "L") and "delins" not in mutation:
				mut = mutation[i+1:]
				for char in mut:
					if char in "ATGCatgc":
						dele += char.upper()
					if char == "i" or char == "I":
						break
			if mutation[i] == ("s" or "S"):
				mut = mutation[i+1:]
				for char in mut:
					if char in "ATGCatgc":
						ins += char.upper()
					if char == "d" or char == "D":
						break
			i+=1
		dele = ''.join(dele)
		ins = ''.join(ins)
		if len(dele)!=0:
			total = len(ins)-len(dele)
			statement = "Total = ", str(total)
			statement = ''.join(statement)
			len_mut = [len(dele), len(ins), statement]
		if len(dele)==0:
			total = len(ins)-len_dele
			statement = "Total = ", str(total)
			statement = ''.join(statement)
			len_mut = [len_dele, len(ins), statement]
		point_change = ''.join(ins)
		pre_change = ''.join(dele)
	if mutation_type == "Point mutation":
		i = 0
		while i<len(mutation):
			if mutation[i] == ">":
				j=1
				while i-j > 0:
					if mutation[i-j] in "ATGCatgc":
						pre_change = mutation[i-j]
						break
					elif mutation[i-j] in "0123456789":
						pre_change = dna[positions[0] + 243]
						print "here: ", pre_change
						break
					j-=1
				mut = mutation[i+1:]
				point_change = ''
				for char in mut:
					if char != ' ' and char not in "0123456789":
						point_change += char.upper()
						len_mut += 1
				point_change = ''.join(point_change)
			if mutation[i].upper() == "O":# and mutation[i-1].upper() != "T":
				j=2
				while i-j > 0:
					if mutation[i-j] in "ATGCatgc":
						pre_change = mutation[i-j]
						break
					elif mutation[i-j] in "0123456789":
						pre_change = dna[positions[0] + 243]
						print "here: ", pre_change
						break
					j+=1   
				mut = mutation[i+1:]
				point_change = ''
				for char in mut:
					if char in "ATGCatgc":
						point_change += char.upper()
						len_mut += 1
				point_change = ''.join(point_change)
				i = len(mutation)
			else:
				i += 1
	if mutation_type == "Deletion" or mutation_type == "Duplication" or mutation_type == "Insertion":
		if len(positions) > 1 and mutation_type != "Insertion":
			len_mut = (max(positions)+1) - min(positions)
			return len_mut, point_change, pre_change
		else:
			if not xon:
				try:
					if mutation_type == "Deletion":
						trunc = mutation.index('del')
					if mutation_type == "Duplication":
						trunc = mutation.index('dup')
					if mutation_type == "Insertion":
						trunc = mutation.index('ins')
					mut = mutation[trunc+3:]
					for char in mut:
						if char in "ATGCatgc":
							len_mut += 1
							point_change += char.upper()
					if len_mut == 0:
						len_mut = 1
				except:
					len_mut = 1
		point_change = ''.join(point_change)
	point_change = point_change.upper()
	pre_change = pre_change.upper()
	return len_mut, point_change, pre_change




def frame(length, mutation_type, nums):
	"""
	Returns boolean (whether the mutated sequence is in frame or out of frame):
		fs: boolean, whether frameshift based on length of mutation.

	Uses the length of the mutated sequence to assign True or False to frame shift.
	Performs the first test of exon skip amenable test: does it correct the reading frame (length mutation + length skipped exon % 3).
	"""
	if mutation_type == "Point mutation":
		fs = False
		return fs
	if mutation_type == "Deletion/Insertion":
		overall = abs(length[0] - length[1])
		if overall%3 == 0:
			fs = False
		else:
			fs = True
	elif length%3 == 0:
		fs = False
	elif length%3 != 0:
		fs = True
	return fs




def exons(positions, ex_bounds, ex_input, mut, part, override=False):
	"""
	Returns list and a boolean:
		result: list of integers, exons directly effected by mutation.
		partial: boolean, does the mutation effect part of the exon

	Finds and stores the exon positions of the mutation from either:
	(1) The user input if they entered the mutation in exon format.
	(2) Comparison between the nucleotide positions (positions variable) and the exon position list (ex_bounds).
	First, direct matches are found to exon boundaries, then ranges are used.
	If ranges are used, the list is eliminated if all positions are in the same exon.
	"""
	partial = False
	if ex_input and override == False:
		result = re.findall("\d+", mut)
		return result, part
	exs = []
	result = []
	count_partial = 0
	for item in positions:
		for index, pos in enumerate(ex_bounds):
			if item <= pos[2] and item >= pos[1]:
				ind = index + 1
				exs += [ind]
	for char in exs:
		if char != 0 and char not in result:
			result += [char]
	for item in positions:
		i=0
		if partial:
			return result, partial 
		if not partial:
			partial = True
			while i<len(ex_bounds):
				if item in ex_bounds[i]:
					partial = False
				i += 1
	return result, partial


	
	
def splice_check(gen_pos, length_mutation, mut_genomic_dna, genomic_dna):
	"""
	Returns a list (len=7):
		affected: list of the significantly changed splice site consensus sequences.
			[0:wildtype fragment (str), 1:wildtype maxent score (float), 2:mutated fragment (str), 3:mutated score (float), 4:acceptor/donor (str), 5:change in score (str), 6:percent change (str)]

	Get fragments of genomic sequence based on cDNA mutations, send to maxent for scoring
	for 5' mutations, take 3 exonic and 6 intronic nucleotides, for 3' (double-check this)
	"""
	affected = []
	cryptic = []
	side = "Neither"
	if length_mutation > 8:
		affected += ["Splice predictions are not made for deletions of this size"]
		return affected
	if length_mutation <= 8:
		donor_mfrag = mut_genomic_dna[min(gen_pos)-8:min(gen_pos)+9]
		acceptor_mfrag = mut_genomic_dna[min(gen_pos)-22:min(gen_pos)+23]
		donor_frag = genomic_dna[min(gen_pos)-8:min(gen_pos)+9]
		acceptor_frag = genomic_dna[min(gen_pos)-22:min(gen_pos)+23]
		i = 0
		donor_frags = []
		donor_mfrags = []
		while i < len(donor_frag) - 8:
			donor_frags += [donor_frag[i:i+9]]
			donor_mfrags += [donor_mfrag[i:i+9]]
			i+=1
		side = 'Donor'
		ref_score = maxent(donor_frags, side)
		mut_score = maxent(donor_mfrags, side)
		for index, ref_item in enumerate(ref_score):
			if float(ref_item) > 3 or float(mut_score[index]) > 3:
				change = float(mut_score[index]) - float(ref_item)
				percent = (change/abs(float(ref_item))) *100
				if percent > 0:
					percent = "<c style='color:red'>+", str("%.2f" % percent), "% </c>"
					percent = ''.join(percent)
				if percent <= 0:
					percent = "<c style='color:red'>",str("%.2f" % percent),"% </c>"
					percent = ''.join(percent)
				if change > 0:
					change = "+", str(change)
					change = ''.join(change)
				if change != 0:
					affected += [donor_frags[index], ref_score[index], donor_mfrags[index], mut_score[index], "<b> Donor </b>", str(change), percent]
		i = 0
		acceptor_frags = []
		acceptor_mfrags = []
		while i < len(acceptor_frag) - 22:
			acceptor_frags += [acceptor_frag[i:i+23]]
			acceptor_mfrags += [acceptor_mfrag[i:i+23]]
			i+=1
		side = 'Acceptor'
		ref_score = maxent(acceptor_frags, side)
		mut_score = maxent(acceptor_mfrags, side)
		for index, ref_item in enumerate(ref_score):
			if float(ref_item) > 3 or float(mut_score[index]) > 3:
				change = float(mut_score[index]) - float(ref_item)
				percent = (change/abs(float(ref_item))) *100
				if percent > 0:
					percent = "<c style='color:red'>+", str("%.2f" % percent), "% </c>"
					percent = ''.join(percent)
				if percent <= 0:
					percent = "<c style='color:red'>",str("%.2f" % percent),"% </c>"
					percent = ''.join(percent)
				if change > 0:
					change = "+", str(change)
					change = ''.join(change)
				if change != 0:
					affected += [acceptor_frags[index], ref_score[index], acceptor_mfrags[index], mut_score[index], "<b> Acceptor </b>", str(change), percent]
	if len(affected)==0:
		affected += ["False"]
	return affected




def cdna2gen(num_list, exon_positions, genomic_positions1, genomic_positions2, intron_list):
	"""
	Returns a list:
		gen_positions: list of ints, genomic positions corresponding to those in the num_list (cDNA positions).
	
	The length of the list will be the same as that of num_list.
	"""
	pos_gen1 = []
	pos_gen2 = []
	for num in num_list:
		for line in exon_positions:
			if num >= line[1] and num <= line[2]:
				exon = line[0]
				difference_3 = num - (line[1])
				gen3 = genomic_positions1[exon][1]
				gen_position = gen3 + difference_3
				if len(intron_list[0]) != 0:
					if intron_list[0][0] == num and intron_list[0][1] != '?':
						gen_position = gen_position + intron_list[0][1]
				if len(intron_list[1]) != 0:
					if intron_list[1][0] == num and intron_list[1][1] != '?':
						gen_position = gen_position + intron_list[1][1]
				pos_gen1 += [gen_position]
	for num in num_list:
		for line in exon_positions:
			if num >= line[1] and num <= line[2]:
				exon = line[0]
				difference_3 = (line[2]) - num +1
				gen3 = genomic_positions2[exon][2]
				gen_position = gen3 + difference_3
				if len(intron_list[0]) != 0:
					if intron_list[0][0] == num and intron_list[0][1] != '?':
						gen_position = gen_position - intron_list[0][1]
				if len(intron_list[1]) != 0:
					if intron_list[1][0] == num and intron_list[1][1] != '?':
						gen_position = gen_position - intron_list[1][1]
				pos_gen2 += [gen_position]
	return pos_gen1, pos_gen2  


def gen_point(genomic_position, pm, intron_only, length_mutation, NSFP):
	"""
	Returns a list and an int:
		result: list of strings (len=9), all scoring algorithm results.
			[0:genomic positions of assessed mutation, 1:sift, 2:polyphen, 3:lrt, 4:mutationtaster, 5:mutationaccessor, 6:fathmm, 7:provean, 8:svm]
				Each list item (except [0]): string of votes, followed by description of vote categories.
		pathogenic_prediction: integer, total number of scoring algorithms who predicted a significantly negative impact.

	Uses the genomic positions of the mutation, and the change at that location (single changes assessed) to search the latest NSFP database.
	"""
	hit = []
	result = []
	empty = ["0","Not missense","Not missense","Not missense","Not missense","Not missense","Not missense","Not missense","Not missense"], 0
	if intron_only or length_mutation > 1:
		return empty
	gen_position = genomic_position[0]
	for index, line in enumerate(NSFP):
		if int(line[1]) == int(gen_position) and line[3] == pm:
			hit = NSFP[index]
	if len(hit) == 0:
		return empty
	else:
		result = [str(gen_position)]
		pathogenic_prediction = 0
		ds = 0
		ts = 0
		for char in hit[4]:
			if char == 'D':
				ds += 1
			if char == 'T':
				ts += 1
		if ds >= ts and ds+ts!=0:
			pathogenic_prediction += 1
		sift = "<b><c style='color:red'>",str(ds),"</c>/<c style='color:green'>",str(ts),'</c></b> (Damaging/Tolerated)'
		sift = ''.join(sift)
		result += [sift]
		dp = 0
		pp = 0
		tp = 0
		for char in hit[5]:
			if char == 'D':
				dp += 1
			if char == 'P':
				pp += 1
			if char == 'B':
				tp += 1
		if dp + pp >= tp and dp+pp+tp!=0:
			pathogenic_prediction += 1
		polyphen = "<b><c style='color:red'>",str(dp),"</c>/<c style='color:orange'>",str(pp),"</c>/<c style='color:green'>",str(tp), '</c></b> (Probably damaging/Possibly damaging/Benign)'
		polyphen= ''.join(polyphen)
		result += [polyphen]
		dl = 0
		pl = 0
		tl = 0
		for char in hit[6]:
			if char == 'D':
				dl += 1
			if char == 'U':
				pl += 1
			if char == 'N':
				tl += 1
		if dl >= pl + tl and dl+pl+tl != 0:
			pathogenic_prediction += 1
		lrt = "<b><c style='color:red'>",str(dl),'</c>/',str(pl),"/<c style='color:green'>",str(tl), '</c></b> (Deleterious/Unknown/Neutral)'
		lrt= ''.join(lrt)
		result += [lrt]
		ddmt = 0
		dmt = 0
		pmt = 0
		tmt = 0
		for char in hit[7]:
			if char == 'A':
				ddmt += 1
			if char == 'D':
				dmt += 1
			if char == 'N':
				tmt += 1
			if char == 'P':
				pmt += 1
		if (ddmt + dmt)>=(tmt+pmt) and (ddmt+dmt+tmt+pmt)!=0:
			pathogenic_prediction += 1
		mt = "<b><c style='color:red'>",str(ddmt),'/',str(dmt),"</c>/<c style='color:green'>",str(tmt),'/',str(pmt),'</c></b> (Disease causing automatic/Disease causing/Polymorphism/Polymorphism automatic)'
		mt= ''.join(mt)
		result += [mt]
		ddma = 0
		dma = 0
		pma = 0
		tma = 0
		for char in hit[8]:
			if char == 'H':
				ddma += 1
			if char == 'M':
				dma += 1
			if char == 'L':
				pma += 1
			if char == 'N':
				tma += 1
		if (ddma+dma)>(pma+tma) and (ddma+dma+pma+tma)!=0:
			pathogenic_prediction += 1
		ma = "<b><c style='color:red'>",str(ddma),"</c>/<c style='color:orange'>",str(dma),"</c>/<c style='color:orange'>",str(pma),"</c>/<c style='color:green'>",str(tma),'</c></b> (High impact/Medium/Low/Neutral)'
		ma= ''.join(ma)
		result += [ma]
		df = 0
		tf = 0
		for char in hit[9]:
			if char == 'D':
				df += 1
			if char == 'T':
				tf += 1
		if df>tf and df+df != 0:
			pathogenic_prediction += 1
		fat = "<b><c style='color:red'>",str(df),"</c>/<c style='color:green'>",str(tf),'</c></b> (Damaging/Tolerated)'
		fat = ''.join(fat)
		result += [fat]
		dpr = 0
		tpr = 0
		for char in hit[10]:
			if char == 'D':
				dpr += 1
			if char == 'N':
				tpr += 1
		if dpr>tpr and dpr+tpr!=0:
			pathogenic_prediction += 1
		pro = "<b><c style='color:red'>",str(dpr),"</c>/<c style='color:green'>",str(tpr),'</c></b> (Damaging/Neutral)'
		pro = ''.join(pro)
		result += [pro]
		dsvm = 0
		tsvm = 0
		for char in hit[11]:
			if char == 'D':
				dsvm += 1
			if char == 'T':
				tsvm += 1
		if dsvm>tsvm and dsvm+tsvm!=0:
			pathogenic_prediction += 1
		svm = "<b><c style='color:red'>",str(dsvm),"</c>/<c style='color:green'>",str(tsvm),'</c></b> (Damaging/Tolerated)'
		svm = ''.join(svm)
		result += [svm]
		return result, pathogenic_prediction


def maxent(sequence, side):
	"""
	Returns a single string:
		nums: string, maxent scores as line-separated strings (ex. "3.145\n5.678\n1.456")

	Sequences can be sent in as a list, nums will be sent as line-separated scores.
	Uses a Perl script subprocess. The script is directly from the authors/publishers of the Maxent algorithm.
	"""
	if side == "Neither":
		return None
	f = open("./interpreter/temp.txt",'w')
	for item in sequence:
		f.write(item)
		f.write('\n')
	f.close()
	var = "./interpreter/temp.txt"
	if side == "Acceptor":
		args_test = ['perl', './interpreter/score3.pl', var]
	if side == "Donor":
		args_test = ['perl', './interpreter/score5.pl', var]
	pipe = subprocess.Popen(args_test, stdout=subprocess.PIPE)
	out = pipe.communicate()
	nums = re.findall(r"[-+]?\d*\.\d+|\d+",str(out))
	return nums




def alter_cDNA(ref, positions, gen_pos, length, typ, pm_change, genomic_dna, intron_only):
	"""
	Returns two strings:
		seq: string, cDNA sequence.
		gen_seq: string, genomic sequence.

	Takes the cDNA reference sequence, and alters the sequence according to the position, length and type specified from the input and interpreted in earlier steps.
	The altered sequence is returned.
	"""
	start = min(positions) + 243
	start_gen = min(gen_pos)
	end = max(positions) + 243
	end_gen = max(gen_pos)
	if typ == 'Deletion':
		seq = ref[0:start] +  ref[(start+length):]
		seq = ''.join(seq)
		gen_seq = genomic_dna[0:start_gen] + genomic_dna[(end_gen + 1):]
		gen_seq = ''.join(gen_seq)
	if typ == 'Duplication':
		seq = ref[0:start] + ref[start:(start+length)] + ref[start:]
		seq = ''.join(seq)
		gen_seq = genomic_dna[0:start_gen] + genomic_dna[start_gen:end_gen+1] + genomic_dna[start_gen:]
		gen_seq = ''.join(gen_seq)
	if typ == 'Point mutation':
		seq = ref[0:start] + pm_change + ref[(start+length):]
		gen_seq = genomic_dna[0:start_gen] + pm_change + genomic_dna[(start_gen+1):]
		seq = ''.join(seq)
		gen_seq = ''.join(gen_seq)
	if typ == 'Insertion':
		seq = ref[0:start+1] + pm_change + ref[start+1:]
		gen_seq = genomic_dna[0:start_gen+1] + pm_change + genomic_dna[start_gen+1:]
		seq = ''.join(seq)
		gen_seq = ''.join(gen_seq)
	if typ == 'Deletion/Insertion':
		ins = pm_change
		seq = ref[0:start] + ins + ref[start + length[0]:]
		gen_seq = ref[0:start_gen] + ins + ref[start_gen + length[0]:]
		seq = ''.join(seq)
		gen_seq = ''.join(gen_seq)
	if intron_only:
		seq = ref
	return seq, gen_seq




def translate(sequence):
	"""
	Returns a single string:
		aa: string, amino acid sequence.

	The sequence will end with "*STOP*" -- this means that the stop index will be the first '*'
		The only exception to this is in frame duplications, or rare mutations near the end of the sequence which make the reading frame extend.
	Translates nucleotide sequence to corresponding amino acids.
	"""
	aa = ''
	i = 0
	if len(sequence) == 0:
		return aa
	while i<len(sequence)-3:
		codon = sequence[i:i+3]
		if codon[0] == "T":
			if codon[1] == "T":
				if codon[2] == "T" or codon[2] == "C":
					aa += "F"
				if codon[2] == "A" or codon[2] == "G":
					aa += "L"
			if codon[1] == "C":
				aa += "S"
			if codon[1] == "A":
				if codon[2] == "T" or codon[2] == "C":
					aa += "Y"
				else:
					aa += "*STOP*"
					i = len(sequence)
			if codon[1] == "G":
				if codon[2] == "T" or codon[2] == "C":
					aa += "C"
				if codon[2] == "G":
					aa += "W"
				if codon[2] == "A":
					aa += "*STOP*"
					i = len(sequence)
		if codon[0] == "C":
			if codon[1] == "T":
				aa += "L"
			if codon[1] == "C":
				aa += "P"
			if codon[1] == "A":
				if codon[2] == "A" or codon[2] == "G":
					aa += "Q"
				if codon[2] == "T" or codon[2] == "C":
					aa += "H"
			if codon[1] == "G":
				aa += "R"
		if codon[0] == "A":
			if codon[1] == "T":
				if codon[2] == "T" or codon[2] == "C" or codon[2] == "A":
					aa += "I"
				if codon[2] == "G":
					aa += "M"
			if codon[1] == "C":
				aa += "T"
			if codon[1] == "A":
				if codon[2] == "A" or codon[2] == "G":
					aa += "K"
				if codon[2] == "T" or codon[2] == "C":
					aa += "N"
			if codon[1] == "G":
				if codon[2] == "T" or codon[2] == "C":
					aa += "S"
				if codon[2] == "A" or codon[2] == "G":
					aa += "R"
		if codon[0] == "G":
			if codon[1] == "T":
				aa += "V"
			if codon[1] == "C":
				aa += "A"
			if codon[1] == "A":
				if codon[2] == "T" or codon[2] == "C":
					aa += "D"
				if codon[2] == "A" or codon[2] == "G":
					aa += "E"
			if codon[1] == "G":
				aa += "G"
		i += 3
	aa = ''.join(aa)
	return aa




def easy_align(ref, mut, mutype, frame_shift, single=False):
	"""
	Returns a list (len=variable), and two strings:
		changes: list of lists:
			each item (len=2): [str:new symbol, int:index]
		alignment: string, the aligned sequences with '-' indicating a misalignment and '*' indicating the first stop
		ptc: the proper notation for a premature stop codon (SymNum*)

	Directly aligns two sequences (ref and mut) position by position.
	Any discrepancies are stored in changes[].
	The sequence disagreements are returned (changes[]), as well as a sequence with only the agreeing members, unequal positions are marked with '-'.
	The alignment is terminated by the first encountered stop codon (represented by '*').
	"""
	changes = []
	ptc = ''
	lst = ['A','Ala','C','Cys','D','Asp','E','Glu','F','Phe','G','Gly','H','His','I','Ile','K','Lys','L','Leu','M','Met','N','Asn','P','Pro','Q','Gln','R','Arg','S','Ser','T','Thr','V','Val','W','Trp','Y','Tyr','-','-']
	alignment = ''
	i = 0
	lengths = [len(ref), len(mut)]
	if (mutype == "Duplication" or mutype == "Insertion") and not frame_shift:
		return changes, alignment, ptc
	while i < max(lengths):
		if i < len(ref):
			refsym = ref[i]
		if i >= len(ref):
			refsym = '-'
		if i < len(mut):
			mutsym = mut[i]
		if i >= len(mut):
			mutsym = '-'
		if mutsym == "*" and refsym != "*" and len(ref) > len(mut):
			changes += ['STOP', i+1]
			alignment += '*'
			alignment = ''.join(alignment)
			index_minus1 = lst.index(refsym)
			refsym = lst[index_minus1+1]
			ptc = refsym, str(i+1), "*"
			ptc = ''.join(ptc)
			return changes, alignment, ptc
		if refsym == "*" and mutsym != "*" and len(ref) < len(mut): #switched conditional to correct logic
			alignment += 'ext'
			alignment = ''.join(alignment)
			index = len(mut) - 3685
			changes += ['STOP', index]
			start_index = re.findall("\d+", changes[0])
			index_minus1 = lst.index(changes[0][0])
			start_aa1 = lst[index_minus1 + 1]
			index_minus1 = lst.index(changes[0][len(changes[0])-1])
			start_aa2 = lst[index_minus1 + 1]
			ext = "p.(",start_aa1, start_index[0], start_aa2, "ext*", str(index), ")"
			ext = ''.join(ext)
			return changes, alignment, ext
		if mutsym != refsym and i < 3685:
			change = refsym, str(i+1), mutsym
			change = ''.join(change)
			changes += [change]
			alignment+= "-"
			if single and mutsym != '-':
				index_minus1 = lst.index(refsym)
				refsym = lst[index_minus1+1]
				index_minus1 = lst.index(mutsym)
				mutsym = lst[index_minus1+1]
				ptc = refsym, str(i+1), mutsym
				ptc = ''.join(ptc)
				return changes, alignment, ptc
		if mutsym == refsym:
			alignment += mutsym
		i+=1
	if changes == []:
		changes = ['STOP', 3686]
	alignment = ''.join(alignment)
	return changes, alignment, ptc




def possible_skips(ex_ints, length, refcdna, frame_shift, exon_positions, mutype, num_list, aa_ref):
	"""
	Returns a list (len=variable) of strings:
		skips: list, each item is a single string containining possible skip or skip combination and the change at the boundary if the skip is performed.

	Tests whether:
	(1) The reading frame is corrected.
	(2) Changes (including nonsense) created at the border.
	by skipping either:
	(1) The exon bordering the left of the mutation.
	(2) The exon bordering the right of the mutation.
	(3) Both the exon bordering the rigth and left of the mutation.
	"""
	if length < 32:
		return [""]
	if not frame_shift:
		return ["""There is no frame shift predicted, thus exon skipping to correct the reading frame is not predicted to have a therapeutic benefit."""]
	if 1 in ex_ints:
		return ["""This includes the first exon of the gene - exon skipping will not have an effect. An internal ribosome entry site has been described downstream of this deletion (<a target="_blank" href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2663021/">PubMed</a>)."""]
	if 79 in ex_ints:
		return ["""This includes the last exon of the <i>DMD</i> gene. Exon skipping to restore the reading frame is not predicted to have an effect."""]
	left_ind = min(ex_ints)-2
	right_ind = max(ex_ints)
	left2_ind = min(ex_ints)-3
	right2_ind = max(ex_ints)+1
	left2_length = 0
	right2_length = 0
	left_length = 0
	right_length = 0
	if left_ind > 0:
		left_length = exon_positions[left_ind][2]-exon_positions[left_ind][1]+1
	if right_ind != 79:
		right_length = exon_positions[right_ind][2]-exon_positions[right_ind][1]+1
	if left2_ind > 0:
		left2_length = exon_positions[left2_ind][2]-exon_positions[left2_ind][1]+1
	if right2_ind < 79:
		right2_length = exon_positions[right2_ind][2]-exon_positions[right2_ind][1]+1
	skips = []
	delta = []
	delta2 = []
	if mutype == "Deletion":
		if (length + left_length) %3 == 0 and left_length != 0:
			delta = []
			delta2 = []				
			cdnaskip = refcdna[0:(exon_positions[left_ind][1])+243] + refcdna[(exon_positions[left_ind][1]) + length + left_length + 243:]
			cdnaskip = ''.join(cdnaskip)
			aa_skip = translate(cdnaskip[244:])
			i = 0
			filler = ''
			while i < (length + left_length)/3:
				filler += '-'
				i+=1
			end = (exon_positions[left_ind][1])
			if end%3 == 1 or end%3 == 0:
				aa_skipgap = aa_skip[0:(end/3)] + filler + aa_skip[(end/3):]
			if end%3 == 2:
				aa_skipgap = aa_skip[0:(end/3)+1] + filler + aa_skip[(end/3)+1:]	            	
			aa_skipgap = ''.join(aa_skipgap)
			delta, align, PTC = easy_align(aa_ref, aa_skipgap, mutype, frame_shift)
			for index, item in enumerate(delta):
				if type(item) != int and '-' not in item:
					if item != 'STOP':
						delta2+=[item]
				try:
					if item == 'STOP' and delta[index+1] != 3686:
						delta2+=["<c style='color:red'><b>STOP CODON CREATED</b></c>"]
				except:
					pass
			if len(delta2) == 0:
				delta2+=['No amino acid change predicted']	            
			temp = "<c style='color:green'><b>Exon ", str(left_ind +1), "</b></c><br> (Change at boundary: ", delta2[0], ")"
			temp = ''.join(temp)
			skips += [temp]
		if (length + right_length) % 3 == 0 and right_length != 0:
			delta = []
			delta2 = []				
			cdnaskip = refcdna[0:min(num_list) + 243] + refcdna[min(num_list) + length + right_length + 243:]
			cdnaskip = ''.join(cdnaskip)
			aa_skip = translate(cdnaskip[244:])
			i = 0
			filler = ''
			while i < (length + right_length)/3:
				filler += '-'
				i += 1
			aa_skipgap = aa_skip[0:min(num_list)/3] + filler + aa_skip[(min(num_list))/3:]
			end = min(num_list)
			if end%3 == 1:
				aa_skipgap = aa_skip[0:(end/3)] + filler + aa_skip[(end/3):]
			if end%3 == 2:
				aa_skipgap = aa_skip[0:(end/3)+1] + filler + aa_skip[(end/3)+1:]
			aa_skipgap = ''.join(aa_skipgap)
			delta, align, PTC = easy_align(aa_ref, aa_skipgap, mutype, frame_shift)
			for index, item in enumerate(delta):
				if type(item) != int and '-' not in item:
					if item != 'STOP':
						delta2+=[item]
				try:
					if item == 'STOP' and delta[index+1] != 3686:
						delta2+=["<c style='color:red'><b>STOP CODON CREATED</b></c>"]
				except:
					pass
			if len(delta2) == 0:
				delta2+=['No amino acid change predicted']	
			temp = "<c style='color:green'><b>Exon ", str(right_ind+1), "</b></c><br> (Change at boundary: ", delta2[0], ")"
			temp = ''.join(temp)
			skips += [temp]
		if (length + right_length + left_length)%3 == 0 and right_length != 0 and left_length != 0:
			delta = []
			delta2 = []				
			cdnaskip = refcdna[0:(exon_positions[left_ind][1])+243] + refcdna[exon_positions[right_ind][2]+1+243:]
			cdnaskip = ''.join(cdnaskip)
			aa_skip = translate(cdnaskip[244:])
			i = 0
			filler = ''
			while i < (length + right_length + left_length)/3:
				filler += '-'
				i+=1
			end = (exon_positions[left_ind][1])
			if end%3 == 1 or end%3 == 0:
				aa_skipgap = aa_skip[0:(end/3)] + filler + aa_skip[(end/3):]
			if end%3 == 2:
				aa_skipgap = aa_skip[0:(end/3)+1] + filler + aa_skip[(end/3)+1:]		  
			aa_skipgap = ''.join(aa_skipgap)
			delta, align, PTC = easy_align(aa_ref, aa_skipgap, mutype, frame_shift)
			for index, item in enumerate(delta):
				if type(item) != int and '-' not in item:
					if item != 'STOP':
						delta2+=[item]
				try:
					if item == 'STOP' and delta[index+1] != 3686:
						delta2+=["<c style='color:red'><b>STOP CODON CREATED</b></c>"]
				except:
					pass
			if len(delta2) == 0:
				delta2+=['No amino acid change predicted']	
			temp = "<c style='color:green'><b>Exons ", str(left_ind+1), ' and ', str(right_ind + 1), "</b></c><br> (Change at boundary: ", delta2[0], ")"
			skips += [''.join(temp)]
		if (length + left_length + left2_length) % 3 == 0 and left2_length != 0:
			delta = []
			delta2 = []				
			cdnaskip = refcdna[0:(exon_positions[left2_ind][1])+243] + refcdna[max(num_list) + 1 + 243:]
			cdnaskip = ''.join(cdnaskip)
			aa_skip = translate(cdnaskip[244:])
			i = 0
			filler = ''
			while i < (length + left2_length + left_length)/3:
				filler += '-'
				i+=1
			end = (exon_positions[left2_ind][1])
			if end%3 == 2:
				aa_skipgap = aa_skip[0:(end/3)] + filler + aa_skip[(end/3):]
			if end%3 == 1:
				aa_skipgap = aa_skip[0:(end/3)+1] + filler + aa_skip[(end/3)+1:]		  
			aa_skipgap = ''.join(aa_skipgap)
			delta, align, PTC = easy_align(aa_ref, aa_skipgap, mutype, frame_shift)
			for index, item in enumerate(delta):
				print index, item
				if type(item) != int and '-' not in item:
					if item != 'STOP':
						delta2+=[item]
				try:
					if item == 'STOP' and delta[index+1] != 3686:
						delta2+=["<c style='color:red'><b>STOP CODON CREATED</b></c>"]
				except:
					pass
			if len(delta2) == 0:
				delta2+=['No amino acid change predicted']	           
			temp = "<c style='color:green'><b>Exon ", str(left2_ind +1), ' and ', str(left_ind +1), "</b></c><br> (Change at boundary: ", delta2[0], ")"
			temp = ''.join(temp)
			skips += [temp]
		if (length + right_length + right2_length) % 3 == 0 and right2_length != 0:
			delta = []
			delta2 = []				
			cdnaskip = refcdna[0:min(num_list)+243] + refcdna[min(num_list) + right_length + length + right2_length + 243:]
			cdnaskip = ''.join(cdnaskip)
			aa_skip = translate(cdnaskip[244:])
			i = 0
			filler = ''
			while i < (length + right2_length + right_length)/3:
				filler += '-'
				i+=1
			end = min(num_list)
			if end%3 == 1 or end%3 == 0:
				aa_skipgap = aa_skip[0:(end/3)] + filler + aa_skip[(end/3):]
			if end%3 == 2:
				aa_skipgap = aa_skip[0:(end/3)+1] + filler + aa_skip[(end/3)+1:]		  
			aa_skipgap = ''.join(aa_skipgap)
			delta, align, PTC = easy_align(aa_ref, aa_skipgap, mutype, frame_shift)
			for index, item in enumerate(delta):
				if type(item) != int and '-' not in item:
					if item != 'STOP':
						delta2+=[item]
				try:
					if item == 'STOP' and delta[index+1] != 3686:
						delta2+=["<c style='color:red'><b>STOP CODON CREATED</b></c>"]
				except:
					pass
			if len(delta2) == 0:
				delta2+=['No amino acid change predicted']	            
			temp = "<c style='color:green'><b>Exon ", str(right_ind +1), ' and ', str(right2_ind +1), "</b></c><br> (Change at boundary: ", delta2[0], ")"
			temp = ''.join(temp)
			skips += [temp]		               
	if mutype == "Duplication":
		if len(ex_ints)==1:
			temp = "<c style='color:green'><b>Exon ", str(ex_ints[0]), "</b></c> (Only tested <i>in vitro</i>)\n"
			temp = ''.join(temp)
			skips += [temp]
		if (length - left_length) % 3 == 0 and left_ind != 0:
			delta = []
			delta2 = []				
			cdnaskip = refcdna[0:(exon_positions[left_ind][1])+243] + refcdna[(exon_positions[left_ind][1]) + length + left_length + 243:]
			cdnaskip = ''.join(cdnaskip)
			aa_skip = translate(cdnaskip[244:])
			i = 0
			filler = ''
			while i < (length + left_length)/3:
					filler += '-'
					i+=1
			end = (exon_positions[left_ind][1])
			if end%3 == 2:
				aa_skipgap = aa_skip[0:(end/3)] + filler + aa_skip[(end/3):]
			if end%3 == 1:
				aa_skipgap = aa_skip[0:(end/3)+1] + filler + aa_skip[(end/3)+1:]		  
			aa_skipgap = ''.join(aa_skipgap)
			delta, align, PTC = easy_align(aa_ref, aa_skipgap, mutype, frame_shift)
			for index, item in enumerate(delta):
				if type(item) != int and '-' not in item:
					if item != 'STOP':
						delta2+=[item]
				try:
					if item == 'STOP' and delta[index+1] != 3686:
						delta2+=["<c style='color:red'><b>STOP CODON CREATED</b></c>"]
				except:
					pass
			if len(delta2) == 0:
				delta2+=['No amino acid change predicted']	
			temp = "<c style='color:green'><b>Exon ", str(left_ind +1), "</b></c><br> (Change at boundary: ", delta2[0], ")"
			temp = ''.join(temp)
			skips += [temp]
		if (length - right_length) % 3 == 0 and right_ind != 78:
			delta = []
			delta2 = []				
			cdnaskip = refcdna[0:min(num_list) + 243] + refcdna[min(num_list) + length + right_length + 243:]
			cdnaskip = ''.join(cdnaskip)
			aa_skip = translate(cdnaskip[244:])
			i = 0
			filler = ''
			while i < (length + right_length)/3:
				filler += '-'
				i += 1
			end = min(num_list)
			if end%3 == 1:
				aa_skipgap = aa_skip[0:(end/3)] + filler + aa_skip[(end/3):]
			if end%3 == 2:
				aa_skipgap = aa_skip[0:(end/3)+1] + filler + aa_skip[(end/3)+1:]		  
			aa_skipgap = ''.join(aa_skipgap)
			delta, align, PTC = easy_align(aa_ref, aa_skipgap, mutype, frame_shift)
			for index, item in enumerate(delta):
				if type(item) != int and '-' not in item:
					if item != 'STOP':
						delta2+=[item]
				try:
					if item == 'STOP' and delta[index+1] != 3686:
						delta2+=["<c style='color:red'><b>STOP CODON CREATED</b></c>"]
				except:
					pass
			if len(delta2) == 0:
				delta2+=['No amino acid change predicted']	
			temp = "<c style='color:green'><b>Exon ", str(right_ind+1), "</b></c><br> (Change at boundary: ", delta2[0], ")"
			temp = ''.join(temp)
			skips += [temp]
		if (length - (right_length + left_length))%3 == 0 and left_ind != 0 and right_ind != 78:
			delta = []
			delta2 = []				
			cdnaskip = refcdna[0:(exon_positions[left_ind][1])+243] + refcdna[exon_positions[right_ind][2]+1+243:]
			cdnaskip = ''.join(cdnaskip)
			aa_skip = translate(cdnaskip[244:])
			i = 0
			filler = ''
			while i < (length + right_length + left_length)/3:
				filler += '-'
				i+=1
			end = min(num_list)
			if end%3 == 1:
				aa_skipgap = aa_skip[0:(end/3)] + filler + aa_skip[(end/3):]
			if end%3 == 2:
				aa_skipgap = aa_skip[0:(end/3)+1] + filler + aa_skip[(end/3)+1:]		  
			aa_skipgap = ''.join(aa_skipgap)
			delta, align, PTC = easy_align(aa_ref, aa_skipgap, mutype, frame_shift)
			for index, item in enumerate(delta):
				if type(item) != int and '-' not in item:
					if item != 'STOP':
						delta2+=[item]
				try:
					if item == 'STOP' and delta[index+1] != 3686:
						delta2+=["<c style='color:red'><b>STOP CODON CREATED</b></c>"]
				except:
					pass
			if len(delta2) == 0:
				delta2+=['No amino acid change predicted']		
			temp = "<c style='color:green'><b>Exons ", str(left_ind+1), ' and ', str(right_ind + 1), "</b></c><br> (Change at boundary: ", delta2[0], ")"
			skips += [''.join(temp)]
	if max(ex_ints) < 6:
		skips += ["""An internal ribosome entry site has been described downstream of this variant (<a target="_blank" href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2663021/">PubMed</a>)."""]
	return skips




def ese_find(cdna, refcdna, ese_list, ess_list, gen_pos, genomic_dna, length_mutation, mut_genomic_dna):
	"""
	Returns a list: rese, ress, esef, ese_types, ess_types
		rese: boolean
		ress: boolean
		esef: boolean
		ese_types: list, ex. ["IgM_BRCA1", str(mx), perc]
		ess_types: list, EMPTY at this point.

	Checks for ESEs via two methods:
	(1) RESCUE-ESE:
		The mutated sequence (along with 5 positions flanking on either side since most consensus sequences are 6 bp long) is checked against the RESCUE-ESE list before and after mutation simulation.
		If there is an ESE before mutation, but not afterwards, the RESCUE-ESE is considered broken.
	(2) ESEFinder:
		The mutated sequence +/- 5 bp goes through scoring via 5 different matrices provided by ESEFinder.
		If there is an ESE above threshold (determined by ESEFinder) before mutation, which drops below threshold after mutation, it is considered broken.
	A boolean variable corresponding to T-Broken, F-Uneffected is returned for: RESCUE-ESE and ESEFinder repectively.
	ESEFinder variable is a list including the boolean results, as well as the specific type of consensus sequence effected.
	"""
	mut_min = min(gen_pos) - 5
	mut_max = max(gen_pos) + 5
	ref_seg = genomic_dna[mut_min:mut_max]
	ress = False
	rese = False
	esef = False
	ese_types = []
	ess_types = []
	if length_mutation > 30:
		return rese, ress, esef, ese_types, ess_types
	mut_seg = mut_genomic_dna[mut_min:mut_max]
	count_ref = 0
	count_mut = 0
	for item in ese_list:
		if item in ref_seg:
			ind = ese_list.index(item)
			if ese_list[ind] not in mut_seg:
				rese = True
		count_ref +=1
	for item in ess_list:
		if item in mut_seg:
			ind = ess_list.index(item)
			if ess_list[ind] not in ref_seg and item != '':
				ress = True

	ref = sf2_asf(ref_seg)
	mut = sf2_asf(mut_seg)
	i=1
	change=[]
	percent=[]
	while i<len(ref) and i<len(mut):
		if ref[i]!=mut[i]:
			if ref[i]>=1.956 or mut[i]>=1.956:
				change+=[mut[i]-ref[i]]
				percent+=[float(((mut[i]-ref[i])/ref[i])*100)]
		i+=2
	if len(change)!=0:
		abs_change = []
		for item in change:
			abs_change += [abs(item)]
		mx = max(abs_change)
		try:
			index = change.index(mx)
		except:
			temp = mx - mx*2
			index = change.index(temp)
		perc = percent[index]
		if perc > 0:
			perc_sign = "+"
		if perc < 0:
			perc_sign = ""
		perc = str("%.2f" % perc)
		if perc_sign == "+":
			perc = "<c style='color:green'>",perc_sign, perc,"%</c>"
		if perc_sign == "":
			perc = "<c style='color:red'>",perc_sign,perc,"%</c>"
		perc = ''.join(perc)
		if mx not in change:
			mx = mx-2*mx
		if mx > 0:
			mx = '+', str(mx)
			mx = ''.join(mx)
		esef = True
		temp = ["SF2/ASF", str(mx), perc]
		ese_types += [temp]

	ref = igm_brca1(ref_seg)
	mut = igm_brca1(mut_seg)
	i=1
	change=[]
	percent=[]
	while i<len(ref) and i<len(mut):
		if ref[i]!=mut[i]:
			if ref[i]>=1.867 or mut[i]>=1.867:
				change+=[mut[i]-ref[i]]
				percent+=[float(((mut[i]-ref[i])/ref[i])*100)]
		i+=2
	if len(change)!=0:
		abs_change = []
		for item in change:
			abs_change += [abs(item)]
		mx = max(abs_change)
		try:
			index = change.index(mx)
		except:
			temp = mx - mx*2
			index = change.index(temp)
		perc = percent[index]
		if perc > 0:
			perc_sign = "+"
		if perc < 0:
			perc_sign = ""
		perc = str("%.2f" % perc)
		if perc_sign == "+":
			perc = "<c style='color:green'>",perc_sign, perc,"%</c>"
		if perc_sign == "":
			perc = "<c style='color:red'>",perc_sign,perc,"%</c>"
		perc = ''.join(perc)
		if mx not in change:
			mx = mx-2*mx
		if mx > 0:
			mx = '+', str(mx)
			mx = ''.join(mx)
		esef = True
		temp = ["IgM_BRCA1", str(mx), perc]
		ese_types += [temp]

	ref = sc35(ref_seg)
	mut = sc35(mut_seg)
	i=1
	change=[]
	percent=[]
	while i<len(ref) and i<len(mut):
		if ref[i]!=mut[i]:
			if ref[i]>=2.383 or mut[i]>=2.383:
				change+=[mut[i]-ref[i]]
				percent+=[float(((mut[i]-ref[i])/ref[i])*100)]
		i+=2
	if len(change)!=0:
		abs_change = []
		for item in change:
			abs_change += [abs(item)]
		mx = max(abs_change)
		try:
			index = change.index(mx)
		except:
			temp = mx - mx*2
			index = change.index(temp)
		perc = percent[index]
		if perc > 0:
			perc_sign = "+"
		if perc < 0:
			perc_sign = ""
		perc = str("%.2f" % perc)
		if perc_sign == "+":
			perc = "<c style='color:green'>",perc_sign, perc,"%</c>"
		if perc_sign == "":
			perc = "<c style='color:red'>",perc_sign,perc,"%</c>"
		perc = ''.join(perc)
		if mx not in change:
			mx = mx-2*mx
		if mx > 0:
			mx = '+', str(mx)
			mx = ''.join(mx)
		esef = True
		temp = ["SC35", str(mx), perc]
		ese_types += [temp]

	ref = srp40(ref_seg)
	mut = srp40(mut_seg)
	i=1
	change=[]
	percent=[]
	while i<len(ref) and i<len(mut):
		if ref[i]!=mut[i]:
			if ref[i]>=2.670 or mut[i]>=2.670:
				change+=[mut[i]-ref[i]]
				percent+=[float(((mut[i]-ref[i])/ref[i])*100)]
		i+=2
	if len(change)!=0:
		abs_change = []
		for item in change:
			abs_change += [abs(item)]
		mx = max(abs_change)
		try:
			index = change.index(mx)
		except:
			temp = mx - mx*2
			index = change.index(temp)
		perc = percent[index]
		if perc > 0:
			perc_sign = "+"
		if perc < 0:
			perc_sign = ""
		perc = str("%.2f" % perc)
		if perc_sign == "+":
			perc = "<c style='color:green'>",perc_sign, perc,"%</c>"
		if perc_sign == "":
			perc = "<c style='color:red'>",perc_sign,perc,"%</c>"
		perc = ''.join(perc)
		if mx not in change:
			mx = mx-2*mx
		if mx > 0:
			mx = '+', str(mx)
			mx = ''.join(mx)
		esef = True
		temp = ["SRp40", str(mx), perc]
		ese_types += [temp]

	ref = srp55(ref_seg)
	mut = srp55(mut_seg)
	i=1
	change=[]
	percent=[]
	while i<len(ref) and i<len(mut):
		if ref[i]!=mut[i]:
			if ref[i]>=2.676 or mut[i]>=2.676:
				change+=[mut[i]-ref[i]]
				percent+=[float(((mut[i]-ref[i])/ref[i])*100)]
		i+=2
	if len(change)!=0:
		abs_change = []
		for item in change:
			abs_change += [abs(item)]
		mx = max(abs_change)
		try:
			index = change.index(mx)
		except:
			temp = mx - mx*2
			index = change.index(temp)
		perc = percent[index]
		if perc > 0:
			perc_sign = "+"
		if perc < 0:
			perc_sign = ""
		perc = str("%.2f" % perc)
		if perc_sign == "+":
			perc = "<c style='color:green'>",perc_sign, perc,"%</c>"
		if perc_sign == "":
			perc = "<c style='color:red'>",perc_sign,perc,"%</c>"
		perc = ''.join(perc)
		if mx not in change:
			mx = mx-2*mx
		if mx > 0:
			mx = '+', str(mx)
			mx = ''.join(mx)
		esef = True	    	
		temp = ["SRp55", str(mx), perc]
		ese_types += [temp]
	return rese, ress, esef, ese_types, ess_types
			



def sf2_asf(frag):
	"""
	Returns a list of lists:
		positives: list of lists [test, score]
			test: string, DNA fragment
			score: float, matrix score

	ESEFinder positional scoring matrix.
	"""
	i = 0
	positives = []
	while i < len(frag) - 6:
		test = frag[i:i+7]
		score = 0
		scores = [[-1.14, 1.37, -0.21, -1.58], [0.62, -1.1, 0.17, -0.5],[-1.58,0.73,0.48, -1.58],[1.32,0.33,-1.58,-1.13],[-1.58,0.94,0.33,-1.58],[-1.58,-1.58,0.99,-1.13],[0.62,-1.58,-0.11,0.27]]
		position = 0
		while position < 7:
			if test[position] == "A":
				score += scores[position][0]
			if test[position] == "C":
				score += scores[position][1]
			if test[position] == "G":
				score += scores[position][2]
			if test[position] == "T":
				score += scores[position][3]
			position += 1
		positives += [test, score]
		i += 1
	return positives



def igm_brca1(frag):
	"""
	Returns a list of lists:
		positives: list of lists [test, score]
			test: string, DNA fragment
			score: float, matrix score

	ESEFinder positional scoring matrix.
	"""
	i = 0
	positives = []
	while i < len(frag) - 6:
		test = frag[i:i+7]
		score = 0
		scores = [[-1.58,1.55,-1.35,-1.55],[0.15,-0.53,0.44, -0.28],[-0.97,0.79,0.41,-1.28],[0.74,0.33,-0.98,-0.92],[-1.19,0.72,0.51,-1.09],[-0.75,-0.62,1.03,-0.52],[0.43,-0.99,0.00,0.20]]
		position = 0
		while position < 7:
			if test[position] == "A":
				score += scores[position][0]
			if test[position] == "C":
				score += scores[position][1]
			if test[position] == "G":
				score += scores[position][2]
			if test[position] == "T":
				score += scores[position][3]
			position += 1
		positives += [test, score]
		i+=1
	return positives




def sc35(frag):
	"""
	Returns a list of lists:
		positives: list of lists [test, score]
			test: string, DNA fragment
			score: float, matrix score

	ESEFinder positional scoring matrix.
	"""
	i = 0
	positives = []
	while i < len(frag) - 7:
		test = frag[i:i+8]
		score = 0
		scores = [[-0.88,-1.16,0.87,-1.18],[0.09,-1.58,0.45,-0.2],[-0.06,0.95,-1.36,0.38],[-1.58,1.11,-1.58,0.88],[0.09,0.56,-0.33,-0.2],[-0.41,0.86,-0.05,-0.86],[-0.06,0.32,-1.36,0.96],[0.23,-1.58,0.68,-1.58]]
		position = 0
		while position < 7:
			if test[position] == "A":
				score += scores[position][0]
			if test[position] == "C":
				score += scores[position][1]
			if test[position] == "G":
				score += scores[position][2]
			if test[position] == "T":
				score += scores[position][3]
			position += 1
		positives += [test, score]
		i+=1
	return positives



def srp40(frag):
	"""
	Returns a list of lists:
		positives: list of lists [test, score]
			test: string, DNA fragment
			score: float, matrix score

	ESEFinder positional scoring matrix.
	"""
	i = 0
	positives = []
	while i < len(frag) - 6:
		test = frag[i:i+7]
		score = 0
		scores = [[-0.13,0.56,-1.58,0.92],[-1.58,0.68,-0.14,0.37],[1.28,-1.12,-1.33,0.23],[-0.33,1.24,-0.48,-1.14],[0.97,-0.77,-1.58,0.72],[-0.13,0.13,0.44,-1.58],[-1.58,-0.05,0.8,-1.58]]
		position = 0
		while position < 7:
			if test[position] == "A":
				score += scores[position][0]
			if test[position] == "C":
				score += scores[position][1]
			if test[position] == "G":
				score += scores[position][2]
			if test[position] == "T":
				score += scores[position][3]
			position += 1
		positives += [test, score]
		i+=1
	return positives




def srp55(frag):
	"""
	Returns a list of lists:
		positives: list of lists [test, score]
			test: string, DNA fragment
			score: float, matrix score

	ESEFinder positional scoring matrix.
	"""
	i = 0
	positives = []
	while i < len(frag) - 5:
		test = frag[i:i+6]
		score = 0
		scores = [[-0.66,0.39,-1.58,1.22],[0.11,-1.58,0.72,-1.58],[-0.66,1.48,-1.58,-0.07],[0.11,-1.58,0.72,-1.58],[-1.58,-1.58,0.21,1.02],[0.61,0.98,-0.79,-1.58]]
		position = 0
		while position < 6:
			if test[position] == "A":
				score += scores[position][0]
			if test[position] == "C":
				score += scores[position][1]
			if test[position] == "G":
				score += scores[position][2]
			if test[position] == "T":
				score += scores[position][3]
			position += 1
		positives += [test, score]
		i+=1
	return positives




def aa_longform(aa, mutype, frame_shift, num_list, length_mutation, aa_ref):
	"""
	Returns a list (len=variable):
		aa_lchange: list of strings, changes ready to be printed individually (each item is a formatted string).

	Converts old notation (single-letter amino acid) to new notation (three-letter amino acid).
	Adds detail to nonsense mutations (the different accepted notations; Ter and *)
	Creates HGVS for amino acid deletions (i.e in-frame deletions)
	"""
	lst = ['A','Ala','C','Cys','D','Asp','E','Glu','F','Phe','G','Gly','H','His','I','Ile','K','Lys','L','Leu','M','Met','N','Asn','P','Pro','Q','Gln','R','Arg','S','Ser','T','Thr','V','Val','W','Trp','Y','Tyr','*','*']	                                                                                         
	if mutype == "Deletion" and not frame_shift:
		first_position = min(num_list)/3
		if min(num_list)%3==2:
			first_position += 1
		last_position = first_position + length_mutation/3
		first_aa = aa_ref[first_position]
		last_aa = aa_ref[last_position]
		first_aa_index = lst.index(first_aa)
		last_aa_index = lst.index(last_aa)
		first_aa_long = lst[first_aa_index + 1]
		last_aa_long = lst[last_aa_index + 1]
		if length_mutation > 3:
			temp = "p.(",first_aa_long,str(first_position),"_",last_aa_long,str(last_position),"del)"
		if length_mutation <= 3:
			temp = "p.(",first_aa_long,str(first_position),"del)"
		temp = ''.join(temp)
		return temp
	if len(aa) == 0:
		return ["Insertion/Duplication without frameshift"]
	if len(aa) == 1 and re.findall("\d+", aa[0]) == []:
		temp_ind = lst.index(aa[0])
		aa_lchange = lst[temp_ind+1]
		return aa_lchange
	if aa[0] == 'STOP':
		if int(aa[1])!= 3686:
			aa1 = aa_ref[int(aa[1])-1]
			temp_ind = lst.index(aa1)
			aa1 = lst[temp_ind+1]
			location = aa[1]
			temp = "<b>p.",aa1, str(location), "*</b><br>\nAlso known as p.", aa1, str(location), "Ter (", "<i>*Note that the p.", aa1, str(location),"X terminology is outdated*</i>)"
			temp =''.join(temp)
			return temp
		if int(aa[1]) == 3686:
			silent = "(=)"
			return silent
	nums = []
	sym_list = []
	for item in aa:
		if item != "STOP" and type(item)!= int and '-' not in item:
			num = re.findall("\d+", item)
			for char in num:
				nums += [char]
			ind1 = lst.index(item[0])
			ind2 = lst.index(item[len(item)-1])
			sym_list += [lst[ind1+1]]
			sym_list += [lst[ind2+1]]
	aa_lchange = []
	i = 0
	n = 0
	while i<len(sym_list)-1 and n<len(num_list):
		temp = sym_list[i],nums[n],sym_list[i+1]
		temp = ''.join(temp)
		aa_lchange += [temp]
		i+=2
		n+=1
	return aa_lchange




def drna(dna):
	"""
	Returns a single string:
		rna: string, a sequence with RNA symbols.
			If the sequence length is 3 and matches a stop codon, the codon, along with its name will be returned as a single string.
	
	Whenever a T is found in the input string (DNA), a U is added to the new sequence (RNA).
	"""
	rna = ''
	for char in dna:
		if char == "T":
			rna += 'U'
		else:
			rna += char
	if rna == 'UAG':
		rna = rna, ', Amber'
	if rna == 'UAA':
		rna = rna, ', Ochre'
	if rna == 'UGA':
		rna = rna, ', Opal/Umber'
	rna = ''.join(rna)
	return rna




def domain(num_list, domains):
	"""
	Returns a list (len=variable) of strings:
		dmns: list of strings, each domain directly effected by the mutation.

	First checks the domain positions against the amino acid positions taken by (cDNA position/3).
	Then the list of domains is built and check for redundancy (each domain should be shown once, followed by the subdomains within).
	"""
	temp_indexes = []
	temp_domains = []
	check_repeat = []
	indexes = []
	for item in num_list:
		aa_item = item/3
		i = 2
		while i<len(domains):
			if aa_item == domains[i] or aa_item == domains[i+1]:
				temp_indexes += [i]
				i+=6
			if aa_item >= domains[i] and aa_item <= domains[i+1]:
				temp_indexes += [i]
				i+=6
			else:
				i+=6
	i = 0
	if len(temp_indexes) ==0:
		return dmns
	while i<len(temp_indexes):
		if temp_indexes[i] not in indexes:
			indexes += [temp_indexes[i]]
		i+=1
	first_domain_index = min(indexes)
	last_domain_index = max(indexes)
	if last_domain_index != first_domain_index and last_domain_index-first_domain_index>6:
		number_of_domains = (last_domain_index-first_domain_index)/6
		indexes = [first_domain_index]
		i = 0
		counter = first_domain_index
		while i<number_of_domains:
			counter += 6
			indexes += [counter]
			i+=1
	for item in indexes:
		check_repeat = domains[item-2], domains[item-1]
		temp_domains += [check_repeat]
	i = 0
	temp = ''
	while i<len(temp_domains):
		tracker = 0
		if temp_domains[i][1] == '':
			if i == len(temp_domains)-1 or temp_domains[i][0] != temp_domains[i+1][0]:
				if temp_domains[i][0] not in temp:
					temp += str(temp_domains[i][0])
					temp +=  ".<br>\n"
			tracker += 1
		if temp_domains[i][1] != '' and temp_domains[i][0]:
			if temp_domains[i][0] not in temp:
				temp += temp_domains[i][0]
				temp += ': '
				tracker += 1
				temp += temp_domains[i][1]
			j=i+1
			while j < len(temp_domains):
				if temp_domains[i][0] == temp_domains[j][0]:
					if temp_domains[j][1] != '':
						temp += ", "
						tracker +=1							
						temp += temp_domains[j][1]
				j+=1 
		if tracker > 1:
			i+=tracker
		if tracker < 2:
			i+=1
	dmns = temp
	return dmns



def complement(a):
	"""
	Returns a single string: 
		result: string, the complement of the entered DNA sequence.
	"""
	result = ''
	for char in a:
		if char == "A":
			result += "T"
		if char == "T":
			result += "A"
		if char == "G":
			result += "C"
		if char == "C":
			result += "G"
	return result



def output_statements(mut, mutype, length, aa_change, exon_numbers, part, frame_shift, long_aa_change, aa_seq, aa_ref, mut_cdna):
	"""
	Returns 3 boolean, followed by 3 strings:
		silent: boolean, silent mutation?
		missense: boolean, missense mutation?
		nonsense: boolean, nonsense mutation?
		statement1: string, generally the primary, verbal description of the mutation.
		statement2: string, generally a verbal description of PTC/nonsense.
		statement3: string, only applicable for frameshift -- the HGVS notation for protein level change.

	Compilation of all of the possible output sequences.
	Difference combinations of output statements are chosen based on alignment results, mutation type, length etc.
	Mutation attributes (missense, nonsense, silent) [boolean] and output statements [string] are returned.
	"""
	statement1 = ''
	statement2 = ''
	statement3 = ''
	silent = False
	missense = False
	nonsense = False
	if mutype == "Deletion":
		if len(exon_numbers)>0:
			end_del = max(exon_numbers)
			start_del = min(exon_numbers)
			temp = 'empty'
			if len(exon_numbers)>1 or part == False:
				temp = str(start_del), " to ", str(end_del)
			temp = ''.join(temp)
			statement1 = temp
		if "STOP" in aa_change:
			stop_index = aa_change.index("STOP")
			stop_position = aa_change[stop_index + 1]
			if stop_position != 3761 - 75 and len(aa_change) == 2:
				nonsense = True
				nonsense_position = stop_position
				temp = "There is a premature stop codon (PTC) predicted at residue (no treatment) ", str(nonsense_position), ": ", long_aa_change,"\n<br><br>Stop codon created: (",drna(mut_cdna[(stop_position*3)-3+244:(stop_position*3)+244]),")"
				statement2 = ''.join(temp)
		if len(aa_change) > 2:
			missense = True
	if mutype == "Duplication":
		if len(exon_numbers)>0:
			end_del = max(exon_numbers)
			start_del = min(exon_numbers)
			temp = "The duplication is from exon ", str(start_del), " to ", str(end_del), "."
			temp = ''.join(temp)
			statement1 = temp
		if "STOP" in aa_change:
			stop_index = aa_change.index("STOP")
			stop_position = aa_change[stop_index + 1]
			if stop_position < 3761 - 75 and len(aa_change) == 2:
				nonsense = True
				nonsense_position = stop_position
				temp = "There is a premature stop codon (PTC) predicted at residue ", str(nonsense_position), " in exon ", str(exon_numbers[0]), ": ", long_aa_change, "\n<br><br>Stop codon created: (", drna(mut_cdna[(stop_position*3)-3+244:(stop_position*3)+244]),")"
				statement2 = ''.join(temp)
		if len(aa_change) > 2:
			missense = True
	if mutype == ("Point mutation") or mutype == ("Deletion/Insertion") or mutype == "Insertion":
		if "STOP" in aa_change:
			stop_index = aa_change.index("STOP")
			stop_position = aa_change[stop_index + 1]
			if stop_position < 3686 and len(aa_change) == 2 and missense == False and not frame_shift:
				nonsense = True
				nonsense_position = stop_position
				temp = long_aa_change, "\n<br><br>Stop codon created: (",drna(mut_cdna[(stop_position*3)-3+244:(stop_position*3)+244]),")" #"There is a premature stop codon (PTC) predicted at residue ", str(nonsense_position), " in exon ", str(exon_numbers[0]),": ", 
				statement1 = ''.join(temp)
			elif stop_position == 3761 - 75 and len(aa_change) == 2: #DOUBLE CHECK THIS
				silent = True
				temp = "The ", mut, " mutation is predicted to be a silent mutation in exon ", str(exon_numbers[0]), ": p.(=)"
				statement1 = ''.join(temp)
		if len(aa_change) == 1 and frame_shift == False:
			missense = True
			missense_change = aa_change[0]
			if length < 50 or mutype == "Deletion/Insertion":
				chngs = ''
				for item in long_aa_change:
					if item == 'STOP':
						break
					else:
						chng = "p.(",str(item),")"
						chng = ''.join(chng)
						chngs += chng
				temp =  chngs
				statement1 = ''.join(temp)
		   
	if frame_shift and not nonsense:
		change = []
		for item in long_aa_change:
		   if type(item) == str:
				if '-' not in item:
					change += [item]
		if len(change)>0 and len(aa_seq) < len(aa_ref):
			change1 = change[0]
			start = re.findall("\d+", change1)
			start = int(start[0])
			end = int(aa_change[stop_index +1])
			shift_stop = (end - start) + 1
			temp = "p.(", change1, "fs*", str(shift_stop), ")"
			statement3 = ''.join(temp)
		# if len(change)>0 and len(aa_seq) > len(aa_ref):
		# 	statement3 = PTC
	return silent, missense, nonsense, statement1, statement2, statement3



def myVariantSearch(hg):
	import myvariant
	mv = myvariant.MyVariantInfo()
	gene_symbol = "DMD"
	query = "'%s' AND '%s'" % (hg, gene_symbol)
	tempres = mv.query(query)
	mvres = []
	for item in tempres['hits']:
		if hg in str(item):
			mvres += [item]
		else:
			pass
	return mvres


def testcv(entry):
	try:
		cs = entry['clinical_significance']
	except:
		cs = "Not available"
	try:
		le = entry['last_evaluated']
	except:
		le = "Not available"
	try:
		ac = entry['accession']
	except:
		ac = "Not available"
	return cs, le, ac


def ClinVar(mvr, hgvs):
	"""
	Returns a single string:
		result: string, the classification of the variant in clinvar and the referenced file's date of release.

	Uses the HGVS notation of mutation to search a static clinvar export.
	"""
	hits = 0
	cv = ''
	for item in mvr:
		try:
			cv = item['clinvar']
		except:
			pass
	# if hits != 0:
	# 	cv_table = "Number of submitters: " + str(hits)
	# 	cv_table += """
	# 	<table style="width:100%">
	# 		<tr>
	# 			<th>Clinical Significance</th>
	# 			<th>Last Evaluated</th>
	# 			<th>Link to accession</th>
	# 		</tr>
	# 		"""
	# 	count = 0
	# 	while count<hits:
	# 		row = ''
	# 		if hits > 1:
	# 			try:
	# 				row = "<tr><td>%s</td><td>%s</td><td><a href='https://www.ncbi.nlm.nih.gov/clinvar/%s/' target='_blank'>Click</a></td></tr>" % (cv['hits'][count]['clinvar']['rcv'][0]['clinical_significance'],cv['hits'][count]['clinvar']['rcv'][0]['last_evaluated'], cv['hits'][count]['clinvar']['rcv'][0]['accession'])
	# 			except:
	# 				pass
	# 		if hits == 1:
	# 			row = "<tr><td>%s</td><td>%s</td><td><a href='https://www.ncbi.nlm.nih.gov/clinvar/%s/' target='_blank'>Click</a></td></tr>" % (cv['hits'][index]['clinvar']['rcv'][0]['clinical_significance'],cv['hits'][index]['clinvar']['rcv'][0]['last_evaluated'], cv['hits'][index]['clinvar']['rcv'][0]['accession'])				
	# 		cv_table += row
	# 		count += 1
	if len(cv)>0:
		cv_table = """
		<table style="width:100%">
			<tr>
				<th>Clinical Significance</th>
				<th>Last Evaluated</th>
				<th>Link to accession</th>
			</tr>
			"""
		if type(cv['rcv']) == list:
			for item in cv['rcv']:
				clinsig, last_eval, access = testcv(item)
				row = "<tr><td>%s</td><td>%s</td><td><a href='https://www.ncbi.nlm.nih.gov/clinvar/%s/' target='_blank'>Click</a></td></tr>" % (clinsig, last_eval, access)				
				cv_table += row
				hits += 1
		else:
			clinsig, last_eval, access = testcv(cv['rcv'])
			row = "<tr><td>%s</td><td>%s</td><td><a href='https://www.ncbi.nlm.nih.gov/clinvar/%s/' target='_blank'>Click</a></td></tr>" % (clinsig, last_eval, access)
			cv_table += row
			hits += 1
		cv_table += "</table>"
		cv_table = "Number of submitters: " + str(hits) + cv_table
	else:
		cv_table = "Variant not found in ClinVar via myVariant.info, however, please <a href='https://www.ncbi.nlm.nih.gov/clinvar/?term=%s' target='_blank'> click here </a>to search ClinVar directly" % hgvs
	return cv_table


def map_c2g(clist, exStart, exEnd, cStart, cEnd, strnd, chrom):
    """
    does this provide help?
    """
    #cnums for only numbers from cHGVS
    #other for everything after numbers in hgvs
    cnums = []
    other = []
    for item in clist:
        temp = re.findall("\d+", item)
        if len(temp) != 0:
            cnums += [temp[0]]
            ind = item.index(temp[len(temp)-1])
            therest = item[ind + len(temp[len(temp)-1]):]
            other += [therest]
        else:
            cnums += [-1]
    #gnums for only numbers for genomic
    gnums = []
    cum = []
    cumCount = 0
    cStart = int(cStart)
    cEnd = int(cEnd)
    if strnd == "-":
        temp = str(exStart)
        temp2 = str(exEnd)
        exStart = (temp2.split(","))[::-1]
        exEnd = (temp.split(","))[::-1]
        temp = cStart
        temp2 = cEnd
        cStart = temp2
        cEnd = temp
    else:
        exStart = (str(exStart).split(","))
        exEnd = (str(exEnd).split(","))
    while '' in exStart:
        exStart.remove('')
    while '' in exEnd:
        exEnd.remove('')
    i = 0
    while i<len(exStart) and i <len(exEnd):
        exStart[i] = int(exStart[i])
        exEnd[i] = int(exEnd[i])
        i += 1
    lengths = []
    if strnd == "+":
        i = 0
        if exEnd[0] > cStart: 
            lengths += [abs(exEnd[i] - cStart)]
            i += 1
        else:
            while exEnd[0] < cStart:
                exEnd.remove(exEnd[0])
                exStart.remove(exStart[0])
        while i <len(exEnd):
            lengths += [abs(exEnd[i] - exStart[i])]
            i += 1
    if strnd == "-":
        i = 0
        if exEnd[0] < cStart:
            lengths += [abs(exEnd[i] - cStart)]
            i += 1
        else:
            while exEnd[0] > cStart:
                exEnd.remove(exEnd[0])
                exStart.remove(exStart[0])
        while i <len(exEnd):
            lengths += [abs(exEnd[i] - exStart[i])]
            i += 1
    for item in lengths:
        cum += [cumCount]
        cumCount += item
    for item in cnums:
        found = False
##        if type(item) == list:
##            item = item[0]
        item = int(item)
        i = 0
        while i<len(cum)-1 and not found:
            if item > cum[i] and item < cum[i+1]:
                diff = abs(item-cum[i])
                if strnd == "-":
                    temp = exStart[i] - diff + 1
                if strnd == "+":
                    temp = exStart[i] + diff
                gnums += [temp]
                found = True
            i += 1
        if not found:
            if item > cum[len(cum)-1]:
                diff = abs(item-cum[len(cum)-1])
                if strnd == "-":
                    temp = exStart[i] - diff + 1
                if strnd == "+":
                    temp = exStart[i] + diff
                gnums += [temp]
                found = True
            else:
                gnums += [-1]
    gnote = []
    gother = []
    for item in other:
        gother += [rev_comp(item, False, False)]
    i = 0
    while i<len(gnums) and i <len(gother):
        gnote += ["%s:g.%s%s" % (chrom, gnums[i], gother[i])]
        i += 1
    return gnums, gnote

def get_exons(inf, sequence):
	starts = inf[9]
	stops = inf[10]
	codingStart = inf[6]
	codingStop = inf[7]
	strnd = inf[3]
	tStart = inf[4]
	tEnd =inf[5]
	temp = []
	temp2 = []
	starts = starts.split(",")
	stops = stops.split(",")
	for item in starts:
		if item != '':
			temp += [int(item)]
	for item in stops:
		if item != '':
			temp2 += [int(item)]
	if strnd == '-':
		starts = temp2[::-1]
		stops = temp[::-1]
		temp3 = codingStart
		temp4 = codingStop
		codingStart = temp4
		codingStop = temp3
		temp4 = tStart
		temp3 = tEnd
		tEnd = temp4
		tStart = temp3
	else:
		starts = temp
		stops = temp2
	codingStop = int(codingStop)
	codingStart = int(codingStart)
	codingSeq = sequence[abs(codingStart-tStart):len(sequence)-abs(tEnd-codingStop)]
	if strnd == "+":
		codingSeq = codingSeq[1:]
	lengths = []
	exons = []
	exons_styled = []
	i = 0
	length_cum = 0
	while stops[0] < codingStart and strnd == "+":
		lengths += [abs(int(starts[0]) - int(stops[0]))]
		starts.remove(starts[0])
		stops.remove(stops[0])
	while i<len(starts) and i<len(stops):
		if starts[i] < codingStart and strnd == "+":
			starts[i] = starts[i] + abs(codingStart - starts[i])
		if i<len(starts)-1:
			length_gap = abs(stops[i] - starts[i+1])
		temp_length = abs(stops[i]-starts[i])
		if strnd == "-" and i == 0:
			temp_length -= abs(codingStart-tStart)
		if temp_length < 0:
			temp_length = abs(stops[i]-starts[i])
 		lengths += [temp_length]
 		temp_start = length_cum
 		temp_end = length_cum + temp_length
 		# if strnd == "-":
 		# 	temp_end +=1
 		if i<len(starts)-1:
			exons += [codingSeq[temp_start: temp_end]]
			exons_styled += ["<c class='exon' id='exon%s' style='color:%s'>%s</c>" % (str(i+1),'gray', codingSeq[temp_start: temp_end])]
		else:
			exons += [codingSeq[temp_start: ]]
			exons_styled += ["<c class='exon' id='exon%s' style='color:%s'>%s</c>" % (str(i+1),'gray', codingSeq[temp_start: ])]
		length_cum += (temp_length + length_gap)
		i += 1
	total_length = abs(sum(lengths))
	exons_styled = ''.join(exons_styled)
	return lengths, exons, exons_styled