from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect, HttpResponse
from django.template import RequestContext, loader
from django.urls import reverse
from .forms import IndexForm, ACMGForm

from .interpreter_functions import *
import os,sys
import subprocess
import csv
import re
import pprint
from math import floor
from urllib.parse import quote

interpreter_dir = os.path.dirname(__file__)

# Create your views here.
def index(request):
	if request.method == "POST":
		user = request.POST.get('user', None)
	return render(request, 'index.html')


def results(request):
	#Get the user input from HTML form
	if request.method == 'POST':
		mut = request.POST.get('mutation', None)

	ucsc_info = [12, 'NM_004006', 'chrX', '-', 31137344, 33229673, 31140035, 33229429, 79, '31137344,31144758,31152218,31164407,31165391,31187559,31190464,31191655,31196048,31196785,31198486,31200854,31222077,31224698,31227614,31241163,31279071,31341714,31366672,31462597,31496222,31497099,31514904,31525397,31645789,31676106,31697491,31747747,31792076,31838091,31854834,31893307,31947712,31950196,31986455,32235032,32305645,32328198,32360216,32361250,32364059,32366522,32380904,32382698,32383136,32398626,32404426,32407617,32408187,32429868,32456357,32459296,32466572,32472778,32481555,32482702,32486614,32490280,32503035,32509393,32519871,32536124,32563275,32583818,32591646,32591861,32613873,32632419,32662248,32663080,32715986,32717228,32827609,32834584,32841411,32862899,32867844,33038255,33229398,', '31140047,31144790,31152311,31164531,31165635,31187718,31190530,31191721,31196087,31196922,31198598,31201021,31222235,31224784,31227816,31241238,31279133,31341775,31366751,31462744,31496491,31497220,31515061,31525570,31645979,31676261,31697703,31747865,31792309,31838200,31854939,31893490,31947862,31950344,31986631,32235180,32305818,32328393,32360399,32361403,32364197,32366645,32381075,32382827,32383316,32398797,32404582,32407791,32408298,32430030,32456507,32459431,32466755,32472949,32481711,32482816,32486827,32490426,32503216,32509635,32519959,32536248,32563451,32583998,32591754,32591963,32613993,32632570,32662430,32663269,32716115,32717410,32827728,32834757,32841504,32862977,32867937,33038317,33229673,', 0, 'DMD', 'cmpl', 'cmpl', '0,1,1,0,2,2,2,2,2,0,2,0,1,2,1,1,2,1,0,0,1,0,2,0,2,0,1,0,1,0,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,2,0,0,0,0,0,2,0,0,0,1,2,0,0,0,0,1,0,']


	#Read in static reference files

	#Rescue ESE sites http://genes.mit.edu/burgelab/rescue-ese/ESE.txt
	f = open(os.path.join(interpreter_dir, "rescue_esesites.txt"))
	ese_sites = []
	for line in f:
	    temp = ''
	    for char in line:
	        if char != "\n":
	            temp += char
	        temp = ''.join(temp)
	    ese_sites += [temp]
	f.close()

	#Fas-ESS hexamers *citation needed*
	f = open(os.path.join(interpreter_dir, "esssites.txt"))
	ess_sites = []
	for line in f:
	    temp = ''
	    for char in line:
	        if char != "\n" and char not in "1234567890E-.":
	            if char != " " and char != ".":
	                temp += char
	        temp = ''.join(temp)
	    ess_sites += [temp]
	f.close()

	#Coding region location (c.###) to exon number location conversion
	with open(os.path.join(interpreter_dir, "exon_positions.csv"), 'r', encoding='utf-8') as csvfile:
	    spamreader = csv.reader(csvfile)
	    temp = []
	    for row in spamreader:
	            temp += [row]
	    i = 0
	    exon_positions = []
	    while i<len(temp):
	            temp2 = []
	            if "#" in temp[i][0]:
	                    i+=1
	            else:
	                    for char in temp[i]:
	                            if char!='':
	                                    temp2 += [int(char)]
	                    exon_positions += [temp2]
	                    i+=1

	#Domain positions
	with open(os.path.join(interpreter_dir, "domains.txt"),'r') as csvfile:
		spamreader = csv.reader(csvfile)
		temp = []
		for row in spamreader:
			temp += [row]
		domains = []
		for item in temp:
			i=0
			while i<len(item):
				if i==0 or i==1:
					domains+=[item[i]]
					i+=1
				else:
					domains+=[int(item[i])]
					i+=1


	#Exon positions in genomic DNA (GRCh38, hg38)
	with open(os.path.join(interpreter_dir, "exons_in_genomic_old.txt"), 'r') as csvfile:
		spamreader = csv.reader(csvfile)
		temp = []
		for row in spamreader:
			temp += [row]
		genomic_positions1 = []
		for item in temp:
			intermediate = []
			for subitem in item:
				intermediate+=[int(subitem)]
			genomic_positions1 += [intermediate]

	with open(os.path.join(interpreter_dir, "exons_in_genomic.txt"), 'r') as csvfile:
		spamreader = csv.reader(csvfile)
		temp = []
		for row in spamreader:
			temp += [row]
		genomic_positions2 = []
		for item in temp:
			intermediate = []
			for subitem in item:
				intermediate+=[int(subitem)]
			genomic_positions2 += [intermediate]

	with open("./interpreter/exons_in_genomic_hg19.txt", 'r') as csvfile:
		spamreader = csv.reader(csvfile)
		temp = []
		for row in spamreader:
			temp += [row]
		genomic_positions3 = []
		for item in temp:
			intermediate = []
			for subitem in item:
				intermediate+=[int(subitem)]
			genomic_positions3 += [intermediate]

	#Reference genomic DNA (NG_12232.1)
	f = open(os.path.join(interpreter_dir, "NG_012232.1.txt"))
	genomic_dna = ''
	for char in f:
		if char != '\n' and char != '\t':
			genomic_dna += char
	f.close()

	#Reference coding sequence (NM_004006.2)
	f = open(os.path.join(interpreter_dir, "DMD_cDNA.txt"))
	reference_cdna = ''
	for line in f:
	    reference_cdna += line
	f.close()

	#dbNSFP
	NSFP_headers = ['Genomic Position','Negative strand (i.e. NG_012232.1) position','WT Nucleotide','Mutant Nucleotide','SIFT', 'PolyPhen2','LRT','MutationTaster','MutationAssessor','FATHMM','PROVEAN','MetaSVM']
	NSFP = []
	with open(os.path.join(interpreter_dir, "filtered_X.txt")) as csvfile:
		spamreader = csv.reader(csvfile, delimiter=',')
		for row in spamreader:
			negative_strand = (int(row[1]) - 33344609)*(-1)
			NSFP += [[row[1],negative_strand,row[2],row[3],row[25],row[31],row[37],row[41],row[48],row[51],row[54],row[57]]]

	ex_lengths, ex_seqs, ex_styled = get_exons(ucsc_info, genomic_dna)


	#Initialize variables to defaults
	posskip = [""]
	boundary_change = []
	past_max = False
	ex_input = False
	part = False
	large_deletion = False
	triage = False
	path_pred = 0
	len_ref = len(reference_cdna)
	sift = 'Not missense'
	polyphen = 'Not missense'
	lrt = 'Not missense'
	mutationtaster = 'Not missense'
	mutationassessor = 'Not missense'
	fathmm = 'Not missense'
	provean = 'Not missense'
	metasvm = 'Not missense'
	nsfp_score = ''
	nsfp_message = ''
	cv = ''
	consequence = ''
	consequence_statement = ''
	ese_message = ''
	esefinder = ''
	rescueese = ''
	ess_message = ''
	splice_message = ''
	splice_message2 = ''
	motif_message = ''
	ds = ''
	splice_result = []
	long_aa_change = []
	readthrough_elig = "<c style='color:red'><b>Not eligible</b></c>"

	#Look for numbers and key phrases in user input.
	ex_input = exoninput(mut)
	mutype, search = type_inp(mut, ex_input)
	num_list, intron_list = parse(mut, mutype, exon_positions, ex_input)


	#If no numbers, or key phrases are detected in user input, redirect to input page with error message.
	if mutype == "Invalid mutation" and len(num_list) == 0:
		return render(request, 'index.html',{'success':False})
	else:
		pass


	#Check if numbers extend past the length of the coding region.
	for num in num_list:
		if int(num) > 3685*3 and not ex_input: #AA length *3
			past_max = True
		elif int(num) > len_ref:
			past_max = True
		else:
			pass


	#If there is only one cDNA position and the intron list is not empty, mutation type is intron_only OR if the num_list is empty
	intron_only = intron_count(intron_list, num_list)


	#Check to see if the user meant exon if the numbers entered were below 80 (the number of exons in DMD = 79).
	if (mutype == "Deletion" or mutype == "Duplication") and "c." not in mut and len(num_list) > 0 and not ex_input and "Nucleotide" not in mut:
		if max(num_list) < 80 and not ex_input:
			with open(os.path.join(interpreter_dir, 'templates/triage.html'), 'w') as g:
				num_list = re.findall("\d+", mut)
				g.write("""<div class='container' style='text-align:center'><h1>You entered: '<c style='color:red'>""")
				g.write(mut)
				g.write("""</c>'</h1><br><h2>Did you mean...<br><br>'<i>""")
				temp = ''
				temp += mutype
				if len(num_list) == 1:
					temp += " of exon "
					temp += str(num_list[0])
				if len(num_list) > 1:
					temp += " of exons "
					temp += str(min(num_list))
					temp += "-"
					temp += str(max(num_list))
				g.write(temp)
				g.write("""</i>'?</h2><br><div class='row'><form action="{% url 'interpreter:results' %}" method="post"> {% csrf_token %}{% csrf_token %}	<button name='mutation' id='mutation' value=' """)
				g.write(temp)
				g.write("""' class='btn-success btn-lg'>Yes, I meant exon.<br></button><button name='mutation' id='mutation' value='(Nucleotide) """)
				g.write(mut)
				g.write("""'class='btn-danger btn-lg'>No, I meant nucleotide.<br></button></form></div></div></body>""")
			g.closed
			return render(request, 'index.html', {'triage':True})


#Triage other common input errors (3 cases):

	if past_max or mutype == "Invalid mutation":
		with open(os.path.join(interpreter_dir, 'templates/triage.html'), 'w') as g:
			g.write("""<div class='container' style='text-align:center'><h1>You entered: '<c style='color:red'>""")
			g.write(mut)
			g.write("""</c>'<br><br>Did you mean...<br>""")

	#Numbers past the length of the cDNA:
			if past_max:
				g.write(("""</h1><h2>A smaller number?</h2><br>You entered: <c style='color:red'>%s""" % str(max(num_list))) + """</c><br>Which is larger than the length of the <i>DMD, Dp427m</i> coding sequence (11055 nt)<br><br>Please enter your revised mutation below, or <a href="{% url 'interpreter:index' %}">return to the previous page</a> to try the "Exon deletion selector" feature.<br><br><form action="{% url 'interpreter:results' %}" method="post"> {% csrf_token %}<input type = "text" id='userAnswer' name = "mutation" id='mutation' style='height:70px;width:100%;font-size:24pt;'/><br><br><button type = "submit" style = 'height:70px;width:150px;font-size:24px;align:center' class = 'btn btn-lg btn-primary'>Interpret</button></form></div></body>""")

	#Numbers, but no key words:
			if mutype == "Invalid mutation" and not past_max and len(num_list)>0:# and not ex_input:
				tridel = "Deletion of Nucleotide(s) "
				tridelex = "Deletion of Exon(s) "
				tridup = "Duplication of Nucleotide(s) "
				tridupex = "Duplication of Exon(s) "
				triins = "Insertion at Nucleotide position(s) "
				tridelins = "Deletion/Insertion at Nucleotide position(s) "
				tripoint = "Point mutation at Nucleotide position(s) "
				if len(num_list) == 1:
					tridel += str(num_list[0])
					tridup += str(num_list[0])
					tridelex += str(num_list[0])
					tridupex += str(num_list[0])
					triins += str(num_list[0])
					tridelins += str(num_list[0])
					tripoint += str(num_list[0])
				if len(num_list) > 1:
					temp1 = str(min(num_list)),'_',str(max(num_list))
					temp1 = ''.join(temp1)
					temp2 = str(min(num_list)),'-',str(max(num_list))
					temp2 = ''.join(temp2)
					tridup += temp1
					tridel += temp1
					triins += temp1
					tridelins += temp1
					tripoint += temp1
					tridelex += temp2
					tridupex += temp2

				g.write("""</b></h1><br><form action="{% url 'interpreter:results' %}" method="post"> {% csrf_token %}<div class='row'><button name='mutation' id='mutation' value='""")
				g.write(tridel)
				g.write("' class='btn-success btn-lg'>")
				g.write(tridel)

				if max(num_list) < 80:
					g.write("</button><button name='mutation' id='mutation' value='")
					g.write(tridelex)
					g.write("' class='btn-info btn-lg'>")
					g.write(tridelex)

				g.write("</div><div class='row'>")

				g.write("</button><button name='mutation' id='mutation' value='")
				g.write(tridup)
				g.write("' class='btn-success btn-lg'>")
				g.write(tridup)

				if max(num_list) < 80:
					g.write("</button><button name='mutation' id='mutation' value='")
					g.write(tridupex)
					g.write("' class='btn-info btn-lg'>")
					g.write(tridupex)

				g.write("</div><div class='row'>")

				g.write("</button><button name='mutation' id='mutation' value='")
				g.write(tripoint)
				g.write("' class='btn-primary btn-lg'>")
				g.write(tripoint)

				g.write("</div><div class='row'>")

				g.write("</button><button name='mutation' id='mutation' value='")
				g.write(triins)
				g.write("' class='btn-primary btn-lg'>")
				g.write(triins)

				g.write("</div><div class='row'>")

				g.write("</button><button name='mutation' id='mutation' value='")
				g.write(tridelins)
				g.write("' class='btn-primary btn-lg'>")
				g.write(tridelins)

				g.write("""</div></form></div>""")

		g.closed
		return render(request, 'index.html', {'triage':True})

	#Key word found in input, but no numbers.
	if mutype != "Invalid mutation" and len(num_list) == 0 and not intron_only:
		with open(os.path.join(interpreter_dir, 'templates/triage.html'), 'w') as g:
			g.write("""<div class='container' style='text-align:center'><h1>You entered: '<c style='color:red'>""")
			g.write(mut)
			g.write("</c>'</h1><br><br><h2>This mutation type (")
			g.write(mutype)
			g.write(""") requires a position to be interpreted.</h2><br><br>Please complete the mutation in the box below, or <a href="{% url 'interpreter:index' %}">return to the previous page</a> to try the "Exon deletion selector" feature.<br><br><form action="{% url 'interpreter:results' %}" method="post"> {% csrf_token %}<input type = "text" id='userAnswer' name = "mutation" id='mutation' value='""")
			try:
				g.write(pm_pre)
			except:
				pass
			g.write(search)
			try:
				g.write(pm_change)
			except:
				pass
			g.write(""" of ' style='height:70px;width:100%;font-size:24pt;'/> <br> *Note that the symbol '>' denotes point mutation*<br><br><button type = "submit" style = 'height:70px;width:150px;font-size:24px;align:center' class = 'btn btn-lg btn-primary'>Interpret</button></form></div>""")
		g.closed
		return render(request, 'index.html', {'triage':True})

	#No input errors.
	if not past_max and mutype != "Invalid mutation": # and len(num_list) != 0:
		if intron_only:
			missense = False
			nonsense = False
			silent = True
			frame_shift = False
			length_mutation, pm_change, pm_pre = get_length(mut, mutype, num_list, mut, ex_input, reference_cdna)
			positions_in_genomic1, positions_in_genomic2, positions_in_genomic3 = cdna2gen(num_list, exon_positions, genomic_positions1, genomic_positions2, genomic_positions3, intron_list)
			ghgvs, ghgvs19 = HGVS(positions_in_genomic2, positions_in_genomic3, mutype, search, length_mutation, ex_input, [[],[]], pm_change, pm_pre, True)


			if mutype == "Point mutation" and len(pm_pre) == 0:
				if mutype == "Point mutation" or mutype == "Deletion/Insertion":
					if len(num_list) == 1:
						pm_pre = genomic_dna[min(positions_in_genomic1)]
					if len(num_list)>1:
						pm_pre = genomic_dna[min(positions_in_genomic1):max(positions_in_genomic1)]
				else:
					if len(num_list) == 1:
						pm_change = genomic_dna[min(positions_in_genomic1)]
					if len(num_list)>1 and length_mutation < 4:
						pm_change = genomic_dna[min(positions_in_genomic1):max(positions_in_genomic1)]
			cDNA, mut_genomic_dna =  alter_cDNA(reference_cdna, num_list, positions_in_genomic1, length_mutation, mutype, pm_change, genomic_dna, intron_only)
			ese, ess, ese_finder, type_ese, type_ess = ese_find(cDNA, reference_cdna, ese_sites, ess_sites, positions_in_genomic1, genomic_dna, length_mutation, mut_genomic_dna)
			splice_test = splice_check(positions_in_genomic1, length_mutation, mut_genomic_dna, genomic_dna)
			nsfp_results, path_pred = gen_point(positions_in_genomic1, pm_change, intron_only, length_mutation, NSFP)
			standard_hgvs, catcher = HGVS(num_list, num_list, mutype, search, length_mutation, ex_input, intron_list, pm_change, pm_pre)
			mv_results = myVariantSearch(standard_hgvs)
			CV = ClinVar(mv_results, standard_hgvs)
			mutype = "Intronic"
			exon_numbers = "Intronic"
			domains = "Intron only"
			consequence_statement = "See <i>In Silico</i> predictions - splice site and splice motif changes"
			therapy = "If splicing changes are predicted, RNA analysis is required to determine eligibility for therapy."
			if len(splice_test) > 1:
			    splice_message = "<b><c style='color:red'>Splice site alteration predicted</c></b>"
			if len(splice_test)>=7:
					i=0
					while i<len(splice_test)-6:
						temp = """<td>%s</td><td>%s</td><td>%s</td><td>%s</td>""" % (str(splice_test[i+4]), str(splice_test[i+1]), str(splice_test[i+3]), str(splice_test[i+6]))
						splice_result += [temp]
						i += 7
			if (ese or ese_finder or ess) and (length_mutation < 30 or mutype=="Deletion/Insertion"):
			    ese_message = "<c style='color:red'><b>Changes to splice regulatory element(s) predicted</b></c>"
			if path_pred > 0 and length_mutation < 30:
			    nsfp_score = ''
			    if path_pred > 2 and path_pred < 6:
			        nsfp_score += "<c style='color:orange'><b>"
			    if path_pred <= 2:
			        nsfp_score += "<c style='color:green'><b>"
			    if path_pred >=6:
			        nsfp_score += "<c style='color:red'><b>"
			    nsfp_message += "%s %s/8 </b></c> <b>prediction algorithms indicate negative impact on protein function.</b><br>" % str(nsfp_score, path_pred)
			if length_mutation < 30 or mutype=="Deletion/Insertion":
			    insilico_message = "<i>Please check the 'In Silico Predictions' for more details</i>"
			if length_mutation >= 30 and mutype!="Deletion/Insertion":
			    insilico_message = "<i>In Silico predictions are not made for large mutations.<br><br>Please see </i><b>'Predicted consequence'</b><i> for reading frame changes, and </i><b>'The Reading-frame Rule'</b><i> in resources for more information</i>"
			if 'Pathogenic' in CV:
			    cv += "<c style='color:red'>"
			if 'Benign' in CV:
			    cv += "<c style='color:green'>"
			cv += CV

            #Start output body for all parseable mutations (non-intron only, valid input)
		elif mutype != "Invalid mutation" and not intron_only:
			length_mutation = 0
			missense = False
			nonsense = False
			silent = False
			frame_shift = False
			if len(num_list) != 0:
				length_mutation, pm_change, pm_pre = get_length(mut, mutype, num_list, mut, ex_input, reference_cdna)
				positions_in_genomic1, positions_in_genomic2, positions_in_genomic3 = cdna2gen(num_list, exon_positions, genomic_positions1, genomic_positions2, genomic_positions3, intron_list)
				ghgvs, ghgvs19 = HGVS(positions_in_genomic2, positions_in_genomic3, mutype, search, length_mutation, ex_input, [[],[]], pm_change, pm_pre, True)
				if len(pm_pre) == 0:
					if mutype == "Point mutation" or "Insertion" in mutype:
						if len(num_list) == 1:
							pm_pre = genomic_dna[min(positions_in_genomic1)]
						if len(num_list)>1:
							pm_pre = genomic_dna[min(positions_in_genomic1)-1:max(positions_in_genomic1)]
					if mutype == "Deletion" or mutype == "Duplication":
						if len(num_list) == 1:
							pm_change = genomic_dna[min(positions_in_genomic1)]
						if len(num_list)>1 and length_mutation < 4:
							pm_change = genomic_dna[min(positions_in_genomic1)-1:max(positions_in_genomic1)]
				cDNA, mut_genomic_dna = alter_cDNA(reference_cdna, num_list, positions_in_genomic1, length_mutation, mutype, pm_change, genomic_dna, intron_only)
				if ("Insertion" in mutype or mutype == "Point mutation") and len(pm_change) == 0:
					with open(os.path.join(interpreter_dir, 'templates/triage.html'), 'w') as g:
						g.write("""<div class='container' style='text-align:center'><h1>You entered: '<c style='color:red'>""")
						g.write(mut)
						g.write("""</c>'<h1><h2>This mutation type requires that you input the inserted nucleotides (i.e. an A, T, G, C or combination thereof).</h2><br><br>Please add the inserted bases in the box below, or <a href="{% url 'interpreter:index' %}">return to the previous page</a> to try the "Exon deletion selector" feature.<br><br><form action="{% url 'interpreter:results' %}" method="post"> {% csrf_token %}<input type = "text" id='userAnswer' name = "mutation" id='mutation' value='c.""")
						if len(num_list) == 1:
							g.write(str(num_list[0]))
						if len(num_list) > 1:
							g.write(str(min(num_list)))
							g.write("_")
							g.write(str(max(num_list)))
						if mutype == "Point mutation":
							try:
								g.write(pm_pre)
							except:
								g.write(" ")
						g.write(search)
						g.write("""' style='height:70px;width:100%;font-size:24pt;'/><br><br><button type = "submit" style = 'height:70px;width:150px;font-size:24px;align:center' class = 'btn btn-lg btn-primary'>Interpret</button></form></div>""")
					g.closed
					return render(request, 'index.html', {'triage':True})
				frame_shift = frame(length_mutation, mutype, num_list)
				exon_numbers = []
				exon_ints = []
				if (mutype == "Deletion" or mutype == "Duplication"):
					exon_numbers, part = exons(num_list, exon_positions, ex_input, mut, part)
					if len(exon_numbers) != 0:
						#ex_input = True
						for item in exon_numbers:
							exon_ints += [int(item)]
							temp = str(exon_numbers[0]), "-", str(exon_numbers[len(exon_numbers)-1])
							temp = ''.join(temp)
				standard_hgvs, catcher = HGVS(num_list, num_list, mutype, search, length_mutation, ex_input, intron_list, pm_change, pm_pre)
				mv_results = myVariantSearch(standard_hgvs)
				CV = ClinVar(mv_results, standard_hgvs)
				aa_ref = translate(reference_cdna[244:])
				nsfp_results, path_pred = gen_point(positions_in_genomic1, complement(pm_change), intron_only, length_mutation, NSFP)

				if mutype == "Deletion" and len(num_list) != 0:
					cDNA, mut_genomic_dna = alter_cDNA(reference_cdna, num_list, positions_in_genomic1, length_mutation, mutype, pm_change, genomic_dna, intron_only)
					aa_seq = translate(cDNA[244:])
					if frame_shift == False:
						aa_length = (length_mutation)/3
						i = 0
						filler = ''
						while i < aa_length:
							filler += '-'
							i+=1
						start1 = 0
						end1 = num_list[0] #if there's no frame shift, the filler will line up perfectly.
						aa_seqgapped = aa_seq[start1:floor(end1/3)] + filler + aa_seq[floor(end1/3):]
						aa_seqgapped = ''.join(aa_seqgapped)
						boundary_change, useless, PTC = easy_align(aa_ref, aa_seqgapped, mutype, frame_shift, True)
					if frame_shift:
						aa_seqgapped = aa_seq
					aa_change, align, PTC = easy_align(aa_ref, aa_seqgapped, mutype, frame_shift)
					if len(exon_ints)!=0 and length_mutation > 31:
						posskip = possible_skips(exon_ints, length_mutation, reference_cdna, frame_shift, exon_positions, mutype, num_list, aa_ref)

				if mutype == "Duplication":
					cDNA, mut_genomic_dna = alter_cDNA(reference_cdna, num_list, positions_in_genomic1, length_mutation, mutype, pm_change, genomic_dna, intron_only)
					aa_seq = translate(cDNA[244:])
					aa_change, align, PTC = easy_align(aa_ref, aa_seq, mutype, frame_shift)
					posskip = possible_skips(exon_ints, length_mutation, reference_cdna, frame_shift, exon_positions, mutype, num_list, aa_ref)

				if mutype == "Deletion/Insertion" or mutype == "Point mutation" or mutype == "Insertion":
					cDNA, mut_genomic_dna = alter_cDNA(reference_cdna, num_list, positions_in_genomic1, length_mutation, mutype, pm_change, genomic_dna, intron_only)
					aa_seq = translate(cDNA[244:])
					aa_change, align, PTC = easy_align(aa_ref, aa_seq, mutype, frame_shift)
					exon_numbers, part = exons(num_list, exon_positions, ex_input, mut, part)
				if len(num_list) != 0:
					ese, ess, ese_finder, type_ese, type_ess = ese_find(cDNA, reference_cdna, ese_sites, ess_sites, positions_in_genomic1, genomic_dna, length_mutation, mut_genomic_dna)

				if mutype != "Duplication":
					long_aa_change = aa_longform(aa_change, mutype, frame_shift, num_list, length_mutation, aa_ref)

				if mutype == "Duplication":
					if frame_shift:
						long_aa_change = aa_longform(aa_change, mutype, frame_shift, num_list, length_mutation, aa_ref)

				if len(num_list) != 0:
					silent, missense, nonsense, statement1, statement2, statement3 = output_statements(mut, mutype, length_mutation, aa_change, exon_numbers, part, frame_shift, long_aa_change, aa_seq, aa_ref, cDNA)

				if len(num_list) == 0:
					frame_shift = False
					silent = False
					missense = False
					nonsense = False
					ese = False
					ese_finder = False
					ess = False
					type_ese = ''
					statement1 = "This mutation is solely in the intronic region of the DMD gene."
					statement2 = ''
				if (frame_shift == True and len(num_list) != 0) or mutype != "Point mutation":
					if mutype != "Deletion/Insertion":
						silent = False
						missense = False
						nonsense = False
				splice_test = splice_check(positions_in_genomic1, length_mutation, mut_genomic_dna, genomic_dna)
				if mutype == "Deletion" and length_mutation == 1:
					mutype = "Point mutation (deletion)"
				if mutype == "Insertion" and length_mutation == 1:
					mutype = "Point mutation (insertion)"
				ds = domain(num_list, domains)

				if "Deletion of exon" in mut:
					if len(exon_ints) == 1:
						if mut[len(mut)-2] == ",":
							mut = mut[0:len(mut)-2]
						else:
							mut = mut
					if len(exon_ints) > 1:
						mut = "Deletion of exons ", str(min(exon_ints)),"-",str(max(exon_ints))
						mut = ''.join(mut)
				if frame_shift and aa_change[0] != "STOP":
					consequence = "<c style='color:red'>Frame Shift</c><br>"
					consequence_statement = statement3
				if silent:
					consequence = "<c style='color:orange'>Silent</c>, p.(=)"
				if missense:
					consequence = "<c style='color:orange'>Missense</c><br>"
					consequence_statement = statement1
				if nonsense or (frame_shift and aa_change[0] == "STOP"):
					consequence = "<c style='color:red'>Nonsense</c><br>"
					consequence_statement = statement1
				if (mutype == "Duplication" or mutype == "Insertion") and not frame_shift:
					consequence_statement = "Insertion without frameshift. New amino acids incorporated, with remaining sequence intact."
				if mutype == "Deletion" and not frame_shift:
					temp = 'Deletion without frame shift: '
					temp += aa_longform("", mutype, frame_shift, num_list, length_mutation, aa_ref)
					i=0
					if PTC != '':
						temp2 = "\n<br><c style='color:orange'>(Change at boundary: ", PTC, ")</c>"
						temp2 = ''.join(temp2)
						temp += temp2
					if length_mutation > 31:
						if PTC == '':
							temp2 = "\n<br><c style='color:green'>(No change predicted at boundary of exon ",str(min(exon_ints)-1), " and exon ", str(max(exon_ints)+1), ")</c>"
							temp2 = ''.join(temp2)
							temp += temp2
					consequence_statement = temp
				if (len(posskip)>0 and frame_shift and length_mutation>30) or nonsense:
					therapy = "<c style='color:green'>Possibly</c> - Please see 'Therapies' tab"
				elif (len(posskip)==0 and not nonsense) or (len(posskip)>0 and not frame_shift) or (length_mutation <= 30 and not nonsense):
					therapy = "<c style='color:orange'>Not currently</c> - Please see 'Therapies' tab"
				if len(splice_test) > 1:
					splice_message = "<c style='color:red'><b>Splice site alteration predicted</b></c>"
				if (ese or ese_finder or ess) and (length_mutation < 30 or mutype=="Deletion/Insertion"):
					motif_message = "<c style='color:red'><b>Changes to splice regulatory element(s) predicted</b></c>"
				if missense:
					nsfp_message =''
					if path_pred > 2 and path_pred < 6:
						nsfp_message += "<c style='color:orange'><b>"
					if path_pred <= 2:
						nsfp_message += "<c style='color:green'><b>"
					if path_pred >=6:
						nsfp_message += "<c style='color:red'><b>"
					nsfp_message += str(path_pred)
					nsfp_message += "/8 </c></b><b> prediction algorithms indicate negative impact on protein function.</b><br>"
				if length_mutation < 30 or mutype=="Deletion/Insertion":
					insilico_message = "<i>Please check the 'In Silico Predictions' for more details</i>"
				if length_mutation >= 30 and mutype!="Deletion/Insertion":
					insilico_message = "<i>In Silico predictions are not made for this mutation type.<br><br>Please see </i><b>'Mutation type'</b> and <b>'Predicted consequence'</b><i> for reading frame changes, and </i><b>'Reading-frame rule'</b><i> in resources for more information</i>"
				if 'Pathogenic' in CV:
					cv += "<c style='color:red'>"
				if 'Benign' in CV:
					cv += "<c style='color:green'>"
				cv += CV

		            #Splicing
			if length_mutation < 30 or mutype == "Deletion/Insertion":
				if missense:
					if path_pred > 2 and path_pred < 6:
						nsfp_score += "<c style='color:orange'>"
					if path_pred <= 2:
						nsfp_score += "<c style='color:green'>"
					if path_pred >=6:
						nsfp_score += "<c style='color:red'>"
					nsfp_score += str(path_pred)
					nsfp_score += "/8 </c>"
				sift = nsfp_results[1]
				polyphen = nsfp_results[2]
				lrt = nsfp_results[3]
				mutationtaster = nsfp_results[4]
				mutationassessor = nsfp_results[5]
				fathmm = nsfp_results[6]
				provean = nsfp_results[7]
				metasvm = nsfp_results[8]
				if ese:
					rescueese = "<c style='color:red'>Disruption of ESE</c>"
				if not ese:
					rescueese = "<c style='color:green'>No change predicted</c>"
				if not ese_finder:
					esefinder += """
					<td>None</td>
					<td><c style='color:green'>No motifs significantly changed</c></td>"""
				if ese_finder:
					i = 0
					esefinder += """
				    <td>"""
					while i<len(type_ese):
						esefinder += type_ese[i][0]
						esefinder += "<br>"
						i+=1
					esefinder += """
				    </td>
				    <td>"""
					i=0
					while i<len(type_ese):
						esefinder += type_ese[i][2]
						esefinder += "<br>"
						i+=1
					esefinder += "</td>"
				if ess:
					ess_message = "<c style='color:red'>Mutation creates a novel ESS motif</c>"
				if not ess:
					ess_message = "<c style='color:green'>No ESS motifs created</c>"
				if len(splice_test)<7:
					splice_message2 = "<c style='color:green'>No changes predicted</c>"
				if len(splice_test)>=7:
					i=0
					while i<len(splice_test)-6:
						temp = """<td>%s</td><td>%s</td><td>%s</td><td>%s</td>""" % (str(splice_test[i+4]), str(splice_test[i+1]), str(splice_test[i+3]), str(splice_test[i+6]))
						splice_result += [temp]
						i += 7

			if length_mutation >= 30 and mutype != "Deletion/Insertion":
				large_deletion = True
				splice_message2 = "<c style='color:orange'><b>Predictions are not made for large mutations</b></c>"

		# ESEFinder predictions are based on changes in consensus scores after mutation according to matrices at the Cold Spring Harbor Laboratory: <a href='http://rulai.cshl.edu/cgi-bin/tools/ESE3/esefinder.cgi?process=matrices'>http://rulai.cshl.edu/cgi-bin/tools/ESE3/esefinder.cgi?process=matrices</a> <br>Specific binding sites predicted to be effected are listed in brackets following ESEFinder.<br><br>
		# Questions/Issues: <a href='mailto:dove.dmd.interpreter@gmail.com?Subject=DOVE%20Web%20App'>Contact Us</a>.<br>


		##EXON SKIPPING
			if len(posskip)>0 and frame_shift:
				if part == True:
					posskip = []
			if part and ex_input:
				print(part)
				posskip = []
				posskip += ["This variant includes a partial exon deletion."]
			elif not frame_shift:
				posskip = []
				posskip += ["<c style='color:gray'><b>There is no frameshift predicted.</b></c>"]
			elif length_mutation < 32:
				posskip = []
				posskip += ["This variant type (missense; nonsense; small insertion, deletion, indel; or splice-affecting) has not been clinically tested in <i>DMD</i> with exon skip therapy."]
			elif len(posskip) == 0 and length_mutation >= 32:
				posskip += ["There are no theoretical exon skips predicted to apply to this mutation."]

			if nonsense and not frame_shift:
				readthrough_elig = "<c style='color:green'><b>Eligible</b></c>"
			else:
				pass


		else:
			return render(request, 'index.html',{'success':False})


		if mutype != "Deletion" and mutype != "Deletion/Insertion" and mutype != "Duplication" and not intron_only:
			temp_cdna = "<em class='gray'>",cDNA[0:min(num_list)+243],"</em><em class='red'>",cDNA[min(num_list)+243:min(num_list)+243+length_mutation],"</em><em class='gray'>",cDNA[max(num_list)+243+length_mutation:len(cDNA)],"</em>"
		if mutype == "Point mutation (insertion)":
			temp_cdna = "<em class='gray'>",cDNA[0:max(num_list)+244],"</em><em class='red'>",cDNA[max(num_list)+244:max(num_list)+244+length_mutation],"</em><em class='gray'>",cDNA[max(num_list)+244+length_mutation:len(cDNA)],"</em>"
		if mutype == "Duplication" and not intron_only:
			temp_cdna = "<em class='gray'>",cDNA[0:min(num_list)+244],"</em><em class='green'>",cDNA[min(num_list)+244:min(num_list)+244+(length_mutation)],"</em><em class='red'>",cDNA[min(num_list)+244+(length_mutation):min(num_list)+244+(length_mutation*2)],"</em><em class='gray'>",cDNA[min(num_list)+244+(length_mutation*2):len(cDNA)],"</em>"
		if mutype == "Insertion" and not intron_only:
			temp_cdna = "<em class='gray'>",cDNA[0:min(num_list)+243],"</em><em class='red'>",cDNA[min(num_list)+244:min(num_list)+244+length_mutation],"</em><em class='gray'>",cDNA[max(num_list)+244+length_mutation:len(cDNA)],"</em>"
		if mutype == "Deletion" and not intron_only:
			fill = ''
			i = 0
			while i<length_mutation:
				fill += 'x'
				i+=1
			temp_cdna = "<em class='gray'>",cDNA[0:min(num_list)+244],"</em><em class='red'>",fill,"</em><em class='gray'>",cDNA[min(num_list)+244:len(cDNA)],"</em>"
		if intron_only:
			temp_refcdna = "<em class='color:gray'>",reference_cdna,"</em>"
			temp_cdna = temp_refcdna
			refcdna_for_print = ''.join(temp_refcdna)
		if mutype == "Deletion/Insertion" and not intron_only:
			i = 0
			fill = ''
			overall = length_mutation[0] - length_mutation[1]
			if overall > 0:
				while i<overall:
					fill += 'x'
					i += 1
				temp_cdna = "<em class='gray'>",cDNA[0:min(num_list)+243],"</em><em class='red'>",fill,"</em><em class='gray'>",cDNA[max(num_list)+243+1:len(cDNA)],"</em>"
			if overall <= 0:
				temp_cdna = "<em class='gray'>",cDNA[0:min(num_list)+243],"</em><em class='red'>",cDNA[min(num_list)+243:min(num_list)+243+max(length_mutation[0],length_mutation[1])],"</em><em class='gray'>",cDNA[min(num_list)+243+max(length_mutation[0],length_mutation[1]):len(cDNA)],"</em>"
			temp_refcdna = "<em class='gray'>",reference_cdna[0:min(num_list)+243],"</em><em class='green'>",reference_cdna[min(num_list)+243:min(num_list)+243+max(length_mutation[0],length_mutation[1])],"</em><em class='gray'>",reference_cdna[min(num_list)+243+max(length_mutation[0],length_mutation[1]):len(reference_cdna)],"</em>"
			refcdna_for_print = ''.join(temp_refcdna)

		cdna_for_print = ''.join(temp_cdna)

		if mutype != "Deletion/Insertion" and not intron_only:
			if mutype == "Point mutation (insertion)" or mutype == "Insertion":
				temp_refcdna = "<em class='color:gray'>",reference_cdna,"</em>"
				temp_cdna = temp_refcdna
				refcdna_for_print = ''.join(temp_refcdna)
			elif mutype == "Deletion":
				temp_refcdna = "<em class='gray'>",reference_cdna[0:min(num_list)+244]+"</em><em class='green'>"+reference_cdna[min(num_list)+244:min(num_list)+244+length_mutation]+"</em><em class='gray'>"+reference_cdna[min(num_list)+244+length_mutation:len(reference_cdna)]+"</em>"
			elif mutype == "Duplication":
				temp_refcdna = "<em class='gray'>",reference_cdna[0:min(num_list)+244]+"</em><em class='green'>"+reference_cdna[min(num_list)+244:min(num_list)+244+length_mutation]+"</em><em class='gray'>"+reference_cdna[min(num_list)+244+length_mutation:len(reference_cdna)]+"</em>"
			else:
				temp_refcdna = "<em class='gray'>",reference_cdna[0:min(num_list)+243],"</em><em class='green'>",reference_cdna[min(num_list)+243:min(num_list)+243+length_mutation],"</em><em class='gray'>",reference_cdna[min(num_list)+243+length_mutation:len(reference_cdna)],"</em>"
			refcdna_for_print = ''.join(temp_refcdna)
			try:
				refcdna_print_preview = refcdna_for_print[min(num_list)+17+243-100:min(num_list)+1+17+243+23+22+(length_mutation*2)+100]
				cdna_print_preview = cdna_for_print[min(num_list)+17+243-100:min(num_list)+1+17+243+22+21+(length_mutation*2)+100]
			except:
				refcdna_print_preview = refcdna_for_print[min(num_list)+17+243:min(num_list)+1+17+243+23+22+length_mutation]
				cdna_print_preview = cdna_for_print[min(num_list)+17+243:min(num_list)+1+17+243+22+21+length_mutation]
		if mutype == "Deletion/Insertion" and not intron_only:
			try:
				refcdna_print_preview = refcdna_for_print[min(num_list)+17+243-100:min(num_list)+1+17+243+23+22+max(length_mutation[0],length_mutation[1])+100]
				cdna_print_preview = cdna_for_print[min(num_list)+17+243-100:min(num_list)+1+17+243+22+21+max(length_mutation[0],length_mutation[1])+100]
			except:
				refcdna_print_preview = refcdna_for_print[min(num_list)+17+243:min(num_list)+1+17+243+23+22+max(length_mutation[0],length_mutation[1])]
				cdna_print_preview = cdna_for_print[min(num_list)+17+243:min(num_list)+1+17+243+22+21+max(length_mutation[0],length_mutation[1])]
		if intron_only:
			refcdna_print_preview = refcdna_for_print[0:200] + "</em>"
			cdna_print_preview = cdna_for_print[0:200] + "</em>"
			exon_numbers = "NA"

		if mutype == "Point mutation (deletion)":
			standard_hgvs_search = standard_hgvs[0:len(standard_hgvs)-1]
		else:
			standard_hgvs_search = standard_hgvs



		search_nums = ''
		search_nums += str(min(num_list)) + " "
		if min(num_list) != max(num_list):
			search_nums += str(max(num_list)) + " "
		if len(intron_list[0])>1:
			search_nums += str(intron_list[0][1]) + " "
		if len(intron_list[1])>1:
			search_nums += str(intron_list[1][1]) + " "

		# Example: https://databases.lovd.nl/shared/variants/DMD/unique?search_VariantOnTranscript/Exon=6&search_VariantOnTranscript/DNA=358%20530%20del
		leiden_base_link ="https://databases.lovd.nl/shared/variants/DMD/unique?"
		leiden_exon = f"search_VariantOnTranscript/Exon={' '.join([str(x) for x in exon_numbers])}"
		temp_nums = re.findall('\d+', standard_hgvs.split('c.')[-1])
		leiden_positions_type = f"search_VariantOnTranscript/DNA={'%20'.join([str(x) for x in temp_nums])}"
		if '>' in standard_hgvs:
			leiden_positions_type += f" {standard_hgvs.split('>')[-1]}"
		for kwd in ['del', 'ins', 'dup']:
			if kwd in standard_hgvs:
				leiden_positions_type += f"%20{kwd}"
		leiden_link = '&'.join([leiden_base_link, leiden_exon, leiden_positions_type])
		iframe = "<iframe src='"+leiden_link+ "'align = 'center' width = '100%' height = '800' style = 'background-color: white'></iframe>"

		if nsfp_score == '':
			nsfp_score = 'N/A'
		else:
			pass

		try:
			mv_results = mv_results[0]
			mv_results_formatted = pprint.pformat(mv_results, indent=4)
			mv_disable = False
		except:
			mv_results_formatted = "No myVariantInfo results available"
			mv_disable = True
	return render(request, 'main.html',
		{'user_inp':mut,
		'mutype':mutype,
		'hgvs':standard_hgvs,
		'ghgvs':ghgvs,
		'ghgvs19':ghgvs19,
		'length_mutation':length_mutation,
		'exons': exon_numbers,
		'domains':ds,
		'therapy':therapy,
		'insilico_message':insilico_message,
		'splice_message':splice_message,
		'ese_message':ese_message,
		'consequence':consequence,
		'consequence_message':consequence_statement,
		'cv':cv,
		'leiden_link':leiden_link,
		'leiden_frame':iframe,
		'cdna_for_print':cdna_for_print, #full mutated sequence
		'cdna_print_preview':cdna_print_preview, #partial mutated
		'refcdna_for_print':refcdna_for_print, #full wt cdna
		'refcdna_print_preview':refcdna_print_preview, #wt preview
		'readthrough_elig':readthrough_elig,
		'posskip':posskip, #list of skips (empty if there are none)
		'sift':sift,
		'polyphen':polyphen,
		'lrt':lrt,
		'mutationtaster':mutationtaster,
		'mutationassessor':mutationassessor,
		'fathmm':fathmm,
		'provean':provean,
		'metasvm':metasvm,
		'esefinder':esefinder,
		'ess_message':ess_message,
		'rescueese':rescueese,
		'large_deletion':large_deletion,
		'nsfp_score':nsfp_score,
		'nsfp_message':nsfp_message,
		'splice_result':splice_result,
		'splice_message2':splice_message2,
		'motif_message':motif_message,
		'mv_results': mv_results,
		'mv_format': mv_results_formatted,
		'mv_disable': mv_disable
		}
		)
