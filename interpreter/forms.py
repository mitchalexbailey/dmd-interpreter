from django import forms

choices = [(True,'Yes'),(False,'No')]

class IndexForm(forms.Form):
	mutation = forms.CharField(label = 'Mutation', max_length = 100)

class ACMGForm(forms.Form):
	pvs1 = forms.TypedChoiceField(label = 'Does the variant cause a premature stop codon (nonsense)?', choices=choices, widget=forms.RadioSelect)
	pvs2 = forms.TypedChoiceField(label = 'Does the variant cause a frameshift?', choices=choices, widget=forms.RadioSelect)
	pvs3 = forms.TypedChoiceField(label = 'Does the variant cause a multiexon deletion involving key functional domains?', choices=choices, widget=forms.RadioSelect)
	pvs4 = forms.TypedChoiceField(label = 'Does the variant cause a change in splice site?', choices=choices, widget=forms.RadioSelect)
	pvs5 = forms.TypedChoiceField(label = 'Does the variant cause a change in initiation codon?', choices=choices, widget=forms.RadioSelect)
	ps1 = forms.TypedChoiceField(label = 'Is there a reported pathogenic variant causing the same amino acid change?', choices=choices, widget=forms.RadioSelect)
	ps2 = forms.TypedChoiceField(label = 'Is the variant <i>de novo</i> (confirmed to not be present in either parent)?', choices=choices, widget=forms.RadioSelect)
	ps3 = forms.TypedChoiceField(label = 'Are there well-established <i>in vitro</i> studies predicting a damaging effect of this variant on the gene or gene product?', choices=choices, widget=forms.RadioSelect)
	ps4 = forms.TypedChoiceField(label = 'Is this variant more prevalence in affected individuals versus controls? (OR > 5.0)', choices=choices, widget=forms.RadioSelect)
	pm1 = forms.TypedChoiceField(label = 'Is this variant located in a mutational hot spot and/or critical functional domain?', choices=choices, widget=forms.RadioSelect)
	pm2 = forms.TypedChoiceField(label = 'Is the variant absent from controls (autosomal dominant), or found at an extremely low frequency (autosomal recessive)? (Ex. in ExAC or 1000 genomes)', choices=choices, widget=forms.RadioSelect)
	pm3 = forms.TypedChoiceField(label = 'Is this variant in a gene linked to an autosomal recessive condition and in <i>trans</i> with a pathogenic variant?', choices=choices, widget=forms.RadioSelect)
	pm4 = forms.TypedChoiceField(label = 'Does the variant change the protein length (while preserving reading frame; deletion, insertion, stop-loss)?', choices=choices, widget=forms.RadioSelect)
	pm5 = forms.TypedChoiceField(label = 'Does the variant cause a missense change at a residue where a different change is known to be pathogenic?', choices=choices, widget=forms.RadioSelect)
	pm6 = forms.TypedChoiceField(label = 'Do you think the variant is <i>de novo</i> but there has not been confirmation (sequencing of parents)?', choices=choices, widget=forms.RadioSelect)
	pp1 = forms.TypedChoiceField(label = 'Is the variant in a known disease-causing gene and has it co-segregated with affected family members?', choices=choices, widget=forms.RadioSelect)
	pp2 = forms.TypedChoiceField(label = 'Is the variant in a gene where disease-causing variants are not commonly missense, and in which missense variants are a common mechanism of disease?', choices=choices, widget=forms.RadioSelect)
	pp3 = forms.TypedChoiceField(label = 'Do multiple <i>in silico</i> functional predication tools support a deleterious effect on the gene or gene product?', choices=choices, widget=forms.RadioSelect)
	pp4 = forms.TypedChoiceField(label = 'Does the patient\'s phenotype and/or family history strongly indicate a disease with a single genetic ontology?', choices=choices, widget=forms.RadioSelect)
	pp5 = forms.TypedChoiceField(label = 'Does a reputable source report the variant as pathogenic (but the evidence is not available)?', choices=choices, widget=forms.RadioSelect)