import datetime
from django import forms
from .models import Plasmid
from .models import Primer
from .models import GlycerolStock
from .custom.standards import ligation_standards
from organization.views import get_projects_where_member_can
from organization.views import get_projects_where_member_can_any


class PlasmidNameInput(forms.Form):
    plasmid_name = forms.CharField(max_length=50, min_length=1)


class GlycerolQRInput(forms.Form):
    glycerol_qr_id = forms.CharField()


class BlastSequenceInput(forms.Form):
    def __init__(self, project_choices, *args, **kwargs):
        super(BlastSequenceInput, self).__init__(*args, **kwargs)
        self.fields['project'].choices = project_choices

    project = forms.ChoiceField(label="Project to search in", choices=())
    fasta_sequence = forms.CharField(label="Fasta Text Input (Preferred)", widget=forms.Textarea(attrs={}), required=False)
    fasta_file = forms.FileField(label="Fasta File Input", required=False)
    short_blast = forms.BooleanField(label="Use short input BLAST parameters?", required=False)


class L0SequenceInput(forms.Form):
    l0_sequence_input = forms.CharField(widget=forms.Textarea, label="Sequence input")
    # Todo append None
    ligation_standard_slug = forms.ChoiceField(choices=tuple([(index, ligation_standard['name']) for index, ligation_standard in ligation_standards.items()]),
                                required=True, label="Ligation standard")
    l0_oh_5 = forms.ChoiceField(choices=tuple([(oh_slug, oh['name'] + " [" + oh['oh'] + "]") for index, ligation_standard in ligation_standards.items() for oh_slug, oh in ligation_standard['ohs']['l0'].items()]),
                                required=True, label="L0 OH 5'")
    l0_oh_3 = forms.ChoiceField(choices=tuple([(oh_slug, oh['name'] + " [" + oh['oh'] + "]") for index, ligation_standard in ligation_standards.items() for oh_slug, oh in ligation_standard['ohs']['l0'].items()]),
                                required=True, label="L0 OH 3'")
    enzyme = forms.CharField(required=True, widget=forms.HiddenInput())


class FastaAlignForm(forms.Form):
    fasta_sequence = forms.CharField(widget=forms.Textarea(attrs={}), required=False)
    save_clustal_file = forms.BooleanField(required=False)
    is_reversed = forms.BooleanField(label="Is reversed?", required=False)
    fasta_file = forms.FileField(label="Fasta File", required=False)


class SangerAlignForm(forms.Form):
    save_clustal_file = forms.BooleanField(required=False)
    is_reversed = forms.BooleanField(label="Is reversed?", required=False)
    ab1 = forms.FileField(label="AB1 File", required=True)


class DateInput(forms.DateInput):
    input_type = 'date'


class GstockCreateForm(forms.ModelForm):
    class Meta:
        model = GlycerolStock
        fields = '__all__'
        exclude = ('project',)
        widgets = {
            'created_on': DateInput()
        }


class GstockEditForm(forms.ModelForm):
    class Meta:
        model = GlycerolStock
        fields = '__all__'
        widgets = {
            'created_on': DateInput()
        }


class PlasmidCreateForm(forms.ModelForm):
    def __init__(self, *args, **kwargs):
        member = kwargs.pop('user')
        super(PlasmidCreateForm, self).__init__(*args, **kwargs)
        user_visible_plasmids = Plasmid.objects.filter(project__in=get_projects_where_member_can_any(member)).order_by('name')
        # self.fields['backbone'].queryset = user_visible_plasmids.filter(type=1)
        # self.fields['inserts'].queryset = user_visible_plasmids.filter(type=0)

    class Meta:
        model = Plasmid
        fields = ['name', 'selectable_markers', 'sequence', 'backbone', 'inserts', 'intended_use', 'type', 'level',
                  'description', 'created_on']


class PlasmidEditForm(forms.ModelForm):
    def __init__(self, *args, **kwargs):
        member = kwargs.pop('user')
        super(PlasmidEditForm, self).__init__(*args, **kwargs)
        user_visible_plasmids = Plasmid.objects.filter(project__in=get_projects_where_member_can_any(member)).order_by('name')
        self.fields['project'].queryset = get_projects_where_member_can(member, ['a', 'w'])
        # self.fields['backbone'].queryset = user_visible_plasmids.filter(type=1)
        # self.fields['inserts'].queryset = user_visible_plasmids.filter(type=0)

    class Meta:
        model = Plasmid
        fields = ['name', 'selectable_markers', 'sequence', 'backbone', 'inserts', 'intended_use', 'type', 'level',
                  'description', 'created_on', 'project', 'reference_sequence', 'under_construction']


class PlasmidValidationForm(forms.ModelForm):
    class Meta:
        model = Plasmid
        fields = ['working_colony',
                  'colonypcr_state', 'colonypcr_date', 'colonypcr_observations',
                  'digestion_state', 'digestion_date', 'digestion_observations',
                  'sequencing_state', 'sequencing_date', 'sequencing_observations',
                  'sequencing_observations', 'sequencing_clustal_file']
        widgets = {
            'colonypcr_date': DateInput(),
            'digestion_date': DateInput(),
            'sequencing_date': DateInput()
        }


class DigestForm(forms.Form):
    enzymes = forms.CharField(widget=forms.HiddenInput(attrs={'id': 'digest_enzymes'}))


class PCRForm(forms.Form):
    primer_f = forms.ModelChoiceField(queryset=Primer.objects.all(), to_field_name="id", label="Primer F",
                                      required=False)
    primer_r = forms.ModelChoiceField(queryset=Primer.objects.all(), to_field_name="id", label="Primer R",
                                      required=False)
    primer_f_seq = forms.CharField(label="Primer F sequence", required=False)
    primer_r_seq = forms.CharField(label="Primer R sequence", required=False)


class MsaUploadAb1FilesForm(forms.Form):
    ab1_file_1 = forms.FileField(label="AB1 File 1")
    ab1_file_2 = forms.FileField(label="AB1 File 2", required=False)


class MsaChromatosStep2Form(forms.Form):
    from1 = forms.IntegerField(label="From Chromato 1")
    to1 = forms.IntegerField(label="To Chromato 1")
    from2 = forms.IntegerField(label="From Chromato 2", required=False)
    to2 = forms.IntegerField(label="To Chromato 2", required=False)
    sequence1 = forms.CharField(widget=forms.HiddenInput)
    sequence2 = forms.CharField(widget=forms.HiddenInput, required=False)
    target = forms.ModelChoiceField(queryset=Plasmid.objects.all(), to_field_name="id")


class MsaUploadFastaFileForm(forms.Form):
    fasta_text = forms.CharField(widget=forms.Textarea)
    fasta_file = forms.FileField(label="Fasta file")


class PlasmidLabel(forms.Form):
    date = forms.DateField(widget=DateInput(), initial=datetime.date.today)
    colony = forms.CharField()
    concentration = forms.CharField()
