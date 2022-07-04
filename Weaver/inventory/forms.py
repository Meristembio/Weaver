from django import forms
from .models import Plasmid
from .models import Primer
from .models import GlycerolStock
from .custom.standards import CURRENT_ASSEMBLY_STANDARD

LO_OHS = CURRENT_ASSEMBLY_STANDARD['odd_custom']


class GlycerolQRInput(forms.Form):
    glycerol_id = forms.CharField()


class BlastSequenceInput(forms.Form):
    sequence_input = forms.CharField(widget=forms.Textarea)


class L0SequenceInput(forms.Form):
    l0_sequence_input = forms.CharField(widget=forms.Textarea, label="Sequence input")
    l0_oh_5 = forms.ChoiceField(choices=tuple([(name, name + " [" + value + "]") for name, value in LO_OHS]), required=True, label="L0 OH 5'")
    l0_oh_3 = forms.ChoiceField(choices=tuple([(name, name + " [" + value + "]") for name, value in LO_OHS]), initial="B", required=True, label="L0 OH 3'")
    enzyme = forms.CharField(required=True, widget=forms.HiddenInput())


class SangerForms(forms.Form):
    ab1 = forms.FileField(label="AB1 File", required=True)
    is_reverse = forms.BooleanField(label="Is reverse?", required=False)


class DateInput(forms.DateInput):
    input_type = 'date'


class GstockEditForm(forms.ModelForm):
    class Meta:
        model = GlycerolStock
        fields = '__all__'
        widgets = {
            'created_on': DateInput()
        }


class PlasmidValidationForm(forms.ModelForm):
    class Meta:
        model = Plasmid
        fields = ['working_colony', 'check_state', 'check_method', 'check_date', 'digestion_check_enzymes', 'check_observations', 'sequencing_state', 'sequencing_date', 'sequencing_observations', 'sequencing_ab1_file_1', 'sequencing_ab1_file_2', 'sequencing_observations']
        widgets = {
            'check_date': DateInput(),
            'sequencing_date': DateInput()
        }


class DigestForm(forms.Form):
    enzymes = forms.CharField(widget=forms.HiddenInput(attrs={'id': 'digest_enzymes'}))


class PCRForm(forms.Form):
    primer_f = forms.ModelChoiceField(queryset=Primer.objects.all(), to_field_name="id", label="Primer F", required=False)
    primer_r = forms.ModelChoiceField(queryset=Primer.objects.all(), to_field_name="id", label="Primer R", required=False)
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
