import datetime
from django import forms
from django.urls import reverse
from .models import Plasmid
from .models import Primer
from .models import GlycerolStock
from .models import Experiment
from .models import RestrictionBuffer
from .models import RestrictionEnzyme
from .models import RestrictionEnzymeBuffer
from .custom.primer_access import visible_primers_for_user
from .custom.standards import assembly_standards
from organization.models import Project
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
    assembly_standard_slug = forms.ChoiceField(choices=tuple([(index, assembly_standard['name']) for index, assembly_standard in assembly_standards.items()]),
                                required=True, label="Ligation standard")
    l0_oh_5 = forms.ChoiceField(choices=tuple([(oh_slug, oh['name'] + " [" + oh['oh'] + "]") for index, assembly_standard in assembly_standards.items() for oh_slug, oh in assembly_standard['ohs']['l0'].items()]),
                                required=True, label="L0 OH 5'")
    l0_oh_3 = forms.ChoiceField(choices=tuple([(oh_slug, oh['name'] + " [" + oh['oh'] + "]") for index, assembly_standard in assembly_standards.items() for oh_slug, oh in assembly_standard['ohs']['l0'].items()]),
                                required=True, label="L0 OH 3'")
    enzyme = forms.CharField(required=True, widget=forms.HiddenInput())


class FastaFileInput(forms.ClearableFileInput):
    allow_multiple_selected = True


class FastaFileField(forms.FileField):
    def clean(self, data, initial=None):
        single_file_clean = super().clean
        if isinstance(data, (list, tuple)):
            return [single_file_clean(d, initial) for d in data]
        return single_file_clean(data, initial)


class FastaAlignForm(forms.Form):
    def __init__(self, *args, **kwargs):
        super(FastaAlignForm, self).__init__(*args, **kwargs)
        self.fields['fasta_sequence'].widget.attrs.update({'class': 'form-control', 'rows': 10, 'placeholder': '>query\nACGT...'})
        self.fields['fasta_file'].widget.attrs.update({'class': 'form-control', 'accept': '.fa,.fas,.fasta,.txt'})
        self.fields['alignment_view_mode'].widget.attrs.update({'class': 'form-select'})
        self.fields['save_clustal_file'].widget.attrs.update({'class': 'form-check-input'})

    fasta_sequence = forms.CharField(widget=forms.Textarea(attrs={}), required=False)
    alignment_view_mode = forms.ChoiceField(
        label="View mode",
        choices=(("combined", "Together"), ("individual", "One at a time")),
        initial="combined",
    )
    save_clustal_file = forms.BooleanField(required=False)
    fasta_file = FastaFileField(label="Fasta File", required=False, widget=FastaFileInput(attrs={"multiple": True}))


class MultipleFileInput(forms.ClearableFileInput):
    allow_multiple_selected = True


class MultipleFileField(forms.FileField):
    def clean(self, data, initial=None):
        single_file_clean = super().clean
        if isinstance(data, (list, tuple)):
            return [single_file_clean(d, initial) for d in data]
        return single_file_clean(data, initial)


class SangerAlignForm(forms.Form):
    save_clustal_file = forms.BooleanField(required=False)
    label = forms.CharField(label="Run label", required=False, max_length=200)
    notes = forms.CharField(label="Notes", required=False, widget=forms.Textarea(attrs={"rows": 2}))
    sanger_files = MultipleFileField(label="Sanger files", required=False, widget=MultipleFileInput(attrs={
        "multiple": True,
        "accept": ".ab1,.phd.1,.seq,.fa,.fas,.fasta",
        "id": "id_sanger_files",
    }))
    ab1 = forms.FileField(label="AB1 File", required=False)

    def clean(self):
        cleaned_data = super().clean()
        if not self.files.getlist("sanger_files") and not self.files.get("ab1"):
            raise forms.ValidationError("Upload at least one .ab1, .phd.1, or .seq file.")
        return cleaned_data


class SangerBatchUploadForm(forms.Form):
    ab1_files = MultipleFileField(
        label="AB1 files",
        widget=MultipleFileInput(attrs={
            "multiple": True,
            "accept": ".ab1",
            "id": "id_ab1_files",
        }),
    )
    mapping_csv = forms.FileField(label="CSV mapping file")
    replace_existing = forms.BooleanField(
        label="Replace existing files",
        required=False,
        help_text="Process a new version when an AB1 with the same name already exists for that plasmid.",
    )

    def clean(self):
        cleaned_data = super().clean()
        files = self.files.getlist("ab1_files")
        invalid_files = [file.name for file in files if not file.name.lower().endswith(".ab1")]
        if invalid_files:
            self.add_error("ab1_files", "Only .ab1 files are accepted: {}".format(", ".join(invalid_files)))
        if not files:
            self.add_error("ab1_files", "Upload at least one .ab1 file.")
        return cleaned_data


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


class ExperimentPlasmidField(forms.ModelMultipleChoiceField):
    def label_from_instance(self, obj):
        return f'{obj.name} [{obj.idx}]'


class ExperimentForm(forms.ModelForm):
    plasmids = ExperimentPlasmidField(queryset=Plasmid.objects.none())

    class Meta:
        model = Experiment
        fields = ('name', 'description', 'project', 'plasmids')
        widgets = {
            'description': forms.Textarea(attrs={'rows': 3}),
            'plasmids': forms.SelectMultiple(attrs={'size': 12}),
        }

    def __init__(self, *args, **kwargs):
        user = kwargs.pop('user', None)
        super().__init__(*args, **kwargs)

        self.fields['plasmids'].label = 'Target plasmids'
        self.fields['plasmids'].help_text = (
            'Select only the final plasmids. Their backbone and inserts are included automatically.'
        )

        writable_projects = (
            get_projects_where_member_can(user, ['a', 'w']).order_by('name')
            if user else Project.objects.none()
        )
        self.fields['project'].queryset = writable_projects
        self.fields['project'].required = False
        self.fields['project'].empty_label = 'Unassigned project'
        visible_projects = get_projects_where_member_can_any(user) if user else Project.objects.none()
        plasmid_queryset = Plasmid.objects.filter(
            project__in=visible_projects
        ).select_related('project').prefetch_related('inserts').order_by(
            'project__name', 'name', 'idx'
        )
        self.fields['plasmids'].queryset = plasmid_queryset

        plasmids = list(plasmid_queryset)
        self.has_level_3 = any(plasmid.level == 3 for plasmid in plasmids)
        self.has_level_4 = any(plasmid.level == 4 for plasmid in plasmids)
        self.has_unidentified_level = any(plasmid.level is None for plasmid in plasmids)
        self.plasmid_dependencies = {}
        self.plasmid_labels = {}
        self.plasmid_levels = {}
        self.plasmid_urls = {}
        for plasmid in plasmids:
            dependencies = list(plasmid.inserts.all())
            if plasmid.backbone:
                dependencies.insert(0, plasmid.backbone)
            self.plasmid_dependencies[str(plasmid.id)] = [str(item.id) for item in dependencies]
            for item in [plasmid, *dependencies]:
                self.plasmid_labels[str(item.id)] = f'{item.name} [{item.idx}]'
                self.plasmid_levels[str(item.id)] = (
                    f'L{item.level}' if item.level is not None else 'Unidentified Level'
                )
                self.plasmid_urls[str(item.id)] = reverse(
                    'plasmid', kwargs={'plasmid_id': item.id}
                )

        for field in self.fields.values():
            field.widget.attrs.setdefault('class', 'form-control')
        self.fields['project'].widget.attrs['class'] = 'form-select'
        self.fields['plasmids'].widget.attrs.update({
            'class': 'form-select',
            'data-plasmid-selector': 'true',
        })

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
                  'description', 'project', 'ligation_state', 'created_on']


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
                  'description', 'created_on', 'project', 'ligation_state', 'reference_sequence']


class PlasmidValidationForm(forms.ModelForm):
    def clean(self):
        cleaned_data = super().clean()
        if cleaned_data.get('no_colony'):
            cleaned_data['working_colony'] = None
        return cleaned_data

    class Meta:
        model = Plasmid
        fields = ['ligation_state', 'working_colony',
                  'no_colony',
                  'colonypcr_state', 'colonypcr_date', 'colonypcr_observations',
                  'digestion_state', 'digestion_date', 'digestion_observations',
                  'sequencing_state', 'sequencing_date', 'sequencing_observations']
        widgets = {
            'colonypcr_date': DateInput(),
            'digestion_date': DateInput(),
            'sequencing_date': DateInput()
        }
        labels = {
            'colonypcr_observations': 'cPCR description',
            'digestion_observations': 'Digestion description',
            'sequencing_observations': 'Sequencing description',
        }
        help_texts = {
            'colonypcr_observations': 'Primers used, phenotype, or other interpretation.',
            'digestion_observations': 'Enzymes used and expected or observed band pattern.',
            'sequencing_observations': 'Reads or primers used, coverage, or interpretation.',
        }


class DigestForm(forms.Form):
    enzymes = forms.CharField(widget=forms.HiddenInput(attrs={'id': 'digest_enzymes'}))


class PCRForm(forms.Form):
    def __init__(self, *args, **kwargs):
        user = kwargs.pop('user', None)
        super(PCRForm, self).__init__(*args, **kwargs)
        if user is not None:
            queryset = visible_primers_for_user(user).order_by('name')
            self.fields['primer_f'].queryset = queryset
            self.fields['primer_r'].queryset = queryset

    primer_f = forms.ModelChoiceField(queryset=Primer.objects.all(), to_field_name="id", label="Primer F",
                                      required=False)
    primer_r = forms.ModelChoiceField(queryset=Primer.objects.all(), to_field_name="id", label="Primer R",
                                      required=False)
    primer_f_seq = forms.CharField(label="Primer F sequence", required=False)
    primer_r_seq = forms.CharField(label="Primer R sequence", required=False)


class PrimerBatchUploadForm(forms.Form):
    def __init__(self, *args, **kwargs):
        super(PrimerBatchUploadForm, self).__init__(*args, **kwargs)
        self.fields['fasta_file'].widget.attrs.update({'class': 'form-control', 'accept': '.fa,.fas,.fasta,.txt'})
        self.fields['name_source'].widget.attrs.update({'class': 'form-select'})
        self.fields['default_direction'].widget.attrs.update({'class': 'form-select'})
        self.fields['update_existing'].widget.attrs.update({'class': 'form-check-input'})

    fasta_file = forms.FileField(label="FASTA file")
    update_existing = forms.BooleanField(label="Update existing primers with the same name", required=False)
    name_source = forms.ChoiceField(
        label="Primer name source",
        choices=(("id", "FASTA ID"), ("description", "Full FASTA description")),
        initial="id",
    )
    default_direction = forms.ChoiceField(
        label="Direction when F/R is not in the name",
        choices=(("f", "Forward"), ("r", "Reverse"), ("", "Leave empty")),
        initial="f",
        required=False,
    )


class RestrictionEnzymeCreateForm(forms.ModelForm):
    def __init__(self, *args, **kwargs):
        super(RestrictionEnzymeCreateForm, self).__init__(*args, **kwargs)
        self.fields['name'].widget.attrs.update({'class': 'form-select'})
        self.fields['hf_version'].widget.attrs.update({'class': 'form-check-input'})
        self.fields['link_datasheet'].widget.attrs.update({
            'class': 'form-control',
            'placeholder': 'Optional datasheet URL',
        })
        self.fields['description'].widget.attrs.update({
            'class': 'form-control',
            'placeholder': 'Optional lab note',
        })

    class Meta:
        model = RestrictionEnzyme
        fields = [
            'name',
            'hf_version',
            'link_datasheet',
            'description',
        ]

    def clean(self):
        cleaned_data = super().clean()
        name = cleaned_data.get("name")
        hf_version = bool(cleaned_data.get("hf_version"))
        duplicate_qs = RestrictionEnzyme.objects.filter(name=name, hf_version=hf_version)
        if self.instance.pk:
            duplicate_qs = duplicate_qs.exclude(pk=self.instance.pk)
        if name and duplicate_qs.exists():
            display_name = f"{name}-HF" if hf_version else name
            self.add_error("name", f"{display_name} is already loaded in Weaver.")
        return cleaned_data


class RestrictionEnzymeBufferForm(forms.Form):
    NEW_BUFFER_CHOICE = "__new__"

    existing_buffer = forms.ChoiceField(
        choices=(),
        required=False,
        label="Buffer",
    )
    new_buffer_name = forms.CharField(
        required=False,
        max_length=100,
        label="New buffer",
    )
    activity_percent = forms.IntegerField(
        required=False,
        min_value=0,
        max_value=100,
        label="Activity %",
    )
    DELETE = forms.BooleanField(required=False, label="Remove")

    def __init__(self, *args, **kwargs):
        super(RestrictionEnzymeBufferForm, self).__init__(*args, **kwargs)
        self.fields['existing_buffer'].choices = [
            ("", "Select a buffer"),
            *[
                (str(buffer.id), buffer.name)
                for buffer in RestrictionBuffer.objects.order_by('name')
            ],
            (self.NEW_BUFFER_CHOICE, "+ New buffer..."),
        ]
        self.fields['existing_buffer'].widget.attrs.update({
            'class': 'form-select buffer-choice-select',
            'data-new-buffer-value': self.NEW_BUFFER_CHOICE,
        })
        self.fields['new_buffer_name'].widget.attrs.update({
            'class': 'form-control buffer-new-name',
            'placeholder': 'e.g. NEB CutSmart',
        })
        self.fields['activity_percent'].widget.attrs.update({
            'class': 'form-control',
            'min': 0,
            'max': 100,
            'placeholder': '0-100',
        })
        self.fields['DELETE'].widget.attrs.update({'class': 'form-check-input'})

    def clean(self):
        cleaned_data = super().clean()
        buffer_choice = cleaned_data.get("existing_buffer")
        new_buffer_name = (cleaned_data.get("new_buffer_name") or "").strip()
        activity_percent = cleaned_data.get("activity_percent")
        delete_row = cleaned_data.get("DELETE")

        existing_buffer = None
        if buffer_choice and buffer_choice != self.NEW_BUFFER_CHOICE:
            existing_buffer = RestrictionBuffer.objects.filter(pk=buffer_choice).first()
            if existing_buffer is None:
                self.add_error("existing_buffer", "Choose a valid buffer.")

        is_empty = not existing_buffer and not new_buffer_name and activity_percent in (None, "")
        cleaned_data["new_buffer_name"] = new_buffer_name
        cleaned_data["existing_buffer"] = existing_buffer
        cleaned_data["buffer_choice"] = buffer_choice
        cleaned_data["is_empty"] = is_empty

        if is_empty or delete_row:
            return cleaned_data

        if buffer_choice == self.NEW_BUFFER_CHOICE:
            if not new_buffer_name:
                self.add_error("new_buffer_name", "Enter the name for the new buffer.")
        else:
            cleaned_data["new_buffer_name"] = ""
            if not existing_buffer:
                raise forms.ValidationError("Choose a buffer or create a new one.")
        if activity_percent is None:
            self.add_error("activity_percent", "Enter the activity percentage for this buffer.")

        return cleaned_data


class BaseRestrictionEnzymeBufferFormSet(forms.BaseFormSet):
    def clean(self):
        if any(self.errors):
            return

        seen_buffers = set()
        for form in self.forms:
            cleaned_data = getattr(form, "cleaned_data", None)
            if not cleaned_data or cleaned_data.get("DELETE") or cleaned_data.get("is_empty"):
                continue
            existing_buffer = cleaned_data.get("existing_buffer")
            buffer_name = existing_buffer.name if existing_buffer else cleaned_data.get("new_buffer_name", "")
            normalized_name = buffer_name.strip().lower()
            if normalized_name in seen_buffers:
                raise forms.ValidationError("Each buffer can only be added once per enzyme.")
            seen_buffers.add(normalized_name)


RestrictionEnzymeBufferFormSet = forms.formset_factory(
    RestrictionEnzymeBufferForm,
    formset=BaseRestrictionEnzymeBufferFormSet,
    extra=1,
    can_delete=True,
)


def save_restriction_enzyme_buffers(restriction_enzyme, buffer_formset):
    for cleaned_data in buffer_formset.cleaned_data:
        if not cleaned_data or cleaned_data.get("DELETE") or cleaned_data.get("is_empty"):
            continue

        buffer = cleaned_data.get("existing_buffer")
        if buffer is None:
            buffer_name = cleaned_data["new_buffer_name"]
            buffer = RestrictionBuffer.objects.filter(name__iexact=buffer_name).first()
            if buffer is None:
                buffer = RestrictionBuffer.objects.create(name=buffer_name)

        RestrictionEnzymeBuffer.objects.create(
            restriction_enzyme=restriction_enzyme,
            buffer=buffer,
            activity_percent=cleaned_data["activity_percent"],
        )


class GenBankBatchUploadForm(forms.Form):
    PROJECT_TARGET_CHOICES = (
        ("existing", "Upload to an existing project"),
        ("new", "Create a new project and upload there"),
    )

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user", None)
        current_project = kwargs.pop("current_project", None)
        super(GenBankBatchUploadForm, self).__init__(*args, **kwargs)
        writable_projects = get_projects_where_member_can(user, ['a', 'w']).order_by('name') if user else Project.objects.none()
        self.fields['target_mode'].widget.attrs.update({'class': 'form-select'})
        self.fields['project'].widget.attrs.update({'class': 'form-select'})
        self.fields['genbank_files'].widget.attrs.update({
            'class': 'form-control',
            'accept': '.gb,.gbk,.genbank',
        })
        self.fields['name_source'].widget.attrs.update({'class': 'form-select'})
        self.fields['new_project_name'].widget.attrs.update({'class': 'form-control'})
        self.fields['new_project_assembly_standard'].widget.attrs.update({'class': 'form-select'})
        self.fields['update_existing'].widget.attrs.update({'class': 'form-check-input'})
        self.fields['public_visibility'].widget.attrs.update({'class': 'form-check-input'})
        self.fields['reference_sequence'].widget.attrs.update({'class': 'form-check-input'})
        self.fields['infer_ytk_metadata'].widget.attrs.update({'class': 'form-check-input'})
        self.fields['new_project_public'].widget.attrs.update({'class': 'form-check-input'})

        self.fields['project'].choices = [
            (str(project.id), project.name)
            for project in writable_projects
        ]

        if current_project and any(str(current_project.id) == choice[0] for choice in self.fields['project'].choices):
            self.fields['project'].initial = str(current_project.id)
            self.fields['target_mode'].initial = "existing"
        elif self.fields['project'].choices:
            self.fields['project'].initial = self.fields['project'].choices[0][0]
            self.fields['target_mode'].initial = "existing"
        else:
            self.fields['target_mode'].initial = "new"

        self.fields['new_project_assembly_standard'].initial = (
            current_project.assembly_standard
            if current_project and current_project.assembly_standard
            else "ytk"
        )

    target_mode = forms.ChoiceField(
        label="Upload destination",
        choices=PROJECT_TARGET_CHOICES,
        initial="existing",
    )
    project = forms.ChoiceField(label="Existing project", required=False, choices=())
    genbank_files = MultipleFileField(
        label="GenBank files",
        required=False,
        widget=MultipleFileInput(attrs={"multiple": True}),
    )
    new_project_name = forms.CharField(label="New project name", required=False, max_length=128)
    new_project_public = forms.BooleanField(label="Make the new project public", required=False)
    new_project_assembly_standard = forms.ChoiceField(
        label="Assembly standard",
        choices=tuple((index, assembly_standard['name']) for index, assembly_standard in assembly_standards.items()),
        required=False,
    )
    update_existing = forms.BooleanField(label="Update existing plasmids with the same name", required=False)
    public_visibility = forms.BooleanField(label="Mark imported plasmids as public", required=False)
    reference_sequence = forms.BooleanField(label="Import as reference sequences", required=False, initial=False)
    infer_ytk_metadata = forms.BooleanField(
        label="Infer assembly type, level, and resistance",
        required=False,
        initial=True,
    )
    name_source = forms.ChoiceField(
        label="Plasmid name source",
        choices=(("filename", "Filename stem"), ("record", "GenBank LOCUS / record name")),
        initial="filename",
    )

    def clean(self):
        cleaned_data = super().clean()
        files = self.files.getlist("genbank_files")
        if not files:
            raise forms.ValidationError("Upload at least one .gb, .gbk, or .genbank file.")
        if cleaned_data.get("target_mode") == "existing":
            if not cleaned_data.get("project"):
                self.add_error("project", "Choose the destination project.")
        else:
            project_name = str(cleaned_data.get("new_project_name") or "").strip()
            if not project_name:
                self.add_error("new_project_name", "Provide a name for the new project.")
            elif Project.objects.filter(name__iexact=project_name).exists():
                self.add_error("new_project_name", "A project with this name already exists.")
            cleaned_data["new_project_name"] = project_name
        return cleaned_data


class ServicesPCRForm(forms.Form):
    def __init__(self, *args, **kwargs):
        user = kwargs.pop('user', None)
        super(ServicesPCRForm, self).__init__(*args, **kwargs)
        if user is not None:
            primers = visible_primers_for_user(user).order_by('name')
            self.fields['primer_f'].queryset = primers.filter(fwd_or_rev='f')
            self.fields['primer_r'].queryset = primers.filter(fwd_or_rev='r')
        for field in self.fields.values():
            field.widget.attrs.update({'class': 'form-control'})
        self.fields['primer_f'].widget.attrs.update({'class': 'form-select'})
        self.fields['primer_r'].widget.attrs.update({'class': 'form-select'})

    primer_f = forms.ModelChoiceField(queryset=Primer.objects.none(), to_field_name="id", label="Forward primer")
    primer_r = forms.ModelChoiceField(queryset=Primer.objects.none(), to_field_name="id", label="Reverse primer")
    min_product_size = forms.IntegerField(label="Minimum product size", min_value=1, initial=100)
    max_product_size = forms.IntegerField(label="Maximum product size", min_value=1, required=False)

    def clean(self):
        cleaned_data = super().clean()
        min_product_size = cleaned_data.get('min_product_size')
        max_product_size = cleaned_data.get('max_product_size')
        if min_product_size and max_product_size and max_product_size < min_product_size:
            raise forms.ValidationError("Maximum product size must be greater than or equal to minimum product size.")
        return cleaned_data


class BatchPrintsForm(forms.Form):
    LABEL_TYPES = (
        ("plasmids", "Plasmids"),
        ("glycerolstocks", "Glycerol stocks"),
    )

    def __init__(self, *args, **kwargs):
        super(BatchPrintsForm, self).__init__(*args, **kwargs)
        self.fields["label_type"].widget.attrs.update({"class": "form-select"})
        self.fields["identifiers"].widget.attrs.update({
            "class": "form-control",
            "rows": 8,
            "placeholder": "One ID, code, QR, or name per line",
        })
        self.fields["date"].widget.attrs.update({"class": "form-control"})
        self.fields["concentration"].widget.attrs.update({"class": "form-control"})

    label_type = forms.ChoiceField(label="Label type", choices=LABEL_TYPES)
    identifiers = forms.CharField(label="Items", widget=forms.Textarea)
    date = forms.DateField(widget=DateInput(), initial=datetime.date.today)
    concentration = forms.FloatField(label="Concentration (ng/ul)", min_value=0.01, initial=1)


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
