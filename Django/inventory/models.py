import os
import uuid
from django.db import models
import datetime
from django.conf import settings
from django.core.validators import MaxValueValidator
from django.core.validators import MinValueValidator
from shortuuidfield import ShortUUIDField
import shortuuid
from .custom.box import BOX_ROWS
from .custom.box import BOX_COLUMNS
from .custom.general import FWD_OR_REV
from .custom.general import CHECK_STATES
from .custom.general import LIGATION_STATES
from .custom.general import COLORS
from Bio.Restriction.Restriction_Dictionary import rest_dict, suppliers
from organization.models import Project
from django.dispatch import receiver
from .custom.standards import assembly_standards


from .validators import clustal_validate

RE_Choices = []
for key in rest_dict:
    if not key.startswith("_"):
        RE_Choices.append((key, key))

LEGACY_RESTRICTION_BUFFERS = (
    ("activity_buffer_1_1", "NEB 1.1"),
    ("activity_buffer_2_1", "NEB 2.1"),
    ("activity_buffer_3_1", "NEB 3.1"),
    ("activity_buffer_CS", "NEB CutSmart"),
    ("activity_buffer_aari", "Thermo AarI"),
)


def generate_shortuuid():
    return shortuuid.uuid()


class Resistance(models.Model):
    id = models.AutoField(primary_key=True)
    three_letter_code = models.CharField(max_length=3, blank=True)
    name = models.CharField(max_length=20, blank=True)

    def __str__(self):
        return self.name


class TableFilter(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=20, blank=True)
    color = models.CharField(choices=COLORS, max_length=10, blank=True, null=True)
    options = models.CharField(max_length=200, blank=True, help_text="Format: \'x|X,y|Y\'. Left side is the name of the filter and right side is the start-with text")

    def __str__(self):
        return self.name


class PlasmidType(models.Model):
    id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=10, blank=True)

    def __str__(self):
        return self.name


class RestrictionBuffer(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=100, unique=True)

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['name']


class RestrictionEnzyme(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    # List https://github.com/biopython/biopython/blob/master/Bio/Restriction/Restriction_Dictionary.py
    name = models.CharField(choices=RE_Choices, max_length=20)
    buffers = models.ManyToManyField(
        RestrictionBuffer,
        through='RestrictionEnzymeBuffer',
        blank=True,
        related_name='restriction_enzymes',
    )
    activity_buffer_1_1 = models.IntegerField(blank=True, null=True)
    activity_buffer_2_1 = models.IntegerField(blank=True, null=True)
    activity_buffer_3_1 = models.IntegerField(blank=True, null=True)
    activity_buffer_CS = models.IntegerField(blank=True, null=True)
    activity_buffer_aari = models.IntegerField(blank=True, null=True)
    hf_version = models.BooleanField(blank=True)
    link_datasheet = models.CharField(max_length=200, blank=True)
    description = models.CharField(max_length=200, blank=True)

    @property
    def max_activity_temperature(self):
        bio_python_re = rest_dict[self.name]
        if bio_python_re:
            return bio_python_re['opt_temp']
        return None

    @property
    def inactivation_temperature(self):
        bio_python_re = rest_dict[self.name]
        if bio_python_re:
            return bio_python_re['inact_temp']
        return None

    @property
    def fcut(self):
        bio_python_re = rest_dict[self.name]
        if bio_python_re:
            return bio_python_re['fst5']
        return None

    @property
    def rcut(self):
        bio_python_re = rest_dict[self.name]
        if bio_python_re:
            return bio_python_re['size'] + bio_python_re['fst3']
        return None

    @property
    def recognition_site(self):
        bio_python_re = rest_dict[self.name]
        if bio_python_re:
            return bio_python_re['site']
        return None

    @property
    def suppliers(self):
        bio_python_re = rest_dict[self.name]
        if bio_python_re:
            suppliers_list = []
            for sup_code in bio_python_re['suppl']:
                suppliers_list.append(suppliers[sup_code][0])
            return suppliers_list
        return None

    @property
    def buffer_activities(self):
        prefetched = getattr(self, '_prefetched_objects_cache', {})
        if 'buffer_links' in prefetched:
            links = prefetched['buffer_links']
        else:
            links = list(self.buffer_links.select_related('buffer').all())
        return sorted(links, key=lambda link: link.buffer.name.lower())

    @property
    def buffer_activity_entries(self):
        if self.buffer_activities:
            return [
                {
                    'name': link.buffer.name,
                    'activity_percent': link.activity_percent,
                }
                for link in self.buffer_activities
            ]
        legacy_entries = []
        for legacy_field, buffer_name in LEGACY_RESTRICTION_BUFFERS:
            activity = getattr(self, legacy_field, None)
            if activity is not None:
                legacy_entries.append({
                    'name': buffer_name,
                    'activity_percent': activity,
                })
        return legacy_entries

    def buffer_activity_map(self):
        return {
            entry['name']: entry['activity_percent']
            for entry in self.buffer_activity_entries
        }

    def __str__(self):
        if self.hf_version:
            return self.name + "-HF"
        else:
            return self.name

    class Meta:
        ordering = ['name']


class RestrictionEnzymeBuffer(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    restriction_enzyme = models.ForeignKey(
        RestrictionEnzyme,
        on_delete=models.CASCADE,
        related_name='buffer_links',
    )
    buffer = models.ForeignKey(
        RestrictionBuffer,
        on_delete=models.CASCADE,
        related_name='enzyme_links',
    )
    activity_percent = models.PositiveSmallIntegerField(
        validators=[MinValueValidator(0), MaxValueValidator(100)],
    )

    def __str__(self):
        return f"{self.restriction_enzyme} / {self.buffer} / {self.activity_percent}%"

    class Meta:
        ordering = ['buffer__name']
        constraints = [
            models.UniqueConstraint(
                fields=['restriction_enzyme', 'buffer'],
                name='unique_restriction_enzyme_buffer',
            ),
        ]


class Location(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=200, blank=True)

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['name']


class Box(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=200, blank=True)
    location = models.ForeignKey(Location, on_delete=models.CASCADE, blank=True)

    def __str__(self):
        return self.name + " / " + self.location.__str__()

    class Meta:
        ordering = ['name']


class Plasmid(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    idx = models.IntegerField(blank=True, null=True, editable=False)
    name = models.CharField(max_length=50)
    selectable_markers = models.ManyToManyField(Resistance, blank=True, symmetrical=False, related_name='+', help_text='Use CTRL for multiple select')
    sequence = models.FileField(upload_to='uploads/plasmids/', blank=True)
    backbone = models.ForeignKey('self', models.SET_NULL, blank=True, null=True)
    computed_size = models.IntegerField(blank=True, null=True, editable=False)
    insert_computed_size = models.IntegerField(blank=True, null=True, editable=False)
    inserts = models.ManyToManyField('self', blank=True, symmetrical=False, related_name='+')
    intended_use = models.CharField(max_length=200)
    type = models.ForeignKey(PlasmidType, models.SET_NULL, blank=True, null=True)
    level = models.IntegerField(blank=True, null=True, choices=(
        (None, 'Not defined'),
        (0, 'Level 0'),
        (1, 'Level 1'),
        (2, 'Level 2'),
        (3, 'Level 3'),
        (4, 'Level 4')
    ))
    description = models.CharField(max_length=1000, blank=True, help_text="Allows markdown")
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    created_on = models.DateField(auto_now_add=False, default=datetime.date.today)
    assembly_metadata = models.JSONField(default=dict, blank=True)

    reference_sequence = models.BooleanField(blank=True, default=0)
    public_visibility = models.BooleanField(blank=True, default=0)

    ligation_state = models.IntegerField(choices=LIGATION_STATES, default=2)

    qr_id = ShortUUIDField(default=generate_shortuuid, editable=False)

    # validation

    working_colony = models.IntegerField(blank=True, null=True)
    no_colony = models.BooleanField(
        default=False,
        help_text="Use when the construct is used directly or comes from synthesis.",
    )

    colonypcr_state = models.IntegerField(choices=CHECK_STATES, blank=True, default=1)
    colonypcr_observations = models.CharField(max_length=1000, blank=True, null=True)
    colonypcr_date = models.DateField(blank=True, null=True)

    digestion_state = models.IntegerField(choices=CHECK_STATES, blank=True, default=1)
    digestion_observations = models.CharField(max_length=1000, blank=True, null=True)
    digestion_date = models.DateField(blank=True, null=True)

    sequencing_state = models.IntegerField(choices=CHECK_STATES, blank=True, default=0)
    sequencing_date = models.DateField(blank=True, null=True)
    sequencing_observations = models.CharField(max_length=1000, blank=True, null=True)

    sequencing_clustal_file = models.FileField(upload_to="uploads/sequencing_clustal", blank=True, null=True, validators=[clustal_validate])

    def __str__(self):
        return self.name + " | " + str(self.idx)

    class Meta:
        ordering = ['name']

    def working_colony_text_short(self):
        if self.reference_sequence:
            return "RS"
        elif self.no_colony:
            return "NC"
        elif self.ligation_state != 1:
            return "UC"
        elif self.is_validated():
            if self.working_colony:
                return "c"+str(self.working_colony) + "-V"
            else:
                return "NS"
        else:
            if self.working_colony:
                return "c"+str(self.working_colony) + "-NV"
            else:
                return "NS"

    def working_colony_text(self):
        if self.reference_sequence:
            return "Reference sequence"
        elif self.no_colony:
            return "No colony"
        elif self.ligation_state != 1:
            return "Under construction"
        elif self.is_validated():
            if self.working_colony:
                return str(self.working_colony) + " (Validated)"
            else:
                return "Not set"
        else:
            if self.working_colony:
                return str(self.working_colony) + " (Not validated)"
            else:
                return "Not set"

    def colony_source_text(self):
        """Return the explicit colony status without validation state details."""
        if self.no_colony:
            return "No colony"
        if self.working_colony is not None:
            return "#" + str(self.working_colony)
        return "Not set"

    def is_validated(self):
        if self.reference_sequence:
            return None
        if self.ligation_state != 1:
            return False
        if self.colonypcr_state != 1 and self.digestion_state != 1 and self.sequencing_state != 1:
            return True
        return False

    def get_check_state(self):
        if self.reference_sequence:
            return 'r'
        if self.is_validated():
            return 'v'
        return 'c'

    def save(self, *args, **kwargs):
        if not self.idx:
            try:
                last_plasmid_idx = Plasmid.objects.order_by("idx").last().idx
                if last_plasmid_idx:
                    self.idx = last_plasmid_idx + 1
                else:
                    self.idx = 1
            except:
                self.idx = 1
        super(Plasmid, self).save(*args, **kwargs)

    def get_insert_of(self):
        return Plasmid.objects.filter(inserts__in=[self])

    def get_backbone_of(self):
        return Plasmid.objects.filter(backbone=self)

    @property
    def detected_assembly(self):
        return (self.assembly_metadata or {}).get("detected", {})

    @property
    def confirmed_assembly(self):
        return (self.assembly_metadata or {}).get("confirmed", {})

    @property
    def assembly_detection_confidence_percent(self):
        confidence = self.detected_assembly.get("confidence")
        if confidence is None:
            return None
        return round(float(confidence) * 100, 1)

    def ligation_concentration(self, units=True):
        if self.computed_size:
            if self.type:
                if str(self.type) == "Insert":
                    result = str(round(self.computed_size / 100, 1))
                    if units:
                        result += " ng / ul"
                    return result
                elif str(self.type) == "Receiver":
                    result = str(round(self.computed_size / 300, 1))
                    if units:
                        result += " ng / ul"
                    return result
                else:
                    return "Plasmid type no formula"
            else:
                return "No plasmid type set"
        else:
            return "No plasmid computed size"

    def ligation_concentration_no_units(self):
        return self.ligation_concentration(units=False)

    def recommended_enzyme_for_create(self, return_name=False):
        try:
            if return_name:
                re = RestrictionEnzyme.objects.filter(
                    name__iexact=assembly_standards[self.project.assembly_standard]['enzymes'][self.level]
                ).order_by('hf_version', 'id').first()
                return re.name
            else:
                return assembly_standards[self.project.assembly_standard]['enzymes'][self.level]
        except:
            return "No level set"

    def getPlasmidResistanceForLigation(self):
        if self.selectable_markers.count() == 1:
            return self.selectable_markers.all()[0].three_letter_code
        elif self.selectable_markers.count() == 0:
            return 'No selectable marker set'
        else:
            res_txt=[]
            for res in self.selectable_markers.all():
                res_txt.append(res.three_letter_code)
            return 'More than one selectable marker set: ' + ' / '.join(res_txt)

    def ligation_raw(self):
        tab = "	"
        ligation_raw = self.__str__() + tab
        if self.backbone:
            ligation_raw += self.backbone.__str__() + " [" + self.backbone.working_colony_text_short() + "]"

        ligation_raw += tab

        inserts = []
        for plasmid in self.inserts.all():
            inserts.append(plasmid.__str__() + " [" + plasmid.working_colony_text_short() + "]")

        if self.level:
            ligation_raw = ligation_raw + " + ".join(inserts) + tab + tab +\
                           self.recommended_enzyme_for_create(return_name=True) + tab +\
                           self.getPlasmidResistanceForLigation().upper()
        else:
            if self.level == 0:
                ligation_raw = "Level 0 ligation is not supported"
            else:
                ligation_raw = "Level not set"

        return ligation_raw

    def dependencies_validated(self):
        if self.backbone:
            if not self.backbone.is_validated():
                return False
        if self.dependencies_validated:
            for pi in self.inserts.all():
                if not pi.is_validated():
                    return False
        return True


@receiver(models.signals.post_delete, sender=Plasmid)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    """
    Deletes file from filesystem
    when corresponding `Plasmid` object is deleted.
    """
    try:
        if instance.sequencing_clustal_file:
            if os.path.isfile(instance.sequencing_clustal_file.path):
                os.remove(instance.sequencing_clustal_file.path)
        if instance.sequence:
            if os.path.isfile(instance.sequence.path):
                os.remove(instance.sequence.path)
    except:
        return False


@receiver(models.signals.pre_save, sender=Plasmid)
def auto_delete_file_on_change(sender, instance, **kwargs):
    """
    Deletes old file from filesystem
    when corresponding `Plasmid` object is updated
    with new file.
    """
    if not instance.pk:
        return False

    try:
        old_file_sequencing_clustal_file = Plasmid.objects.get(pk=instance.pk).sequencing_clustal_file
        old_file_sequence = Plasmid.objects.get(pk=instance.pk).sequence
    except Plasmid.DoesNotExist:
        return False

    try:
        new_file_sequencing_clustal_file = instance.sequencing_clustal_file
        new_file_sequence = instance.sequence

        if not old_file_sequencing_clustal_file == new_file_sequencing_clustal_file:
            if os.path.isfile(old_file_sequencing_clustal_file.path):
                os.remove(old_file_sequencing_clustal_file.path)

        if not old_file_sequence == new_file_sequence:
            if os.path.isfile(old_file_sequence.path):
                os.remove(old_file_sequence.path)
    except:
        return False


SANGER_AUTOMATED_STATES = (
    ("PASS", "Verifica"),
    ("REVIEW", "Requiere revisión"),
    ("FAIL", "No verifica"),
    ("NO_DATA", "Sin datos utilizables"),
)

SANGER_MANUAL_DECISIONS = (
    ("", "Pending"),
    ("VERIFIED", "Verified"),
    ("REJECTED", "Not verified"),
    ("INCONCLUSIVE", "Inconclusive"),
    ("PENDING", "Pending"),
)

SANGER_READ_ORIENTATIONS = (
    ("forward", "Forward"),
    ("reverse", "Reverse complement"),
    ("ambiguous", "Ambiguous"),
    ("unmapped", "Unmapped"),
)

SANGER_FILE_FORMATS = (
    ("ab1", "AB1"),
    ("phd1", "PHD.1"),
    ("seq", "SEQ"),
)


def sanger_read_file_upload_to(instance, filename):
    return "uploads/sanger/{}/{}/{}".format(instance.read.run.plasmid_id, instance.read.run_id, filename)


class SangerVerificationRun(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    plasmid = models.ForeignKey(Plasmid, on_delete=models.CASCADE, related_name="sanger_verification_runs")
    created_by = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.SET_NULL, blank=True, null=True, related_name="sanger_verification_runs")
    created_at = models.DateTimeField(auto_now_add=True)
    sequencing_datetime = models.DateTimeField(blank=True, null=True)
    label = models.CharField(max_length=200, blank=True)
    colony = models.CharField(max_length=100, blank=True)
    sample = models.CharField(max_length=100, blank=True)
    notes = models.TextField(blank=True)
    parameters = models.JSONField(default=dict, blank=True)
    automated_state = models.CharField(max_length=20, choices=SANGER_AUTOMATED_STATES, default="NO_DATA")
    automated_reasons = models.JSONField(default=list, blank=True)
    combined_metrics = models.JSONField(default=dict, blank=True)
    manual_decision = models.CharField(max_length=20, choices=SANGER_MANUAL_DECISIONS, blank=True, default="")
    manual_decision_by = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.SET_NULL, blank=True, null=True, related_name="sanger_manual_decisions")
    manual_decision_at = models.DateTimeField(blank=True, null=True)
    manual_decision_effective_date = models.DateField(blank=True, null=True)
    manual_decision_comment = models.TextField(blank=True)
    clustal_file = models.FileField(upload_to="uploads/sequencing_clustal", blank=True, null=True, max_length=500, validators=[clustal_validate])

    class Meta:
        ordering = ["-created_at"]

    def __str__(self):
        when = self.created_at.strftime("%Y-%m-%d %H:%M") if self.created_at else ""
        return "{} Sanger {}".format(self.plasmid.name, when)


class SangerRead(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    run = models.ForeignKey(SangerVerificationRun, on_delete=models.CASCADE, related_name="reads")
    name = models.CharField(max_length=255)
    detected_orientation = models.CharField(max_length=20, choices=SANGER_READ_ORIENTATIONS, default="unmapped")
    forced_orientation = models.CharField(max_length=20, choices=SANGER_READ_ORIENTATIONS, blank=True, default="")
    raw_sequence = models.TextField(blank=True)
    trimmed_sequence = models.TextField(blank=True)
    selected_source = models.CharField(max_length=20, blank=True)
    parsing_result = models.JSONField(default=dict, blank=True)
    quality_metrics = models.JSONField(default=dict, blank=True)
    alignment_metrics = models.JSONField(default=dict, blank=True)
    warnings = models.JSONField(default=list, blank=True)
    is_usable = models.BooleanField(default=False)

    class Meta:
        ordering = ["name"]

    def __str__(self):
        return self.name


class SangerReadFile(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    read = models.ForeignKey(SangerRead, on_delete=models.CASCADE, related_name="files")
    primer = models.ForeignKey("Primer", on_delete=models.SET_NULL, blank=True, null=True, related_name="sanger_read_files")
    format = models.CharField(max_length=10, choices=SANGER_FILE_FORMATS)
    original_name = models.CharField(max_length=255)
    file = models.FileField(upload_to=sanger_read_file_upload_to, blank=True, max_length=500)
    sha256 = models.CharField(max_length=64)
    size = models.PositiveIntegerField(default=0)
    metadata = models.JSONField(default=dict, blank=True)

    class Meta:
        ordering = ["format", "original_name"]


class SangerVariant(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    run = models.ForeignKey(SangerVerificationRun, on_delete=models.CASCADE, related_name="variants")
    read = models.ForeignKey(SangerRead, on_delete=models.CASCADE, related_name="variants", blank=True, null=True)
    coordinate = models.PositiveIntegerField(help_text="0-based plasmid coordinate")
    variant_type = models.CharField(max_length=20)
    expected_base = models.CharField(max_length=20, blank=True)
    observed_base = models.CharField(max_length=20, blank=True)
    quality = models.IntegerField(blank=True, null=True)
    evidence = models.JSONField(default=dict, blank=True)
    flags = models.JSONField(default=list, blank=True)

    class Meta:
        ordering = ["coordinate", "variant_type"]


class Strain(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=200, blank=True)
    selectable_markers = models.ManyToManyField(Resistance, blank=True, symmetrical=False, related_name='+', help_text='Use CTRL for multiple select')
    description = models.CharField(max_length=1000, blank=True)
    for_primary_gs = models.BooleanField(default=False)

    def __str__(self):
        result = self.name
        if self.for_primary_gs:
            result += ' (p)'
        return result

    class Meta:
        ordering = ['name']


class GlycerolStock(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    strain = models.ForeignKey(Strain, on_delete=models.CASCADE, blank=True)
    parent = models.ForeignKey("self", on_delete=models.CASCADE, blank=True, null=True)
    plasmid = models.ForeignKey(Plasmid, on_delete=models.CASCADE, blank=True, null=True)
    created_on = models.DateField(auto_now_add=False, default=datetime.date.today)
    box_row = models.CharField(max_length=1, choices=BOX_ROWS, help_text="Click box position (below) to autocomplete")
    box_column = models.IntegerField(choices=BOX_COLUMNS, help_text="Click box position (below) to autocomplete")
    box = models.ForeignKey(Box, on_delete=models.CASCADE, help_text="Click box position (below) to autocomplete")
    qr_id = ShortUUIDField(default=generate_shortuuid, editable=False)
    details = models.CharField(max_length=1000, blank=True)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)

    def __str__(self):
        if self.plasmid is None:
            return self.strain.name
        else:
            return self.strain.name + " / " + self.plasmid.__str__()

    class Meta:
        ordering = ['strain', 'plasmid']


class Primer(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=100, blank=True)
    sequence_3 = models.CharField(max_length=200, blank=True, verbose_name="Sequence (3' end)", help_text="5' → 3' direction")
    sequence_5 = models.CharField(max_length=200, blank=True, verbose_name="Sequence (5' end / overhang)",
                                  help_text="5' → 3' direction")
    fwd_or_rev = models.CharField(choices=FWD_OR_REV, max_length=1, blank=True)
    intended_use = models.CharField(max_length=1000, blank=True)
    qr_id = ShortUUIDField(default=generate_shortuuid, editable=False)

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['name']


class Stats(models.Model):
    id = models.AutoField(primary_key=True)
    plasmid_count = models.CharField(max_length=200, blank=True)
    plasmids_by_month = models.JSONField(null=True)
    plasmids_with_sequence = models.JSONField(null=True)
    plasmids_with_gs = models.JSONField(null=True)
    plasmids_by_type = models.JSONField(null=True)
    plasmids_by_level = models.JSONField(null=True)
    gs_box_fill = models.JSONField(null=True)
    last_update = models.DateField(auto_now_add=False, default=datetime.date.today)


class Experiment(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, help_text="Experiment name")
    description = models.CharField(max_length=500, help_text="Experiment description", null=True, blank=True)
    project = models.ForeignKey(Project, on_delete=models.CASCADE, blank=True, null=True)
    plasmids = models.ManyToManyField(Plasmid, blank=True, related_name='+')

    def __str__(self):
        return self.name
