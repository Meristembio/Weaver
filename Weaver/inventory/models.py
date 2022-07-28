import uuid
from django.db import models
import datetime
from shortuuidfield import ShortUUIDField
import shortuuid
from .custom.box import BOX_ROWS
from .custom.box import BOX_COLUMNS
from .custom.general import FWD_OR_REV
from .custom.general import CHECK_STATES
from .custom.general import CHECK_METHODS
from .custom.general import SEQUENCING_STATES
from .custom.general import COLORS
from Bio.Restriction.Restriction_Dictionary import rest_dict, suppliers
from organization.models import Project

RE_Choices = []
for key in rest_dict:
    if not key.startswith("_"):
        RE_Choices.append((key, key))


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


class RestrictionEnzyme(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    # List https://github.com/biopython/biopython/blob/master/Bio/Restriction/Restriction_Dictionary.py
    name = models.CharField(choices=RE_Choices, max_length=20)
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

    def __str__(self):
        if self.hf_version:
            return self.name + "-HF"
        else:
            return self.name

    class Meta:
        ordering = ['name']


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
    name = models.CharField(max_length=50)
    selectable_markers = models.ManyToManyField(Resistance, blank=True, symmetrical=False, related_name='+', help_text='Use CTRL for multiple select')
    sequence = models.FileField(upload_to='uploads/plasmids/', blank=True)
    backbone = models.ForeignKey('self', models.SET_NULL, blank=True, null=True)
    computed_size = models.IntegerField(blank=True, null=True, editable=False)
    insert_computed_size = models.IntegerField(blank=True, null=True, editable=False)
    inserts = models.ManyToManyField('self', blank=True, symmetrical=False, related_name='+')
    intended_use = models.CharField(max_length=200)
    type = models.ForeignKey(PlasmidType, models.SET_NULL, blank=True, null=True)
    level = models.IntegerField(blank=True, null=True)
    description = models.CharField(max_length=1000, blank=True)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    created_on = models.DateField(auto_now_add=False, default=datetime.date.today)
    qr_id = ShortUUIDField(default=shortuuid.uuid(), editable=False)

    # validation

    working_colony = models.IntegerField(blank=True, null=True)
    check_state = models.IntegerField(choices=CHECK_STATES, blank=True, default=1)
    check_method = models.IntegerField(choices=CHECK_METHODS, blank=True, default=0, null=True)
    check_date = models.DateField(blank=True, null=True)
    digestion_check_enzymes = models.CharField(max_length=100, blank=True, null=True)
    check_observations = models.CharField(max_length=1000, blank=True, null=True)
    sequencing_state = models.IntegerField(choices=SEQUENCING_STATES, blank=True, default=0)
    sequencing_date = models.DateField(blank=True, null=True)
    sequencing_observations = models.CharField(max_length=1000, blank=True, null=True)

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['name']


class Strain(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=200, blank=True)
    selectable_markers = models.ManyToManyField(Resistance, blank=True, symmetrical=False, related_name='+', help_text='Use CTRL for multiple select')
    description = models.CharField(max_length=1000, blank=True)

    def __str__(self):
        return self.name

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
    qr_id = ShortUUIDField(default=shortuuid.uuid(), editable=False)
    details = models.CharField(max_length=1000, blank=True)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)

    def __str__(self):
        if self.plasmid is None:
            return self.strain.name
        else:
            return self.strain.name + " / " + self.plasmid.name

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
    qr_id = ShortUUIDField(default=shortuuid.uuid(), editable=False)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)

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
