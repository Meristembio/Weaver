import os
import uuid
from django.db import models
import datetime
from shortuuidfield import ShortUUIDField
import shortuuid
from .custom.box import BOX_ROWS
from .custom.box import BOX_COLUMNS
from .custom.general import FWD_OR_REV
from .custom.general import CHECK_STATES
from .custom.general import COLORS
from Bio.Restriction.Restriction_Dictionary import rest_dict, suppliers
from organization.models import Project
from django.dispatch import receiver


from .validators import clustal_validate

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
    level = models.IntegerField(blank=True, null=True)
    description = models.CharField(max_length=1000, blank=True)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    created_on = models.DateField(auto_now_add=False, default=datetime.date.today)

    reference_sequence = models.BooleanField(blank=True, default=0)
    under_construction = models.BooleanField(blank=True, default=0)
    public_visibility = models.BooleanField(blank=True, default=0)

    qr_id = ShortUUIDField(default=shortuuid.uuid(), editable=False)

    # validation

    working_colony = models.IntegerField(blank=True, null=True)

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

    def working_colony_text_for_ligation(self):
        if self.under_construction:
            return "UC"
        elif self.reference_sequence:
            return "RS"
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
        if self.under_construction:
            return "Under construction"
        elif self.reference_sequence:
            return "Reference sequence"
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

    def is_validated(self):
        if self.colonypcr_state != 1 and self.digestion_state != 1 and self.sequencing_state != 1:
            return True
        return False

    def save(self, *args, **kwargs):
        if not self.idx:
            last_plasmid_idx = Plasmid.objects.order_by("idx").last().idx
            if last_plasmid_idx:
                new_idx = last_plasmid_idx + 1
            else:
                new_idx = 1
            self.idx = new_idx
        super(Plasmid, self).save(*args, **kwargs)

    def get_insert_of(self):
        # ToDo optimize query
        insert_of = []
        for plasmid in Plasmid.objects.all():
            for insert in plasmid.inserts.all():
                if insert == self:
                    insert_of.append(plasmid)

        return insert_of
    def get_backbone_of(self):
        # ToDo optimize query
        backbone_of = []
        for plasmid in Plasmid.objects.all():
            if plasmid.backbone == self:
                backbone_of.append(plasmid)

        return backbone_of

    def ligation_concentration(self):
        if self.computed_size:
            if self.type:
                if str(self.type) == "Insert":
                    return str(round(self.computed_size / 300, 1)) + " ng / ul"
                elif str(self.type) == "Receiver":
                    return str(round(self.computed_size / 600, 1)) + " ng / ul"
                else:
                    return "Plasmid type no formula"
            else:
                return "No plasmid type set"
        else:
            return "No plasmid computed size"

    def recommended_enzyme_for_create(self):
        # ToDo generalize
        output = "No level set"
        if self.level is not None:
            output = "SapI"
            if self.level % 2:
                output = "BsaI"
        return output

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
            ligation_raw += self.backbone.__str__() + " [" + self.backbone.working_colony_text_for_ligation() + "]" + tab

        inserts = []
        for plasmid in self.inserts.all():
            inserts.append(plasmid.__str__() + " [" + plasmid.working_colony_text_for_ligation() + "]")

        if self.level:
            ligation_raw = ligation_raw + " + ".join(inserts) + tab + tab +\
                           self.recommended_enzyme_for_create() + tab +\
                           self.getPlasmidResistanceForLigation().capitalize()
        else:
            if self.level == 0:
                ligation_raw = "Level 0 ligation is not supported"
            else:
                ligation_raw = "Level not set"

        return ligation_raw


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
    qr_id = ShortUUIDField(default=shortuuid.uuid(), editable=False)
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