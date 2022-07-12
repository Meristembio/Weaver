from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt

from .models import Plasmid
from .models import GlycerolStock
from .models import RestrictionEnzyme
from .models import Primer
from .models import Box
from .models import Location
from .custom.general import SEQUENCING_STATES
from .custom.general import CHECK_STATES
from .custom.general import CHECK_METHODS
from .models import Stats
from .models import PlasmidType
from .models import TableFilter
from .custom.standards import CURRENT_ASSEMBLY_STANDARD
from .custom.standards import recommended_enzyme_for_create

from .forms import PlasmidValidationForm

from django.http import HttpResponseRedirect

import html
from django.conf import settings
from django.http import HttpResponse, Http404
from django.core.exceptions import ObjectDoesNotExist
from .custom.box import BOX_ROWS
from .custom.box import BOX_COLUMNS
from django.views.generic.edit import UpdateView
from django.views.generic.edit import CreateView
from django.views.generic.edit import DeleteView
from django.urls import reverse

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Restriction import RestrictionBatch
from Bio.Restriction.Restriction_Dictionary import rest_dict

from .forms import DigestForm
from .forms import PCRForm
import json
from .forms import SangerForms
from .forms import L0SequenceInput
from .forms import BlastSequenceInput
from .forms import GstockEditForm
from .forms import GlycerolQRInput

import os
import tempfile
import re
import Bio
import django
from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove
from django.core.files.base import ContentFile
from datetime import datetime
from datetime import date
from django.http import JsonResponse

import plotly.express as px
import pandas as pd

LO_OHS = CURRENT_ASSEMBLY_STANDARD['odd_custom']


def get_table_filters(level_from_table_filters, level_to_table_filters):
    pt_table_filters = []
    for pt in PlasmidType.objects.all():
        pt_table_filters.append((pt.name, 't' + str(pt.id), 'success'))

    level_table_filters = []
    for level in range(level_from_table_filters, level_to_table_filters + 1):
        level_table_filters.append(('L' + str(level), 'l' + str(level), 'warning'))

    sw_table_filters = []
    for tf in TableFilter.objects.all():
        filters = []
        color = 'info'
        if tf.color:
            color = tf.color
        for the_filter in tf.options.split(","):
            name, search = the_filter.split("|")
            filters.append((name, search, color))
        sw_table_filters.append(['startswith', tf.name, filters])

    table_filters = [
                        ['all', 'All', [
                            ('All', 'all', 'primary'),
                        ]],
                        ['type', 'Type', pt_table_filters],
                        ['level', 'Level', level_table_filters]
                    ] + sw_table_filters
    return table_filters


def json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""
    if isinstance(obj, (datetime, date)):
        return obj.isoformat()
    raise TypeError("Type %s not serializable" % type(obj))


def restrictionenzyme(request, restrictionenzyme_id):
    try:
        restrictionenzyme_to_detail = RestrictionEnzyme.objects.get(id=restrictionenzyme_id)
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'restrictionenzyme': restrictionenzyme_to_detail,
    }
    return render(request, 'inventory/restrictionenzyme.html', context)


class RestrictionenzymeEdit(UpdateView):
    model = RestrictionEnzyme
    fields = '__all__'
    template_name_suffix = '_update_form'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['model_name'] = "Restriction Enzyme"
        return context

    def get_success_url(self, **kwargs):
        return reverse('restrictionenzyme', args=(self.object.id,)) + '?form_result_restrictionenzyme_edit_success=true'


class RestrictionenzymeCreate(CreateView):
    model = RestrictionEnzyme
    fields = '__all__'
    template_name_suffix = '_create_form'

    def get_success_url(self, **kwargs):
        return reverse('restrictionenzyme',
                       args=(self.object.id,)) + '?form_result_restrictionenzyme_create_success=true'


def restrictionenzymes(request):
    context = {
        'restrictionenzymes': RestrictionEnzyme.objects.all(),
    }
    return render(request, 'inventory/restrictionenzymes.html', context)


def glycerolstocks(request):
    glycerolstocks = GlycerolStock.objects.all()

    level_from_table_filters = 0
    level_to_table_filters = 0
    for glycerolstock in glycerolstocks:
        if glycerolstock.plasmid:
            if glycerolstock.plasmid.level:
                if glycerolstock.plasmid.level > level_to_table_filters:
                    level_to_table_filters = glycerolstock.plasmid.level
                if glycerolstock.plasmid.level < level_from_table_filters:
                    level_from_table_filters = glycerolstock.plasmid.level
    context = {
        'glycerolstocks': glycerolstocks,
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
    }
    return render(request, 'inventory/glycerolstocks.html', context)


def glycerolstock(request, glycerolstock_id):
    try:
        glycerolstock_to_detail = GlycerolStock.objects.get(id=glycerolstock_id)
    except ObjectDoesNotExist:
        raise Http404

    resistantes_human_context = "None"
    if glycerolstock_to_detail.plasmid:
        resistantes_human_context = resistantes_human(glycerolstock_to_detail.plasmid.resistances)

    context = {
        'glycerolstock': glycerolstock_to_detail,
        'resistantes_human': resistantes_human_context,
        'resistantes_strain_human': resistantes_human(glycerolstock_to_detail.strain.resistances),
    }
    return render(request, 'inventory/glycerolstock.html', context)


def glycerolstock_qr(request):
    if request.method == 'POST' and 'glycerol_id' in request.POST:
        form = GlycerolQRInput(request.POST)
        if form.is_valid():
            return glycerolstock_from_qr(request, form.cleaned_data['glycerol_id'])
    else:
        context = {
            'glycerol_qr_nput_form': GlycerolQRInput()
        }
    return render(request, 'inventory/glycerolstock_qr.html', context)


def glycerolstock_from_qr(request, glycerolstock_id):
    try:
        glycerolstock_to_detail = GlycerolStock.objects.filter(qr_id=glycerolstock_id)[0]
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'glycerolstock': glycerolstock_to_detail,
    }
    return render(request, 'inventory/glycerolstock.html', context)


class GstockEdit(UpdateView):
    model = GlycerolStock
    template_name_suffix = '_update_form'
    form_class = GstockEditForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["collection"] = build_boxes()
        context["render_mod"] = 'p'
        return context

    def get_success_url(self, **kwargs):
        return reverse('glycerolstock', args=(self.object.id,)) + '?form_result_glycerolstock_edit_success=true'


def build_boxes():
    output = {
        'BOX_ROWS': BOX_ROWS,
        'BOX_COLUMNS': BOX_COLUMNS,
        'locations': []
    }
    glycerolstocks = GlycerolStock.objects.all()
    for location in Location.objects.all():
        boxes = []
        boxes_at_location = Box.objects.filter(location=location).order_by('name')

        for box in boxes_at_location:
            box_output = {
                'name': box.name,
                'id': box.id
            }
            for glycerolstock in glycerolstocks:
                if glycerolstock.box and glycerolstock.box == box and glycerolstock.box_row and glycerolstock.box_column:
                    box_output[str(glycerolstock.box_row) + str(glycerolstock.box_column)] = glycerolstock

            boxes.append(box_output)

        output['locations'].append({
            'name': location.name,
            'boxes': boxes
        })
    return output


class GstockCreate(CreateView):
    model = GlycerolStock
    template_name_suffix = '_create_form'
    form_class = GstockEditForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["collection"] = build_boxes()
        context["render_mod"] = 'p'
        return context

    def get_success_url(self, **kwargs):
        return reverse('glycerolstock', args=(self.object.id,)) + '?form_result_glycerolstock_create_success=true'


class GstockCreatePlasmidDefined(CreateView):
    model = GlycerolStock
    form_class = GstockEditForm
    template_name_suffix = '_create_form'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["collection"] = build_boxes()
        context["pid"] = self.kwargs['pid']
        context["render_mod"] = 'p'
        return context

    def get_success_url(self, **kwargs):
        return reverse('glycerolstock', args=(self.object.id,)) + '?form_result_glycerolstock_create_success=true'


def glycerolstock_label(request, glycerolstock_id):
    try:
        glycerolstock_to_label = GlycerolStock.objects.get(id=glycerolstock_id)
    except ObjectDoesNotExist:
        raise Http404

    resistantes_human_context = "None"
    if glycerolstock_to_label.plasmid:
        resistantes_human_context = resistantes_human(glycerolstock_to_label.plasmid.resistances, True)

    context = {
        'glycerolstock': glycerolstock_to_label,
        'resistantes_human': resistantes_human_context,
        'resistantes_strain_human': resistantes_human(glycerolstock_to_label.strain.resistances, True),
    }
    return render(request, 'inventory/glycerolstock_label.html', context)


def glycerolstock_boxes(request):
    context = {
        'collection': build_boxes(),
        'render_mod': 'n'
    }
    return render(request, 'inventory/glycerolstock_boxes.html', context)


def plasmids(request):
    plasmids = Plasmid.objects.all()
    level_from_table_filters = 0
    level_to_table_filters = 0
    for plasmid in plasmids:
        plasmid.refc = recommended_enzyme_for_create(plasmid.level)
        if plasmid.level:
            if plasmid.level > level_to_table_filters:
                level_to_table_filters = plasmid.level
            if plasmid.level < level_from_table_filters:
                level_from_table_filters = plasmid.level
    context = {
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
        'plasmids': plasmids,
        'RESTRICTION_ENZYMES': RestrictionEnzyme.objects.all,
    }
    return render(request, 'inventory/plasmids.html', context)


def getPlasmidResistanceForLigation(plasmid):
    if plasmid.level:
        if plasmid.level % 2 == 0:
            is_spe = False
            for resistance in plasmid.resistances.all():
                if resistance.three_letter_code == 'SPE':
                    is_spe = True
                    break
            if is_spe:
                return 'SPE'
            else:
                return 'Error: SPE not in resistance list'
        else:
            is_kan = False
            for resistance in plasmid.resistances.all():
                if resistance.three_letter_code == 'KAN':
                    is_kan = True
                    break
            if is_kan:
                return 'KAN'
            else:
                return 'Error: KAN not in resistance list'
    if len(plasmid.resistances):
        return plasmid.resistances[0]
    else:
        return 'More than one'


def resistantes_human(resistances, short=False):
    resistantes_human_return = []
    if resistances:
        for resistance in resistances.all():
            if short:
                resistantes_human_return.append(str(resistance.three_letter_code))
            else:
                resistantes_human_return.append(resistance.name + " (" + str(resistance.three_letter_code) + ")")
            continue
        return " / ".join(resistantes_human_return)
    else:
        return "None"


def plasmid(request, plasmid_id):
    try:
        plasmid_to_detail = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    insert_of = []
    backbone_of = []
    tab = "	"
    ligation_raw = plasmid_to_detail.name + tab
    if plasmid_to_detail.backbone:
        ligation_raw += plasmid_to_detail.backbone.name + tab
    for plasmid in Plasmid.objects.all():
        if plasmid.backbone == plasmid_to_detail:
            backbone_of.append(plasmid)
        for insert in plasmid.inserts.all():
            if insert == plasmid_to_detail:
                insert_of.append(plasmid)

    inserts = []
    for plasmid in plasmid_to_detail.inserts.all():
        inserts.append(plasmid.name)

    if plasmid_to_detail.level:
        ligation_raw = ligation_raw + " + ".join(inserts) + tab + recommended_enzyme_for_create(
            plasmid_to_detail.level) + tab + getPlasmidResistanceForLigation(plasmid_to_detail).capitalize()
    else:
        if plasmid_to_detail.level == 0:
            ligation_raw = "Level 0 ligation is not supported"
        else:
            ligation_raw = "Level not set"

    plasmid_to_detail.refc = recommended_enzyme_for_create(plasmid_to_detail.level)

    context = {
        'plasmid': plasmid_to_detail,
        'resistantes_human': resistantes_human(plasmid_to_detail.resistances),
        'insert_of': insert_of,
        'backbone_of': backbone_of,
        'ligation_raw': ligation_raw,
        'SEQUENCING_STATES': SEQUENCING_STATES,
        'CHECK_STATES': CHECK_STATES,
        'CHECK_METHODS': CHECK_METHODS,
        'RESTRICTION_ENZYMES': RestrictionEnzyme.objects.all()
    }

    if request.method == 'POST' and 'l0_sequence_input' in request.POST:
        form = L0SequenceInput(request.POST)
        if form.is_valid():
            plasmid_create_from_inserts(plasmid_to_detail, context, insert=form.cleaned_data['l0_sequence_input'],
                                        oh_5=form.cleaned_data['l0_oh_5'], oh_3=form.cleaned_data['l0_oh_3'],
                                        the_re=RestrictionEnzyme.objects.get(name=request.POST.get('enzyme')))
        else:
            context['plasmid_create_result'] = ("Bad inputs", "danger")

    if request.method == 'POST' and 'create_from_parts' in request.POST and 'enzyme' in request.POST:
        if plasmid_to_detail.level == 0:
            form = L0SequenceInput(initial={'enzyme': request.POST.get('enzyme')})
            return render(request, 'inventory/plasmid.html',
                          {'L0SequenceInputForm': form, 'plasmid': plasmid_to_detail})
        elif plasmid_to_detail.level == -1:
            form = L0SequenceInput(initial={'enzyme': request.POST.get('enzyme')})
            return render(request, 'inventory/plasmid.html',
                          {'L_1SequenceInputForm': form, 'plasmid': plasmid_to_detail})
        else:
            if 'enzyme' in request.POST:
                plasmid_create_from_inserts(plasmid_to_detail, context)
            else:
                context['plasmid_create_result'] = ("No enzyme selected", "danger")

    if request.method == 'GET' and 'ac' in request.GET:
        plasmid_create_from_inserts(plasmid_to_detail, context)

    if request.method == 'POST' and 'params' in request.POST:
        context['plasmid_create_result'] = ("Plasmid create wizard is complete.", "success")

    # in case of update or never computed
    plasmid_update_computed_size(plasmid_to_detail)

    return render(request, 'inventory/plasmid.html', context)


def plasmid_create_from_inserts(plasmid_to_build, context, insert=None, oh_5=None, oh_3=None, the_re=None):
    if not the_re:
        the_re = RestrictionEnzyme.objects.get(name=recommended_enzyme_for_create(plasmid_to_build.level))
    plasmid_record = plasmid_record_from_inserts(plasmid_to_build, insert, oh_5, oh_3, the_re)
    if plasmid_record[0]:
        plasmid_record_final = plasmid_record[1]
        plasmid_record_final.name = plasmid_record_final.name.replace(" ", "_")
        plasmid_to_build.sequence.save(plasmid_to_build.name + ".gb", ContentFile(plasmid_record[1].format("gb")))
        context['plasmid_create_result'] = ("Plasmid sequence built from backbone / insert data", "success")
    else:
        context['plasmid_create_result'] = (plasmid_record[1], "danger")


def plasmid_record_from_inserts(plasmid_to_build, insert, oh_5, oh_3, the_re):
    if plasmid_to_build.level != 0 and plasmid_to_build.level != -1 and len(plasmid_to_build.inserts.all()) == 0:
        return False, "No inserts defined"
    if not plasmid_to_build.backbone:
        return False, "No backbone set"

    # take the backbone
    plasmid_backbone_seq_result = seqio_get(plasmid_to_build.backbone)
    oh_length = abs(the_re.rcut - the_re.fcut)

    if plasmid_backbone_seq_result[0]:
        backbone_record = plasmid_backbone_seq_result[1]
        hits = re_find_cut_positions(backbone_record.seq, the_re, True, True)
        if len(hits) == 2 and hits[1] > hits[0]:

            # go for the inserts
            if plasmid_to_build.level == 0 or plasmid_to_build.level == -1:
                final_record = backbone_record[0:hits[0] + oh_length - 1]

                for name, value in LO_OHS:
                    if oh_5 == name:
                        oh_5_tuple = (name, value)
                    if oh_3 == name:
                        oh_3_tuple = (name, value)

                if oh_5_tuple[1]:
                    rec_oh_5 = SeqRecord(
                        Seq(oh_5_tuple[1].upper()),
                        id="oh_5",
                        annotations={"molecule_type": "DNA"}
                    )
                    rec_oh_5.features.append(SeqFeature(
                        FeatureLocation(0, len(oh_5_tuple[1])),
                        type="misc_feature",
                        strand=1,
                        qualifiers={
                            'ApEinfo_label': oh_5_tuple[0],
                            'label': oh_5_tuple[0],
                            'locus_tag': oh_5_tuple[0],
                        }
                    ))
                    final_record = final_record + rec_oh_5
                rec_insert = SeqRecord(
                    Seq(insert.lower()),
                    id="oh_h",
                    annotations={"molecule_type": "DNA"}
                )
                rec_insert.features.append(SeqFeature(
                    FeatureLocation(0, len(insert)),
                    type="misc_feature",
                    strand=1,
                    qualifiers={
                        'ApEinfo_label': plasmid_to_build.name,
                        'label': plasmid_to_build.name,
                        'locus_tag': plasmid_to_build.name,
                    }
                ))
                final_record = final_record + rec_insert

                if oh_3_tuple[1]:
                    rec_oh_3 = SeqRecord(
                        Seq(oh_3_tuple[1].upper()),
                        id="oh_3",
                        annotations={"molecule_type": "DNA"}
                    )
                    rec_oh_3.features.append(SeqFeature(
                        FeatureLocation(0, len(oh_3_tuple[1])),
                        type="misc_feature",
                        strand=1,
                        qualifiers={
                            'ApEinfo_label': oh_3_tuple[0],
                            'label': oh_3_tuple[0],
                            'locus_tag': oh_3_tuple[0],
                        }
                    ))
                    final_record = final_record + rec_oh_3

                final_record = final_record + backbone_record[hits[1] - 1:]
            else:
                final_record = backbone_record[0:hits[0] - 1]
                first_oh = backbone_record[hits[0] - 1:hits[0] + oh_length - 1]
                last_oh = backbone_record[hits[1] - 1: hits[1] + oh_length - 1]

                inserts = []
                for insert in plasmid_to_build.inserts.all():
                    if seqio_get(insert)[0]:
                        insert_record = seqio_get(insert)[1]
                        insert_hits = re_find_cut_positions(insert_record.seq, the_re, True, True)

                        if len(insert_hits) == 0:
                            return False, "No " + re.name + " sites found at " + insert.name
                        if len(insert_hits) != 2:
                            return False, "!= 2 " + re.name + " sites found at " + insert.name + ". Found sites #: " \
                                   + len(insert_hits)
                        inserts.append((
                            insert_record[insert_hits[0] - 1:insert_hits[1] - 1],
                            str(insert_record.seq[insert_hits[0] - 1:insert_hits[0] + oh_length - 1]),
                            str(insert_record.seq[insert_hits[1] - 1:insert_hits[1] + oh_length - 1])
                        ))
                    else:
                        return False, "Error reading insert sequence file [" + insert.name + "]"

                last_oh_added = first_oh.seq
                joined = []
                while last_oh_added.lower() != last_oh.seq.lower():
                    init_last_oh_added = last_oh_added
                    for insert in inserts:
                        # print("Comparing: " + insert[1].lower() + " == " + last_oh_added.lower())
                        if insert[1].lower() == last_oh_added.lower():
                            final_record = final_record + insert[0]
                            last_oh_added = insert[2]
                            joined.append(insert[0].name + "/" + insert[1] + "/" + insert[2])
                            break

                    # if couldnt find a new part, something is missing in part list
                    if init_last_oh_added.lower() == last_oh_added.lower():
                        return False, "Inserts are not concatenated from " + first_oh.seq + " to " + \
                               last_oh.seq + ". Joined = " + " + ".join(joined) + "."

                final_record = final_record + backbone_record[hits[1] - 1:]

            final_record.id = str(plasmid_to_build.id)
            final_record.name = plasmid_to_build.name
            final_record.description = plasmid_to_build.description
            final_record.annotations = {"molecule_type": "DNA", "topology": "circular"}

            return True, final_record
        else:
            return False, "No restriction sites found at backbone (" + the_re.name + ")"
    else:
        return False, "Error reading backbone sequence file"


def plasmid_from_qr(request, plasmid_id):
    try:
        plasmid_to_detail = Plasmid.objects.filter(qr_id=plasmid_id)[0]
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'plasmid': plasmid_to_detail,
    }
    return render(request, 'inventory/plasmid.html', context)


class PlasmidEdit(UpdateView):
    model = Plasmid
    fields = ['name', 'resistances', 'sequence', 'backbone', 'inserts', 'intended_use', 'type', 'level', 'author',
              'description', 'created_on']
    template_name_suffix = '_update_form'

    def get_success_url(self, **kwargs):
        return reverse('plasmid', args=(self.object.id,)) + '?form_result_plasmid_edit_success=true'


class PlasmidCreate(CreateView):
    model = Plasmid
    fields = '__all__'
    template_name_suffix = '_create_form'

    def get_success_url(self, **kwargs):
        auto_create = ""
        if self.request.GET.get('b'):
            # comes from wizard --> auto assemble
            auto_create = "&ac"
        return reverse('plasmid', args=(self.object.id,)) + '?form_result_plasmid_create_success=true' + auto_create


class PlasmidCreateWizard(CreateView):
    model = Plasmid
    fields = '__all__'
    template_name_suffix = '_create_form_wizard'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["create_url"] = reverse('plasmid_create')
        context["next_url"] = reverse('plasmid_create_wizard_end')
        return context

    def get_success_url(self, **kwargs):
        return reverse('plasmid', args=(self.object.id,))


def plasmid_create_wizard_end(request):
    context = {}
    if request.method == 'POST' and 'params' in request.POST:
        params = {}
        for param in request.POST.get('params').split('&'):
            name, values = param.split('=')
            if len(values.split('+')) > 1:
                values = values.split('+')
            params[name] = values
        if 'n' in params and 'b' in params and 'i' in params and params['n'] and params['b'] and len(params['i']):
            backbone = Plasmid.objects.get(id=params['b'])
            if backbone:
                description = ""
                if 'd' in params:
                    description = params['d']
                plasmid_created = Plasmid.objects.create(
                    name=params['n'],
                    description=description,
                    backbone=backbone,
                    type=PlasmidType.objects.get(id=0),
                    level=backbone.level,
                )
                for i in params['i']:
                    plasmid_created.inserts.add(i)
                for r in backbone.resistances.all():
                    plasmid_created.resistances.add(r.id)
                plasmid_create_from_inserts(plasmid_created, context)
                return plasmid(request, plasmid_created.id)
            else:
                context['wizard_error'] = 'Backbone not found.'
        else:
            context['wizard_error'] = 'Name & backbone & inserts are required fields.'
    else:
        context['wizard_error'] = "Plasmid can\'t be created. No parameters set."
    return render(request, 'inventory/plasmid.html', context)


@csrf_exempt
def plasmid_view_edit(request, plasmid_id):
    try:
        plasmid_to_detail = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    warnings = []

    if request.user.has_perm('inventory.change_plasmid'):
        if request.method == 'POST' and 'saveOve' in request.POST:
            if 'gbContent' in request.POST:
                # saving from OVE
                plasmid_to_detail.sequence.save(plasmid_to_detail.name + ".gb", ContentFile(request.POST['gbContent']))
                result = 'File saved'
            else:
                result = 'Error: no gbContent'
            return HttpResponse(json.dumps({
                'result': result
            }), content_type="application/json")

        if request.method == 'POST' and 'create' in request.POST:
            if not plasmid_to_detail.sequence:
                plasmid_to_detail.sequence.save(plasmid_to_detail.name + ".gb", ContentFile('''LOCUS       ''' + plasmid_to_detail.name.replace(" ", "_") + '''                   0 bp    DNA     circular  31-AUG-2021
    ORIGIN      
    //'''))
            else:
                warnings.append('Can\'t create empty sequence on this plasmid, already has one')

    else:
        warnings.append('Current user can\'t make modifications')

    context = {
        'plasmid': plasmid_to_detail,
        'sequence_file_contents': plasmid_sequence_file_contents(plasmid_to_detail),
        'warnings': warnings
    }
    return render(request, 'inventory/plasmid_view_edit.html', context)


def plasmid_sequence_file_contents(plasmid):
    with open(plasmid.sequence.path, 'r') as file:
        sequence_file_contents = html.unescape(file.read())

    return re.sub(r'[\'"]', '', sequence_file_contents)


def plasmid_download(request, plasmid_id):
    try:
        plasmid_to_download = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    if request.method == 'GET' and 'format' in request.GET:
        record_response = seqio_get(plasmid_to_download)
        if record_response[0]:
            response = HttpResponse(record_response[1].format(request.GET['format']), content_type="plain/text")
            response['Content-Disposition'] = 'inline; filename=' + plasmid_to_download.name + '.' + request.GET[
                'format']
            return response

    # return original file if no format is specified
    file_path = os.path.join(settings.MEDIA_ROOT, plasmid_to_download.sequence.name)
    if os.path.exists(file_path):
        with open(file_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="plain/text")
            response['Content-Disposition'] = 'inline; filename=' + plasmid_to_download.name + \
                                              os.path.splitext(os.path.basename(file_path))[1]
            return response
    raise Http404


def plasmid_label(request, plasmid_id):
    try:
        plasmid_to_label = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'plasmid': plasmid_to_label,
        'today': date.today(),
    }
    return render(request, 'inventory/plasmid_label.html', context)


def plasmid_digest(request, plasmid_id):
    try:
        plasmid_to_digest = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'plasmid': plasmid_to_digest,
    }

    sequence = grab_seq(plasmid_to_digest)

    if sequence[0]:
        context['digest_form'] = DigestForm()
        res = RestrictionEnzyme.objects.all()
        for the_re in res:
            the_re.hits, the_re.fragments = re_find_cut_fragments(sequence[1], the_re, True)

        context['res'] = res

        if request.method == 'POST':
            selected_res = []
            post_enzymes = json.loads(request.POST['enzymes'])
            for the_re in res:
                if the_re.name in post_enzymes:
                    selected_res.append(the_re)
            context['selected_res'] = selected_res
            context['fragments'] = re_digestion_fragments(sequence[1], selected_res, True)

    else:
        context['error'] = sequence[1]

    return render(request, 'inventory/plasmid_digest.html', context)


def plasmid_pcr(request, plasmid_id):
    try:
        plasmid_to_pcr = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'plasmid': plasmid_to_pcr,
        'show_new_PCR': False,
    }

    sequence = grab_seq(plasmid_to_pcr)

    if sequence[0]:
        if request.method == 'POST':
            if request.POST['primer_f'] != "":
                primer_f = Primer.objects.get(id=request.POST['primer_f'])
            else:
                if request.POST['primer_f_seq'] != "":
                    primer_f = Primer(
                        id='custom_f',
                        name='Custom F',
                        sequence_3=request.POST['primer_f_seq'],
                        sequence_5='',
                        fwd_or_rev='f',
                        intended_use='Custom sequence for PCR prediction'
                    )
                else:
                    context['error'] = "No forward primer set"
            if request.POST['primer_r'] != "":
                primer_r = Primer.objects.get(id=request.POST['primer_r'])
            else:
                if request.POST['primer_r_seq'] != "":
                    primer_r = Primer(
                        id='custom_r',
                        name='Custom R',
                        sequence_3=request.POST['primer_r_seq'],
                        sequence_5='',
                        fwd_or_rev='f',
                        intended_use='Custom sequence for PCR prediction'
                    )
                else:
                    context['error'] = "No forward primer set"
            if not 'error' in context:
                context['primer_f'] = primer_f
                context['primer_r'] = primer_r
                if primer_r.sequence_5:
                    context['primer_r_5_rc'] = str(Seq(primer_r.sequence_5).reverse_complement())
                if primer_r.sequence_3:
                    context['primer_r_3_rc'] = str(Seq(primer_r.sequence_3).reverse_complement())
                double_seq = str(sequence[1]) + str(sequence[1])
                pos_f = re.search(primer_f.sequence_3, double_seq, re.IGNORECASE)
                if pos_f:
                    start = pos_f.end()
                    seq_from_f = double_seq[start:]
                    pos_r = re.search(context['primer_r_3_rc'], seq_from_f, re.IGNORECASE)
                    if pos_r:
                        end = pos_r.start()
                        context['amplicon'] = seq_from_f[:end].lower()
                        context['size'] = len(
                            primer_f.sequence_5 + primer_f.sequence_3 + context['amplicon'] + primer_r.sequence_3 +
                            primer_r.sequence_5)
                    else:
                        context['error'] = "REV primer does not hit template"
                else:
                    context['error'] = "FWD primer does not hit template"
            context['show_new_PCR'] = True
        else:
            context['pcr_form'] = PCRForm()

    else:
        context['error'] = sequence[1]

    return render(request, 'inventory/plasmid_pcr.html', context)


def plasmid_sanger(request, plasmid_id):
    try:
        plasmid_to_align = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'plasmid': plasmid_to_align,
    }

    if request.method == 'POST':
        form = SangerForms(request.POST, request.FILES)
        if form.is_valid():
            ab1_chromatos = []
            ab1_chromatos.append({})  # no chromato data for template

            ab1_file = request.FILES["ab1"]
            f = tempfile.NamedTemporaryFile(delete=False)
            for chunk in ab1_file.chunks():
                f.write(chunk)
            f.close()
            record = SeqIO.read(f.name, "abi")

            plasmid_seq = grab_seq(plasmid_to_align)[1]
            ab1_chromatos.append({
                'aTrace': record.annotations['abif_raw']['DATA9'],
                'tTrace': record.annotations['abif_raw']['DATA10'],
                'gTrace': record.annotations['abif_raw']['DATA11'],
                'cTrace': record.annotations['abif_raw']['DATA12'],
                'basePos': record.annotations['abif_raw']['PLOC2'],
                'baseCalls': list(record.annotations['abif_raw']['PBAS2'].decode()),
                'qualNums': []
            })
            if request.POST.get('is_reverse'):
                alignment = pairwise2.align.localxx(plasmid_seq, reverse_complement(record.seq))[0]
            else:
                alignment = pairwise2.align.localxx(plasmid_seq, record.seq)[0]

            context['align_data'] = json.dumps([
                [plasmid_to_align.name, ab1_file.name],
                [str(plasmid_seq), str(record.seq)],
                [alignment[0], alignment[1]],
                ab1_chromatos
            ])
            context['plasmid_sequence_file_contents'] = plasmid_sequence_file_contents(plasmid_to_align)
    else:
        context['upload_form'] = SangerForms()
        context['show_upload_form'] = True
    return render(request, 'inventory/plasmid_sanger.html', context)


def plasmid_update_computed_size(plasmid_to_update):
    sequence = grab_seq(plasmid_to_update)

    if sequence[0]:
        if not plasmid_to_update.level is None:
            if plasmid_to_update.level % 2:
                re = RestrictionEnzyme.objects.filter(name="SapI")[0]
            else:
                re = RestrictionEnzyme.objects.filter(name="BsaI")[0]

            hits = re_find_cut_positions(sequence[1], re, True, True)

            if len(hits) == 2 and hits[1] > hits[0]:
                plasmid_to_update.insert_computed_size = hits[1] - hits[0] + 4

        plasmid_to_update.computed_size = len(sequence[1])
        plasmid_to_update.save()
        return True
    return False


def grab_seq(plasmid_to_grab_from):
    result, gb_record = seqio_get(plasmid_to_grab_from)
    if result:
        return True, gb_record.seq
    return result, gb_record


def grab_features(plasmid_to_grab_from):
    result, gb_record = seqio_get(plasmid_to_grab_from)
    if result:
        return True, gb_record.features
    return result, gb_record


def grab_features_json(plasmid_to_grab_from):
    result, gb_features = grab_features(plasmid_to_grab_from)
    features = []
    if result:
        for gb_feature in gb_features:
            start = int(gb_feature.location.start)
            end = int(gb_feature.location.end)
            if type(gb_feature.location) is Bio.SeqFeature.CompoundLocation:
                if len(gb_feature.location.parts):
                    # asume partes contiguas
                    start = gb_feature.location.parts[1].start
                    end = gb_feature.location.parts[0].end
            forward = False
            if gb_feature.location.strand:
                forward = True
            features.append({
                'name': gb_feature.qualifiers['label'][0],
                'type': gb_feature.type,
                'start': start,
                'end': end,
                'forward': forward,
            })
    return json.dumps(features)


def seqio_get(plasmid_to_grab_from):
    name, extension = os.path.splitext(plasmid_to_grab_from.sequence.name)
    format_name = ''
    if extension == '.gb' or extension == '.gbk':
        format_name = "genbank"
    if extension == '.fasta':
        format_name = "fasta"
    if format_name:
        try:
            for gb_record in SeqIO.parse(plasmid_to_grab_from.sequence.path, format_name):
                return True, gb_record
        except ValueError as e:
            try:
                # Create temp file
                fh, abs_path = mkstemp()
                file_path = plasmid_to_grab_from.sequence.path
                with fdopen(fh, 'w') as new_file:
                    with open(file_path) as old_file:
                        for line in old_file:
                            new_line = line
                            if line.startswith("LOCUS"):
                                line_split = []
                                for idx, val in enumerate(line.split()):
                                    if idx == 3:
                                        continue
                                    if idx == 2:
                                        line_split.append(val + " bp")
                                    else:
                                        line_split.append(val)
                                spaces = [12, 13, 11, 16, 10, 12]
                                new_line = ""
                                for idx, val in enumerate(line_split):
                                    if idx >= len(spaces):
                                        new_line = new_line + val + " "
                                    else:
                                        if spaces[idx] > len(val):
                                            new_line = new_line + val + " " * (spaces[idx] - len(val))
                                        else:
                                            new_line = new_line + val[:spaces[idx] - 1] + " "
                                new_line = new_line + "\n"
                            new_file.write(line.replace(line, new_line))
                # Copy the file permissions from the old file to the new file
                copymode(file_path, abs_path)
                # Remove original file
                remove(file_path)
                # Move new file
                move(abs_path, file_path)
                # Ready
                for gb_record in SeqIO.parse(file_path, format_name):
                    return True, gb_record
            except ValueError as e:
                return False, 'File bad format: ' + e.__str__()
            except FileNotFoundError as e:
                return False, 'File not found: ' + e.__str__()
    return False, 'Unsupported file extension'


def re_digestion_fragments(sequence, the_res, is_circular):
    ordered_results = []
    for the_re in the_res:
        cut_positions = re_find_cut_positions(sequence, the_re, is_circular, True)
        for cp in cut_positions:
            ordered_results.append((the_re, cp))
    ordered_results.sort(key=lambda tup: tup[1])

    fragments = []
    prev_fragment = None
    last_element = ordered_results[len(ordered_results) - 1]
    for ordered_result in ordered_results:
        if len(fragments) == 0:
            # first item
            if is_circular:
                prev_fragment = {
                    'end': last_element[1],
                    'right': last_element[0]
                }
            else:
                prev_fragment = {
                    'end': 0,
                    'right': 'None'
                }
        fragment = {
            'start': prev_fragment['end'],
            'end': ordered_result[1],
            'left': prev_fragment['right'],
            'right': ordered_result[0],
            'length': ordered_result[1] - prev_fragment['end'],
        }
        if len(fragments) == 0 and is_circular:
            fragment['length'] = ordered_result[1] + len(sequence) - last_element[1]
        fragments.append(fragment)
        prev_fragment = fragment

    return fragments


def re_find_cut_fragments(sequence, the_re, is_circular):
    cut_positions = re_find_cut_positions(sequence, the_re, is_circular, True)
    fragments = []
    if cut_positions:
        prev_cut_pos = 0
        for cp in cut_positions:
            fragments.append(cp - prev_cut_pos)
            prev_cut_pos = cp

        last_frag = len(sequence) - cut_positions[len(cut_positions) - 1]
        if is_circular:
            fragments[0] += last_frag
        else:
            fragments.append(last_frag)
    return cut_positions, sorted(fragments)


def re_find_cut_positions(sequence, the_re, is_circular, sort):
    search_results = RestrictionBatch([the_re.name]).search(Seq(sequence), linear=not is_circular)
    found_hits = []
    for key in search_results:
        if str(key) == the_re.name:
            found_hits = search_results[key]
    if sort:
        found_hits = sorted(found_hits)
    return found_hits


def primer(request, primer_id):
    try:
        primer_to_detail = Primer.objects.get(id=primer_id)
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'primer': primer_to_detail,
    }
    return render(request, 'inventory/primer.html', context)


def primers(request):
    context = {
        'primers': Primer.objects.all()
    }
    return render(request, 'inventory/primers.html', context)


class PrimerEdit(UpdateView):
    model = Primer
    fields = '__all__'
    template_name_suffix = '_update_form'

    def get_success_url(self, **kwargs):
        return reverse('primer', args=(self.object.id,)) + '?form_result_primer_edit_success=true'


class PrimerCreate(CreateView):
    model = Primer
    fields = '__all__'
    template_name_suffix = '_create_form'

    def get_success_url(self, **kwargs):
        return reverse('primer', args=(self.object.id,)) + '?form_result_primer_create_success=true'


def primer_label(request, primer_id):
    try:
        primer_to_label = Primer.objects.get(id=primer_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'primer': primer_to_label,
        'date': datetime.now().date(),
    }
    return render(request, 'inventory/primer_label.html', context)


class PlasmidDelete(DeleteView):
    model = Plasmid

    def get_success_url(self, **kwargs):
        return reverse('plasmids') + '?form_result_object_deleted=true'


class GstockDelete(DeleteView):
    model = GlycerolStock

    def get_success_url(self, **kwargs):
        return reverse('glycerolstocks') + '?form_result_object_deleted=true'


class RestrictionenzymeDelete(DeleteView):
    model = RestrictionEnzyme

    def get_success_url(self, **kwargs):
        return reverse('restrictionenzymes') + '?form_result_object_deleted=true'


class PrimerDelete(DeleteView):
    model = Primer

    def get_success_url(self, **kwargs):
        return reverse('primers') + '?form_result_object_deleted=true'


def PlasmidValidationEdit(request, plasmid_id):
    try:
        plasmid_to_validate = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    form = PlasmidValidationForm(request.POST or None, instance=plasmid_to_validate)
    if form.is_valid():
        form.save()
        return HttpResponseRedirect(reverse('plasmid', args=(
            plasmid_to_validate.id,)) + '?form_result_plasmidvalidation_edit_success=true')

    return render(request, 'inventory/plasmidvalidation_update_form.html',
                  {'form': form, 'plasmid': plasmid_to_validate})


def PlasmidValidations(request):
    plasmidsToCheck = []
    plasmidsToStock = []
    plasmidsWithStockWoCheck = []

    all_plasmids = Plasmid.objects

    # Plasmid w/o GS && Check state == Pending
    for plasmid in all_plasmids.filter(check_state=1):
        if plasmid.glycerolstock_set.all().count() == 0:
            plasmidsToCheck.append(plasmid)

    # Plasmid w/o GS && Check state != Pending && Sequencing state != Required
    for plasmid in all_plasmids.exclude(check_state=1).exclude(sequencing_state=1):
        if plasmid.glycerolstock_set.all().count() == 0:
            plasmidsToStock.append(plasmid)

    # Plasmid w/ GS && Check state == Pending && Sequencing state != Required
    for plasmid in all_plasmids.filter(check_state=1).exclude(sequencing_state=1):
        if plasmid.glycerolstock_set.all().count() != 0:
            plasmidsWithStockWoCheck.append(plasmid)

    context = {
        'plasmidsToCheck': plasmidsToCheck,
        # Sequencing state = Required
        'plasmidsToSequence': all_plasmids.filter(sequencing_state=1),
        'plasmidsToStock': plasmidsToStock,
        'plasmidsWithStockWoCheck': plasmidsWithStockWoCheck,
        'CHECK_METHODS': CHECK_METHODS
    }
    return render(request, 'inventory/plasmidvalidations.html', context)


def ServicesStats(request):
    if Stats.objects.all():
        stats = Stats.objects.all()[0]
    else:
        stats = Stats()
        stats.save()

    context = {
        'error': 'No data'
    }

    if request.method == 'POST' and 'refresh_stats' in request.POST:
        plasmids_by_month = {'date': [], 'plasmid_month_count': []}
        current_year = ''
        current_month = ''
        current_month_count = 0
        plasmid_count = 0
        plasmids_with_sequence = 0
        plasmids_with_gs = 0
        plasmids_by_type = {}
        plasmids_by_level = {}
        plasmids_ordered = Plasmid.objects.order_by('created_on')
        for plasmid in plasmids_ordered:
            plasmid_count += 1
            year = plasmid.created_on.year
            month = plasmid.created_on.month
            if current_year == '' or current_month == '':
                current_year = year
                current_month = month
            if current_year == year and current_month == month:
                current_month_count += 1
            else:
                # save
                plasmids_by_month['date'].append(json_serial(datetime(current_year, current_month, 1)))
                plasmids_by_month['plasmid_month_count'].append(current_month_count)
                # update current
                current_year = year
                current_month = month
                current_month_count = 1
            if plasmid.computed_size:
                plasmids_with_sequence += 1
            if plasmid.glycerolstock_set.all().count():
                plasmids_with_gs += 1

            key = str(plasmid.type)
            if key not in plasmids_by_type:
                plasmids_by_type[key] = 0
            plasmids_by_type[key] += 1

            key = "Level " + str(plasmid.level)
            if key not in plasmids_by_level:
                plasmids_by_level[key] = 0
            plasmids_by_level[key] += 1

        # save last
        plasmids_by_month['date'].append(json_serial(datetime(current_year, current_month, 1)))
        plasmids_by_month['plasmid_month_count'].append(current_month_count)

        stats.plasmids_by_month = plasmids_by_month
        stats.plasmids_with_sequence = {'values': [plasmids_with_sequence, plasmid_count - plasmids_with_sequence],
                                        'names': ['With sequence', 'Without sequence']}
        stats.plasmids_with_gs = {'values': [plasmids_with_gs, plasmid_count - plasmids_with_gs],
                                  'names': ['With GStock', 'Without Gstock']}

        stats.plasmids_by_type = {'values': list(plasmids_by_type.values()), 'names': list(plasmids_by_type.keys())}
        stats.plasmids_by_level = {'values': list(plasmids_by_level.values()), 'names': list(plasmids_by_level.keys())}
        stats.last_update = date.today()
        stats.plasmid_count = plasmid_count
        stats.save()

    if stats.plasmids_by_month:
        df = pd.DataFrame(data=stats.plasmids_by_month)
        fig_plasmid_month_count = px.line(df, x="date", y="plasmid_month_count", text="plasmid_month_count",
                                          title="Plasmid creation",
                                          labels={'date': 'Date', 'plasmid_month_count': '# created plasmids'})
        fig_plasmid_month_count.update_traces(textposition="top center")

        df = pd.DataFrame(data=stats.plasmids_with_sequence)
        fig_plasmids_with_sequence = px.pie(df, values='values', names='names', hole=.3, title="Plasmids Sequence")

        df = pd.DataFrame(data=stats.plasmids_with_gs)
        fig_plasmids_with_gs = px.pie(df, values='values', names='names', hole=.3, title="Plasmids GStock")

        df = pd.DataFrame(data=stats.plasmids_by_type)
        fig_plasmids_by_type = px.pie(df, values='values', names='names', hole=.3, title="Type")

        df = pd.DataFrame(data=stats.plasmids_by_level)
        fig_plasmids_by_level = px.pie(df, values='values', names='names', hole=.3, title="Level")

        context = {
            'fig_plasmid_month_count': fig_plasmid_month_count.to_html(),
            'fig_plasmids_with_sequence': fig_plasmids_with_sequence.to_html(),
            'fig_plasmids_with_gs': fig_plasmids_with_gs.to_html(),
            'fig_plasmids_by_type': fig_plasmids_by_type.to_html(),
            'fig_plasmids_by_level': fig_plasmids_by_level.to_html(),
            'last_update': stats.last_update,
            'plasmid_count': stats.plasmid_count
        }
    return render(request, 'inventory/services/stats/stats.html', context)


def ServicesGtr(request):
    return render(request, 'inventory/services/gtr/gtr.html')


def ServicesL0d(request):
    return render(request, 'inventory/services/l0d/l0d.html')


def ServicesBlast(request):
    # collect fastas

    records = []
    for plasmid in Plasmid.objects.all():
        record_output = seqio_get(plasmid)
        if record_output[0]:
            record_output[1].id = str(len(records) + 1) + "-" + plasmid.name.replace(" ", "_")
            record_output[1].description = str(plasmid.id)
            records.append(record_output[1])

    database_folder_path = os.path.join(settings.BASE_DIR, 'blast')
    database_path = os.path.join(database_folder_path, 'all')
    error = ""
    alignment_output = ""
    database_creation_result = ""
    last_update = "No database yet"

    if not os.path.exists(database_path + ".nhr") or (request.method == 'POST' and 'makeblastdb' in request.POST):
        with tempfile.NamedTemporaryFile() as tmp:
            SeqIO.write(records, tmp.name, "fasta")
            cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=tmp.name, out=database_path, parse_seqids=True)
            stdout, stderr = cline()
            database_creation_result = "New BLAST database created with " + str(len(records)) + " plasmids"

    if os.path.exists(database_path + ".nhr"):
        # db exists
        last_update = datetime.fromtimestamp(os.stat(database_path + ".nhr").st_mtime)
        alignments = []
        if request.method == 'POST' and 'doblast' in request.POST:
            form = BlastSequenceInput(request.POST)
            if form.is_valid():
                query_seq = form.cleaned_data['sequence_input']
                query_record = SeqRecord(Seq(query_seq), id="query_seq")
                with tempfile.NamedTemporaryFile() as tmp:
                    SeqIO.write(query_record, tmp.name, "fasta")
                    with tempfile.NamedTemporaryFile() as xml_output:
                        cline = NcbiblastnCommandline(query=tmp.name, db=database_path, evalue=0.001,
                                                      out=xml_output.name, outfmt=5)
                        stdout, stderr = cline()

                        records = NCBIXML.parse(open(xml_output.name, "r"))
                        item = next(records)
                        chunk_size = 200
                        for alignment in item.alignments:
                            for hsp in alignment.hsps:
                                if hsp.expect < 0.01:
                                    alignment.hsp = hsp
                                    alignment.chunk = []
                                    alignment.plasmid = Plasmid.objects.get(id=alignment.title.split(' ')[1])
                                    for i in range(0, len(hsp.query), chunk_size):
                                        alignment.chunk.append(
                                            [
                                                hsp.query[i:i + chunk_size],
                                                hsp.match[i:i + chunk_size],
                                                hsp.sbjct[i:i + chunk_size]
                                            ]
                                        )
                                    alignments.append(alignment)
                        alignment_output = True
        else:
            form = BlastSequenceInput(request.POST)
    else:
        error = "No database exists yet"

    context = {
        'last_update': last_update,
        'error': error,
        'form': form,
        'alignment_output': alignment_output,
        'alignments': alignments,
        'database_creation_result': database_creation_result
    }
    return render(request, 'inventory/services/blast/blast.html', context)


def get_plasmid_type_id(plasmid):
    if plasmid.type:
        return plasmid.type.id
    return None


def api_plasmid_getfasta_byid(request, name):
    plasmids = Plasmid.objects.all()
    for plasmid in plasmids:
        if plasmid.name == name:
            sequence = ""
            result = grab_seq(plasmid)
            if result[0]:
                sequence = result[1]
            return JsonResponse({
                'name': plasmid.name,
                'id': plasmid.id,
                'seq': str(sequence)
            })
    return JsonResponse({
        'error': 'Plasmid not found'
    })


def api_plasmids(request):
    output = []
    level_from_table_filters = 0
    level_to_table_filters = 0
    for plasmid in Plasmid.objects.all():
        plasmid_grab_seq = grab_seq(plasmid)
        gs_out = []
        for gs in plasmid.glycerolstock_set.all():
            gs_out.append({
                'i': str(gs.id),
                's': str(gs.strain),
                'bc': gs.box_column,
                'br': gs.box_row,
                'b': str(gs.box),
            })

        output.append({
            'n': plasmid.name,
            'l': plasmid.level,
            't': get_plasmid_type_id(plasmid),
            'i': plasmid.id,
            'hs': plasmid_grab_seq[0],
            'c': plasmid.computed_size,
            'ic': plasmid.insert_computed_size,
            'cs': plasmid.check_state,
            'ss': plasmid.sequencing_state,
            'r': recommended_enzyme_for_create(plasmid.level),
            'g': gs_out,
            'rh': " / ".join(list(plasmid.resistances.all().values_list('three_letter_code', flat=True)))
        })
        if plasmid.level:
            if plasmid.level > level_to_table_filters:
                level_to_table_filters = plasmid.level
            if plasmid.level < level_from_table_filters:
                level_from_table_filters = plasmid.level
    context = {
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
        'plasmids': output,
        'has_perm_to_edit': request.user.has_perm('inventory.change_plasmid'),
        'csrf_token': django.middleware.csrf.get_token(request)
    }
    return JsonResponse(context, safe=False)


def api_glycerolstocks(request):
    output = []
    level_from_table_filters = 0
    level_to_table_filters = 0
    for glycerolstock in GlycerolStock.objects.all():
        pi = ""
        pn = ""
        pt = ""
        pl = ""
        if glycerolstock.plasmid:
            if glycerolstock.plasmid.level:
                if glycerolstock.plasmid.level > level_to_table_filters:
                    level_to_table_filters = glycerolstock.plasmid.level
                if glycerolstock.plasmid.level < level_from_table_filters:
                    level_from_table_filters = glycerolstock.plasmid.level
            pi = glycerolstock.plasmid.id
            pn = glycerolstock.plasmid.name
            if glycerolstock.plasmid.type:
                pt = glycerolstock.plasmid.type.id
            pl = glycerolstock.plasmid.level
        bn = ""
        bl = ""
        if glycerolstock.box:
            bn = glycerolstock.box.name
            bl = str(glycerolstock.box.location)
        output.append({
            'n': str(glycerolstock),
            'i': glycerolstock.id,
            'pi': pi,
            'pn': pn,
            'pt': pt,
            'pl': pl,
            's': str(glycerolstock.strain),
            'bc': glycerolstock.box_column,
            'br': glycerolstock.box_row,
            'bn': bn,
            'bl': bl
        })
    context = {
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
        'glycerolstocks': output,
        'has_perm_to_edit': request.user.has_perm('inventory.change_glycerolstock'),
    }
    return JsonResponse(context, safe=False)


def createEnzymeFromName(enzyme_name):
    if enzyme_name in rest_dict:
        return RestrictionEnzyme(name=enzyme_name)
    return None


def api_parts(request, enzyme_name, assembly_standard):
    parts = []
    api_error = ""

    the_re = createEnzymeFromName(enzyme_name)
    if not the_re:
        api_error = 'API Error / Restriction enzyme not found (' + enzyme_name + ')'
    else:
        for plasmid in Plasmid.objects.order_by('name'):
            if plasmid.level is not None and plasmid.type is not None and plasmid.type.id == 1 and assembly_standard == 'loop' and (
                    (enzyme_name == 'BsaI' and plasmid.level % 2 == 0) or (
                    enzyme_name == 'SapI' and plasmid.level % 2 == 1)):
                # make sure use correct enzyme at loop
                continue
            plasmid_grab_seq = grab_seq(plasmid)
            if plasmid.level is not None and plasmid_grab_seq[0]:
                found_cut_positions = re_find_cut_positions(plasmid_grab_seq[1], the_re, True, True)
                if len(found_cut_positions) == 2:
                    length = found_cut_positions[1] - found_cut_positions[0]
                    oh5 = str(plasmid_grab_seq[1])[
                          found_cut_positions[0] - 1:found_cut_positions[0] + abs(the_re.fcut - the_re.rcut) - 1]
                    oh3 = str(plasmid_grab_seq[1])[
                          found_cut_positions[1] - 1:found_cut_positions[1] + abs(the_re.fcut - the_re.rcut) - 1]
                    if get_plasmid_type_id(plasmid) == 1:
                        # receiver
                        ohtmp = oh3
                        oh3 = oh5
                        oh5 = ohtmp
                        length = len(str(plasmid_grab_seq[1])) - length
                    parts.append({
                        'n': plasmid.name,
                        'd': plasmid.description,
                        'l': plasmid.level,
                        't': get_plasmid_type_id(plasmid),
                        'i': plasmid.id,
                        'len': length,
                        'o5': oh5,
                        'o3': oh3,
                    })
    context = {
        'error': api_error,
        'parts': parts,
        'csrf_token': django.middleware.csrf.get_token(request),
    }
    return JsonResponse(context, safe=False)
