from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.utils.decorators import method_decorator
from django.shortcuts import redirect

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
from organization.decorators import require_current_project_set
from organization.decorators import require_member_can_write_or_admin_current_project
from organization.decorators import require_member_can_read_project_of_plasmid
from organization.decorators import require_member_can_read_project_of_gs
from organization.decorators import require_member_can_read_project_of_primer
from organization.decorators import require_member_can_write_or_admin_project_of_plasmid
from organization.decorators import require_member_can_write_or_admin_project_of_gs
from organization.decorators import require_member_can_write_or_admin_project_of_primer
from organization.views import has_current_project
from organization.views import get_current_project_id
from organization.views import get_current_project
from organization.views import on_current_project_member_can_write_or_admin
from organization.views import get_projects_where_member_can_any
from organization.views import member_can_write_or_admin_plasmid
from organization.views import member_can_write_or_admin_gs
from organization.views import get_show_from_all_projects
from organization.views import member_can_write_or_admin_primer

from .forms import PlasmidValidationForm
from .forms import PlasmidCreateForm
from .forms import PlasmidEditForm

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
from django.core.exceptions import PermissionDenied

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
from .forms import GstockCreateForm
from .forms import GstockEditForm
from .forms import GlycerolQRInput
from .forms import PlasmidNameInput

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


def restrictionenzymes(request):
    context = {
        'restrictionenzymes': RestrictionEnzyme.objects.all(),
    }
    return render(request, 'inventory/restrictionenzymes.html', context)


@require_current_project_set
def glycerolstocks(request):
    show_from_all_projects = get_show_from_all_projects(request)
    if show_from_all_projects:
        glycerolstocks = GlycerolStock.objects.filter(project_id__in=get_projects_where_member_can_any(request.user))
    else:
        glycerolstocks = GlycerolStock.objects.filter(project_id=get_current_project_id(request))

    level_from_table_filters = 0
    level_to_table_filters = 0
    hasGlycerolStocks = False
    for glycerolstock in glycerolstocks:
        hasGlycerolStocks = True
        if glycerolstock.plasmid:
            if glycerolstock.plasmid.level:
                if glycerolstock.plasmid.level > level_to_table_filters:
                    level_to_table_filters = glycerolstock.plasmid.level
                if glycerolstock.plasmid.level < level_from_table_filters:
                    level_from_table_filters = glycerolstock.plasmid.level
    context = {
        'on_current_project_member_can_write_or_admin': on_current_project_member_can_write_or_admin(request),
        'has_glycerolstocks': hasGlycerolStocks,
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
        'show_from_all_projects': show_from_all_projects
    }
    return render(request, 'inventory/glycerolstocks.html', context)


@require_member_can_read_project_of_gs
def glycerolstock(request, glycerolstock_id):
    try:
        glycerolstock_to_detail = GlycerolStock.objects.get(id=glycerolstock_id)
    except ObjectDoesNotExist:
        raise Http404

    resistantes_human_context = "None"
    if glycerolstock_to_detail.plasmid:
        resistantes_human_context = resistantes_human(glycerolstock_to_detail.plasmid.selectable_markers)

    context = {
        'glycerolstock': glycerolstock_to_detail,
        'resistantes_human': resistantes_human_context,
        'user_can_edit_gs': member_can_write_or_admin_gs(glycerolstock_to_detail, request.user),
        'resistantes_strain_human': resistantes_human(glycerolstock_to_detail.strain.selectable_markers),
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


@require_member_can_read_project_of_gs
def glycerolstock_from_qr(request, glycerolstock_id):
    try:
        glycerolstock_to_detail = GlycerolStock.objects.filter(qr_id=glycerolstock_id)[0]
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'glycerolstock': glycerolstock_to_detail,
    }
    return render(request, 'inventory/glycerolstock.html', context)


def gstock_check_pos(the_class, the_self, the_form):
    try:
        obj_curr_pos = GlycerolStock.objects.get(box_row=the_form.cleaned_data['box_row'],
                                                    box_column=the_form.cleaned_data['box_column'],
                                                    box=the_form.cleaned_data['box'])
        if obj_curr_pos:
            if obj_curr_pos != the_form.instance:
                the_form.add_error(None, 'Box position not available')
                return super(the_class, the_self).form_invalid(the_form)
        return super(the_class, the_self).form_valid(the_form)
    except GlycerolStock.DoesNotExist:
        return super(the_class, the_self).form_valid(the_form)


def build_boxes(request, mode):
    output = {
        'BOX_ROWS': BOX_ROWS,
        'BOX_COLUMNS': BOX_COLUMNS,
        'locations': []
    }
    for location in Location.objects.all():
        boxes = []
        for box in Box.objects.filter(location=location).order_by('name'):
            box_output = {
                'name': box.name,
                'id': box.id
            }
            if mode == 'p':
                glycerolstocks = box.glycerolstock_set.all()
            else:
                if get_show_from_all_projects(request):
                    glycerolstocks = box.glycerolstock_set.filter(project_id__in=get_projects_where_member_can_any(request.user))
                else:
                    glycerolstocks = box.glycerolstock_set.filter(project_id=get_current_project_id(request))

            anyGs = False
            for glycerolstock in glycerolstocks:
                box_output[str(glycerolstock.box_row) + str(glycerolstock.box_column)] = glycerolstock
                anyGs = True

            if anyGs or mode == 'p':
                boxes.append(box_output)
        if boxes or mode == 'p':
            output['locations'].append({
                'name': location.name,
                'boxes': boxes
            })
    return output


class GstockEdit(UpdateView):
    model = GlycerolStock
    template_name_suffix = '_update_form'
    form_class = GstockEditForm

    @method_decorator(require_member_can_write_or_admin_project_of_gs)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def form_valid(self, form):
        return gstock_check_pos(GstockEdit, self, form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["collection"] = build_boxes(self.request, 'p')
        context["render_mod"] = 'p'
        return context

    def get_success_url(self, **kwargs):
        return reverse('glycerolstock', args=(self.object.id,)) + '?form_result_glycerolstock_edit_success=true'


class GstockCreate(CreateView):
    model = GlycerolStock
    template_name_suffix = '_create_form'
    form_class = GstockCreateForm

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def form_valid(self, form):
        form.instance.project = get_current_project(self.request)
        return gstock_check_pos(GstockCreate, self, form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["collection"] = build_boxes(self.request, 'p')
        context["render_mod"] = 'p'
        return context

    def get_success_url(self, **kwargs):
        return reverse('glycerolstock', args=(self.object.id,)) + '?form_result_glycerolstock_create_success=true'


class GstockCreatePlasmidDefined(CreateView):
    model = GlycerolStock
    form_class = GstockCreateForm
    template_name_suffix = '_create_form'

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def form_valid(self, form):
        form.instance.project = get_current_project(self.request)
        return gstock_check_pos(GstockCreatePlasmidDefined, self, form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["collection"] = build_boxes(self.request, 'p')
        context["plasmid_id"] = self.kwargs['pid']
        context["render_mod"] = 'p'
        return context

    def get_success_url(self, **kwargs):
        return reverse('glycerolstock', args=(self.object.id,)) + '?form_result_glycerolstock_create_success=true'


class GstockDelete(DeleteView):
    model = GlycerolStock

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_success_url(self, **kwargs):
        return reverse('glycerolstocks') + '?form_result_object_deleted=true'


@require_member_can_read_project_of_gs
def glycerolstock_label(request, glycerolstock_id):
    try:
        glycerolstock_to_label = GlycerolStock.objects.get(id=glycerolstock_id)
    except ObjectDoesNotExist:
        raise Http404

    resistantes_human_context = "None"
    if glycerolstock_to_label.plasmid:
        resistantes_human_context = resistantes_human(glycerolstock_to_label.plasmid.selectable_markers, True)

    context = {
        'glycerolstock': glycerolstock_to_label,
        'resistantes_human': resistantes_human_context,
        'resistantes_strain_human': resistantes_human(glycerolstock_to_label.strain.selectable_markers, True),
    }
    return render(request, 'inventory/glycerolstock_label.html', context)


def glycerolstock_boxes(request):
    context = {
        'collection': build_boxes(request, 'n'),
        'render_mod': 'n',
        'show_from_all_projects': get_show_from_all_projects(request)
    }
    return render(request, 'inventory/glycerolstock_boxes.html', context)


@require_current_project_set
def plasmids(request):
    show_from_all_projects = get_show_from_all_projects(request)
    if show_from_all_projects:
        plasmids = Plasmid.objects.filter(project_id__in=get_projects_where_member_can_any(request.user))
    else:
        plasmids = Plasmid.objects.filter(project_id=get_current_project_id(request))
    level_from_table_filters = 0
    level_to_table_filters = 0
    hasPlasmids = False
    for plasmid in plasmids:
        hasPlasmids = True
        plasmid.refc = recommended_enzyme_for_create(plasmid.level)
        if plasmid.level:
            if plasmid.level > level_to_table_filters:
                level_to_table_filters = plasmid.level
            if plasmid.level < level_from_table_filters:
                level_from_table_filters = plasmid.level
    context = {
        'on_current_project_member_can_write_or_admin': on_current_project_member_can_write_or_admin(request),
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
        'has_plasmids': hasPlasmids,
        'RESTRICTION_ENZYMES': RestrictionEnzyme.objects.all,
        'show_from_all_projects': show_from_all_projects
    }
    return render(request, 'inventory/plasmids.html', context)


def getPlasmidResistanceForLigation(plasmid):
    if plasmid.level:
        if plasmid.level % 2 == 0:
            is_spe = False
            for resistance in plasmid.selectable_markers.all():
                if resistance.three_letter_code == 'SPE':
                    is_spe = True
                    break
            if is_spe:
                return 'SPE'
            else:
                return 'Error: SPE not in resistance list'
        else:
            is_kan = False
            for resistance in plasmid.selectable_markers.all():
                if resistance.three_letter_code == 'KAN':
                    is_kan = True
                    break
            if is_kan:
                return 'KAN'
            else:
                return 'Error: KAN not in resistance list'
    if len(plasmid.selectable_markers):
        return plasmid.selectable_markers[0]
    else:
        return 'More than one'


def resistantes_human(selectable_markers, short=False):
    resistantes_human_return = []
    if selectable_markers:
        for resistance in selectable_markers.all():
            if short:
                resistantes_human_return.append(str(resistance.three_letter_code))
            else:
                resistantes_human_return.append(resistance.name + " (" + str(resistance.three_letter_code) + ")")
            continue
        return " / ".join(resistantes_human_return)
    else:
        return "None"


@require_member_can_read_project_of_plasmid
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
        'resistantes_human': resistantes_human(plasmid_to_detail.selectable_markers),
        'insert_of': insert_of,
        'backbone_of': backbone_of,
        'ligation_raw': ligation_raw,
        'SEQUENCING_STATES': SEQUENCING_STATES,
        'CHECK_STATES': CHECK_STATES,
        'CHECK_METHODS': CHECK_METHODS,
        'RESTRICTION_ENZYMES': RestrictionEnzyme.objects.all(),
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_detail, request.user)
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


@require_member_can_read_project_of_plasmid
@require_member_can_write_or_admin_current_project
def plasmid_duplicate(request, plasmid_id):
    try:
        plasmid_to_duplicate = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'plasmid': plasmid_to_duplicate
    }
    if 'plasmid_name' in request.POST:
        form = PlasmidNameInput(request.POST)
        context['form'] = form
        if form.is_valid():
            selectable_markers = plasmid_to_duplicate.selectable_markers.all()
            inserts = plasmid_to_duplicate.inserts.all()
            plasmid_to_duplicate.pk = None
            plasmid_to_duplicate.created_on = datetime.now().date()
            plasmid_to_duplicate.working_colony = None
            plasmid_to_duplicate.check_state = 1
            plasmid_to_duplicate.check_method = None
            plasmid_to_duplicate.check_date = None
            plasmid_to_duplicate.digestion_check_enzymes = None
            plasmid_to_duplicate.check_observations = None
            plasmid_to_duplicate.sequencing_state = 0
            plasmid_to_duplicate.sequencing_date = None
            plasmid_to_duplicate.sequencing_observations = None
            plasmid_to_duplicate.name = form.cleaned_data['plasmid_name']
            plasmid_to_duplicate.save()
            plasmid_to_duplicate.selectable_markers.add(*selectable_markers)
            plasmid_to_duplicate.inserts.add(*inserts)
            return redirect('plasmid', plasmid_id=plasmid_to_duplicate.id)
        return render(request, 'inventory/plasmid_duplicate.html', context)
    else:
        context['form'] = PlasmidNameInput()
        return render(request, 'inventory/plasmid_duplicate.html', context)

@require_member_can_read_project_of_plasmid
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
    form_class = PlasmidEditForm
    template_name_suffix = '_update_form'

    @method_decorator(require_member_can_write_or_admin_project_of_plasmid)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_form_kwargs(self):
        kwargs = super(PlasmidEdit, self).get_form_kwargs()
        kwargs.update({'user': self.request.user})
        return kwargs

    def get_success_url(self, **kwargs):
        return reverse('plasmid', args=(self.object.id,)) + '?form_result_plasmid_edit_success=true'


class PlasmidCreate(CreateView):
    model = Plasmid
    form_class = PlasmidCreateForm
    template_name_suffix = '_create_form'

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_form_kwargs(self):
        kwargs = super(PlasmidCreate, self).get_form_kwargs()
        kwargs.update({'user': self.request.user})
        return kwargs

    def form_valid(self, form):
        form.instance.project = get_current_project(self.request)
        return super().form_valid(form)

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

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["create_url"] = reverse('plasmid_create')
        context["next_url"] = reverse('plasmid_create_wizard_end')
        return context

    def form_valid(self, form):
        form.instance.project = get_current_project(self.request)
        return super().form_valid(form)

    def get_success_url(self, **kwargs):
        return reverse('plasmid', args=(self.object.id,))


@require_member_can_write_or_admin_current_project
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
                    project=get_current_project(request),
                )
                for i in params['i']:
                    plasmid_created.inserts.add(i)
                for r in backbone.selectable_markers.all():
                    plasmid_created.selectable_markers.add(r.id)
                plasmid_create_from_inserts(plasmid_created, context)
                return redirect('plasmid', plasmid_id=plasmid_created.id)
            else:
                context['wizard_error'] = 'Backbone not found.'
        else:
            context['wizard_error'] = 'Name & backbone & inserts are required fields.'
    else:
        context['wizard_error'] = "Plasmid can\'t be created. No parameters set."
    return render(request, 'inventory/plasmid.html', context)



@csrf_exempt
@require_member_can_read_project_of_plasmid
def plasmid_view_edit(request, plasmid_id):
    try:
        plasmid_to_detail = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    warnings = []

    if member_can_write_or_admin_plasmid(plasmid_to_detail, request.user):
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
                plasmid_to_detail.sequence.save(plasmid_to_detail.name + ".gb", ContentFile(
                    '''LOCUS       ''' + plasmid_to_detail.name.replace(" ", "_") + '''                   0 bp    DNA     circular  31-AUG-2021
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


@require_member_can_read_project_of_plasmid
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


@require_member_can_read_project_of_plasmid
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


@require_member_can_read_project_of_plasmid
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


@require_member_can_read_project_of_plasmid
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


@require_member_can_read_project_of_plasmid
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


class PlasmidDelete(DeleteView):
    model = Plasmid

    @method_decorator(require_member_can_write_or_admin_project_of_plasmid)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_object(self, *args, **kwargs):
        obj = super(PlasmidDelete, self).get_object(*args, **kwargs)
        if not member_can_write_or_admin_plasmid(obj, self.request.user):
            raise PermissionDenied()
        return obj

    def get_success_url(self, **kwargs):
        return reverse('plasmids') + '?form_result_object_deleted=true'


@require_current_project_set
@require_member_can_write_or_admin_project_of_plasmid
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


@require_current_project_set
def PlasmidValidations(request):
    plasmidsToCheck = []
    plasmidsToStock = []
    # plasmidsWithStockWoCheck = []

    show_from_all_projects = get_show_from_all_projects(request)
    if show_from_all_projects:
        all_plasmids = Plasmid.objects.filter(project_id__in=get_projects_where_member_can_any(request.user))
    else:
        all_plasmids = Plasmid.objects.filter(project_id=get_current_project_id(request))

    # Plasmid w/o GS && Check state == Pending
    for plasmid in all_plasmids.filter(check_state=1):
        if plasmid.glycerolstock_set.all().count() == 0:
            plasmidsToCheck.append(plasmid)

    # Plasmid w/o GS && Check state != Pending && Sequencing state != Required
    for plasmid in all_plasmids.exclude(check_state=1).exclude(sequencing_state=1):
        if plasmid.glycerolstock_set.all().count() == 0:
            plasmidsToStock.append(plasmid)

    # Plasmid w/ GS && Check state == Pending && Sequencing state != Required
    # for plasmid in all_plasmids.filter(check_state=1).exclude(sequencing_state=1):
    #    if plasmid.glycerolstock_set.all().count() != 0:
    #        plasmidsWithStockWoCheck.append(plasmid)

    context = {
        'plasmidsToCheck': plasmidsToCheck,
        # Sequencing state = Required
        'plasmidsToSequence': all_plasmids.filter(sequencing_state=1),
        'plasmidsToStock': plasmidsToStock,
        # 'plasmidsWithStockWoCheck': plasmidsWithStockWoCheck,
        'CHECK_METHODS': CHECK_METHODS,
        'show_from_all_projects': show_from_all_projects
    }
    return render(request, 'inventory/plasmidvalidations.html', context)



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

@require_member_can_read_project_of_primer
def primer(request, primer_id):
    try:
        primer_to_detail = Primer.objects.get(id=primer_id)
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'primer': primer_to_detail,
        'user_can_edit_primer': member_can_write_or_admin_primer(primer_to_detail, request.user)
    }
    return render(request, 'inventory/primer.html', context)


def primers(request):
    show_from_all_projects = get_show_from_all_projects(request)
    if show_from_all_projects:
        primers = Primer.objects.filter(project_id__in=get_projects_where_member_can_any(request.user))
    else:
        primers = Primer.objects.filter(project_id=get_current_project_id(request))
    for primer in primers:
        primer.can_edit = member_can_write_or_admin_primer(primer, request.user)
    context = {
        'primers': primers,
        'show_from_all_projects': show_from_all_projects,
        'on_current_project_member_can_write_or_admin': on_current_project_member_can_write_or_admin(request)
    }
    return render(request, 'inventory/primers.html', context)


class PrimerCreate(CreateView):
    model = Primer
    fields = '__all__'
    template_name_suffix = '_create_form'

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_success_url(self, **kwargs):
        return reverse('primer', args=(self.object.id,)) + '?form_result_primer_create_success=true'


class PrimerEdit(UpdateView):
    model = Primer
    fields = '__all__'
    template_name_suffix = '_update_form'

    @method_decorator(require_member_can_write_or_admin_project_of_primer)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_success_url(self, **kwargs):
        return reverse('primer', args=(self.object.id,)) + '?form_result_primer_edit_success=true'


class PrimerDelete(DeleteView):
    model = Primer

    @method_decorator(require_member_can_write_or_admin_project_of_primer)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_success_url(self, **kwargs):
        return reverse('primers') + '?form_result_object_deleted=true'


@require_member_can_read_project_of_primer
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
    for plasmid in Plasmid.objects.filter(project__in=get_projects_where_member_can_any(request.user)).order_by('name'):
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
    plasmid = Plasmid.objects.get(name=name)
    if plasmid:
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
    if has_current_project(request):
        output = []
        level_from_table_filters = 0
        level_to_table_filters = 0
        if get_show_from_all_projects(request):
            plasmids = Plasmid.objects.filter(project_id__in=get_projects_where_member_can_any(request.user)).order_by('name')
        else:
            plasmids = Plasmid.objects.filter(project_id=get_current_project_id(request)).order_by('name')
        for plasmid in plasmids:
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
                'sm': " + ".join(list(plasmid.selectable_markers.all().values_list('three_letter_code', flat=True))),
                'p': member_can_write_or_admin_plasmid(plasmid, request.user)
            })
            if plasmid.level:
                if plasmid.level > level_to_table_filters:
                    level_to_table_filters = plasmid.level
                if plasmid.level < level_from_table_filters:
                    level_from_table_filters = plasmid.level
        context = {
            'current_project': True,
            'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
            'plasmids': output,
            'csrf_token': django.middleware.csrf.get_token(request)
        }
        return JsonResponse(context, safe=False)
    return JsonResponse({})


def api_glycerolstocks(request):
    output = []
    level_from_table_filters = 0
    level_to_table_filters = 0
    if get_show_from_all_projects(request):
        glycerolstocks = GlycerolStock.objects.filter(project_id__in=get_projects_where_member_can_any(request.user)).order_by('strain', 'plasmid')
    else:
        glycerolstocks = GlycerolStock.objects.filter(project_id=get_current_project_id(request)).order_by('strain', 'plasmid')
    for glycerolstock in glycerolstocks:
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
            'bl': bl,
            'p': member_can_write_or_admin_gs(glycerolstock, request.user)
        })
    context = {
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
        'glycerolstocks': output,
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
        for plasmid in Plasmid.objects.filter(project__in=get_projects_where_member_can_any(request.user)).order_by(
                'name'):
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
