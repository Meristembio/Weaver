from . import views
from .views import PlasmidCreate
from .views import PlasmidCreateWizard
from .views import plasmid_create_wizard_end
from .views import PlasmidEdit
from .views import PlasmidDelete
from .views import PlasmidCreateL0d

from .views import GstockCreate
from .views import GstockCreatePlasmidDefined
from .views import GlycerolstockBatchPrepare
from .views import GstockEdit
from .views import GstockDelete

from .views import PrimerCreate
from .views import PrimerEdit
from .views import PrimerDelete
from .views import RestrictionEnzymeCreate
from .views import RestrictionEnzymeEdit
from .views import RestrictionEnzymeDelete

from .views import PlasmidValidations
from .views import PlasmidValidationEdit

from .views import ServicesBlast
from .views import ServicesGtr
from .views import ServicesL0d
from .views import ServicesSanger
from .views import ServicesStats
from .views import ServicesBatchPrints
from .views import api_plasmid_get_fasta_by_name
from .views import api_plasmid_get_fasta_by_idx
from .views import api_fidelity_calc


from django.contrib.auth.decorators import login_required
from django.urls import path

urlpatterns = [
    path('plasmids/', login_required(views.plasmids, redirect_field_name='next'), name='plasmids'),
    path('plasmids/validations/', login_required(PlasmidValidations, redirect_field_name='next'), name='plasmid_validations'),
    path('plasmid/<uuid:plasmid_id>/', login_required(views.plasmid, redirect_field_name='next'), name='plasmid'),
    path('qr/i/<str:plasmid_id>/', login_required(views.plasmid_from_qr, redirect_field_name='next'), name='plasmid_from_qr'), # for QR
    path('plasmid/create/', login_required(PlasmidCreate.as_view(), redirect_field_name='next'), name='plasmid_create'),
    path('plasmid/import/', login_required(views.plasmid_import, redirect_field_name='next'), name='plasmid_import'),
    path('plasmid/create/l0d', login_required(PlasmidCreateL0d, redirect_field_name='next'), name='plasmid_create_l0d'),
    path('plasmid/create/wizard', login_required(PlasmidCreateWizard.as_view(), redirect_field_name='next'), name='plasmid_create_wizard'),
    path('plasmid/create/wizard/end', login_required(plasmid_create_wizard_end, redirect_field_name='next'), name='plasmid_create_wizard_end'),
    path('plasmid/edit/<uuid:pk>/', login_required(PlasmidEdit.as_view(), redirect_field_name='next'), name='plasmid_edit'),
    path('plasmid/delete/<uuid:pk>/', login_required(PlasmidDelete.as_view(), redirect_field_name='next'), name='plasmid_delete'),
    path('plasmid/deleted/', login_required(views.plasmid_deleted, redirect_field_name='next'), name='plasmid_deleted'),
    path('plasmid/label/<uuid:plasmid_id>/', login_required(views.plasmid_label, redirect_field_name='next'), name='plasmid_label'),
    path('plasmid/download/<uuid:plasmid_id>', login_required(views.plasmid_download, redirect_field_name='next'), name='plasmid_download'),
    path('plasmid/download/clustal/<uuid:plasmid_id>', login_required(views.plasmid_download_clustal, redirect_field_name='next'), name='plasmid_download_clustal'),
    path('plasmid/view_edit/<uuid:plasmid_id>', login_required(views.plasmid_view_edit, redirect_field_name='next'), name='plasmid_view_edit'),
    path('plasmid/seq_verification_<int:weaver_id>', login_required(views.plasmid_seq_verification_entry, redirect_field_name='next'), name='plasmid_seq_verification_entry'),
    path('plasmid/align/<uuid:plasmid_id>', login_required(views.plasmid_align_sanger, redirect_field_name='next'), name='plasmid_align'),
    path('plasmid/align/sanger/<uuid:plasmid_id>', login_required(views.plasmid_align_sanger, redirect_field_name='next'), name='plasmid_align_sanger'),
    path('plasmid/align/sanger/<uuid:plasmid_id>/run/<uuid:run_id>', login_required(views.sanger_run_detail, redirect_field_name='next'), name='sanger_run_detail'),
    path('plasmid/align/sanger/<uuid:plasmid_id>/run/<uuid:run_id>/read/<uuid:read_id>/chromatogram', login_required(views.sanger_read_chromatogram, redirect_field_name='next'), name='sanger_read_chromatogram'),
    path('plasmid/align/sanger/<uuid:plasmid_id>/run/<uuid:run_id>/read/<uuid:read_id>/file/<uuid:file_id>/download', login_required(views.sanger_read_file_download, redirect_field_name='next'), name='sanger_read_file_download'),
    path('plasmid/align/sanger/<uuid:plasmid_id>/run/<uuid:run_id>/download/<str:kind>', login_required(views.sanger_run_download, redirect_field_name='next'), name='sanger_run_download'),
    path('plasmid/align/sanger/<uuid:plasmid_id>/run/<uuid:run_id>/decision', login_required(views.sanger_run_decision, redirect_field_name='next'), name='sanger_run_decision'),
    path('plasmid/align/sanger/<uuid:plasmid_id>/run/<uuid:run_id>/rerun', login_required(views.sanger_run_rerun, redirect_field_name='next'), name='sanger_run_rerun'),
    path('plasmid/align/sanger/<uuid:plasmid_id>/run/<uuid:run_id>/delete', login_required(views.sanger_run_delete, redirect_field_name='next'), name='sanger_run_delete'),
    path('plasmid/align/fasta/<uuid:plasmid_id>', login_required(views.plasmid_align_fasta, redirect_field_name='next'), name='plasmid_align_fasta'),
    path('plasmid/digest/<uuid:plasmid_id>', login_required(views.plasmid_digest, redirect_field_name='next'), name='plasmid_digest'),
    path('plasmid/pcr/<uuid:plasmid_id>', login_required(views.plasmid_pcr, redirect_field_name='next'), name='plasmid_pcr'),
    path('plasmid/validation/edit/<uuid:plasmid_id>/', login_required(PlasmidValidationEdit, redirect_field_name='next'), name='plasmid_validation_edit'),
    path('plasmid/duplicate/<uuid:plasmid_id>', login_required(views.plasmid_duplicate, redirect_field_name='next'), name='plasmid_duplicate'),

    path('public/plasmid/<uuid:plasmid_id>/', views.plasmid_public, name='plasmid_public'),

    path('glycerolstocks/', login_required(views.glycerolstocks, redirect_field_name='next'), name='glycerolstocks'),
    path('glycerolstock/<uuid:glycerolstock_id>/', login_required(views.glycerolstock, redirect_field_name='next'), name='glycerolstock'),
    path('qr/g/', login_required(views.glycerolstock_qr, redirect_field_name='next'), name='glycerolstock_qr'),
    path('qr/g/<str:glycerolstock_qr_id>/', login_required(views.glycerolstock_from_qr, redirect_field_name='next'), name='glycerolstock_from_qr'), # for QR
    path('glycerolstock/create/', login_required(GstockCreate.as_view(), redirect_field_name='next'), name='glycerolstock_create'),
    path('glycerolstock/create/<uuid:pid>', login_required(GstockCreatePlasmidDefined.as_view(), redirect_field_name='next'), name='glycerolstock_create_plasmid_defined'),
    path('glycerolstock/batch/prepare/', login_required(GlycerolstockBatchPrepare, redirect_field_name='next'), name='glycerolstock_batch_prepare'),
    path('glycerolstock/edit/<uuid:pk>/', login_required(GstockEdit.as_view(), redirect_field_name='next'), name='glycerolstock_edit'),
    path('glycerolstock/delete/<uuid:pk>/', login_required(GstockDelete.as_view(), redirect_field_name='next'), name='glycerolstock_delete'),
    path('glycerolstock/deleted/', login_required(views.glycerolstock_deleted, redirect_field_name='next'), name='glycerolstock_deleted'),
    path('glycerolstock/label/<uuid:glycerolstock_id>/', login_required(views.glycerolstock_label, redirect_field_name='next'), name='glycerolstock_label'),
    path('glycerolstock/boxes/', login_required(views.glycerolstock_boxes, redirect_field_name='next'), name='glycerolstock_boxes'),

    path('restrictionenzymes/', login_required(views.restrictionenzymes, redirect_field_name='next'), name='restrictionenzymes'),
    path('restrictionenzyme/create/', login_required(RestrictionEnzymeCreate.as_view(), redirect_field_name='next'), name='restrictionenzyme_create'),
    path('restrictionenzyme/edit/<uuid:pk>/', login_required(RestrictionEnzymeEdit.as_view(), redirect_field_name='next'), name='restrictionenzyme_edit'),
    path('restrictionenzyme/delete/<uuid:pk>/', login_required(RestrictionEnzymeDelete.as_view(), redirect_field_name='next'), name='restrictionenzyme_delete'),
    path('restrictionenzyme/<uuid:restrictionenzyme_id>/', login_required(views.restrictionenzyme, redirect_field_name='next'), name='restrictionenzyme'),

    path('primers/', login_required(views.primers, redirect_field_name='next'), name='primers'),
    path('primer/<uuid:primer_id>/', login_required(views.primer, redirect_field_name='next'), name='primer'),
    path('primer/create/', login_required(PrimerCreate.as_view(), redirect_field_name='next'), name='primer_create'),
    path('primer/edit/<uuid:pk>/', login_required(PrimerEdit.as_view(), redirect_field_name='next'), name='primer_edit'),
    path('primer/delete/<uuid:pk>/', login_required(PrimerDelete.as_view(), redirect_field_name='next'), name='primer_delete'),

    path('services/blast/', login_required(ServicesBlast, redirect_field_name='next'), name='services-blast'),
    path('services/gtr/', login_required(ServicesGtr, redirect_field_name='next'), name='services-gtr'),
    path('services/l0d/', login_required(ServicesL0d, redirect_field_name='next'), name='services-l0d'),
    path('services/sanger/', login_required(ServicesSanger, redirect_field_name='next'), name='services-sanger'),
    path('services/sanger/rerun-all/', login_required(views.sanger_rerun_all, redirect_field_name='next'), name='sanger_rerun_all'),
    path('services/batch-prints/', login_required(ServicesBatchPrints, redirect_field_name='next'), name='services-batch-prints'),
    path('services/stats/', login_required(ServicesStats, redirect_field_name='next'), name='services-stats'),

    path('experiments', login_required(views.experiments, redirect_field_name='next'), name='experiments'),
    path('experiment/create/', login_required(views.experiment_create, redirect_field_name='next'), name='experiment_create'),
    path('experiment/<int:experiment_id>/edit/', login_required(views.experiment_edit, redirect_field_name='next'), name='experiment_edit'),
    path('experiment/<int:experiment_id>/delete/', login_required(views.experiment_delete, redirect_field_name='next'), name='experiment_delete'),
    path('api/experiments-map/', login_required(views.api_experiments_map, redirect_field_name='next'), name='api-experiments-map'),

    path('api/plasmid/get_fasta/by_name/<str:name>/', api_plasmid_get_fasta_by_name, name='api-plasmid-get_fasta_by_name'),
    path('api/plasmid/get_fasta/by_idx/<int:idx>/', api_plasmid_get_fasta_by_idx, name='api-plasmid-get_fasta_by_idx'),
    path('api/fidelity_calc/<str:enzyme>/<str:ohs>', login_required(views.api_fidelity_calc, redirect_field_name='next'), name='api-fidelity_calc'),
    path('api/plasmids/', login_required(views.api_plasmids, redirect_field_name='next'), name='api-plasmids'),
    # path('api/plasmids/', views.api_plasmids, name='api-plasmids'),
    path('api/parts/<str:enzyme_name>/<str:assembly_standard>/', login_required(views.api_parts, redirect_field_name='next'), name='api-parts'),
    # path('api/parts/<str:enzyme_name>/<str:assembly_standard>/', views.api_parts, name='api-parts'),
    path('api/glycerolstocks/', login_required(views.api_glycerolstocks, redirect_field_name='next'), name='api-glycerolstocks'),
    # path('api/glycerolstocks/', views.api_glycerolstocks, name='api-glycerolstocks'),
]
