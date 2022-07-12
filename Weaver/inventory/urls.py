from . import views
from .views import PlasmidCreate
from .views import PlasmidCreateWizard
from .views import plasmid_create_wizard_end
from .views import PlasmidEdit
from .views import PlasmidDelete

from .views import GstockCreate
from .views import GstockCreatePlasmidDefined
from .views import GstockEdit
from .views import GstockDelete

from .views import RestrictionenzymeCreate
from .views import RestrictionenzymeEdit
from .views import RestrictionenzymeDelete

from .views import PrimerCreate
from .views import PrimerEdit
from .views import PrimerDelete

from .views import PlasmidValidations
from .views import PlasmidValidationEdit

from .views import ServicesBlast
from .views import ServicesGtr
from .views import ServicesL0d
from .views import ServicesStats
from .views import api_plasmid_getfasta_byid


from django.contrib.auth.decorators import login_required
from django.urls import path

urlpatterns = [
    path('plasmids/', login_required(views.plasmids, redirect_field_name='next'), name='plasmids'),
    path('plasmids/validations/', login_required(PlasmidValidations, redirect_field_name='next'), name='plasmid_validations'),
    path('plasmid/<uuid:plasmid_id>/', login_required(views.plasmid, redirect_field_name='next'), name='plasmid'),
    path('qr/i/<str:plasmid_id>/', login_required(views.plasmid_from_qr, redirect_field_name='next'), name='plasmid_from_qr'), # for QR
    path('plasmid/create/', login_required(PlasmidCreate.as_view(), redirect_field_name='next'), name='plasmid_create'),
    path('plasmid/create/wizard', login_required(PlasmidCreateWizard.as_view(), redirect_field_name='next'), name='plasmid_create_wizard'),
    path('plasmid/create/wizard/end', login_required(plasmid_create_wizard_end, redirect_field_name='next'), name='plasmid_create_wizard_end'),
    path('plasmid/edit/<uuid:pk>/', login_required(PlasmidEdit.as_view(), redirect_field_name='next'), name='plasmid_edit'),
    path('plasmid/delete/<uuid:pk>/', login_required(PlasmidDelete.as_view(), redirect_field_name='next'), name='plasmid_delete'),
    path('plasmid/label/<uuid:plasmid_id>/', login_required(views.plasmid_label, redirect_field_name='next'), name='plasmid_label'),
    path('plasmid/download/<uuid:plasmid_id>', login_required(views.plasmid_download, redirect_field_name='next'), name='plasmid_download'),
    path('plasmid/view_edit/<uuid:plasmid_id>', login_required(views.plasmid_view_edit, redirect_field_name='next'), name='plasmid_view_edit'),
    path('plasmid/sanger/<uuid:plasmid_id>', login_required(views.plasmid_sanger, redirect_field_name='next'), name='plasmid_sanger'),
    path('plasmid/digest/<uuid:plasmid_id>', login_required(views.plasmid_digest, redirect_field_name='next'), name='plasmid_digest'),
    path('plasmid/pcr/<uuid:plasmid_id>', login_required(views.plasmid_pcr, redirect_field_name='next'), name='plasmid_pcr'),
    path('plasmid/validation/edit/<uuid:plasmid_id>/', login_required(PlasmidValidationEdit, redirect_field_name='next'), name='plasmid_validation_edit'),

    path('glycerolstocks/', login_required(views.glycerolstocks, redirect_field_name='next'), name='glycerolstocks'),
    path('glycerolstock/<uuid:glycerolstock_id>/', login_required(views.glycerolstock, redirect_field_name='next'), name='glycerolstock'),
    path('qr/g/', login_required(views.glycerolstock_qr, redirect_field_name='next'), name='glycerolstock_qr'),
    path('qr/g/<str:glycerolstock_id>/', login_required(views.glycerolstock_from_qr, redirect_field_name='next'), name='glycerolstock_from_qr'), # for QR
    path('glycerolstock/create/', login_required(GstockCreate.as_view(), redirect_field_name='next'), name='glycerolstock_create'),
    path('glycerolstock/create/<uuid:pid>', login_required(GstockCreatePlasmidDefined.as_view(), redirect_field_name='next'), name='glycerolstock_create_plasmid_defined'),
    path('glycerolstock/edit/<uuid:pk>/', login_required(GstockEdit.as_view(), redirect_field_name='next'), name='glycerolstock_edit'),
    path('glycerolstock/delete/<uuid:pk>/', login_required(GstockDelete.as_view(), redirect_field_name='next'), name='glycerolstock_delete'),
    path('glycerolstock/label/<uuid:glycerolstock_id>/', login_required(views.glycerolstock_label, redirect_field_name='next'), name='glycerolstock_label'),
    path('glycerolstock/boxes/', login_required(views.glycerolstock_boxes, redirect_field_name='next'), name='glycerolstock_boxes'),

    path('restrictionenzymes/', login_required(views.restrictionenzymes, redirect_field_name='next'), name='restrictionenzymes'),
    path('restrictionenzyme/<uuid:restrictionenzyme_id>/', login_required(views.restrictionenzyme, redirect_field_name='next'), name='restrictionenzyme'),
    path('restrictionenzyme/create/', login_required(RestrictionenzymeCreate.as_view(), redirect_field_name='next'), name='restrictionenzyme_create'),
    path('restrictionenzyme/edit/<uuid:pk>/', login_required(RestrictionenzymeEdit.as_view(), redirect_field_name='next'), name='restrictionenzyme_edit'),
    path('restrictionenzyme/delete/<uuid:pk>/', login_required(RestrictionenzymeDelete.as_view(), redirect_field_name='next'), name='restrictionenzyme_delete'),

    path('primers/', login_required(views.primers, redirect_field_name='next'), name='primers'),
    path('primer/<uuid:primer_id>/', login_required(views.primer, redirect_field_name='next'), name='primer'),
    path('primer/create/', login_required(PrimerCreate.as_view(), redirect_field_name='next'), name='primer_create'),
    path('primer/edit/<uuid:pk>/', login_required(PrimerEdit.as_view(), redirect_field_name='next'), name='primer_edit'),
    path('primer/delete/<uuid:pk>/', login_required(PrimerDelete.as_view(), redirect_field_name='next'), name='primer_delete'),

    path('services/blast/', login_required(ServicesBlast, redirect_field_name='next'), name='services-blast'),
    path('services/gtr/', login_required(ServicesGtr, redirect_field_name='next'), name='services-gtr'),
    path('services/l0d/', login_required(ServicesL0d, redirect_field_name='next'), name='services-l0d'),
    path('services/stats/', login_required(ServicesStats, redirect_field_name='next'), name='services-stats'),

    path('api/plasmid/get_fasta/by_name/<str:name>/', api_plasmid_getfasta_byid, name='api-plasmid-getfasta_byid'),
    path('api/plasmids/', login_required(views.api_plasmids, redirect_field_name='next'), name='api-plasmids'),
    # path('api/plasmids/', views.api_plasmids, name='api-plasmids'),
    path('api/parts/<str:enzyme_name>/<str:assembly_standard>/', login_required(views.api_parts, redirect_field_name='next'), name='api-parts'),
    # path('api/parts/<str:enzyme_name>/<str:assembly_standard>/', views.api_parts, name='api-parts'),
    path('api/glycerolstocks/', login_required(views.api_glycerolstocks, redirect_field_name='next'), name='api-glycerolstocks'),
    # path('api/glycerolstocks/', views.api_glycerolstocks, name='api-glycerolstocks'),

    # path('__debug__/', include(debug_toolbar.urls)), #delete
]
