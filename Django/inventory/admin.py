from django.contrib import admin

from .models import GlycerolStock
from .models import Strain
from .models import Plasmid
from .models import PlasmidType
from .models import TableFilter
from .models import Resistance
from .models import Location
from .models import Box
from .models import RestrictionEnzyme
from .models import Primer

admin.site.register(GlycerolStock)
admin.site.register(Strain)
# admin.site.register(Plasmid)
admin.site.register(PlasmidType)
admin.site.register(TableFilter)
admin.site.register(Resistance)
admin.site.register(Location)
admin.site.register(Box)
admin.site.register(RestrictionEnzyme)
admin.site.register(Primer)

admin.site.site_header = 'Weaver. Admin'
admin.site.site_title = 'Weaver'


@admin.register(Plasmid)
class PlasmidAdmin(admin.ModelAdmin):
    def render_change_form(self, request, context, *args, **kwargs):
        context['adminform'].form.fields['backbone'].queryset = Plasmid.objects.filter(type='1')
        return super(PlasmidAdmin, self).render_change_form(request, context, *args, **kwargs)
