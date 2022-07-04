from django.contrib import admin
from .models import Reactive
from .models import Recipe
from .models import Component

# Register your models here.
admin.site.register(Reactive)
admin.site.register(Recipe)
admin.site.register(Component)
