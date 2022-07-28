from django import forms
from .models import Recipe
from .models import Component
from .models import Reactive
from organization.views import get_projects_where_member_can_write_or_admin

RECIPE_UNITS = (
    ('ml', 'ml'),
    ('lt', 'lt'),
)


class ComponentCreateForm(forms.ModelForm):
    class Meta:
        model = Component
        fields = '__all__'
        exclude = ['owner']


def get_member_components(member):
    return Component.objects.filter(owner=member)


def get_member_reactives(member):
    return Reactive.objects.filter(owner=member)


class RecipeCreateEditForm(forms.ModelForm):
    def __init__(self, *args, **kwargs):
        member = kwargs.pop('user')
        super(RecipeCreateEditForm, self).__init__(*args, **kwargs)
        self.fields['components'].queryset = get_member_components(member)
        self.fields['shared_to_project'].queryset = get_projects_where_member_can_write_or_admin(member)

    class Meta:
        model = Recipe
        fields = '__all__'
        exclude = ['owner']


class ComponentCreateEditForm(forms.ModelForm):
    def __init__(self, *args, **kwargs):
        member = kwargs.pop('user')
        super(ComponentCreateEditForm, self).__init__(*args, **kwargs)
        self.fields['reactive'].queryset = get_member_reactives(member)

    class Meta:
        model = Component
        fields = '__all__'
        exclude = ['owner']


class RecipeForm(forms.Form):
    quantity = forms.FloatField(label='Quantity')
    unit = forms.ChoiceField(label='Unit', choices=RECIPE_UNITS)
    concentration = forms.FloatField(label='Concentration', required=False, initial=1)
