from django import forms

RECIPE_UNITS = (
    ('ml', 'ml'),
    ('lt', 'lt'),
)


class RecipeForm(forms.Form):
    quantity = forms.FloatField(label='Quantity')
    unit = forms.ChoiceField(label='Unit', choices=RECIPE_UNITS)
    concentration = forms.FloatField(label='Concentration', required=False, initial=1)
