from django.shortcuts import render
from .models import Recipe
from .models import Component
from .models import Reactive
from .models import RECIPE_CATEGORIES
from django.core.exceptions import ObjectDoesNotExist
from django.http import Http404
from .models import CONCENTRATION_UNITS
from .forms import RecipeForm
from django.views.generic.edit import UpdateView
from django.views.generic.edit import CreateView
from django.urls import reverse
import datetime


def recipe(request, recipe_id):
    try:
        recipe_to_detail = Recipe.objects.get(id=recipe_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'CONCENTRATION_UNITS': CONCENTRATION_UNITS,
        'recipe': recipe_to_detail
    }

    if request.method == 'POST':
        recipe_form = RecipeForm(request.POST)
        if recipe_form.is_valid():
            context['result'] = []
            context['warnings'] = []

            context['quantity'] = recipe_form.cleaned_data['quantity']
            context['quantity_unit'] = recipe_form.cleaned_data['unit']
            context['concentration'] = recipe_form.cleaned_data['concentration']

            try:
                quantity = recipe_form.cleaned_data['quantity'] * context['concentration']
            except:
                context['error'] = "Quantity or concentration not valid"
            # siempre en ml
            if recipe_form.cleaned_data['unit'] == 'ml':
                quantity = quantity / 1000
            for component in recipe_to_detail.components.all():
                # checks
                # if not component.concentration_unit or not component.concentration:
                #     context['error'] = "Component " + component.name + ": Quantity or concentration not set"
                #     break
                # if not component.reactive.concentration_unit or not component.reactive.concentration:
                #     context['error'] = "Reactive " + component.reactive.name + ": Quantity or concentration not set"
                #     break
                # siempre mMol
                if component.concentration_unit == 'mol':
                    component.concentration_unit = 'mmol'
                    component.concentration = component.concentration * 1000
                if component.reactive.concentration_unit == 'mol':
                    component.reactive.concentration_unit = 'mmol'
                    component.reactive.concentration = component.reactive.concentration * 1000
                # calcular componentes
                if component.concentration_unit == 'grlt':
                    if component.reactive.state != 0:
                        # if reactive isn't solid state
                        context['warnings'].append(
                            "Reactive " + component.reactive.name +
                            " is not in solid state. Calculation of components in \"" +
                            CONCENTRATION_UNITS[2][1] + "\" unit assumes solid state reactive.")
                    component_mass_or_volume = component.concentration * quantity
                    component_unit = 'gr'
                    if component_mass_or_volume < 1:
                        component_mass_or_volume = component_mass_or_volume * 1000
                        component_unit = 'mg'
                elif component.concentration_unit == 'vvp':
                    if component.reactive.state != 1:
                        # if reactive isn't liquid state
                        context['warnings'].append(
                            "Reactive " + component.reactive.name +
                            " is not in solid state. Calculation of components in \"" +
                            CONCENTRATION_UNITS[3][1] + "\" unit assumes liquid state reactive.")
                    component_mass_or_volume = component.concentration * quantity / component.reactive.concentration
                    component_unit = 'lt'
                    if component_mass_or_volume < 1:
                        component_mass_or_volume = component_mass_or_volume * 1000
                        component_unit = 'ml'
                elif component.concentration_unit == 'wvp':
                    if component.reactive.state != 0:
                        # if reactive isn't solid state
                        context['warnings'].append(
                            "Reactive " + component.reactive.name +
                            " is not in solid state. Calculation of components in \"" +
                            CONCENTRATION_UNITS[4][1] + "\" unit assumes solid state reactive.")
                    component_mass_or_volume = component.concentration * 10 * quantity
                    component_unit = 'gr'
                    if component_mass_or_volume < 1:
                        component_mass_or_volume = component_mass_or_volume * 1000
                        component_unit = 'mg'
                elif component.concentration_unit == 'mmol':
                    # solid reactive
                    if component.reactive.state == 0:
                        # mass (gr) = concentration (mM) * volume (Lt) * molecular mass (g/mol) / 1000
                        component_mass_or_volume = component.concentration * quantity * component.reactive.mm / 1000
                        component_unit = 'gr'
                        if component_mass_or_volume < 1:
                            component_mass_or_volume = component_mass_or_volume * 1000
                            component_unit = 'mg'
                    else:
                        # liquid
                        if component.reactive.concentration_unit != 'mmol':
                            context[
                                'error'] = 'I can only calculate liquid dilutions from molarity to molarity (' + component.reactive.name + ')'
                            break
                        if component.reactive.concentration < component.concentration:
                            context[
                                'error'] = 'Your stock concentration (' + str(
                                component.reactive.concentration) + ') of ' + component.reactive.name + ' is lower than desired concentration (' + str(
                                component.concentration) + ')'
                            break
                        component_mass_or_volume = component.concentration * quantity / component.reactive.concentration
                        component_unit = 'lt'
                        if component_mass_or_volume < 1:
                            component_mass_or_volume = component_mass_or_volume * 1000
                            component_unit = 'ml'
                        if component_mass_or_volume < 1:
                            component_mass_or_volume = component_mass_or_volume * 1000
                            component_unit = 'ul'
                else:
                    context[
                        'error'] = 'Unit (' + component.concentration_unit + ') of component (' + component.reactive.name + ') is not valid (yet)'
                    break

                context['result'].append(
                    (component.reactive.name, component_mass_or_volume, component_unit))

    context['recipe_form'] = RecipeForm()

    return render(request, 'protocols/recipe.html', context)


def recipes(request):
    table_filters = [
        ['all', 'All', [
            ('All', 'all', 'primary'),
        ]],
        ['category', 'Category', []],
    ]
    for recipe_cat in RECIPE_CATEGORIES:
        table_filters[1][2].append((recipe_cat[1], recipe_cat[0], 'warning'))
    context = {
        'recipes': Recipe.objects.all(),
        'table_filters': table_filters
    }
    return render(request, 'protocols/recipes.html', context)


class RecipeEdit(UpdateView):
    model = Recipe
    fields = '__all__'
    template_name_suffix = '_update_form'

    def get_success_url(self, **kwargs):
        return reverse('recipe', args=(self.object.id,)) + '?form_result_recipe_edit_success=true'


class RecipeCreate(CreateView):
    model = Recipe
    fields = '__all__'
    template_name_suffix = '_create_form'

    def get_success_url(self, **kwargs):
        return reverse('recipes') + '?form_result_recipe_create_success=true'


def recipe_label(request, recipe_id):
    try:
        recipe_to_label = Recipe.objects.get(id=recipe_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'recipe': recipe_to_label,
        'date': datetime.datetime.now().date(),
        'CONCENTRATION_UNITS': CONCENTRATION_UNITS,
    }
    return render(request, 'protocols/recipe_label.html', context)


class ComponentEdit(UpdateView):
    model = Component
    fields = '__all__'
    template_name_suffix = '_update_form'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['CONCENTRATION_UNITS'] = CONCENTRATION_UNITS
        context['recipe_id'] = self.request.GET['recipe_id']
        return context

    def get_success_url(self, **kwargs):
        return reverse('recipe', args=(self.request.GET['recipe_id'],)) + '?form_result_component_edit_success=true'


class ComponentCreate(CreateView):
    model = Component
    fields = '__all__'
    template_name_suffix = '_create_form'

    def get_success_url(self, **kwargs):
        return reverse('recipe_create') + '?form_result_component_create_success=true'


class ReactiveCreate(CreateView):
    model = Reactive
    fields = '__all__'
    template_name_suffix = '_create_form'

    def get_success_url(self, **kwargs):
        return reverse('component_create') + '?form_result_reactive_create_success=true'


class ReactiveEdit(UpdateView):
    model = Reactive
    fields = '__all__'
    template_name_suffix = '_update_form'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['recipe_id'] = self.request.GET['recipe_id']
        return context

    def get_success_url(self, **kwargs):
        return reverse('recipe',
                       args=(self.request.GET['recipe_id'],)) + '?form_result_reactive_edit_success=true'


def index(request):
    context = {
    }
    return render(request, 'protocols/index.html', context)
