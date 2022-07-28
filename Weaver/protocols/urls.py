from django.urls import path

from . import views
from .views import RecipeEdit
from .views import RecipeCreate
from .views import RecipeDelete
from .views import ComponentEdit
from .views import ComponentCreateReturnToRecipe
from .views import ComponentCreate
from .views import ComponentDelete
from .views import ReactiveCreateReturnToRecipe
from .views import ReactiveCreate
from .views import ReactiveEdit
from .views import ReactiveDelete
from django.contrib.auth.decorators import login_required

urlpatterns = [
    path('', login_required(views.index, redirect_field_name='next'), name='index'),
    path('recipes', login_required(views.recipes, redirect_field_name='next'), name='recipes'),
    path('recipe/<uuid:recipe_id>', login_required(views.recipe, redirect_field_name='next'), name='recipe'),
    path('recipe/create/', login_required(RecipeCreate.as_view(), redirect_field_name='next'), name='recipe_create'),
    path('recipe/edit/<uuid:pk>/', login_required(RecipeEdit.as_view(), redirect_field_name='next'), name='recipe_edit'),
    path('recipe/delete/<uuid:pk>/', login_required(RecipeDelete.as_view(), redirect_field_name='next'), name='recipe_delete'),
    path('recipe/label/<uuid:recipe_id>/', login_required(views.recipe_label, redirect_field_name='next'), name='recipe_label'),

    path('component/create/', login_required(ComponentCreate.as_view(), redirect_field_name='next'), name='component_create'),
    path('component/create/return_to_recipe/<uuid:recipe_id>/', login_required(ComponentCreateReturnToRecipe.as_view(), redirect_field_name='next'), name='component_create_return_to_recipe'),
    path('component/edit/<uuid:pk>/return_to_recipe/<uuid:recipe_id>/', login_required(ComponentEdit.as_view(), redirect_field_name='next'), name='component_edit'),
    path('component/delete/<uuid:pk>/return_to_recipe/<uuid:recipe_id>/', login_required(ComponentDelete.as_view(), redirect_field_name='next'), name='component_delete'),

    path('reactive/create/', login_required(ReactiveCreate.as_view(), redirect_field_name='next'), name='reactive_create'),
    path('reactive/create/return_to_recipe/<uuid:recipe_id>/', login_required(ReactiveCreateReturnToRecipe.as_view(), redirect_field_name='next'), name='reactive_create_return_to_recipe'),
    path('reactive/edit/<uuid:pk>/return_to_recipe/<uuid:recipe_id>/', login_required(ReactiveEdit.as_view(), redirect_field_name='next'), name='reactive_edit'),
    path('reactive/delete/<uuid:pk>/return_to_recipe/<uuid:recipe_id>/', login_required(ReactiveDelete.as_view(), redirect_field_name='next'), name='reactive_delete'),
]
