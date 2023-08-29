from django.shortcuts import render
from .models import Recipe
from .models import Reactive
from .models import Component
from organization.views import get_projects_where_member_can_any


def require_member_can_view_recipe(function):
    def wrap(request, *args, **kwargs):
        key = 'recipe_id'
        if key not in kwargs:
            key = 'pk'
        if key in kwargs:
            recipe = Recipe.objects.get(pk=kwargs[key])
            if recipe:
                if recipe.owner == request.user:
                    return function(request, *args, **kwargs)
                user_projects = get_projects_where_member_can_any(request.user)
                for project in recipe.shared_to_project.all():
                    if project in user_projects:
                        return function(request, *args, **kwargs)
                else:
                    return render(request, 'common/no_permission_to_read.html')
            else:
                return render(request, 'common/no_object_found.html')
        else:
            return render(request, 'common/no_element_id_defined.html')
    return wrap


def require_member_own_recipe(function):
    def wrap(request, *args, **kwargs):
        key = 'recipe_id'
        if key not in kwargs:
            key = 'pk'
        if key in kwargs:
            recipe = Recipe.objects.get(pk=kwargs[key])
            if recipe:
                if recipe.owner == request.user:
                    return function(request, *args, **kwargs)
                else:
                    return render(request, 'common/no_permission_to_edit.html')
            else:
                return render(request, 'common/no_object_found.html')
        else:
            return render(request, 'common/no_element_id_defined.html')
    return wrap


def require_member_own_reactive(function):
    def wrap(request, *args, **kwargs):
        key = 'reactive_id'
        if key not in kwargs:
            key = 'pk'
        if key in kwargs:
            reactive = Reactive.objects.get(pk=kwargs[key])
            if reactive:
                if reactive.owner == request.user:
                    return function(request, *args, **kwargs)
                else:
                    return render(request, 'common/no_permission_to_edit.html')
            else:
                return render(request, 'common/no_object_found.html')
        else:
            return render(request, 'common/no_element_id_defined.html')
    return wrap


def require_member_own_component(function):
    def wrap(request, *args, **kwargs):
        key = 'component_id'
        if key not in kwargs:
            key = 'pk'
        if key in kwargs:
            component = Component.objects.get(pk=kwargs[key])
            if component:
                if component.owner == request.user:
                    return function(request, *args, **kwargs)
                else:
                    return render(request, 'common/no_permission_to_edit.html')
            else:
                return render(request, 'common/no_object_found.html')
        else:
            return render(request, 'common/no_element_id_defined.html')
    return wrap
