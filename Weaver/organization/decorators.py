from django.shortcuts import render
from .views import get_current_project
from .views import has_current_project
from .views import on_project_member_can_any
from .views import on_project_member_can_write_or_admin
from organization.models import access_policies_options
from inventory.models import GlycerolStock
from inventory.models import Plasmid


def require_current_project_set(function):
    def wrap(request, *args, **kwargs):
        if has_current_project(request):
            return function(request, *args, **kwargs)
        else:
            context = {
                'current_project': False,
                'apo': access_policies_options
            }
            return render(request, 'common/set_project.html', context)
    return wrap


def require_member_can_write_or_admin_current_project(function):
    def wrap(request, *args, **kwargs):
        if on_project_member_can_write_or_admin(get_current_project(request), request.user):
            return function(request, *args, **kwargs)
        else:
            return render(request, 'common/no_permission_to_edit.html')
    return wrap


def require_member_can_write_or_admin_project_of_plasmid(function):
    def wrap(request, *args, **kwargs):
        key = 'plasmid_id'
        if key not in kwargs:
            key = 'pk'
        if key in kwargs:
            plasmid = Plasmid.objects.get(pk=kwargs[key])
            if plasmid:
                if on_project_member_can_write_or_admin(plasmid.project, request.user):
                    return function(request, *args, **kwargs)
                else:
                    return render(request, 'common/no_object_found.html')
            else:
                return render(request, 'common/no_permission_to_edit.html')
        else:
            return render(request, 'common/no_element_id_defined.html')
    return wrap


def require_member_can_write_or_admin_project_of_gs(function):
    def wrap(request, *args, **kwargs):
        key = 'glycerolstock_id'
        if key not in kwargs:
            key = 'pk'
        if key in kwargs:
            gs = GlycerolStock.objects.get(pk=kwargs[key])
            if gs:
                if on_project_member_can_write_or_admin(gs.project, request.user):
                    return function(request, *args, **kwargs)
                else:
                    return render(request, 'common/no_object_found.html')
            else:
                return render(request, 'common/no_permission_to_edit.html')
        else:
            return render(request, 'common/no_element_id_defined.html')
    return wrap


def require_member_can_read_project_of_plasmid(function):
    def wrap(request, *args, **kwargs):
        pk = None
        key = 'plasmid_id'
        if key not in kwargs:
            key = 'pk'
        if key in kwargs:
            pk = kwargs[key]
        else:
            if args[0]:
                pk = args[0]
        if pk:
            plasmid = Plasmid.objects.get(pk=pk)
            if plasmid:
                if on_project_member_can_any(plasmid.project, request.user):
                    return function(request, *args, **kwargs)
                else:
                    return render(request, 'common/no_object_found.html')
            else:
                return render(request, 'common/no_permission_to_edit.html')
        else:
            return render(request, 'common/no_element_id_defined.html')
    return wrap



def require_member_can_read_project_of_gs(function):
    def wrap(request, *args, **kwargs):
        key = 'glycerolstock_id'
        if key not in kwargs:
            key = 'pk'
        if key in kwargs:
            gs = GlycerolStock.objects.get(pk=kwargs[key])
            if gs:
                if on_project_member_can_any(gs.project, request.user):
                    return function(request, *args, **kwargs)
                else:
                    return render(request, 'common/no_object_found.html')
            else:
                return render(request, 'common/no_permission_to_edit.html')
        else:
            return render(request, 'common/no_element_id_defined.html')
    return wrap
