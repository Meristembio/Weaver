from django import template
from organization.models import Project

register = template.Library()


@register.filter
def get_project_name(project_id):
    return Project.objects.get(pk=project_id).name


@register.filter
def get_element_by_key(h, key):
    return h[key]


@register.filter
def get_element_id_by_key(h, key):
    return h[key].id


@register.filter
def startswith(text, starts):
    if text.lower().startswith(starts):
        return True
    return False

@register.filter
def divide(value, arg):
    try:
        return int(value) / int(arg)
    except (ValueError, ZeroDivisionError):
        return None