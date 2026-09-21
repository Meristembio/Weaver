from django import template
from organization.models import Project
from markdown import markdown as render_markdown

register = template.Library()


@register.filter
def get_project_name(project_id):
    try:
        return Project.objects.get(pk=project_id).name
    except (Project.DoesNotExist, ValueError, TypeError):
        return ""


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


@register.filter
def markdown(value):
    return render_markdown(value)
