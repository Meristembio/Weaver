from django import template

register = template.Library()


@register.filter(name='get_element_by_key')
def get_element_by_key(h, key):
    return h[key]


@register.filter(name='get_element_id_by_key')
def get_element_id_by_key(h, key):
    return h[key].id


@register.filter('startswith')
def startswith(text, starts):
    if text.lower().startswith(starts):
        return True
    return False