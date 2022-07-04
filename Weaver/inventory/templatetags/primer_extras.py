from django import template

register = template.Library()


@register.filter(name='gc_content')
def gc_content(value):
    gc_count = 0
    for b in value.lower():
        if b == 'c' or b == 'g':
            gc_count = gc_count + 1
    return gc_count / len(value) * 100


@register.filter(name='tm_value')
def tm_value(value):
    tm = 0
    a = 0
    c = 0
    t = 0
    g = 0
    for b in value.lower():
        if b == 'a':
            a = a + 1
        if b == 'c':
            c = c + 1
        if b == 't':
            t = t + 1
        if b == 'g':
            g = g + 1
    if len(value) < 14:
        tm = (a + t) * 2 + (c + g) * 4
    else:
        tm = 64.9 + 41 * (g + c - 16.4) / (a + t + c + g)
    return tm
