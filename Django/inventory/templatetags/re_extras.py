from django import template

register = template.Library()


def complement_base(b):
    if b == 'a':
        return 't'
    if b == 'A':
        return 'T'
    if b == 't':
        return 'a'
    if b == 'T':
        return 'A'
    if b == 'c':
        return 'g'
    if b == 'C':
        return 'G'
    if b == 'g':
        return 'c'
    if b == 'G':
        return 'C'
    if b == 'n':
        return 'n'
    if b == 'N':
        return 'N'


def complement(seq):
    compl_seq = ''
    for b in seq:
        compl_seq = compl_seq + complement_base(b)
    return compl_seq


@register.simple_tag
def re_render(sequence, fcut, rcut):
    extra_ns = 0
    if fcut is not None and rcut is not None:
        if len(sequence) < fcut:
            extra_ns = fcut - len(sequence)
        if len(sequence) < int(rcut):
            if rcut - len(sequence) > extra_ns:
                extra_ns = rcut - len(sequence)
    if extra_ns:
        extra_ns = extra_ns + 2
        for i in range(extra_ns):
            sequence = sequence + 'N'
    complement_seq = complement(sequence)
    cut1 = min(fcut, rcut)
    cut2 = max(fcut, rcut)
    re_center_class = "re_left"
    re_center_class_comp = "re_right"
    if fcut < rcut:
        re_center_class = "re_right"
        re_center_class_comp = "re_left"
    sequence = "<span class=\"re_left\">" + sequence[
                                            0:cut1] + "</span><span class=\"re_center " + re_center_class + "\">" + sequence[
                                                                                                                    cut1:cut2] + "</span>" + "</span><span class=\"re_right\">" + sequence[
                                                                                                                                                                                  cut2:] + "</span>"
    complement_seq = "<span class=\"re_left\">" + complement_seq[
                                                  0:cut1] + "</span><span class=\"re_center " + re_center_class_comp + "\">" + complement_seq[
                                                                                                                               cut1:cut2] + "</span>" + "</span><span class=\"re_right\">" + complement_seq[
                                                                                                                                                                                             cut2:] + "</span>"
    return "<pre class=\"re\">" + sequence + "<br/>" + complement_seq + "</pre>"
