from inventory.models import Primer


def visible_primers_for_user(user):
    if not user or not user.is_authenticated:
        return Primer.objects.none()
    return Primer.objects.all()
