from django.shortcuts import render
from organization.models import access_policies_options


def index(request):
    if request.user.is_authenticated:
        context = {'apo': access_policies_options}
        return render(request, 'organization/projects.html', context)
    else:
        return render(request, 'common/index.html')
