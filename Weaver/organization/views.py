from django.shortcuts import render
from .models import Project
from .models import Membership
from .models import access_policies_options
from django.urls import reverse
from django.views.generic.edit import UpdateView
from django.views.generic.edit import CreateView
from django.views.generic.edit import DeleteView
from django.http import Http404
from django.core.exceptions import ObjectDoesNotExist


def projects(request):
    context = {'apo': access_policies_options}
    return render(request, 'organization/projects.html', context)


def project(request, project_id):
    try:
        project_to_detail = Project.objects.get(id=project_id)
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'apo': access_policies_options,
        'project': project_to_detail
    }
    if 'edit_error' in request.GET:
        context['edit_error'] = request.GET.get('edit_error')
    return render(request, 'organization/project.html', context)


def project_has_at_least_one_admin(project, members):
    for member in members:
        for ms in project.membership_set.all():
            if member == ms.member:
                if ms.access_policies == 'a':
                    return True
    return False


def on_membership_are_other_admins(membership):
    for m in membership.project.members.all():
        if m != membership.member and on_project_member_can(membership.project, m, 'a'):
            return True
    return False


def on_membership_member_can(membership, member, perm):
    for m in membership.project.members.all():
        if m == member and on_project_member_can(membership.project, member, perm):
            return True
    return False


def on_project_member_can(project, member, perm):
    try:
        membership = Membership.objects.get(member=member, project=project)
        if membership.access_policies == perm:
            return True
        return False
    except Membership.DoesNotExist:
        return False


class ProjectEdit(UpdateView):
    model = Project
    template_name_suffix = '_update_form'
    fields = ['name', 'members', 'public']

    def get(self, request, *args, **kwargs):
        if on_project_member_can(self.get_object(), self.request.user, 'a'):
            return super(ProjectEdit, self).get(request, *args, **kwargs)
        else:
            if not self.request.GET._mutable:
                self.request.GET._mutable = True
            self.request.GET['edit_error'] = 'Current user does not have \'admin\' permission'
            return project(self.request, self.object.id)

    def form_valid(self, form):
        if project_has_at_least_one_admin(form.instance, form.cleaned_data['members']):
            return super(ProjectEdit, self).form_valid(form)
        else:
            form.add_error(None, "At least one user with \'admin\' access is required")
            return super(ProjectEdit, self).form_invalid(form)

    def get_success_url(self, **kwargs):
        return reverse('project',
                       args=(self.object.id,)) + '?form_result_project_edit_success=true'


class ProjectCreate(CreateView):
    model = Project
    template_name_suffix = '_create_form'
    fields = ['name', 'members', 'public']

    def form_valid(self, form):
        response = super(ProjectCreate, self).form_valid(form)
        try:
            membership = Membership.objects.get(member=self.request.user, project=form.instance)
            membership.access_policies = 'a'
            membership.save()
        except Membership.DoesNotExist:
            Membership.objects.create(member=self.request.user, project=form.instance, access_policies='w')

        return response

    def get_initial(self):
        return {
            'members': [self.request.user],
        }

    def get_success_url(self, **kwargs):
        return reverse('project',
                       args=(self.object.id,)) + '?form_result_project_create_success=true'


class ProjectDelete(DeleteView):
    model = Project

    def get(self, request, *args, **kwargs):
        if on_project_member_can(self.get_object(), self.request.user, 'a'):
            return super(ProjectDelete, self).get(request, *args, **kwargs)
        else:
            if not self.request.GET._mutable:
                self.request.GET._mutable = True
            self.request.GET['edit_error'] = 'Current user does not have \'admin\' permission'
            return project(self.request, self.object.id)

    def get_success_url(self, **kwargs):
        return reverse('projects') + '?form_result_object_deleted=true'


class ProjectMemberEdit(UpdateView):
    model = Membership
    template_name_suffix = '_update_form'
    fields = ['access_policies']

    def get(self, request, *args, **kwargs):
        if not self.request.GET._mutable:
            self.request.GET._mutable = True
        if on_membership_member_can(self.get_object(), self.request.user, 'a'):
            if on_membership_are_other_admins(self.get_object()):
                return super(ProjectMemberEdit, self).get(request, *args, **kwargs)
            else:
                self.request.GET['edit_error'] = 'This membership can not be modified as user is the last with \'admin\' permission'
                return project(self.request, self.get_object().project.id)
        else:
            self.request.GET['edit_error'] = 'Current user does not have \'admin\' permission'
            return project(self.request, self.get_object().project.id)

    def get_success_url(self, **kwargs):
        return reverse('project',
                       args=(self.object.project.id,)) + '?form_result_membership_edit_success=true'


def profile(request):
    context = {'apo': access_policies_options}
    return render(request, 'registration/profile.html', context)
