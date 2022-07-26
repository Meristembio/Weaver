from . import views
from django.urls import path
from .views import ProjectCreate
from .views import ProjectEdit
from .views import ProjectDelete
from .views import ProjectMemberEdit
from django.contrib.auth.decorators import login_required


urlpatterns = [
    path('projects/', login_required(views.projects, redirect_field_name='next'), name='projects'),
    path('project/<int:project_id>/', login_required(views.project, redirect_field_name='next'), name='project'),
    path('project/edit/<int:pk>/', login_required(ProjectEdit.as_view(), redirect_field_name='next'), name='project_edit'),
    path('project/delete/<int:pk>/', login_required(ProjectDelete.as_view(), redirect_field_name='next'), name='project_delete'),
    path('project/set_current/<int:pk>/', login_required(views.project_set_current, redirect_field_name='next'), name='project_set_current'),
    path('project/membership/edit/<int:pk>/', login_required(ProjectMemberEdit.as_view(), redirect_field_name='next'), name='membership_edit'),
    path('project/create', login_required(ProjectCreate.as_view(), redirect_field_name='next'), name='project_create'),

    path('accounts/profile/', login_required(views.profile), name='profile'),
]