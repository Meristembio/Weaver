from django.db import models
from django.contrib.auth.models import User


class Project(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, help_text="Project name")
    description = models.CharField(max_length=500, help_text="Project description", null=True, blank=True)
    members = models.ManyToManyField(User, through='Membership', related_name='project_members',
                                     help_text="User Ctrl to select more than one member")
    public = models.BooleanField(help_text='Is it publicly viewable?')

    def __str__(self):
        return self.name


access_policies_options = (
    ('r', 'read'),
    ('w', 'write'),
    ('a', 'admin')
)


class Membership(models.Model):
    member = models.ForeignKey(User, on_delete=models.CASCADE)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    access_policies = models.CharField(max_length=1, choices=access_policies_options, default='r')

    def __str__(self):
        access_policies_text = ""
        for apo in access_policies_options:
            if apo[0] == self.access_policies:
                access_policies_text = apo[1]
        return str(self.member) + " - " + str(self.project) + ": " + access_policies_text
