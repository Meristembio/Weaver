# Generated by Django 4.0.4 on 2022-07-11 19:52

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('organization', '0005_rename_group_membership_project_and_more'),
    ]

    operations = [
        migrations.RenameField(
            model_name='membership',
            old_name='person',
            new_name='member',
        ),
    ]
