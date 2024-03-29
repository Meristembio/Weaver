# Generated by Django 4.0.4 on 2022-07-27 14:11

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('organization', '0012_project_description'),
        ('protocols', '0003_remove_tablefilter_options'),
    ]

    operations = [
        migrations.AddField(
            model_name='recipe',
            name='shared_to_project',
            field=models.ManyToManyField(blank=True, help_text='Use CTRL for multiple select', to='organization.project'),
        ),
        migrations.AlterField(
            model_name='recipe',
            name='category',
            field=models.ManyToManyField(blank=True, help_text='Use CTRL for multiple select', to='protocols.tablefilter'),
        ),
    ]
