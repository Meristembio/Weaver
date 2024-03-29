# Generated by Django 4.0.4 on 2022-05-24 18:20

from django.db import migrations, models
import django.db.models.deletion
import uuid


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Component',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('concentration', models.FloatField()),
                ('concentration_unit', models.CharField(choices=[('mol', 'M'), ('mmol', 'mM'), ('grlt', 'gr / lt'), ('vvp', 'Volume / Volume %'), ('wvp', 'Weight / Volume %')], max_length=4)),
            ],
            options={
                'ordering': ['reactive'],
            },
        ),
        migrations.CreateModel(
            name='Reactive',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=200)),
                ('state', models.IntegerField(choices=[(0, 'Solid'), (1, 'Liquid')])),
                ('mm', models.FloatField(blank=True, help_text='g/mol', null=True)),
                ('concentration', models.FloatField(blank=True, null=True)),
                ('concentration_unit', models.CharField(blank=True, choices=[('mol', 'M'), ('mmol', 'mM'), ('grlt', 'gr / lt'), ('vvp', 'Volume / Volume %'), ('wvp', 'Weight / Volume %')], max_length=4, null=True)),
            ],
            options={
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='Recipe',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=200)),
                ('ph', models.CharField(blank=True, max_length=10, null=True)),
                ('description', models.CharField(blank=True, max_length=1000, null=True)),
                ('category', models.IntegerField(blank=True, choices=[(0, 'General'), (1, 'DNA'), (2, 'RNA'), (3, 'Protein'), (4, 'Cell'), (5, 'Base editors'), (6, 'Protoplasts')], null=True)),
                ('components', models.ManyToManyField(blank=True, to='protocols.component')),
            ],
            options={
                'ordering': ['name'],
            },
        ),
        migrations.AddField(
            model_name='component',
            name='reactive',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='protocols.reactive'),
        ),
    ]
