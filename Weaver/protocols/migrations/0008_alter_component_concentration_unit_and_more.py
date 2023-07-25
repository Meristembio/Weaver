# Generated by Django 4.0.4 on 2022-11-28 16:12

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('protocols', '0007_alter_component_owner_alter_reactive_owner_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='component',
            name='concentration_unit',
            field=models.CharField(choices=[('mol', 'M'), ('mmol', 'mM'), ('grlt', 'mg / ml (= g/l)'), ('vvp', 'Volume / Volume %'), ('wvp', 'Weight / Volume %')], max_length=4),
        ),
        migrations.AlterField(
            model_name='reactive',
            name='concentration_unit',
            field=models.CharField(blank=True, choices=[('mol', 'M'), ('mmol', 'mM'), ('grlt', 'mg / ml (= g/l)'), ('vvp', 'Volume / Volume %'), ('wvp', 'Weight / Volume %')], max_length=4, null=True),
        ),
    ]
