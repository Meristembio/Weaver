# Generated by Django 4.0.4 on 2022-11-29 15:43

from django.db import migrations, models
import shortuuidfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('inventory', '0030_alter_glycerolstock_qr_id_alter_plasmid_qr_id_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='plasmid',
            name='idx',
            field=models.IntegerField(default=None, editable=False),
        ),
        migrations.AlterField(
            model_name='glycerolstock',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='JVBLk8wQ847pijm6P8Df55', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='plasmid',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='EoWZ7PeExxPQgxihhmzJk3', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='primer',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='atoK6KCXP9qPPtd9xyUp2p', editable=False, max_length=22),
        ),
    ]
