# Generated by Django 4.0.4 on 2022-07-11 17:00

from django.db import migrations
import shortuuidfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('inventory', '0006_alter_glycerolstock_qr_id_alter_plasmid_qr_id_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='glycerolstock',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='2jxaJDLjEMcSYWBbjA7CpB', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='plasmid',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='o3KDWbm4jneimHKJHiuEyb', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='primer',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='ejhwv2BudtYXT9jP53YbtN', editable=False, max_length=22),
        ),
    ]
