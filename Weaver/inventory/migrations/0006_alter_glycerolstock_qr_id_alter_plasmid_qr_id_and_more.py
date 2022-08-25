# Generated by Django 4.0.4 on 2022-07-11 16:59

from django.db import migrations
import shortuuidfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('inventory', '0005_alter_glycerolstock_qr_id_alter_plasmid_qr_id_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='glycerolstock',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='79uNAREGoK2eNxmKeSas9z', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='plasmid',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='76WfxnxPgCmZDh9KkX4uok', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='primer',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='2XaqvUeQYASEdpad3sWshP', editable=False, max_length=22),
        ),
    ]