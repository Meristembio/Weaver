# Generated by Django 4.0.4 on 2022-07-26 13:08

from django.db import migrations
import shortuuidfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('inventory', '0022_alter_glycerolstock_box_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='glycerolstock',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='Yae7qgJncdP4VPGejhiC43', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='plasmid',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='a2a4JYsd7EJbpH8tBfwyKV', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='primer',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='AYNG92p6cgbUMeep6cT2W4', editable=False, max_length=22),
        ),
    ]
