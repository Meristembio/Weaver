# Generated by Django 4.0.4 on 2022-07-19 16:40

from django.db import migrations
import shortuuidfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('inventory', '0016_alter_glycerolstock_qr_id_alter_plasmid_qr_id_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='glycerolstock',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='2HXMQomomVzo98RGhHhUJK', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='plasmid',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='PCdnnsKsMpzWgunA4WFheN', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='primer',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='mQooALoWTgamevbvidPfJH', editable=False, max_length=22),
        ),
    ]
