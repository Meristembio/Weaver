# Generated by Django 4.0.4 on 2023-04-20 11:36

from django.db import migrations, models
import inventory.validators
import shortuuidfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('inventory', '0034_alter_glycerolstock_qr_id_alter_plasmid_idx_and_more'),
    ]

    operations = [
        migrations.RenameField(
            model_name='plasmid',
            old_name='check_date',
            new_name='colonypcr_date',
        ),
        migrations.RenameField(
            model_name='plasmid',
            old_name='check_observations',
            new_name='colonypcr_observations',
        ),
        migrations.RenameField(
            model_name='plasmid',
            old_name='check_state',
            new_name='colonypcr_state',
        ),
        migrations.RemoveField(
            model_name='plasmid',
            name='check_method',
        ),
        migrations.RemoveField(
            model_name='plasmid',
            name='digestion_check_enzymes',
        ),
        migrations.AddField(
            model_name='plasmid',
            name='clustal',
            field=models.FileField(blank=True, null=True, upload_to='uploads/sequencing_clustal', validators=[inventory.validators.clustal_validate]),
        ),
        migrations.AddField(
            model_name='plasmid',
            name='digestion_date',
            field=models.DateField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='plasmid',
            name='digestion_observations',
            field=models.CharField(blank=True, max_length=1000, null=True),
        ),
        migrations.AddField(
            model_name='plasmid',
            name='digestion_state',
            field=models.IntegerField(blank=True, choices=[(0, 'Not required'), (1, 'Pending'), (2, 'Correct')], default=1),
        ),
        migrations.AddField(
            model_name='plasmid',
            name='reference_sequence',
            field=models.BooleanField(blank=True, default=0),
        ),
        migrations.AddField(
            model_name='plasmid',
            name='under_construction',
            field=models.BooleanField(blank=True, default=0),
        ),
        migrations.AlterField(
            model_name='glycerolstock',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='LQBW3rigRxcBoyKEoAw56W', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='plasmid',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='cjSADrRcNpxqfee2QWzsMp', editable=False, max_length=22),
        ),
        migrations.AlterField(
            model_name='plasmid',
            name='sequencing_state',
            field=models.IntegerField(blank=True, choices=[(0, 'Not required'), (1, 'Pending'), (2, 'Correct')], default=0),
        ),
        migrations.AlterField(
            model_name='primer',
            name='qr_id',
            field=shortuuidfield.fields.ShortUUIDField(blank=True, default='Eqf6fs5SVgE32fqWLUaEMH', editable=False, max_length=22),
        ),
    ]