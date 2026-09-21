from django.db import migrations
from django.db import models


class Migration(migrations.Migration):

    dependencies = [
        ("inventory", "0013_seed_ytk_enzyme_metadata"),
    ]

    operations = [
        migrations.AddField(
            model_name="plasmid",
            name="assembly_metadata",
            field=models.JSONField(blank=True, default=dict),
        ),
    ]
