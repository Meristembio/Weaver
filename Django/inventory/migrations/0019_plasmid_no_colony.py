from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("inventory", "0018_sangerverificationrun_sequencing_datetime"),
    ]

    operations = [
        migrations.AddField(
            model_name="plasmid",
            name="no_colony",
            field=models.BooleanField(
                default=False,
                help_text="Use when the construct is used directly or comes from synthesis.",
            ),
        ),
    ]
