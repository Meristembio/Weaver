from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ("inventory", "0016_sangerreadfile_primer"),
    ]

    operations = [
        migrations.RemoveField(
            model_name="primer",
            name="project",
        ),
    ]
