from django.core.management.base import BaseCommand, CommandError

from inventory.custom.primer_import import PrimerImportError
from inventory.custom.primer_import import import_primers_from_fasta

class Command(BaseCommand):
    help = "Import primers from a FASTA file into Weaver's Primer database."

    def add_arguments(self, parser):
        parser.add_argument("fasta_path", help="Path to the FASTA file to import.")
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Parse and report changes without writing primers to the database.",
        )
        parser.add_argument(
            "--update-existing",
            action="store_true",
            help="Update existing primers with the same name.",
        )
        parser.add_argument(
            "--name-source",
            choices=("id", "description"),
            default="id",
            help="Use FASTA record ID or full description as the primer name.",
        )
        parser.add_argument(
            "--require-direction",
            action="store_true",
            help="Skip primers whose F/R orientation cannot be inferred from the name.",
        )
        parser.add_argument(
            "--default-direction",
            choices=("f", "r", ""),
            default="",
            help="Direction to use when F/R cannot be inferred from the name.",
        )
    def handle(self, *args, **options):
        try:
            with open(options["fasta_path"], encoding="utf-8-sig") as fasta_handle:
                result = import_primers_from_fasta(
                    fasta_handle,
                    dry_run=options["dry_run"],
                    update_existing=options["update_existing"],
                    name_source=options["name_source"],
                    require_direction=options["require_direction"],
                    default_direction=options["default_direction"],
                )
        except PrimerImportError as error:
            raise CommandError(str(error))

        for message in result["messages"]:
            if message["level"] == "danger":
                self.stderr.write(message["text"])
            else:
                self.stdout.write(message["text"])

        action = "Would create" if options["dry_run"] else "Created"
        update_action = "would update" if options["dry_run"] else "updated"
        self.stdout.write(
            self.style.SUCCESS(
                f"{action} {result['created']}, {update_action} {result['updated']}, skipped {result['skipped']}, errors {result['errors']}."
            )
        )
