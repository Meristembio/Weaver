from django.core.management.base import BaseCommand

from inventory.custom.pcr import infer_type_iis_overhang
from inventory.models import Primer


class Command(BaseCommand):
    help = "Identify Type IIS cloning overhangs in primers whose full sequence is stored as the 3' sequence."

    def add_arguments(self, parser):
        parser.add_argument(
            "--apply",
            action="store_true",
            help="Write inferred sequence_5 and trimmed sequence_3 back to matching primers.",
        )

    def handle(self, *args, **options):
        primers = Primer.objects.order_by("name")

        scanned = 0
        matched = 0
        updated = 0
        self.stdout.write("\t".join([
            "primer_uuid",
            "name",
            "direction",
            "sequence_5_overhang",
            "hybridizing_primer",
            "ytk_junction_overhang",
            "ytk_key",
            "ytk_name",
            "ytk_canonical_overhang",
            "ytk_orientation",
            "type_iis_site",
        ]))
        for primer in primers:
            scanned += 1
            if primer.sequence_5:
                continue

            inferred = infer_type_iis_overhang(primer.sequence_3)
            if not inferred:
                continue

            ytk = inferred["ytk"]
            matched += 1
            self.stdout.write(
                "\t".join([
                    str(primer.id),
                    primer.name,
                    primer.fwd_or_rev or "",
                    inferred["sequence_5"],
                    inferred["sequence_3"],
                    inferred["cloning_overhang"],
                    ytk["key"],
                    ytk["name"],
                    ytk["canonical_overhang"],
                    ytk["orientation"],
                    inferred["site"],
                ])
            )

            if options["apply"]:
                primer.sequence_5 = inferred["sequence_5"]
                primer.sequence_3 = inferred["sequence_3"]
                primer.save(update_fields=["sequence_5", "sequence_3"])
                updated += 1

        action = "Updated" if options["apply"] else "Would update"
        self.stdout.write(self.style.SUCCESS(
            f"Scanned {scanned} primers. Found {matched} inferred overhangs. {action} {updated if options['apply'] else matched}."
        ))
