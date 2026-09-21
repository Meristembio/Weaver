import os

from django.core.management.base import BaseCommand
from django.core.management.base import CommandError

from inventory.custom.primer_dimers import PrimerDimerConditions
from inventory.custom.primer_dimers import PrimerDimerInputError
from inventory.custom.primer_dimers import PrimerDimerThresholds
from inventory.custom.primer_dimers import analyze_primers
from inventory.custom.primer_dimers import read_primers
from inventory.custom.primer_dimers import write_markdown_report
from inventory.custom.primer_dimers import write_results_csv


class Command(BaseCommand):
    help = "Analyze primer homodimer and heterodimer risk with primer3-py."

    def add_arguments(self, parser):
        parser.add_argument("input_path", help="FASTA, CSV, or TSV file with primers.")
        parser.add_argument("--format", choices=("auto", "fasta", "csv", "tsv", "table"), default="auto")
        parser.add_argument("--output-prefix", default="primer_dimer_analysis")
        parser.add_argument("--annealing-temp-c", type=float, default=None)
        parser.add_argument("--mv-conc-mM", type=float, default=50.0)
        parser.add_argument("--dv-conc-mM", type=float, default=1.5)
        parser.add_argument("--dntp-conc-mM", type=float, default=0.6)
        parser.add_argument("--primer-conc-nM", type=float, default=250.0)
        parser.add_argument("--calculation-temp-c", type=float, default=37.0)
        parser.add_argument("--critical-3prime-bases", type=int, default=8)
        parser.add_argument("--high-dg", type=float, default=-9.0)
        parser.add_argument("--moderate-dg", type=float, default=-6.0)

    def handle(self, *args, **options):
        conditions = PrimerDimerConditions(
            mv_conc_mM=options["mv_conc_mM"],
            dv_conc_mM=options["dv_conc_mM"],
            dntp_conc_mM=options["dntp_conc_mM"],
            primer_conc_nM=options["primer_conc_nM"],
            annealing_temp_c=options["annealing_temp_c"],
            calculation_temp_c=options["calculation_temp_c"],
            critical_3prime_bases=options["critical_3prime_bases"],
        )
        thresholds = PrimerDimerThresholds(
            high_delta_g_kcal_mol=options["high_dg"],
            moderate_delta_g_kcal_mol=options["moderate_dg"],
        )

        try:
            with open(options["input_path"], encoding="utf-8-sig") as handle:
                primers = read_primers(handle, options["format"])
            analysis = analyze_primers(primers, conditions=conditions, thresholds=thresholds)
        except PrimerDimerInputError as error:
            raise CommandError(str(error))

        output_prefix = options["output_prefix"]
        output_dir = os.path.dirname(output_prefix)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        full_csv = f"{output_prefix}.csv"
        summary_csv = f"{output_prefix}.moderate_high.csv"
        markdown = f"{output_prefix}.md"
        write_results_csv(full_csv, analysis["results"])
        write_results_csv(
            summary_csv,
            [result for result in analysis["results"] if result["risk"] in ("MODERATE", "HIGH", "CALCULATION_ERROR")],
        )
        write_markdown_report(markdown, analysis)

        self.stdout.write(self.style.SUCCESS(
            f"Analyzed {len(primers)} primers and {len(analysis['results'])} interactions."
        ))
        self.stdout.write(f"Full CSV: {full_csv}")
        self.stdout.write(f"Moderate/high CSV: {summary_csv}")
        self.stdout.write(f"Markdown report: {markdown}")
        if conditions.annealing_temp_c is None:
            self.stdout.write(self.style.WARNING(
                "No annealing temperature was provided; Tm-margin risk rules were not applied."
            ))
        for warning in analysis["input_warnings"]:
            self.stdout.write(self.style.WARNING(warning))
