import datetime
import tempfile
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from uuid import uuid4

from Bio import SeqIO
from django.core.files.uploadedfile import SimpleUploadedFile
from django.contrib.auth.models import User
from django.test import RequestFactory
from django.test import SimpleTestCase
from django.test import TestCase
from django.test import override_settings
from django.urls import reverse
from django.utils.datastructures import MultiValueDict
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from django.core.files.base import ContentFile

from inventory.custom.genbank_import import import_plasmids_from_genbank_dir
from inventory.custom.genbank_import import import_plasmids_from_uploaded_genbanks
from inventory.custom.glycerolstock_storage import suggest_storage_positions
from inventory.custom.assembly_classification import classify_assembly_record
from inventory.custom.primer_access import visible_primers_for_user
from inventory.custom.primer_dimers import PrimerDimerConditions
from inventory.custom.primer_dimers import PrimerDimerInputError
from inventory.custom.primer_dimers import PrimerDimerThresholds
from inventory.custom.primer_dimers import PrimerInput
from inventory.custom.primer_dimers import analyze_pair
from inventory.custom.primer_dimers import analyze_primers
from inventory.custom.primer_dimers import primer_combinations
from inventory.custom.primer_dimers import read_primers
from inventory.custom.primer_dimers import validate_primers
from inventory.custom.primer_import import import_primers_from_fasta
from inventory.custom.primer_import import primer_entries_from_fasta
from inventory.custom.pcr import classify_ytk_overhang
from inventory.custom.pcr import find_primer_binding_hits
from inventory.custom.pcr import infer_type_iis_overhang
from inventory.custom.pcr import inferred_primer_parts
from inventory.custom.pcr import matching_amplicon_annotations
from inventory.custom.pcr import primer_pair_amplicons
from inventory.custom.pcr import primer_pair_complementarity
from inventory.custom.pcr import select_non_overlapping_amplicons
from inventory.custom.pcr import suggest_pcr_primers
from inventory.custom.restriction_digest import DigestConstraints
from inventory.custom.restriction_digest import LabEnzyme
from inventory.custom.restriction_digest import compatible_buffers
from inventory.custom.restriction_digest import compatible_temperature
from inventory.custom.restriction_digest import digest_candidates
from inventory.custom.restriction_digest import enzymes_with_effective_cuts
from inventory.custom.restriction_digest import effective_cut_sites
from inventory.custom.restriction_digest import evaluate_digest
from inventory.custom.restriction_digest import fragment_sizes
from inventory.custom.restriction_digest import min_band_difference
from inventory.custom.restriction_digest import normalize_regions
from inventory.custom.restriction_digest import region_contains_position
from inventory.custom.restriction_digest import serialize_digest_response
from inventory.custom.sanger import detect_format
from inventory.custom.sanger import ab1_run_datetime_label
from inventory.custom.sanger import filename_run_datetime
from inventory.custom.sanger import detect_confidence_regions
from inventory.custom.sanger import display_trim_range
from inventory.custom.sanger import include_unaligned_display_flanks
from inventory.custom.sanger import normalized_group_name
from inventory.custom.sanger import parse_phd1
from inventory.custom.sanger import parse_seq
from inventory.custom.sanger import process_sanger_files
from inventory.custom.sanger import UploadedSangerFile
from inventory.custom.sanger import preferred_sanger_run
from inventory.custom.sanger import read_is_usable
from inventory.custom.sanger import SangerProcessingParameters
from inventory.custom.sanger import latest_confirmation_pair
from inventory.views import sanger_failed_read_groups
from inventory.custom.sanger import select_nonredundant_read_candidates
from inventory.custom.sanger import select_sanger_review_candidates
from inventory.forms import PlasmidValidationForm
from inventory.forms import SangerAlignForm
from inventory.models import Primer
from inventory.models import Plasmid
from inventory.models import Experiment
from inventory.models import Box
from inventory.models import GlycerolStock
from inventory.models import Location
from inventory.models import RestrictionBuffer
from inventory.models import RestrictionEnzyme
from inventory.models import RestrictionEnzymeBuffer
from inventory.models import SangerRead
from inventory.models import SangerReadFile
from inventory.models import SangerVerificationRun
from inventory.models import Strain
from inventory.views import fasta_alignment_result
from inventory.views import fasta_record_from_text
from inventory.views import fasta_records_from_text
from inventory.views import amplicon_contains_region
from inventory.views import amplicon_matches_primer_binding_regions
from inventory.views import amplicon_matches_any_primer_id
from inventory.views import amplicon_matches_primer_id
from inventory.views import optional_int_query_param
from inventory.views import plasmid_update_computed_size
from inventory.views import plasmid_validation_initial_from_payload
from inventory.views import build_restriction_enzymes
from inventory.views import run_local_blast
from inventory.views import sanger_feature_color
from inventory.views import sanger_browser_data
from organization.models import Membership
from organization.models import Project


class GlycerolstockStorageSuggestionTests(SimpleTestCase):
    def test_suggestions_order_plasmids_by_idx_and_fill_positions(self):
        boxes = [
            SimpleNamespace(id="box-1", name="L0 B1"),
            SimpleNamespace(id="box-2", name="L0 B2"),
        ]
        plasmids = [
            SimpleNamespace(idx=20, level=0),
            SimpleNamespace(idx=10, level=0),
        ]

        suggestions = suggest_storage_positions(plasmids, boxes)

        self.assertEqual([item["plasmid"].idx for item in suggestions], [10, 20])
        self.assertEqual([item["box"].name for item in suggestions], ["L0 B1", "L0 B1"])
        self.assertEqual([item["position"] for item in suggestions], ["A1", "A2"])

    def test_suggestions_skip_occupied_positions_before_assigning(self):
        box = SimpleNamespace(id="box-1", name="L1 B1")
        plasmid = SimpleNamespace(idx=30, level=1)

        suggestions = suggest_storage_positions(
            [plasmid],
            [box],
            occupied_positions=[("box-1", "A", 1)],
        )

        self.assertEqual(suggestions[0]["position"], "A2")


class GlycerolstockBatchFlowTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="gs-batch-user", password="pw")
        self.project = Project.objects.create(name="GS Batch Project", public=False)
        Membership.objects.create(member=self.user, project=self.project, access_policies="w")
        self.location = Location.objects.create(name="Freezer 1")
        self.box = Box.objects.create(name="L0 B1", location=self.location)
        self.primary_strain = Strain.objects.create(name="Top10", for_primary_gs=True)
        self.plasmid = Plasmid.objects.create(
            idx=100,
            name="Batch plasmid",
            intended_use="Test",
            level=0,
            ligation_state=1,
            colonypcr_state=2,
            digestion_state=2,
            sequencing_state=2,
            project=self.project,
        )
        self.client.force_login(self.user)
        self.client.cookies["current_project_id"] = str(self.project.id)

    def test_prepare_page_tracks_created_gs_and_links_to_print_preview(self):
        response = self.client.post(
            reverse("glycerolstock_batch_prepare"),
            {"plasmid_idx": [str(self.plasmid.idx)]},
        )

        self.assertContains(response, "Pending GS")
        self.assertContains(response, "Create selected GS")
        self.assertContains(response, "batch_plasmid_idx")

        glycerolstock = GlycerolStock.objects.create(
            strain=self.primary_strain,
            plasmid=self.plasmid,
            box=self.box,
            box_row="A",
            box_column=1,
            project=self.project,
        )
        response = self.client.get(
            reverse("glycerolstock_batch_prepare"),
            {"plasmid_idx": str(self.plasmid.idx)},
        )

        self.assertContains(response, "GS created")
        self.assertContains(response, "Print GS labels")
        print_url = response.context["print_url"]
        self.assertIn(str(glycerolstock.id), print_url)

        print_response = self.client.get(print_url)
        self.assertContains(print_response, "Print GS labels")
        self.assertNotContains(print_response, "GS labels ready")
        self.assertContains(print_response, self.plasmid.name)

    def test_create_selected_gs_uses_suggested_position_and_refreshes_batch(self):
        response = self.client.post(
            reverse("glycerolstock_batch_prepare"),
            {
                "batch_create_gs": "1",
                "batch_plasmid_idx": [str(self.plasmid.idx)],
                "plasmid_idx": [str(self.plasmid.idx)],
            },
        )

        self.assertEqual(response.status_code, 302)
        glycerolstock = GlycerolStock.objects.get(plasmid=self.plasmid)
        self.assertEqual(glycerolstock.strain, self.primary_strain)
        self.assertEqual(glycerolstock.box, self.box)
        self.assertEqual(glycerolstock.box_row, "A")
        self.assertEqual(glycerolstock.box_column, 1)

        refreshed = self.client.get(response.url)
        self.assertContains(refreshed, "GS created")
        self.assertContains(refreshed, "Print GS labels")

    def test_batch_creation_uses_selected_same_level_box_and_position(self):
        selected_box = Box.objects.create(name="L0 B2", location=self.location)
        Box.objects.create(name="L1 B1", location=self.location)

        response = self.client.get(
            reverse("glycerolstock_batch_prepare"),
            {"batch_plasmid_idx": str(self.plasmid.idx)},
        )
        self.assertContains(response, "L0 B2")
        self.assertNotContains(response, "L1 B1")

        self.client.post(
            reverse("glycerolstock_batch_prepare"),
            {
                "batch_create_gs": "1",
                "batch_plasmid_idx": [str(self.plasmid.idx)],
                "plasmid_idx": [str(self.plasmid.idx)],
                f"box_id-{self.plasmid.idx}": str(selected_box.id),
                f"position-{self.plasmid.idx}": "B2",
            },
        )

        glycerolstock = GlycerolStock.objects.get(plasmid=self.plasmid)
        self.assertEqual(glycerolstock.box, selected_box)
        self.assertEqual(glycerolstock.box_row, "B")
        self.assertEqual(glycerolstock.box_column, 2)

    def test_batch_creation_can_create_a_new_level_box(self):
        response = self.client.post(
            reverse("glycerolstock_batch_prepare"),
            {
                "batch_create_gs": "1",
                "batch_plasmid_idx": [str(self.plasmid.idx)],
                "plasmid_idx": [str(self.plasmid.idx)],
                f"box_id-{self.plasmid.idx}": "new",
                f"new_box_location-{self.plasmid.idx}": str(self.location.id),
                f"position-{self.plasmid.idx}": "C3",
            },
        )

        self.assertEqual(response.status_code, 302)
        glycerolstock = GlycerolStock.objects.get(plasmid=self.plasmid)
        self.assertEqual(glycerolstock.box.name, "L0 B2")
        self.assertEqual(glycerolstock.box.location, self.location)
        self.assertEqual(glycerolstock.box_row, "C")
        self.assertEqual(glycerolstock.box_column, 3)

    def test_batch_creation_fills_a_new_box_after_existing_positions_run_out(self):
        for row in "ABCDEFG":
            for column in range(1, 10):
                if row == "G" and column > 7:
                    break
                GlycerolStock.objects.create(
                    strain=self.primary_strain,
                    plasmid=None,
                    box=self.box,
                    box_row=row,
                    box_column=column,
                    project=self.project,
                )

        plasmids = [self.plasmid]
        for index in range(101, 130):
            plasmids.append(Plasmid.objects.create(
                idx=index,
                name=f"Overflow plasmid {index}",
                intended_use="Test",
                level=0,
                ligation_state=1,
                colonypcr_state=2,
                digestion_state=2,
                sequencing_state=2,
                project=self.project,
            ))
        indices = [str(plasmid.idx) for plasmid in plasmids]

        response = self.client.get(
            reverse("glycerolstock_batch_prepare"),
            {"batch_plasmid_idx": indices},
        )
        suggestions = response.context["suggestions"]
        self.assertEqual(len(suggestions), 30)
        self.assertEqual(
            [suggestion["position"] for suggestion in suggestions[:20]],
            [f"G{column}" for column in range(8, 10)]
            + [f"{row}{column}" for row in "HI" for column in range(1, 10)],
        )
        self.assertEqual(suggestions[20]["default_box_id"], "new")
        self.assertEqual(suggestions[20]["default_position"], "A1")
        self.assertEqual(suggestions[28]["default_position"], "A9")
        self.assertEqual(suggestions[29]["default_position"], "B1")

        post_data = {
            "batch_create_gs": "1",
            "batch_plasmid_idx": indices,
            "plasmid_idx": indices,
        }
        for suggestion in suggestions:
            plasmid_idx = str(suggestion["plasmid"].idx)
            post_data[f"box_id-{plasmid_idx}"] = suggestion["default_box_id"]
            post_data[f"position-{plasmid_idx}"] = suggestion["default_position"]
            if suggestion["default_box_id"] == "new":
                post_data[f"new_box_location-{plasmid_idx}"] = str(self.location.id)

        response = self.client.post(reverse("glycerolstock_batch_prepare"), post_data)

        self.assertEqual(response.status_code, 302)
        created_stocks = list(
            GlycerolStock.objects.filter(plasmid__in=plasmids).select_related("box")
        )
        self.assertEqual(len(created_stocks), 30)
        new_box_stocks = sorted(
            (stock for stock in created_stocks if stock.box != self.box),
            key=lambda stock: ("ABCDEFGHI".index(stock.box_row), stock.box_column),
        )
        self.assertEqual(len(new_box_stocks), 10)
        self.assertEqual({stock.box.name for stock in new_box_stocks}, {"L0 B2"})
        self.assertEqual(
            [(stock.box_row, stock.box_column) for stock in new_box_stocks],
            [("A", column) for column in range(1, 10)] + [("B", 1)],
        )

    def test_deleting_gs_from_box_keeps_back_to_box_navigation(self):
        glycerolstock = GlycerolStock.objects.create(
            strain=self.primary_strain,
            plasmid=self.plasmid,
            box=self.box,
            box_row="A",
            box_column=1,
            project=self.project,
        )

        box_response = self.client.get(
            reverse("glycerolstock_box", kwargs={"box_id": self.box.id})
        )
        self.assertContains(box_response, f"return_to_box={self.box.id}")

        return_query = f"?return_to_box={self.box.id}"
        delete_response = self.client.post(
            reverse("glycerolstock_delete", kwargs={"pk": glycerolstock.id}) + return_query
        )

        self.assertEqual(delete_response.status_code, 302)
        deleted_response = self.client.get(delete_response.url)
        self.assertContains(deleted_response, "Back to box")
        self.assertContains(
            deleted_response,
            reverse("glycerolstock_box", kwargs={"box_id": self.box.id}),
        )

    def test_prepare_page_accepts_selected_plasmids_from_multiple_projects(self):
        other_project = Project.objects.create(name="Other GS Batch Project", public=False)
        Membership.objects.create(member=self.user, project=other_project, access_policies="w")
        other_plasmid = Plasmid.objects.create(
            idx=101,
            name="Other batch plasmid",
            intended_use="Test",
            level=0,
            ligation_state=1,
            colonypcr_state=2,
            digestion_state=2,
            sequencing_state=2,
            project=other_project,
        )

        response = self.client.post(
            reverse("glycerolstock_batch_prepare"),
            {"plasmid_idx": [str(self.plasmid.idx), str(other_plasmid.idx)]},
        )

        self.assertEqual(
            [suggestion["plasmid"].idx for suggestion in response.context["suggestions"]],
            [self.plasmid.idx, other_plasmid.idx],
        )
        self.assertContains(response, self.plasmid.name)
        self.assertContains(response, other_plasmid.name)

    def test_pichia_gs_does_not_count_as_library_gs(self):
        pichia_strain = Strain.objects.create(name="Pichia Producción")
        GlycerolStock.objects.create(
            strain=pichia_strain,
            plasmid=self.plasmid,
            box=self.box,
            box_row="A",
            box_column=1,
            project=self.project,
        )

        response = self.client.get(
            reverse("glycerolstock_batch_prepare"),
            {"batch_plasmid_idx": str(self.plasmid.idx)},
        )

        self.assertContains(response, "Pending GS")
        self.assertEqual(response.context["created_count"], 0)
        self.assertFalse(response.context["all_created"])


class BatchPrintLabelFormatTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="batch-print-user", password="pw")
        self.project = Project.objects.create(name="Batch Print Project", public=False)
        Membership.objects.create(member=self.user, project=self.project, access_policies="w")
        self.plasmid = Plasmid.objects.create(
            idx=200,
            name="Print format plasmid",
            computed_size=1234,
            intended_use="Test",
            project=self.project,
        )
        self.client.force_login(self.user)

    def test_plasmid_batch_uses_individual_label_markup_and_date_format(self):
        response = self.client.post(
            reverse("services-batch-prints"),
            {
                "label_type": "plasmids",
                "identifier": str(self.plasmid.idx),
                "colony": "7",
                "date": "2026-01-09",
                "concentration": "12.5",
            },
        )

        self.assertContains(response, 'class="label-name"')
        self.assertContains(response, 'class="label-created"')
        self.assertContains(response, 'class="label-conc-quantus"')
        self.assertContains(response, 'class="label-size"')
        self.assertContains(response, 'class="batch-print-page"')
        self.assertContains(response, "09.Jan.2026")
        self.assertContains(response, "ID 200 ~ c7")
        self.assertNotContains(response, "batch-label-title")


class GlycerolstockListViewTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="stock-list-user", password="pw")
        self.project = Project.objects.create(name="Stock List Project", public=False)
        Membership.objects.create(member=self.user, project=self.project, access_policies="w")
        self.location = Location.objects.create(name="Stock freezer")
        self.box = Box.objects.create(name="Stock box", location=self.location)
        self.strain = Strain.objects.create(name="Stock strain")
        self.plasmids = [
            Plasmid.objects.create(
                idx=idx,
                name=name,
                intended_use="Test",
                project=self.project,
            )
            for idx, name in ((101, "Stock A"), (102, "Stock B"), (103, "Stock C"))
        ]
        for column, plasmid in enumerate(self.plasmids, start=1):
            GlycerolStock.objects.create(
                strain=self.strain,
                plasmid=plasmid,
                box=self.box,
                box_row="A",
                box_column=column,
                project=self.project,
            )
        self.client.force_login(self.user)
        self.client.cookies["current_project_id"] = str(self.project.id)

    def test_stocks_view_uses_server_search_and_versioned_asset(self):
        response = self.client.get(reverse("glycerolstocks"))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'id="glycerolstocks-list" data-server-search="true"')
        self.assertContains(response, 'placeholder="Search Stocks ..."')
        self.assertContains(response, 'src="/static/js/glycerolstocks_table.js?version=14"')

    def test_api_glycerolstocks_paginates_in_stable_order(self):
        response = self.client.get(reverse("api-glycerolstocks"), {"page_size": 2})

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload["total"], 3)
        self.assertEqual(payload["page"], 1)
        self.assertEqual(payload["page_size"], 2)
        self.assertEqual(payload["num_pages"], 2)
        self.assertEqual([stock["pn"] for stock in payload["glycerolstocks"]], ["Stock A", "Stock B"])

        second_page = self.client.get(reverse("api-glycerolstocks"), {
            "page": 2,
            "page_size": 2,
        }).json()
        self.assertEqual([stock["pn"] for stock in second_page["glycerolstocks"]], ["Stock C"])

    def test_api_glycerolstocks_searches_name_id_and_all_fields(self):
        response = self.client.get(reverse("api-glycerolstocks"), {
            "q": "Stock strain",
            "search": "name",
        })
        self.assertEqual(response.json()["total"], 3)

        response = self.client.get(reverse("api-glycerolstocks"), {
            "q": "102",
            "search": "idx",
        })
        self.assertEqual([stock["pn"] for stock in response.json()["glycerolstocks"]], ["Stock B"])

        response = self.client.get(reverse("api-glycerolstocks"), {
            "q": "Stock freezer",
            "search": "all",
        })
        self.assertEqual(response.json()["total"], 3)


def primer(name, sequence_3, direction):
    return SimpleNamespace(
        id=uuid4(),
        name=name,
        sequence_3=sequence_3,
        sequence_5="",
        fwd_or_rev=direction,
    )


def primer_with_overhang(name, sequence_3, sequence_5, direction):
    return SimpleNamespace(
        id=uuid4(),
        name=name,
        sequence_3=sequence_3,
        sequence_5=sequence_5,
        fwd_or_rev=direction,
    )


class RestrictionEnzymeStub(SimpleNamespace):
    def __str__(self):
        return self.name


def restriction_enzyme(name, **activities):
    defaults = {
        "activity_buffer_1_1": 100,
        "activity_buffer_2_1": 100,
        "activity_buffer_3_1": 100,
        "activity_buffer_CS": 100,
        "activity_buffer_aari": 100,
    }
    defaults.update(activities)
    return RestrictionEnzymeStub(name=name, hf_version=False, **defaults)


def lab_enzyme(name, temperature=37, activities=None):
    return LabEnzyme(
        name=name,
        display_name=name,
        recognition_site="TEST",
        fcut=1,
        rcut=1,
        temperature=temperature,
        activities=activities or {
            "NEB 1.1": 100,
            "NEB 2.1": 100,
            "NEB 3.1": 100,
            "NEB CutSmart": 100,
            "Thermo AarI": 100,
        },
    )


class FastaParsingTests(SimpleTestCase):
    def test_fasta_records_accept_comment_lines_before_header(self):
        records = fasta_records_from_text("; exported from aligner\n# comment\n>read1\nACGT\n")

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, "read1")
        self.assertEqual(str(records[0].seq), "ACGT")

    def test_fasta_record_from_text_wraps_plain_sequence(self):
        record = fasta_record_from_text("ACGTACGT", name="Amplicon")

        self.assertEqual(record.id, "Amplicon")
        self.assertEqual(str(record.seq), "ACGTACGT")

    def test_fasta_alignment_result_uses_sanger_result_shape(self):
        record = fasta_record_from_text(">query\nGTACGT\n")
        result = fasta_alignment_result("AACGTACGTT", record)

        self.assertEqual(result["combined"]["useful_reads"], 1)
        self.assertEqual(result["reads"][0]["formats"], ["FASTA"])
        self.assertEqual(result["reads"][0]["alignment"]["start_display"], 4)
        self.assertIn("reference_projection", result["reads"][0]["alignment"])

    def test_fasta_alignment_result_accepts_multiple_records_and_detects_reverse(self):
        records = fasta_records_from_text(">forward\nACGTTGC\n>reverse\nGCAACGT\n")
        result = fasta_alignment_result("AAAACGTTGCAAAA", records)

        self.assertEqual(result["combined"]["useful_reads"], 2)
        self.assertEqual(result["reads"][0]["alignment"]["best_orientation"], "forward")
        self.assertEqual(result["reads"][1]["alignment"]["best_orientation"], "reverse")


class SangerVerificationServiceTests(SimpleTestCase):
    def test_short_trimmed_sequence_has_user_friendly_reason(self):
        usable, reason = read_is_usable(
            "A" * 12,
            {"trimmed_length": 12},
            [],
            SangerProcessingParameters(),
        )

        self.assertFalse(usable)
        self.assertEqual(reason, "The usable sequence is too short for reliable alignment.")

    def test_display_alignment_exposes_unaligned_flanks_as_low_quality_insertions(self):
        alignment = {
            "start": 20,
            "end": 40,
            "query_start": 2,
            "query_end": 8,
            "best_orientation": "forward",
            "oriented_sequence": "AACCGGTTAA",
            "variants": [],
        }

        display_alignment = include_unaligned_display_flanks(
            alignment,
            [5, 8, 35, 35, 35, 35, 35, 35, 7, 4],
            0,
            10,
            100,
        )

        self.assertEqual(len(display_alignment["variants"]), 4)
        self.assertEqual(
            [variant["base_index"] for variant in display_alignment["variants"]],
            [0, 1, 8, 9],
        )
        self.assertTrue(all(variant["low_quality"] for variant in display_alignment["variants"]))

    def test_display_alignment_keeps_low_quality_flanks_visible(self):
        self.assertEqual(
            display_trim_range(
                100,
                {"alignment_blocks": [(12, 82)], "low_confidence_regions": [{"start": 0, "end": 99}]},
            ),
            (0, 100),
        )

    def test_sanger_browser_coordinates_start_at_plasmid_origin(self):
        reference = "ACGT" * 10
        alignment = {
            "start": 7,
            "end": 15,
            "segments": [],
            "crosses_origin": False,
            "identity": 100.0,
            "reference_projection": reference,
            "reference_projection_base_indices": list(range(len(reference))),
            "variants": [],
        }

        browser_data = sanger_browser_data(
            reference,
            {
                "reads": [{
                    "name": "read-1",
                    "is_usable": True,
                    "alignment": alignment,
                    "display_alignment": alignment,
                    "chromatogram": {},
                    "quality_metrics": {},
                }],
                "combined": {"depth": [], "uncovered_regions": []},
            },
        )

        self.assertEqual(browser_data["displayOrigin"], 0)

    def test_sanger_browser_data_preserves_display_only_insertions(self):
        reference = "ACGT" * 10
        alignment = {
            "start": 0,
            "end": 9,
            "segments": [],
            "crosses_origin": False,
            "identity": 100.0,
            "reference_projection": reference,
            "reference_projection_base_indices": list(range(len(reference))),
            "variants": [{
                "coordinate": 9,
                "type": "insertion",
                "observed": "A",
                "quality": 4,
                "low_quality": True,
                "display_only": True,
                "base_index": 0,
            }],
        }

        browser_data = sanger_browser_data(
            reference,
            {
                "reads": [{
                    "name": "read-1",
                    "is_usable": True,
                    "alignment": alignment,
                    "display_alignment": alignment,
                    "chromatogram": {},
                    "quality_metrics": {},
                }],
                "combined": {"depth": [], "uncovered_regions": []},
            },
        )

        self.assertTrue(browser_data["reads"][0]["variants"][0]["display_only"])

    def test_normalizes_complementary_sanger_files_to_one_group(self):
        base = "sample_C02"

        self.assertEqual(normalized_group_name(base + ".ab1"), base)
        self.assertEqual(normalized_group_name(base + ".phd.1"), base)
        self.assertEqual(normalized_group_name(base + ".seq"), base)
        self.assertEqual(detect_format(base + ".phd.1"), "phd1")

    def test_parse_seq_accepts_plain_or_fasta_sequence(self):
        parsed = parse_seq(b">read one\nacgt n\nRYSW\n")

        self.assertFalse(parsed.errors)
        self.assertEqual(parsed.sequence, "ACGTNRYSW")

    def test_parse_seq_reports_invalid_characters(self):
        parsed = parse_seq(b"ACGTZ\n")

        self.assertTrue(parsed.errors)
        self.assertIn("Z", parsed.errors[0])

    def test_parse_phd1_reads_bases_qualities_and_peaks(self):
        text = b"""BEGIN_SEQUENCE read1\r\nBEGIN_COMMENT\r\nCHROMAT_FILE: read1.ab1\r\nEND_COMMENT\r\nBEGIN_DNA\r\nA 31 10\r\nC 20 22\r\nEND_DNA\r\nEND_SEQUENCE\r\n"""

        parsed = parse_phd1(text)

        self.assertFalse(parsed.errors)
        self.assertEqual(parsed.sequence, "AC")
        self.assertEqual(parsed.qualities, [31, 20])
        self.assertEqual(parsed.peak_positions, [10, 22])
        self.assertEqual(parsed.metadata["CHROMAT_FILE"], "read1.ab1")

    def test_duplicate_same_format_in_group_is_rejected(self):
        files = [
            SimpleUploadedFile("sample.seq", b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
            SimpleUploadedFile("sample.seq", b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        ]

        result = process_sanger_files(files, "ACGT" * 30)

        self.assertEqual(len(result["reads"]), 1)
        self.assertFalse(result["reads"][0]["is_usable"])
        self.assertIn("Multiple SEQ files", result["reads"][0]["errors"][0])

    def test_detects_reverse_complement_orientation(self):
        reference = "ATGCGTACCTAGGATCCGATTAACCGGTTAGCTAGGCTAATCGGATCCGGAATTCGCGGCCGC"
        reverse_read = str(Seq(reference[8:58]).reverse_complement())
        files = [SimpleUploadedFile("reverse.seq", reverse_read.encode())]

        result = process_sanger_files(files, reference)

        self.assertTrue(result["reads"][0]["is_usable"])
        self.assertEqual(result["reads"][0]["alignment"]["best_orientation"], "reverse")

    def test_partial_sanger_coverage_does_not_prevent_pass_by_itself(self):
        reference = "ACGT" * 500

        from inventory.custom.sanger import SangerProcessingParameters
        from inventory.custom.sanger import align_read
        from inventory.custom.sanger import classify_run
        from inventory.custom.sanger import combined_metrics

        alignment = align_read(reference, "NNNNNN" + reference[100:600], [40] * 506, 0, SangerProcessingParameters())
        reads = [{
            "name": "partial",
            "is_usable": True,
            "alignment": alignment,
            "warnings": [],
            "errors": [],
        }]
        combined = combined_metrics(reference, reads, SangerProcessingParameters())
        classification = classify_run(combined, reads, SangerProcessingParameters())

        self.assertLess(combined["combined_coverage"], 95)
        self.assertEqual(classification["state"], "PASS")
        self.assertEqual(alignment["query_start"], 6)
        self.assertEqual(alignment["reference_projection_base_indices"][alignment["start"]], 6)

    def test_seq_without_quality_reduces_confidence_without_high_quality_variants(self):
        reference = "ACGT" * 80
        read_sequence = reference[:40] + "T" + reference[41:80]
        files = [SimpleUploadedFile("no-quality.seq", read_sequence.encode())]

        result = process_sanger_files(files, reference)

        self.assertEqual(result["combined"]["high_quality_variant_count"], 0)
        self.assertEqual(result["classification"]["state"], "REVIEW")
        self.assertIn("Quality scores unavailable", result["classification"]["reasons"][0])

    def test_confidence_regions_detect_progressive_3prime_drop(self):
        sequence = "A" * 120
        qualities = [32] * 80 + [18] * 10 + [9] * 30

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertTrue(metrics["low_confidence_regions"])
        self.assertGreaterEqual(metrics["low_confidence_regions"][0]["start"], 75)
        self.assertEqual(metrics["low_confidence_regions"][-1]["end"], 119)

    def test_terminal_quality_collapse_backtracks_to_first_nearby_warning(self):
        sequence = "A" * 1200
        qualities = [32] * 1200
        qualities[977] = 15
        qualities[1056:] = [8] * (1200 - 1056)

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(metrics["low_confidence_regions"][-1]["start"], 977)
        self.assertEqual(metrics["low_confidence_regions"][-1]["end"], 1199)
        self.assertEqual(metrics["alignment_blocks"][-1][1], 977)
        self.assertIn("Terminal quality collapse", " ".join(metrics["low_confidence_regions"][-1]["reasons"]))

    def test_confidence_regions_detect_bad_5prime_start(self):
        sequence = "A" * 120
        qualities = [8] * 25 + [31] * 95

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(metrics["low_confidence_regions"][0]["start"], 0)
        self.assertGreater(metrics["accepted_blocks"][0][0], 0)

    def test_confidence_regions_ignore_single_bad_base_inside_good_read(self):
        sequence = "A" * 100
        qualities = [32] * 100
        qualities[50] = 3

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(metrics["low_confidence_regions"], [])

    def test_confidence_regions_merge_bad_regions_separated_by_short_gap(self):
        sequence = "A" * 130
        qualities = [32] * 130
        qualities[30:45] = [8] * 15
        qualities[52:67] = [8] * 15

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(len(metrics["low_confidence_regions"]), 1)
        self.assertLessEqual(metrics["low_confidence_regions"][0]["start"], 30)
        self.assertGreaterEqual(metrics["low_confidence_regions"][0]["end"], 66)

    def test_confidence_regions_short_recovery_does_not_split_region(self):
        sequence = "A" * 150
        qualities = [32] * 150
        qualities[30:60] = [8] * 30
        qualities[60:70] = [32] * 10
        qualities[70:100] = [8] * 30

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(len(metrics["low_confidence_regions"]), 1)

    def test_confidence_regions_sustained_recovery_closes_region(self):
        sequence = "A" * 180
        qualities = [32] * 180
        qualities[25:45] = [8] * 20
        qualities[80:100] = [8] * 20

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(len(metrics["low_confidence_regions"]), 2)

    def test_high_phred_low_signal_is_warning_not_low_region(self):
        sequence = "A" * 100
        qualities = [32] * 100
        base_pos = [i * 10 for i in range(100)]
        trace_len = base_pos[-1] + 20
        a_trace = [100] * trace_len
        for index in range(35, 50):
            a_trace[base_pos[index]] = 10
        chromatogram = {
            "basePos": base_pos,
            "aTrace": a_trace,
            "cTrace": [1] * trace_len,
            "gTrace": [1] * trace_len,
            "tTrace": [1] * trace_len,
        }

        metrics = detect_confidence_regions(sequence, qualities, chromatogram)

        self.assertEqual(metrics["low_confidence_regions"], [])
        self.assertTrue(metrics["signal"]["warnings"])

    def test_secondary_peaks_with_low_phred_support_low_region(self):
        sequence = "A" * 100
        qualities = [32] * 100
        qualities[35:55] = [10] * 20
        base_pos = [i * 10 for i in range(100)]
        trace_len = base_pos[-1] + 20
        a_trace = [100] * trace_len
        c_trace = [1] * trace_len
        for index in range(35, 55):
            c_trace[base_pos[index]] = 70
        chromatogram = {
            "basePos": base_pos,
            "aTrace": a_trace,
            "cTrace": c_trace,
            "gTrace": [1] * trace_len,
            "tTrace": [1] * trace_len,
        }

        metrics = detect_confidence_regions(sequence, qualities, chromatogram)

        self.assertTrue(metrics["low_confidence_regions"])
        self.assertIn("Chromatogram signal morphology", metrics["low_confidence_regions"][0]["reasons"])

    def test_clear_range_intersects_weaver_alignment_blocks(self):
        sequence = "A" * 120
        qualities = [32] * 120
        qualities[:20] = [8] * 20
        qualities[100:] = [8] * 20

        metrics = detect_confidence_regions(sequence, qualities, clear_range=(10, 90))

        self.assertEqual(metrics["file_clear_range"], (10, 90))
        self.assertGreaterEqual(metrics["alignment_blocks"][0][0], 10)
        self.assertLessEqual(metrics["alignment_blocks"][-1][1], 90)

    def test_seq_without_quality_reports_quality_unavailable_and_aligns_with_warning(self):
        reference = "ACGT" * 40
        files = [SimpleUploadedFile("plain.seq", reference[:100].encode())]

        result = process_sanger_files(files, reference)

        self.assertFalse(result["reads"][0]["quality_metrics"]["quality_available"])
        self.assertTrue(result["reads"][0]["is_usable"])
        self.assertIn("Quality scores unavailable", result["reads"][0]["warnings"][0])

    def test_internal_low_quality_region_splits_read_into_aligned_blocks(self):
        reference = "ACGT" * 90
        read = reference[:150] + reference[180:330]
        qualities = [35] * 150 + [35] * 150
        sequence = read[:80] + ("N" * 40) + read[120:]
        qualities[80:120] = [3] * 40
        phd_lines = ["BEGIN_SEQUENCE split", "BEGIN_COMMENT", "END_COMMENT", "BEGIN_DNA"]
        phd_lines.extend("{} {} {}".format(base, quality, index * 10) for index, (base, quality) in enumerate(zip(sequence, qualities)))
        phd_lines.extend(["END_DNA", "END_SEQUENCE"])

        result = process_sanger_files([SimpleUploadedFile("split.phd.1", "\n".join(phd_lines).encode())], reference)

        alignment = result["reads"][0]["alignment"]
        self.assertTrue(result["reads"][0]["quality_metrics"]["low_confidence_regions"])
        self.assertGreaterEqual(len(alignment["segments"]), 2)
        self.assertEqual(result["combined"]["high_quality_variant_count"], 0)

    def test_multiple_forward_and_reverse_reads_combine_high_confidence_coverage(self):
        reference = "ATGCGTAC" * 80
        reverse_read = str(Seq(reference[300:520]).reverse_complement())
        files = [
            SimpleUploadedFile("left.seq", reference[:220].encode()),
            SimpleUploadedFile("right.seq", reverse_read.encode()),
        ]

        result = process_sanger_files(files, reference)

        self.assertEqual(result["combined"]["useful_reads"], 2)
        self.assertIn("reverse", {read["alignment"]["best_orientation"] for read in result["reads"]})
        self.assertGreater(result["combined"]["combined_coverage"], 30)

    def test_real_complementary_files_are_grouped_when_available(self):
        base = Path("/home/tproschle-sticta/Downloads/PUC-6.jul.26/Diego Lagos 060726/6096-Diego-495_P6_Stuffer_c1-457-pL0-chech-F_2026-07-07-13-15-37_C02")
        paths = [Path(str(base) + suffix) for suffix in (".ab1", ".phd.1", ".seq")]
        if not all(path.exists() for path in paths):
            self.skipTest("Real Sanger fixture files are not available on this machine")
        files = [SimpleUploadedFile(path.name, path.read_bytes()) for path in paths]
        reference = parse_seq(paths[2].read_bytes()).sequence

        result = process_sanger_files(files, reference)

        self.assertEqual(len(result["reads"]), 1)
        self.assertIn("ab1", result["reads"][0]["formats"])
        self.assertIn("phd1", result["reads"][0]["formats"])
        self.assertIn("seq", result["reads"][0]["formats"])


class SangerFeatureColorTests(SimpleTestCase):
    def test_known_feature_roles_use_stable_palette(self):
        self.assertEqual(sanger_feature_color("CDS", "CamR", {"color": ["#0000ff"]}), "#2fb344")
        self.assertEqual(sanger_feature_color("misc_feature", "pTDH3 promoter", {}), "#f0b429")
        self.assertEqual(sanger_feature_color("terminator", "CYC1 terminator", {}), "#d94841")
        self.assertEqual(sanger_feature_color("restriction_site", "BsaI", {}), "#8f63d9")
        self.assertEqual(sanger_feature_color("protein_bind", "BsaI", {}), "#8f63d9")
        self.assertEqual(sanger_feature_color("misc_feature", "attB1 recombination site", {}), "#e66fb2")
        self.assertEqual(sanger_feature_color("misc_recomb", "Hr1Chr4_Up", {}), "#e66fb2")
        self.assertEqual(sanger_feature_color("misc_recomb", "HR1-Chr4", {}), "#e66fb2")

    def test_unknown_feature_keeps_explicit_genbank_color(self):
        self.assertEqual(sanger_feature_color("misc_feature", "custom tag", {"ApEinfo_fwdcolor": ["#123456"]}), "#123456")


class PcrSuggestionTests(SimpleTestCase):
    def test_infers_type_iis_overhang_from_full_primer_sequence(self):
        inferred = infer_type_iis_overhang("aaCGTCTCtctccTATGcgtaaaggcgaagag")

        self.assertEqual(inferred["sequence_5"], "AACGTCTCTCTCCTATG")
        self.assertEqual(inferred["sequence_3"], "CGTAAAGGCGAAGAG")
        self.assertEqual(inferred["cloning_overhang"], "TATG")
        self.assertEqual(inferred["ytk"]["key"], "2-3")
        self.assertEqual(inferred["ytk"]["orientation"], "forward")

    def test_infers_reverse_primer_type_iis_overhang_from_full_sequence(self):
        inferred = infer_type_iis_overhang("aaCGTCTCtctcgGGATtttgtacagttcatccataccatg")

        self.assertEqual(inferred["sequence_5"], "AACGTCTCTCTCGGGAT")
        self.assertEqual(inferred["sequence_3"], "TTTGTACAGTTCATCCATACCATG")
        self.assertEqual(inferred["cloning_overhang"], "GGAT")
        self.assertEqual(inferred["ytk"]["key"], "3-4")
        self.assertEqual(inferred["ytk"]["canonical_overhang"], "ATCC")
        self.assertEqual(inferred["ytk"]["orientation"], "reverse_complement")

    def test_classifies_non_ytk_overhang_as_empty(self):
        self.assertEqual(classify_ytk_overhang("AAAA")["key"], "")

    def test_amplicon_finder_uses_inferred_hybridizing_sequence(self):
        sequence = (
            "CGTAAAGGCGAAGAG" +
            "AAAAC" +
            str(Seq("TTTGTACAGTTCATCCATACCATG").reverse_complement())
        )
        primers = [
            primer("693-L0-P3-sfGFP-F", "aaCGTCTCtctccTATGcgtaaaggcgaagag", "f"),
            primer("694-L0-P3-sfGFP-R", "aaCGTCTCtctcgGGATtttgtacagttcatccataccatg", "r"),
        ]

        amplicons = matching_amplicon_annotations(
            sequence,
            primers,
            min_product_size=1,
            max_product_size=999,
            max_tm_difference=99,
        )

        self.assertEqual(len(amplicons), 1)
        self.assertEqual(amplicons[0]["notes"]["fwd_primer"], ["L0-P3-sfGFP-F"])
        self.assertEqual(amplicons[0]["notes"]["rev_primer"], ["L0-P3-sfGFP-R"])
        self.assertEqual(
            amplicons[0]["notes"]["amplicon_sequence"],
            ["AACGTCTCTCTCCTATG" + sequence + str(Seq("AACGTCTCTCTCGGGAT").reverse_complement())],
        )

    def test_existing_declared_overhang_is_not_reinferred(self):
        parts = inferred_primer_parts(
            primer_with_overhang("declared-F", "aaCGTCTCtctccTATGcgtaaaggcgaagag", "TT", "f")
        )

        self.assertFalse(parts["inferred"])
        self.assertEqual(parts["sequence_5"], "TT")
        self.assertEqual(parts["sequence_3"], "AACGTCTCTCTCCTATGCGTAAAGGCGAAGAG")

    def test_primer_binding_hits_can_match_from_3prime_end(self):
        hits = find_primer_binding_hits(
            "AAAACCCCGGGGTTTT",
            "GGGGGAAAACCCCGGGGTTTT",
            min_binding_length=15,
        )

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["start"], 0)
        self.assertEqual(hits[0]["binding_sequence"], "AAAACCCCGGGGTTTT")
        self.assertEqual(hits[0]["unmatched_5"], "GGGGG")

    def test_amplicon_finder_allows_3prime_binding_without_full_primer_containment(self):
        fwd_binding = "AAAACCCCGGGGTTTT"
        rev_binding = "TTTTGGGGCCCCAAAA"
        sequence = fwd_binding + "GG" + str(Seq(rev_binding).reverse_complement())
        primers = [
            primer_with_overhang("partial-F", "GGGGG" + fwd_binding, "AA", "f"),
            primer_with_overhang("partial-R", "CCCCC" + rev_binding, "TT", "r"),
        ]

        amplicons = matching_amplicon_annotations(
            sequence,
            primers,
            min_product_size=1,
            max_product_size=999,
            max_tm_difference=99,
        )

        self.assertEqual(len(amplicons), 1)
        self.assertEqual(amplicons[0]["notes"]["partial_primer_binding"], ["true"])
        self.assertEqual(amplicons[0]["notes"]["fwd_partial_binding"], ["true"])
        self.assertEqual(amplicons[0]["notes"]["rev_partial_binding"], ["true"])
        self.assertEqual(amplicons[0]["notes"]["fwd_binding_length"], ["16"])
        self.assertEqual(amplicons[0]["notes"]["rev_binding_length"], ["16"])
        self.assertEqual(amplicons[0]["notes"]["fwd_unmatched_5"], ["GGGGG"])
        self.assertEqual(amplicons[0]["notes"]["rev_unmatched_5"], ["CCCCC"])
        self.assertIn("FWD partial 3' binding: 16 bp aligned, 5 bp 5' unaligned", amplicons[0]["notes"]["warnings"])
        self.assertIn("REV partial 3' binding: 16 bp aligned, 5 bp 5' unaligned", amplicons[0]["notes"]["warnings"])
        self.assertEqual(
            amplicons[0]["notes"]["amplicon_sequence"],
            ["AA" + "GGGGG" + sequence + str(Seq("CCCCC").reverse_complement()) + "AA"],
        )

    def test_primer_pair_amplicons_flags_partial_binding(self):
        fwd_binding = "AAAACCCCGGGGTTTT"
        rev_binding = "TTTTGGGGCCCCAAAA"
        sequence = fwd_binding + "GG" + str(Seq(rev_binding).reverse_complement())
        primer_f = primer("partial-F", "GGGGG" + fwd_binding, "f")
        primer_r = primer("partial-R", "CCCCC" + rev_binding, "r")

        amplicons = primer_pair_amplicons(
            sequence,
            primer_f,
            primer_r,
            max_product_size=999,
        )

        self.assertEqual(len(amplicons), 1)
        self.assertTrue(amplicons[0]["partial_primer_binding"])
        self.assertTrue(amplicons[0]["fwd_partial_binding"])
        self.assertTrue(amplicons[0]["rev_partial_binding"])

    def test_amplicon_region_filter_requires_region_inside_amplicon(self):
        amplicon = {"start": 10, "end": 30}
        inside = SimpleNamespace(start=12, end=20)
        outside = SimpleNamespace(start=5, end=12)

        self.assertTrue(amplicon_contains_region(amplicon, inside, 100))
        self.assertFalse(amplicon_contains_region(amplicon, outside, 100))

    def test_amplicon_region_filter_allows_configured_flank(self):
        amplicon = {"start": 6374, "end": 6844}
        region = SimpleNamespace(start=6369, end=6743)

        self.assertFalse(amplicon_contains_region(amplicon, region, 9171))
        self.assertTrue(amplicon_contains_region(amplicon, region, 9171, flank_bp=30))

    def test_amplicon_region_filter_supports_circular_amplicons(self):
        amplicon = {"start": 90, "end": 10, "overlapsSelf": True}
        inside = SimpleNamespace(start=95, end=5)
        outside = SimpleNamespace(start=20, end=30)

        self.assertTrue(amplicon_contains_region(amplicon, inside, 100))
        self.assertFalse(amplicon_contains_region(amplicon, outside, 100))

    def test_amplicon_binding_region_filter_uses_primer_three_prime_position(self):
        amplicon = {
            "notes": {
                "fwd_3prime_position": ["3"],
                "rev_3prime_position": ["8"],
            },
        }
        fwd_zone = SimpleNamespace(start=3, end=3)
        rev_zone = SimpleNamespace(start=8, end=8)
        wrong_zone = SimpleNamespace(start=4, end=7)

        self.assertTrue(amplicon_matches_primer_binding_regions(amplicon, (fwd_zone,), "fwd", 100))
        self.assertTrue(amplicon_matches_primer_binding_regions(amplicon, (rev_zone,), "rev", 100))
        self.assertFalse(amplicon_matches_primer_binding_regions(amplicon, (wrong_zone,), "fwd", 100))

    def test_amplicon_binding_region_filter_supports_circular_zones_and_alternatives(self):
        amplicon = {
            "notes": {
                "fwd_3prime_position": ["2"],
            },
        }
        circular_zone = SimpleNamespace(start=90, end=5)
        alternative_zone = SimpleNamespace(start=2, end=2)

        self.assertTrue(amplicon_matches_primer_binding_regions(amplicon, (circular_zone,), "fwd", 100))
        self.assertTrue(amplicon_matches_primer_binding_regions(amplicon, (alternative_zone,), "fwd", 100))

    def test_matching_amplicons_record_primer_three_prime_positions(self):
        amplicons = matching_amplicon_annotations(
            "AAAACCCCGGGGTTTT",
            [
                primer("1-left-F", "AAAA", "f"),
                primer("2-right-R", "CCCC", "r"),
            ],
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(amplicons[0]["notes"]["fwd_3prime_position"], ["3"])
        self.assertEqual(amplicons[0]["notes"]["rev_3prime_position"], ["8"])

    def test_amplicon_primer_id_filter_matches_either_primer(self):
        amplicon = {"notes": {"fwd_primer_id": ["693"], "rev_primer_id": ["694"]}}

        self.assertTrue(amplicon_matches_primer_id(amplicon, "693"))
        self.assertTrue(amplicon_matches_primer_id(amplicon, "694"))
        self.assertFalse(amplicon_matches_primer_id(amplicon, "695"))
        self.assertTrue(amplicon_matches_any_primer_id(amplicon, ("695", "694")))
        self.assertFalse(amplicon_matches_any_primer_id(amplicon, ("695", "696")))

    def test_amplicon_size_filter_limits_are_optional(self):
        factory = RequestFactory()

        empty_request = factory.get("/amplicons", {"min_size": "", "max_size": ""})
        min_request = factory.get("/amplicons", {"min_size": "200", "max_size": ""})
        max_request = factory.get("/amplicons", {"min_size": "", "max_size": "1200"})

        self.assertEqual(optional_int_query_param(empty_request, "min_size", 100), 100)
        self.assertIsNone(optional_int_query_param(empty_request, "max_size"))
        self.assertEqual(optional_int_query_param(min_request, "min_size", 100), 200)
        self.assertEqual(optional_int_query_param(max_request, "max_size"), 1200)

    def test_suggests_pair_covering_selection(self):
        sequence = "AAAACCCCGGGGTTTT"
        primers = [
            primer("1-left-F", "AAAA", "f"),
            primer("2-right-R", "CCCC", "r"),
        ]

        suggestions = suggest_pcr_primers(
            sequence,
            primers,
            4,
            7,
            margin=8,
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(len(suggestions), 1)
        self.assertEqual(suggestions[0]["category"], "Full selection")
        self.assertEqual(suggestions[0]["selection_coverage"], 4)
        self.assertEqual(suggestions[0]["product_range"], "1..12")

    def test_suggests_circular_pair_crossing_origin(self):
        sequence = "AAACCCGGGTTT"
        primers = [
            primer("1-origin-F", "TTT", "f"),
            primer("2-origin-R", "GGG", "r"),
        ]

        suggestions = suggest_pcr_primers(
            sequence,
            primers,
            10,
            2,
            margin=5,
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(len(suggestions), 1)
        self.assertIn("Circular", suggestions[0]["category"])
        self.assertEqual(suggestions[0]["selection_coverage"], 5)
        self.assertEqual(suggestions[0]["product_range"], "10..6 (circular)")

    def test_default_minimum_product_size_excludes_short_amplicons(self):
        sequence = "AAAACCCCGGGGTTTT"
        primers = [
            primer("1-left-F", "AAAA", "f"),
            primer("2-right-R", "CCCC", "r"),
        ]

        suggestions = suggest_pcr_primers(sequence, primers, 4, 7, margin=8)

        self.assertEqual(suggestions, [])

    def test_default_tm_filter_excludes_dissimilar_pairs(self):
        sequence = "GGGGAAAACCCCTTTTGGGGCCCC"
        primers = [
            primer("1-low-tm-F", "AAAA", "f"),
            primer("2-high-tm-R", "GGGG", "r"),
        ]

        suggestions = suggest_pcr_primers(sequence, primers, 4, 7, margin=20, min_product_size=1)

        self.assertEqual(suggestions, [])

    def test_matching_amplicons_are_returned_as_parts(self):
        sequence = "AAAACCCCGGGGTTTT"
        primers = [
            primer("1-left-F", "AAAA", "f"),
            primer("2-right-R", "CCCC", "r"),
        ]

        amplicons = matching_amplicon_annotations(
            sequence,
            primers,
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(len(amplicons), 1)
        self.assertEqual(amplicons[0]["annotationTypePlural"], "parts")
        self.assertTrue(amplicons[0]["id"].startswith("weaver-amplicon-"))
        self.assertIn("tm_difference", amplicons[0]["notes"])
        self.assertEqual(amplicons[0]["notes"]["fwd_tm"], ["8.0"])
        self.assertEqual(amplicons[0]["notes"]["rev_tm"], ["16.0"])
        self.assertEqual(amplicons[0]["notes"]["product_tm"], ["40.0"])
        self.assertEqual(amplicons[0]["notes"]["recommended_annealing_tm"], ["15.5"])

    def test_matching_amplicons_include_copyable_amplicon_sequence(self):
        sequence = "AAAACCCCGGGGTTTT"
        primers = [
            primer_with_overhang("1-left-F", "AAAA", "TT", "f"),
            primer_with_overhang("2-right-R", "CCCC", "GG", "r"),
        ]

        amplicons = matching_amplicon_annotations(
            sequence,
            primers,
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(amplicons[0]["notes"]["amplicon_sequence"], ["TTAAAACCCCGGGGCC"])
        self.assertEqual(amplicons[0]["notes"]["product_size"], ["16"])
        self.assertIn("primer_complementarity_alignment_fwd", amplicons[0]["notes"])
        self.assertIn("primer3_dimer_risk", amplicons[0]["notes"])
        self.assertIn("primer3_dimer_dg", amplicons[0]["notes"])

    def test_circular_amplicon_visual_interval_matches_real_product(self):
        sequence = "AAACCCGGGTTT"
        primers = [
            primer("1-origin-F", "TTT", "f"),
            primer("2-origin-R", "GGG", "r"),
        ]

        amplicons = matching_amplicon_annotations(
            sequence,
            primers,
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(len(amplicons), 1)
        self.assertTrue(amplicons[0]["overlapsSelf"])
        self.assertEqual(amplicons[0]["notes"]["template_size"], ["9"])
        self.assertEqual(amplicons[0]["notes"]["visual_start"], ["9"])
        self.assertEqual(amplicons[0]["notes"]["visual_end"], ["5"])
        self.assertEqual(amplicons[0]["notes"]["visual_size"], ["9"])
        self.assertEqual(amplicons[0]["notes"]["visual_is_alternative"], ["false"])

    def test_selects_non_overlapping_amplicons(self):
        amplicons = [
            {"id": "a", "start": 0, "end": 9},
            {"id": "b", "start": 5, "end": 15},
            {"id": "c", "start": 16, "end": 20},
        ]

        selected = select_non_overlapping_amplicons(amplicons, 30)

        self.assertEqual([amplicon["id"] for amplicon in selected], ["a", "c"])

    def test_selects_non_overlapping_amplicons_across_origin(self):
        amplicons = [
            {"id": "a", "start": 18, "end": 3, "overlapsSelf": True},
            {"id": "b", "start": 2, "end": 5},
            {"id": "c", "start": 8, "end": 10},
        ]

        selected = select_non_overlapping_amplicons(amplicons, 20)

        self.assertEqual([amplicon["id"] for amplicon in selected], ["a", "c"])

    def test_primer_pair_amplicons_finds_product(self):
        sequence = "AAAACCCCGGGGTTTT"
        primer_f = primer("1-left-F", "AAAA", "f")
        primer_r = primer("2-right-R", "CCCC", "r")

        amplicons = primer_pair_amplicons(sequence, primer_f, primer_r)

        self.assertEqual(len(amplicons), 1)
        self.assertEqual(amplicons[0]["product_range"], "1..12")
        self.assertEqual(amplicons[0]["product_size"], 12)
        self.assertFalse(amplicons[0]["circular"])

    def test_primer_pair_amplicons_finds_circular_product(self):
        sequence = "AAACCCGGGTTT"
        primer_f = primer("1-origin-F", "TTT", "f")
        primer_r = primer("2-origin-R", "GGG", "r")

        amplicons = primer_pair_amplicons(sequence, primer_f, primer_r)

        self.assertEqual(len(amplicons), 1)
        self.assertEqual(amplicons[0]["product_range"], "10..6 (circular)")
        self.assertEqual(amplicons[0]["product_size"], 9)
        self.assertTrue(amplicons[0]["circular"])
        self.assertIn("primer3_dimer", amplicons[0])

    def test_primer_pair_complementarity_flags_3prime_dimer_risk(self):
        primer_f = primer("1-risk-F", "AAAACCCC", "f")
        primer_r = primer("2-risk-R", "GGGGTTTT", "r")

        complementarity = primer_pair_complementarity(primer_f, primer_r)

        self.assertEqual(complementarity["severity"], "high")
        self.assertEqual(complementarity["max_both_3prime_contiguous"], 8)
        self.assertEqual(set(complementarity["alignment"].keys()), {"fwd", "match", "rev"})
        self.assertIn("|", complementarity["alignment"]["match"])
        self.assertIn("FWD 5'", complementarity["alignment"]["fwd"])
        self.assertIn("REV 3'", complementarity["alignment"]["rev"])
        self.assertIn("AAAACCCC", complementarity["alignment"]["fwd"])
        self.assertIn("TTTTGGGG", complementarity["alignment"]["rev"])

    def test_primer_pair_complementarity_allows_unmatched_pair(self):
        primer_f = primer("1-ok-F", "AAAA", "f")
        primer_r = primer("2-ok-R", "AAAA", "r")

        complementarity = primer_pair_complementarity(primer_f, primer_r)

        self.assertEqual(complementarity["severity"], "none")
        self.assertEqual(complementarity["max_contiguous"], 0)
        self.assertEqual(set(complementarity["alignment"].keys()), {"fwd", "match", "rev"})
        self.assertIn("AAAA", complementarity["alignment"]["fwd"])
        self.assertIn("AAAA", complementarity["alignment"]["rev"])

    def test_primer_pair_complementarity_does_not_warn_for_trivial_matches(self):
        primer_f = primer("1-trivial-F", "GCCTTTTGCTGGCCTTTTGC", "f")
        primer_r = primer("2-trivial-R", "CTGTGTTGACATCTGGTTTG", "r")

        complementarity = primer_pair_complementarity(primer_f, primer_r)

        self.assertEqual(complementarity["severity"], "none")
        self.assertEqual(complementarity["max_contiguous"], 2)
        self.assertEqual(complementarity["max_both_3prime_contiguous"], 1)
        self.assertEqual(complementarity["warnings"], [])
        self.assertIn("||", complementarity["alignment"]["match"])


class RestrictionDigestTests(SimpleTestCase):
    def test_circular_plasmid_without_cuts_has_no_fragments(self):
        result = evaluate_digest(
            "AAAACCCCGGGG",
            [lab_enzyme("EcoRI")],
            DigestConstraints(min_fragment_size_bp=0, min_band_difference_bp=0),
            is_circular=True,
        )

        self.assertEqual(result["cut_count"], 0)
        self.assertEqual(result["fragment_count"], 0)
        self.assertFalse(result["exact"])

    def test_circular_plasmid_with_one_cut_returns_full_length_fragment(self):
        result = evaluate_digest(
            "AAAAGAATTCTTT",
            [lab_enzyme("EcoRI")],
            DigestConstraints(min_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=0),
            is_circular=True,
        )

        self.assertEqual(result["fragments_map_order"], [13])

    def test_single_fragment_does_not_check_band_separation(self):
        result = evaluate_digest(
            "AAAAGAATTCTTT",
            [lab_enzyme("EcoRI")],
            DigestConstraints(min_fragments=1, max_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=500),
            is_circular=True,
        )

        self.assertTrue(result["exact"])
        self.assertFalse(any(item["criterion"] == "minimum band separation bp" for item in result["violations"]))

    def test_circular_plasmid_with_multiple_cuts_includes_origin_fragment(self):
        result = evaluate_digest(
            "AAAAGAATTCTTTGAATTCAAA",
            [lab_enzyme("EcoRI")],
            DigestConstraints(min_fragments=2, min_fragment_size_bp=1, min_band_difference_bp=1),
            is_circular=True,
        )

        self.assertEqual([site["position"] for site in result["cut_sites"]], [5, 14])
        self.assertEqual(result["fragments_map_order"], [9, 13])

    def test_linear_fragment_generation_zero_one_and_many_cuts(self):
        self.assertEqual(fragment_sizes(10, [], is_circular=False), [10])
        self.assertEqual(fragment_sizes(10, [4], is_circular=False), [4, 6])
        self.assertEqual(fragment_sizes(10, [2, 7], is_circular=False), [2, 5, 3])

    def test_coincident_cuts_keep_enzyme_identity(self):
        result = evaluate_digest(
            "AAAAGAATTCTTT",
            [lab_enzyme("EcoRI"), lab_enzyme("EcoRI")],
            DigestConstraints(min_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=0),
            is_circular=True,
        )

        self.assertEqual(result["cut_count"], 1)
        self.assertEqual(result["cut_sites"][0]["enzymes"], ["EcoRI"])

    def test_enzyme_with_multiple_sites(self):
        sites = effective_cut_sites("AAAAGAATTCTTTGAATTCAAA", lab_enzyme("EcoRI"), is_circular=False)

        self.assertEqual([site["position"] for site in sites], [5, 14])

    def test_iupac_motif_site_is_found(self):
        sites = effective_cut_sites("AAAGAATCAAA", lab_enzyme("HinfI"), is_circular=False)

        self.assertEqual([site["position"] for site in sites], [4])

    def test_type_iis_cut_can_be_outside_recognition_site(self):
        sites = effective_cut_sites("AAAAGGTCTCAAAAATTTT", lab_enzyme("BsaI"), is_circular=False)

        self.assertEqual([site["position"] for site in sites], [11])

    def test_reverse_strand_type_iis_site_on_circular_sequence(self):
        sites = effective_cut_sites("AAAAGAGACCTTTTTAAAA", lab_enzyme("BsaI"), is_circular=True)

        self.assertEqual([site["position"] for site in sites], [18])

    def test_regions_include_normal_and_origin_crossing_intervals(self):
        normal, circular = normalize_regions([
            {"start": 2, "end": 5},
            {"start": 8, "end": 1},
        ], 10)

        self.assertTrue(region_contains_position(normal, 3, 10))
        self.assertFalse(region_contains_position(normal, 7, 10))
        self.assertTrue(region_contains_position(circular, 9, 10))
        self.assertTrue(region_contains_position(circular, 0, 10))
        self.assertFalse(region_contains_position(circular, 4, 10))

    def test_minimum_band_difference_at_threshold_is_accepted(self):
        self.assertEqual(min_band_difference([500, 1000, 1500]), 500)

    def test_equal_fragment_sizes_have_zero_band_separation(self):
        self.assertEqual(min_band_difference([1000, 1000]), 0)

    def test_fragment_count_above_range_is_not_exact(self):
        result = evaluate_digest(
            "AAAAGAATTCTTTGAATTCAAAGAATTCTTT",
            [lab_enzyme("EcoRI")],
            DigestConstraints(min_fragments=1, max_fragments=2, min_fragment_size_bp=1, min_band_difference_bp=0),
            is_circular=True,
        )

        self.assertEqual(result["fragment_count"], 3)
        self.assertFalse(result["exact"])
        self.assertIn("above range 1-2", result["violations"][0]["message"])

    def test_best_buffer_maximizes_minimum_then_average_activity(self):
        left = lab_enzyme("Left", activities={
            "NEB 1.1": 100,
            "NEB 2.1": 90,
            "NEB 3.1": 80,
            "NEB CutSmart": 80,
            "Thermo AarI": 10,
        })
        right = lab_enzyme("Right", activities={
            "NEB 1.1": 80,
            "NEB 2.1": 90,
            "NEB 3.1": 100,
            "NEB CutSmart": 80,
            "Thermo AarI": 10,
        })

        buffers = compatible_buffers([left, right], 75)

        self.assertEqual(buffers[0]["name"], "NEB 2.1")
        self.assertEqual(buffers[0]["min_activity"], 90)

    def test_incompatible_temperatures_are_not_exact(self):
        temperature, is_compatible = compatible_temperature([
            lab_enzyme("EcoRI", temperature=37),
            lab_enzyme("BamHI", temperature=42),
        ])

        self.assertEqual(temperature, 37)
        self.assertFalse(is_compatible)

    def test_unknown_activity_or_temperature_is_not_exact(self):
        enzyme = lab_enzyme("EcoRI", temperature=None, activities={
            "NEB 1.1": None,
            "NEB 2.1": None,
            "NEB 3.1": None,
            "NEB CutSmart": None,
            "Thermo AarI": None,
        })

        result = evaluate_digest(
            "AAAAGAATTCTTT",
            [enzyme],
            DigestConstraints(min_fragments=0, min_fragment_size_bp=0, min_band_difference_bp=0),
        )

        self.assertEqual(result["temperature_status"], "unknown")
        self.assertFalse(result["exact"])
        self.assertTrue(any(item["criterion"] == "shared buffer" for item in result["violations"]))

    def test_exact_results_rank_before_closest_matches_deterministically(self):
        enzymes = [
            restriction_enzyme("EcoRI"),
            restriction_enzyme("HinfI"),
        ]

        results = digest_candidates(
            "AAAGAATCAAAAGAATTCAAA",
            enzymes,
            DigestConstraints(min_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=0, limit=10),
            is_circular=True,
        )

        self.assertTrue(results[0]["exact"])
        self.assertEqual(results, sorted(results, key=lambda result: (not result["exact"], len(result["enzymes"]), "+".join(result["enzymes"]))))

    def test_closest_match_reports_concrete_violations(self):
        results = digest_candidates(
            "AAAAGAATTCTTT",
            [restriction_enzyme("EcoRI")],
            DigestConstraints(min_fragments=2, min_fragment_size_bp=500, min_band_difference_bp=500),
            is_circular=True,
        )

        self.assertFalse(results[0]["exact"])
        self.assertTrue(results[0]["violations"])
        self.assertIn("missing", results[0]["violations"][0])

    def test_digest_candidates_skip_enzymes_without_effective_cuts(self):
        results = digest_candidates(
            "AAAAGAATTCTTT",
            [restriction_enzyme("EcoRI"), restriction_enzyme("BamHI")],
            DigestConstraints(min_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=0, limit=10),
            is_circular=True,
        )

        self.assertTrue(results)
        self.assertTrue(all("BamHI" not in result["enzyme_names"] for result in results))

    def test_enzymes_with_effective_cuts_filters_non_cutters(self):
        enzymes = enzymes_with_effective_cuts(
            "AAAAGAATTCTTT",
            [restriction_enzyme("EcoRI"), restriction_enzyme("BamHI")],
            is_circular=True,
        )

        self.assertEqual([enzyme.name for enzyme in enzymes], ["EcoRI"])

    def test_digest_candidates_can_require_an_enzyme(self):
        results = digest_candidates(
            "AAAGAATCAAAAGAATTCAAA",
            [restriction_enzyme("EcoRI"), restriction_enzyme("HinfI")],
            DigestConstraints(
                min_fragments=1,
                min_fragment_size_bp=1,
                min_band_difference_bp=0,
                limit=10,
                required_enzymes=("HinfI",),
            ),
            is_circular=True,
        )

        self.assertTrue(results)
        self.assertTrue(all("HinfI" in result["enzyme_names"] for result in results))

    def test_serialized_contract_contains_frontend_card_fields(self):
        response = serialize_digest_response(
            "AAAAGAATTCTTTGAATTCAAA",
            [restriction_enzyme("EcoRI")],
            DigestConstraints(min_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=0),
            is_circular=True,
        )
        result = response["results"][0]

        for key in [
            "enzymes",
            "best_buffer",
            "cut_sites",
            "fragment_count",
            "fragments_map_order",
            "fragments_by_size",
            "min_band_difference",
            "regions",
            "status",
            "violations",
        ]:
            self.assertIn(key, result)


class PrimerBatchImportTests(TestCase):
    def test_defaults_missing_direction_to_forward(self):
        fasta = StringIO(">123-no-direction\nAAAACCCC\n")

        result = import_primers_from_fasta(fasta, default_direction="f")

        imported = Primer.objects.get(name="123-no-direction")
        self.assertEqual(result["created"], 1)
        self.assertEqual(imported.name, "123-no-direction")
        self.assertEqual(imported.sequence_3, "AAAACCCC")
        self.assertEqual(imported.fwd_or_rev, "f")

    def test_parser_infers_reverse_direction_from_name(self):
        entries = primer_entries_from_fasta(StringIO(">826_Primer-R\nAAAACCCC\n"), default_direction="f")

        self.assertEqual(entries[0]["direction"], "r")

    def test_fasta_header_overhang_is_split_from_full_sequence(self):
        entries = primer_entries_from_fasta(
            StringIO(">1001-L0-P3a-ATR-F overhang=aaCGTCTCtctcc\naaCGTCTCtctccTATGACTTCTGCTTTGTATGCATCAG\n")
        )

        self.assertEqual(entries[0]["name"], "1001-L0-P3a-ATR-F")
        self.assertEqual(entries[0]["direction"], "f")
        self.assertEqual(entries[0]["sequence_5"], "aaCGTCTCtctcc")
        self.assertEqual(entries[0]["sequence"], "TATGACTTCTGCTTTGTATGCATCAG")

    def test_primer_import_view_redirects_to_primer_list_with_summary_counts(self):
        user = User.objects.create_user(username="primer-import-user", password="pw")
        project = Project.objects.create(name="Primer Import Project", public=False)
        Membership.objects.create(member=user, project=project, access_policies="w")
        self.client.force_login(user)
        self.client.cookies["current_project_id"] = str(project.id)

        upload = SimpleUploadedFile(
            "primers.fasta",
            b">1001-Test-F\nAAAACCCC\n>1002-Test-R\nGGGGTTTT\n",
            content_type="text/plain",
        )

        response = self.client.post(reverse("primer_import"), {
            "fasta_file": upload,
            "name_source": "id",
            "default_direction": "f",
            "update_existing": "",
        })

        self.assertEqual(response.status_code, 302)
        self.assertEqual(
            response.url,
            reverse("primers") + "?form_result_primer_import_success=true&primer_import_changed=2&primer_import_created=2&primer_import_updated=0&primer_import_skipped=0&primer_import_errors=0",
        )
        self.assertTrue(Primer.objects.filter(name="1001-Test-F").exists())
        self.assertTrue(Primer.objects.filter(name="1002-Test-R").exists())

        summary_response = self.client.get(response.url)
        self.assertContains(summary_response, "Primer batch import complete! 2 primers loaded into Weaver.")
        self.assertContains(summary_response, "Created 2, updated 0, skipped 0, errors 0.")

    def test_primer_delete_view_renders_confirmation_page(self):
        user = User.objects.create_user(username="primer-delete-user", password="pw")
        project = Project.objects.create(name="Primer Delete Project", public=False)
        Membership.objects.create(member=user, project=project, access_policies="w")
        primer = Primer.objects.create(name="1001-Test-F", sequence_3="AAAACCCC", fwd_or_rev="f")
        self.client.force_login(user)

        response = self.client.get(reverse("primer_delete", args=(primer.id,)))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Are you sure you want to delete")
        self.assertContains(response, "1001-Test-F")


class PrimerAccessTests(TestCase):
    def test_visible_primers_are_global_for_authenticated_users(self):
        user = User.objects.create_user(username="weaver-user")
        visible_project_a = Project.objects.create(name="Visible A", public=False)
        visible_project_b = Project.objects.create(name="Visible B", public=False)
        hidden_project = Project.objects.create(name="Hidden", public=False)
        Primer.objects.create(name="1-global-a-F", sequence_3="AAAA", fwd_or_rev="f")
        Primer.objects.create(name="2-global-b-R", sequence_3="CCCC", fwd_or_rev="r")
        Primer.objects.create(name="3-global-c-F", sequence_3="GGGG", fwd_or_rev="f")

        primer_names = set(visible_primers_for_user(user).values_list("name", flat=True))

        self.assertEqual(primer_names, {"1-global-a-F", "2-global-b-R", "3-global-c-F"})

    def test_fasta_header_direction_can_replace_name_suffix(self):
        entries = primer_entries_from_fasta(
            StringIO(">1001-L0-P3a-ATR overhang=aaCGTCTCtctcc direction=F\naaCGTCTCtctccTATGACTTCTGCTTTGTATGCATCAG\n"),
            require_direction=True,
        )

        self.assertEqual(entries[0]["name"], "1001-L0-P3a-ATR")
        self.assertEqual(entries[0]["direction"], "f")
        self.assertEqual(entries[0]["sequence_5"], "aaCGTCTCtctcc")
        self.assertEqual(entries[0]["sequence"], "TATGACTTCTGCTTTGTATGCATCAG")


class PrimerListViewTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="primer-list-user", password="pw")
        self.project_a = Project.objects.create(name="Primer Project A", public=False)
        self.project_b = Project.objects.create(name="Primer Project B", public=False)
        Membership.objects.create(member=self.user, project=self.project_a, access_policies="w")
        Membership.objects.create(member=self.user, project=self.project_b, access_policies="w")
        Primer.objects.create(name="1001-Project-A-F", sequence_3="AAAACCCC", fwd_or_rev="f")
        Primer.objects.create(name="2001-Project-B-R", sequence_3="GGGGTTTT", fwd_or_rev="r")
        self.client.force_login(self.user)
        self.client.cookies["current_project_id"] = str(self.project_a.id)

    def test_primers_view_shows_all_global_primers(self):
        response = self.client.get(reverse("primers"))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'id="root-primers"')
        self.assertContains(response, 'src="/static/js/primers_table.js?version=13"')

    def test_api_primers_paginates_in_descending_display_sort(self):
        Primer.objects.create(name="3001-Project-C-F", sequence_3="CCCC", fwd_or_rev="f")

        response = self.client.get(reverse("api-primers"), {"page_size": 2})

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload["total"], 3)
        self.assertEqual(payload["page"], 1)
        self.assertEqual(payload["page_size"], 2)
        self.assertEqual(payload["num_pages"], 2)
        self.assertEqual(
            [primer["display_idx"] for primer in payload["primers"]],
            [3001, 2001],
        )

        second_page = self.client.get(reverse("api-primers"), {
            "page": 2,
            "page_size": 2,
        }).json()
        self.assertEqual([primer["display_idx"] for primer in second_page["primers"]], [1001])

    def test_api_primers_searches_id_name_and_all_fields(self):
        response = self.client.get(reverse("api-primers"), {
            "q": "2001",
            "search": "idx",
        })
        self.assertEqual([primer["display_name"] for primer in response.json()["primers"]], ["Project-B-R"])

        response = self.client.get(reverse("api-primers"), {
            "q": "Project-B",
            "search": "name",
        })
        self.assertEqual([primer["display_name"] for primer in response.json()["primers"]], ["Project-B-R"])

        response = self.client.get(reverse("api-primers"), {
            "q": "GGGG",
            "search": "all",
        })
        self.assertEqual([primer["display_name"] for primer in response.json()["primers"]], ["Project-B-R"])

    def test_primers_view_ignores_project_toggle_cookie(self):
        self.client.cookies["show_from_all_projects"] = "True"

        response = self.client.get(reverse("primers"))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'id="root-primers"')
        self.assertNotContains(response, 'id="show_from_all_projects"')
        self.assertContains(response, 'title="Batch upload primers"')
        self.assertContains(response, 'title="Create primer"')

    def test_primers_view_does_not_render_project_column(self):
        response = self.client.get(reverse("primers"))

        self.assertEqual(response.status_code, 200)
        self.assertNotContains(response, "<th scope=\"col\">Project</th>", html=True)


class LocalBlastFormattingTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="blast-user", password="pw")
        self.project = Project.objects.create(name="Blast Project", public=False)
        Membership.objects.create(member=self.user, project=self.project, access_policies="w")
        self.plasmid = Plasmid.objects.create(
            idx=96,
            name="pYTK096",
            intended_use="test",
            project=self.project,
        )

    def test_run_local_blast_enriches_subject_display_fields(self):
        request = RequestFactory().get("/inventory/services/blast/")
        request.user = self.user
        context = {}
        record = fasta_record_from_text("ACGTACGT", name="FWD + REV")
        blast_result = [{
            "query": {"name": "FWD + REV"},
            "subject": {"name": "Subject fallback", "origin_record_id": str(self.plasmid.id)},
            "meta": {"query seq": "ACGT", "subject seq": "ACGT"},
        }]

        original_seqio_get = run_local_blast.__globals__["seqio_get"]
        original_make_circular = run_local_blast.__globals__["make_circular"]
        original_make_linear = run_local_blast.__globals__["make_linear"]
        original_bio_blast = run_local_blast.__globals__["BioBlast"]
        original_run_pyblast_compat = run_local_blast.__globals__["run_pyblast_compat"]

        class FakeBioBlast:
            def __init__(self, subjects, queries):
                self.subjects = subjects
                self.queries = queries

            def blastn(self):
                return blast_result

            def blastn_short(self):
                return blast_result

        try:
            run_local_blast.__globals__["seqio_get"] = lambda plasmid: (True, SeqRecord(Seq("ACGTACGT"), id=str(plasmid.id), name=plasmid.name))
            run_local_blast.__globals__["make_circular"] = lambda records: records
            run_local_blast.__globals__["make_linear"] = lambda records: records
            run_local_blast.__globals__["BioBlast"] = FakeBioBlast
            run_local_blast.__globals__["run_pyblast_compat"] = lambda callback: callback()

            run_local_blast(request, context, record, project_id='a', short_blast=False)
        finally:
            run_local_blast.__globals__["seqio_get"] = original_seqio_get
            run_local_blast.__globals__["make_circular"] = original_make_circular
            run_local_blast.__globals__["make_linear"] = original_make_linear
            run_local_blast.__globals__["BioBlast"] = original_bio_blast
            run_local_blast.__globals__["run_pyblast_compat"] = original_run_pyblast_compat

        self.assertEqual(context["results"][0]["subject_plasmid_id"], self.plasmid.id)
        self.assertEqual(context["results"][0]["subject_display_name"], "pYTK096")
        self.assertEqual(context["results"][0]["subject_display_idx"], 96)
        self.assertIn("FWD_+_REV", context["results"][0]["alignment"])


class SangerVerificationEntryTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="seq-user", password="pw")
        self.project = Project.objects.create(name="Sticta", public=False)
        Membership.objects.create(member=self.user, project=self.project, access_policies="r")
        self.sanger_primer = Primer.objects.create(
            name="7100-Sanger-Sequence-F",
            sequence_3="ACGTACGTACGT",
            fwd_or_rev="f",
        )
        self.plasmid = Plasmid.objects.create(
            idx=499,
            name="L0-P8b-Hr1Chr4_Up",
            intended_use="test",
            project=self.project,
        )
        self.client.force_login(self.user)

    def test_entry_redirects_to_upload_when_no_sanger_run_exists(self):
        response = self.client.get(reverse("plasmid_seq_verification_entry", kwargs={"weaver_id": self.plasmid.idx}))

        self.assertRedirects(response, reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}))

    def test_entry_with_trailing_slash_redirects_to_upload(self):
        response = self.client.get(reverse("plasmid_seq_verification_entry_slash", kwargs={"weaver_id": self.plasmid.idx}))

        self.assertRedirects(response, reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}))

    def test_uuid_entry_redirects_to_upload_when_no_idx_link_is_available(self):
        response = self.client.get(reverse("plasmid_seq_verification_entry_uuid", kwargs={"plasmid_id": self.plasmid.id}))

        self.assertRedirects(response, reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}))

    def test_entry_keeps_upload_page_when_a_verified_run_exists(self):
        SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user)
        SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user, manual_decision="VERIFIED")

        response = self.client.get(reverse("plasmid_seq_verification_entry", kwargs={"weaver_id": self.plasmid.idx}))

        self.assertRedirects(response, reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}), fetch_redirect_response=False)

    def test_entry_keeps_upload_page_when_previous_run_exists(self):
        SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user)
        SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user)

        response = self.client.get(reverse("plasmid_seq_verification_entry", kwargs={"weaver_id": self.plasmid.idx}))

        self.assertRedirects(response, reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}), fetch_redirect_response=False)

    def test_upload_page_lists_previous_files_without_loading_alignment(self):
        run = SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user)
        read = SangerRead.objects.create(run=run, name="forward")
        SangerReadFile.objects.create(
            read=read,
            format="ab1",
            original_name="forward.ab1",
            sha256="a" * 64,
            size=7,
        )

        response = self.client.get(reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Saved Sanger runs")
        self.assertContains(response, "forward.ab1")
        self.assertContains(response, "Return to run tables")
        self.assertNotContains(response, 'id="sanger-run-selector-list"')
        self.assertNotContains(response, "Select one or more Sanger trace files.")
        self.assertContains(response, 'id="sanger-upload-form"')
        self.assertContains(response, 'accept=".ab1,.phd.1,.seq,.fa,.fas,.fasta"')
        self.assertContains(response, "multiple")
        self.assertNotContains(response, 'id="sanger-map"')
        self.assertLess(
            response.content.index(b"Saved Sanger runs"),
            response.content.index(b'id="sanger-upload-form"'),
        )
        self.assertContains(response, "Select primer")

    def test_upload_page_orders_runs_by_sequencing_date(self):
        latest_sequencing_run = SangerVerificationRun.objects.create(
            plasmid=self.plasmid,
            created_by=self.user,
        )
        latest_read = SangerRead.objects.create(run=latest_sequencing_run, name="latest")
        SangerReadFile.objects.create(
            read=latest_read,
            format="ab1",
            original_name="latest_2025-02-01-12-00-00.ab1",
            sha256="b" * 64,
            size=7,
        )
        older_sequencing_run = SangerVerificationRun.objects.create(
            plasmid=self.plasmid,
            created_by=self.user,
        )
        older_read = SangerRead.objects.create(run=older_sequencing_run, name="older")
        SangerReadFile.objects.create(
            read=older_read,
            format="ab1",
            original_name="older_2025-01-01-12-00-00.ab1",
            sha256="c" * 64,
            size=7,
        )

        response = self.client.get(reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}))

        self.assertEqual(response.status_code, 200)
        self.assertLess(
            response.content.index(b"latest_2025-02-01-12-00-00.ab1"),
            response.content.index(b"older_2025-01-01-12-00-00.ab1"),
        )

    def test_upload_requires_a_primer_for_every_sanger_file(self):
        response = self.client.post(
            reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}),
            {"sanger_files": SimpleUploadedFile("forward.ab1", b"trace-data")},
        )

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Select a valid primer for every sequencing file.")
        self.assertFalse(SangerVerificationRun.objects.filter(plasmid=self.plasmid).exists())

    def test_successful_upload_redirects_to_saved_run(self):
        service_result = {
            "parameters": {},
            "reads": [],
            "combined": {"variants": []},
            "classification": {"state": "NO_DATA", "reasons": []},
        }
        reference = SeqRecord(Seq("ACGT" * 20), id="redirect-reference")

        with patch("inventory.views.process_sanger_files", return_value=service_result), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))), \
                patch("inventory.views.plasmid_seqrecord", return_value=reference), \
                patch("inventory.views.plasmid_sequence_file_contents", return_value=""):
            response = self.client.post(
                reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}),
                {
                    "label": "first upload",
                    "sanger_files": SimpleUploadedFile("forward.ab1", b"trace-data"),
                    "primer_id": str(self.sanger_primer.id),
                },
            )

            run = SangerVerificationRun.objects.get(plasmid=self.plasmid)
            detail_url = reverse("sanger_run_detail", kwargs={
                "plasmid_id": self.plasmid.id,
                "run_id": run.id,
            })
            self.assertRedirects(response, detail_url + "?uploaded=1", fetch_redirect_response=False)

            refreshed = self.client.get(response.url)
            self.assertEqual(refreshed.status_code, 200)
            self.assertContains(refreshed, "Sanger sequencing files uploaded and saved.")
            self.assertEqual(SangerVerificationRun.objects.filter(plasmid=self.plasmid).count(), 1)

    def test_rerun_reprocesses_saved_files_into_new_run(self):
        membership = Membership.objects.get(member=self.user, project=self.project)
        membership.access_policies = "w"
        membership.save(update_fields=["access_policies"])
        original_run = SangerVerificationRun.objects.create(
            plasmid=self.plasmid,
            created_by=self.user,
            label="Original run",
            notes="Keep this note",
            manual_decision="VERIFIED",
        )
        original_read = SangerRead.objects.create(
            run=original_run,
            name="forward",
            selected_source="ab1",
        )
        source_file = SangerReadFile.objects.create(
            read=original_read,
            primer=self.sanger_primer,
            format="ab1",
            original_name="forward.ab1",
            sha256="a" * 64,
            size=10,
        )
        source_file.file.save("forward.ab1", ContentFile(b"trace-data"), save=True)
        service_result = {
            "parameters": {},
            "uploaded_files": [],
            "reads": [{
                "name": "forward",
                "formats": ["ab1"],
                "selected_source": "ab1",
                "files": [UploadedSangerFile(
                    original_name="forward.ab1",
                    data=b"trace-data",
                    size=10,
                    sha256="a" * 64,
                    format="ab1",
                    group_name="forward",
                )],
                "alignment": {},
                "display_alignment": {},
                "raw_sequence": "",
                "trimmed_sequence": "",
                "quality_metrics": {},
                "chromatogram": {},
                "warnings": [],
                "errors": [],
                "is_usable": False,
            }],
            "combined": {"variants": [], "read_count": 1, "useful_reads": 0},
            "classification": {"state": "NO_DATA", "reasons": []},
        }
        with patch("inventory.views.process_sanger_files", return_value=service_result), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))), \
                patch("inventory.views.plasmid_seqrecord", return_value=SeqRecord(Seq("ACGT" * 20), id="rerun-reference")):
            response = self.client.post(reverse("sanger_run_rerun", kwargs={
                "plasmid_id": self.plasmid.id,
                "run_id": original_run.id,
            }))

        self.assertEqual(response.status_code, 302)
        self.assertIn("?rerun=1", response.url)
        with patch("inventory.views.plasmid_sequence_file_contents", return_value=""):
            detail_response = self.client.get(response.url)
        self.assertEqual(detail_response.status_code, 200)
        self.assertContains(detail_response, "Sanger run reprocessed successfully.")
        self.assertContains(detail_response, "Re-run analysis")
        self.assertContains(detail_response, 'btn btn-outline-warning btn-sm')
        self.assertContains(detail_response, 'btn btn-outline-danger btn-sm sanger-delete-button')
        self.assertEqual(SangerVerificationRun.objects.filter(plasmid=self.plasmid).count(), 1)
        original_run.refresh_from_db()
        self.assertEqual(original_run.label, "Original run")
        self.assertEqual(original_run.notes, "Keep this note")
        self.assertEqual(original_run.manual_decision, "")
        self.assertTrue(SangerReadFile.objects.filter(
            original_name="forward.ab1",
            read__run=original_run,
            primer=self.sanger_primer,
        ).exists())

class SangerUploadFormTests(SimpleTestCase):
    def test_accepts_multiple_ab1_files_in_one_submission(self):
        files = MultiValueDict({
            "sanger_files": [
                SimpleUploadedFile("forward.ab1", b"trace-a"),
                SimpleUploadedFile("reverse.ab1", b"trace-b"),
            ],
        })

        form = SangerAlignForm({}, files)

        self.assertTrue(form.is_valid())
        self.assertEqual([file.name for file in form.cleaned_data["sanger_files"]], ["forward.ab1", "reverse.ab1"])


class SangerBatchUploadTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="sanger-batch-user", password="pw")
        self.project = Project.objects.create(name="Sanger batch project", public=False)
        Membership.objects.create(member=self.user, project=self.project, access_policies="w")
        self.batch_primer = Primer.objects.create(
            name="7000-Batch-Sequence-F",
            sequence_3="ACGTACGTACGT",
            fwd_or_rev="f",
        )
        self.first_plasmid = Plasmid.objects.create(
            idx=601,
            name="Batch plasmid one",
            intended_use="test",
            project=self.project,
        )
        self.second_plasmid = Plasmid.objects.create(
            idx=602,
            name="Batch plasmid two",
            intended_use="test",
            project=self.project,
        )
        self.client.force_login(self.user)

    @staticmethod
    def service_result():
        return {
            "parameters": {},
            "reads": [],
            "combined": {"variants": []},
            "classification": {"state": "NO_DATA", "reasons": []},
        }

    def test_batch_upload_groups_ab1_files_by_plasmid_id(self):
        files = {
            "ab1_files": [
                SimpleUploadedFile("one-forward.ab1", b"one-forward"),
                SimpleUploadedFile("one-reverse.ab1", b"one-reverse"),
                SimpleUploadedFile("two-forward.ab1", b"two-forward"),
            ],
            "mapping_csv": SimpleUploadedFile(
                "mapping.csv",
                b"ab1_file,plasmid_id,primer_id\none-forward.ab1,601,7000\none-reverse.ab1,601,7000\ntwo-forward.ab1,602,7000\n",
                content_type="text/csv",
            ),
        }
        calls = []

        def fake_process(uploaded_files, reference_sequence):
            calls.append([file.name for file in uploaded_files])
            return self.service_result()

        with patch("inventory.views.process_sanger_files", side_effect=fake_process), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))):
            response = self.client.post(reverse("sanger_batch_upload"), files)

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Uploaded 3 AB1 files into 2 Sanger runs.")
        self.assertEqual(sorted(calls), [["one-forward.ab1", "one-reverse.ab1"], ["two-forward.ab1"]])
        self.assertEqual(SangerVerificationRun.objects.filter(plasmid=self.first_plasmid).count(), 1)
        self.assertEqual(SangerVerificationRun.objects.filter(plasmid=self.second_plasmid).count(), 1)

    def test_batch_upload_uses_primer_import_layout(self):
        response = self.client.get(reverse("sanger_batch_upload"))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'class="col-lg-7"')
        self.assertContains(response, 'class="col-lg-5"')
        self.assertContains(response, 'class="alert alert-light border"')
        self.assertContains(response, "ab1_file,plasmid_id,primer_id")
        self.assertContains(response, 'name="replace_existing"')
        self.assertNotContains(response, 'name="replace_existing" checked')

    def test_batch_upload_allows_empty_primer_id(self):
        with patch("inventory.views.process_sanger_files", return_value=self.service_result()), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))):
            response = self.client.post(
                reverse("sanger_batch_upload"),
                {
                    "ab1_files": [SimpleUploadedFile("trace.ab1", b"trace")],
                    "mapping_csv": SimpleUploadedFile(
                        "mapping.csv",
                        b"ab1_file,plasmid_id,primer_id\ntrace.ab1,601,\n",
                        content_type="text/csv",
                    ),
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Uploaded 1 AB1 files into 1 Sanger runs.")
        self.assertEqual(SangerVerificationRun.objects.filter(plasmid=self.first_plasmid).count(), 1)

    def test_batch_upload_allows_omitted_primer_column(self):
        with patch("inventory.views.process_sanger_files", return_value=self.service_result()), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))):
            response = self.client.post(
                reverse("sanger_batch_upload"),
                {
                    "ab1_files": [SimpleUploadedFile("trace.ab1", b"trace")],
                    "mapping_csv": SimpleUploadedFile(
                        "mapping.csv",
                        b"ab1_file,plasmid_id\ntrace.ab1,601\n",
                        content_type="text/csv",
                    ),
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Uploaded 1 AB1 files into 1 Sanger runs.")
        self.assertEqual(SangerVerificationRun.objects.filter(plasmid=self.first_plasmid).count(), 1)

    def test_batch_upload_normalizes_single_uploaded_file_result(self):
        uploaded = UploadedSangerFile(
            original_name="trace.ab1",
            data=b"trace",
            size=5,
            sha256="a" * 64,
            format="ab1",
            group_name="trace",
        )
        service_result = {
            "parameters": {},
            "uploaded_files": uploaded,
            "reads": [{
                "name": "trace",
                "files": uploaded,
                "alignment": {},
                "raw_sequence": "",
                "trimmed_sequence": "",
                "quality_metrics": {},
                "warnings": [],
                "is_usable": False,
            }],
            "combined": {"variants": []},
            "classification": {"state": "NO_DATA", "reasons": []},
        }

        with patch("inventory.views.process_sanger_files", return_value=service_result), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))):
            response = self.client.post(
                reverse("sanger_batch_upload"),
                {
                    "ab1_files": [SimpleUploadedFile("trace.ab1", b"trace")],
                    "mapping_csv": SimpleUploadedFile(
                        "mapping.csv",
                        b"ab1_file,plasmid_id\ntrace.ab1,601\n",
                        content_type="text/csv",
                    ),
                },
            )

        self.assertContains(response, "Uploaded 1 AB1 files into 1 Sanger runs.")
        self.assertEqual(SangerReadFile.objects.get(original_name="trace.ab1").size, 5)

    def test_batch_upload_skips_existing_file_and_uploads_the_rest(self):
        existing_run = SangerVerificationRun.objects.create(plasmid=self.first_plasmid, created_by=self.user)
        existing_read = SangerRead.objects.create(run=existing_run, name="already")
        SangerReadFile.objects.create(
            read=existing_read,
            format="ab1",
            original_name="already.ab1",
            sha256="b" * 64,
            size=5,
        )
        files = {
            "ab1_files": [
                SimpleUploadedFile("already.ab1", b"new-trace"),
                SimpleUploadedFile("new.ab1", b"new-trace"),
            ],
            "mapping_csv": SimpleUploadedFile(
                "mapping.csv",
                b"ab1_file,plasmid_id,primer_id\nalready.ab1,601,7000\nnew.ab1,601,7000\n",
                content_type="text/csv",
            ),
        }
        calls = []

        def fake_process(uploaded_files, reference_sequence):
            calls.append([file.name for file in uploaded_files])
            return self.service_result()

        with patch("inventory.views.process_sanger_files", side_effect=fake_process), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))):
            response = self.client.post(reverse("sanger_batch_upload"), files)

        self.assertContains(response, "Skipped existing AB1 file(s): already.ab1.")
        self.assertContains(response, "Uploaded 1 AB1 files into 1 Sanger runs (1 existing file(s) skipped).")
        self.assertEqual(calls, [["new.ab1"]])

    def test_batch_upload_can_process_existing_file_when_replacement_is_checked(self):
        existing_run = SangerVerificationRun.objects.create(plasmid=self.first_plasmid, created_by=self.user)
        existing_read = SangerRead.objects.create(run=existing_run, name="already")
        SangerReadFile.objects.create(
            read=existing_read,
            format="ab1",
            original_name="already.ab1",
            sha256="b" * 64,
            size=5,
        )
        files = {
            "ab1_files": [SimpleUploadedFile("already.ab1", b"replacement")],
            "mapping_csv": SimpleUploadedFile(
                "mapping.csv",
                b"ab1_file,plasmid_id,primer_id\nalready.ab1,601,7000\n",
                content_type="text/csv",
            ),
            "replace_existing": "on",
        }
        calls = []

        def fake_process(uploaded_files, reference_sequence):
            calls.append([file.name for file in uploaded_files])
            return self.service_result()

        with patch("inventory.views.process_sanger_files", side_effect=fake_process), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))):
            response = self.client.post(reverse("sanger_batch_upload"), files)

        self.assertContains(response, "Uploaded 1 AB1 files into 1 Sanger runs.")
        self.assertEqual(calls, [["already.ab1"]])

    def test_batch_upload_warns_and_skips_missing_plasmid_without_blocking_valid_files(self):
        files = {
            "ab1_files": [
                SimpleUploadedFile("deleted-plasmid.ab1", b"deleted"),
                SimpleUploadedFile("valid-plasmid.ab1", b"valid"),
            ],
            "mapping_csv": SimpleUploadedFile(
                "mapping.csv",
                b"ab1_file,plasmid_id,primer_id\ndeleted-plasmid.ab1,155,7000\nvalid-plasmid.ab1,601,7000\n",
                content_type="text/csv",
            ),
        }
        calls = []

        def fake_process(uploaded_files, reference_sequence):
            calls.append([file.name for file in uploaded_files])
            return self.service_result()

        with patch("inventory.views.process_sanger_files", side_effect=fake_process), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))):
            response = self.client.post(reverse("sanger_batch_upload"), files)

        self.assertContains(response, "Row 2 ignored: unavailable plasmid ID 155.")
        self.assertContains(response, "Uploaded 1 AB1 files into 1 Sanger runs.")
        self.assertEqual(calls, [["valid-plasmid.ab1"]])
        self.assertEqual(SangerVerificationRun.objects.filter(plasmid=self.first_plasmid).count(), 1)

    def test_batch_upload_rejects_unmapped_file_without_creating_runs(self):
        files = {
            "ab1_files": [SimpleUploadedFile("unmapped.ab1", b"unmapped")],
            "mapping_csv": SimpleUploadedFile(
                "mapping.csv",
                b"ab1_file,plasmid_id,primer_id\nother.ab1,601,7000\n",
                content_type="text/csv",
            ),
        }

        response = self.client.post(reverse("sanger_batch_upload"), files)

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "The CSV references other.ab1, but that file was not uploaded.")
        self.assertContains(response, "The uploaded AB1 file unmapped.ab1 is missing from the CSV.")
        self.assertEqual(SangerVerificationRun.objects.count(), 0)

    def test_batch_upload_links_ab1_to_primer_in_same_project(self):
        primer = Primer.objects.create(
            name="7001-Sequence-F",
            sequence_3="ACGTACGTACGT",
            fwd_or_rev="f",
        )
        files = {
            "ab1_files": [SimpleUploadedFile("primer-read.ab1", b"trace")],
            "mapping_csv": SimpleUploadedFile(
                "mapping.csv",
                b"ab1_file,plasmid_id,primer_id\nprimer-read.ab1,601,7001\n",
                content_type="text/csv",
            ),
        }

        def fake_process(uploaded_files, reference_sequence):
            uploaded = uploaded_files[0]
            return {
                "parameters": {},
                "reads": [{
                    "name": "primer-read",
                    "files": [SimpleNamespace(
                        format="ab1",
                        original_name=uploaded.name,
                        sha256="a" * 64,
                        size=5,
                        data=b"trace",
                    )],
                    "alignment": {},
                    "raw_sequence": "",
                    "trimmed_sequence": "",
                    "quality_metrics": {},
                    "warnings": [],
                    "is_usable": False,
                }],
                "combined": {"variants": []},
                "classification": {"state": "NO_DATA", "reasons": []},
            }

        with patch("inventory.views.process_sanger_files", side_effect=fake_process), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))):
            response = self.client.post(reverse("sanger_batch_upload"), files)

        self.assertContains(response, "Uploaded 1 AB1 files into 1 Sanger runs.")
        source_file = SangerReadFile.objects.get(original_name="primer-read.ab1")
        self.assertEqual(source_file.primer, primer)

    def test_batch_upload_links_global_primer(self):
        primer = Primer.objects.create(
            name="7002-Shared-Sequence-R",
            sequence_3="TGCATGCATGCA",
            fwd_or_rev="r",
        )
        files = {
            "ab1_files": [SimpleUploadedFile("shared-primer-read.ab1", b"trace")],
            "mapping_csv": SimpleUploadedFile(
                "mapping.csv",
                b"ab1_file,plasmid_id,primer_id\nshared-primer-read.ab1,601,7002\n",
                content_type="text/csv",
            ),
        }

        def fake_process(uploaded_files, reference_sequence):
            uploaded = uploaded_files[0]
            return {
                "parameters": {},
                "reads": [{
                    "name": "shared-primer-read",
                    "files": [SimpleNamespace(
                        format="ab1",
                        original_name=uploaded.name,
                        sha256="c" * 64,
                        size=5,
                        data=b"trace",
                    )],
                    "alignment": {},
                    "raw_sequence": "",
                    "trimmed_sequence": "",
                    "quality_metrics": {},
                    "warnings": [],
                    "is_usable": False,
                }],
                "combined": {"variants": []},
                "classification": {"state": "NO_DATA", "reasons": []},
            }

        with patch("inventory.views.process_sanger_files", side_effect=fake_process), \
                patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))):
            response = self.client.post(reverse("sanger_batch_upload"), files)

        self.assertContains(response, "Uploaded 1 AB1 files into 1 Sanger runs.")
        source_file = SangerReadFile.objects.get(original_name="shared-primer-read.ab1")
        self.assertEqual(source_file.primer, primer)


class SangerFailedReadReviewTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="failed-read-user", password="pw")
        self.project = Project.objects.create(name="Failed read project", public=False)
        Membership.objects.create(member=self.user, project=self.project, access_policies="r")
        self.plasmid = Plasmid.objects.create(
            idx=504,
            name="Failed read plasmid",
            intended_use="test",
            project=self.project,
            ligation_state=1,
        )
        record = SeqRecord(Seq("ACGT" * 40), id="failed-read-reference", name="failed-read-reference", description=".")
        record.annotations["molecule_type"] = "DNA"
        self.plasmid.sequence.save("failed-read-reference.gb", ContentFile(record.format("genbank")), save=True)
        self.run = SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user)
        chromatogram = {
            "aTrace": [10, 12, 15, 11] * 20,
            "cTrace": [8, 9, 10, 7] * 20,
            "gTrace": [7, 8, 9, 6] * 20,
            "tTrace": [6, 7, 8, 5] * 20,
            "basePos": list(range(0, 80)),
            "baseCalls": list("A" * 80),
            "qualNums": [10] * 80,
        }
        self.read = SangerRead.objects.create(
            run=self.run,
            name="poor-quality-read",
            selected_source="ab1",
            parsing_result={"chromatogram": chromatogram},
            is_usable=False,
        )
        self.client.force_login(self.user)

    def test_saved_run_detail_shows_and_updates_primer_assignment(self):
        assigned_primer = Primer.objects.create(
            name="7101-Assigned-Seq-F",
            sequence_3="ACGTACGTACGT",
            fwd_or_rev="f",
        )
        replacement_primer = Primer.objects.create(
            name="7102-Replacement-Seq-F",
            sequence_3="TTTTCCCCAAAA",
            fwd_or_rev="f",
        )
        source_file = SangerReadFile.objects.create(
            read=self.read,
            primer=assigned_primer,
            format="ab1",
            original_name="poor-quality-read.ab1",
            sha256="a" * 64,
            size=7,
        )
        membership = Membership.objects.get(member=self.user, project=self.project)
        membership.access_policies = "w"
        membership.save(update_fields=["access_policies"])

        response = self.client.get(reverse("sanger_run_detail", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
        }))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Primer")
        self.assertContains(response, "Return to run tables")
        self.assertContains(response, reverse("services-sanger"))
        self.assertContains(response, assigned_primer.name)
        self.assertContains(response, "No primer")

        response = self.client.post(reverse("sanger_read_file_primer_update", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
            "read_id": self.read.id,
            "file_id": source_file.id,
        }), {"primer_id": str(replacement_primer.id)})

        self.assertRedirects(response, reverse("sanger_run_detail", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
        }) + "?primer_updated=1", fetch_redirect_response=False)
        source_file.refresh_from_db()
        self.assertEqual(source_file.primer, replacement_primer)

        self.client.post(reverse("sanger_read_file_primer_update", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
            "read_id": self.read.id,
            "file_id": source_file.id,
        }), {"primer_id": ""})
        source_file.refresh_from_db()
        self.assertIsNone(source_file.primer)

    def test_failed_read_offers_chromatogram_and_local_blast_separately(self):
        response = self.client.get(reverse("sanger_run_detail", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
        }))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "poor-quality-read did not produce a reliable alignment.")
        self.assertContains(response, "Open chromatogram window")
        self.assertContains(response, "Local BLAST")
        self.assertContains(response, reverse("services-blast"))
        self.assertContains(response, 'const plasmidHeader = document.getElementById("header-container")')

    def test_mixed_run_hides_failed_read_from_active_review(self):
        SangerRead.objects.create(
            run=self.run,
            name="usable-forward-read",
            detected_orientation="forward",
            alignment_metrics={
                "orientation": "forward",
                "covered_positions": [0, 1, 2],
                "identity": 100.0,
                "aligned_length": 3,
                "variants": [],
            },
            is_usable=True,
        )

        response = self.client.get(reverse("sanger_run_detail", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
        }))

        self.assertEqual(response.status_code, 200)
        self.assertNotContains(response, "poor-quality-read did not produce a reliable alignment.")

    def test_saved_run_detail_handles_legacy_alignment_without_covered_positions(self):
        SangerRead.objects.create(
            run=self.run,
            name="legacy-aligned-read",
            detected_orientation="forward",
            alignment_metrics={
                "orientation": "forward",
                "identity": 100.0,
                "aligned_length": 3,
                "reference_projection_base_indices": [0, 1, 2],
                "variants": [],
            },
            is_usable=True,
        )

        response = self.client.get(reverse("sanger_run_detail", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
        }))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "legacy-aligned-read")

    def test_failed_read_chromatogram_can_be_opened_without_alignment(self):
        response = self.client.get(reverse("sanger_read_chromatogram", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
            "read_id": self.read.id,
        }))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Sanger chromatogram")
        self.assertContains(response, 'const plasmidHeader = document.getElementById("header-container")')

    def test_short_read_warning_uses_icon_with_tooltip(self):
        self.read.warnings = ["The usable sequence is too short for reliable alignment."]
        self.read.save(update_fields=["warnings"])

        response = self.client.get(reverse("sanger_run_detail", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
        }))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "bi-exclamation-triangle-fill")
        self.assertContains(response, 'data-bs-title="The usable sequence is too short for reliable alignment."')
        self.assertNotContains(response, "trimmed sequence shorter than minimum")

    def test_saved_run_detail_hides_previous_runs(self):
        SangerVerificationRun.objects.create(
            plasmid=self.plasmid,
            created_by=self.user,
            label="previous run",
        )

        response = self.client.get(reverse("sanger_run_detail", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
        }))

        self.assertEqual(response.status_code, 200)
        self.assertNotContains(response, "Saved Sanger runs")
        self.assertContains(response, "Reads")
        self.assertContains(response, "sanger-inline-autoadjust")
        self.assertContains(response, "Auto-adjust")
        self.assertContains(response, "enableHoldToRepeat")
        self.assertContains(response, "inlineDragging")
        self.assertContains(response, "isSkippedReferenceGap")
        self.assertContains(response, "ctx.lineTo(gapStartX, baseline)")


class PlasmidValidationSequencingFilesTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="validation-editor", password="pw")
        self.project = Project.objects.create(name="Validation project", public=False)
        Membership.objects.create(member=self.user, project=self.project, access_policies="w")
        self.plasmid = Plasmid.objects.create(
            idx=503,
            name="Validation plasmid",
            intended_use="test",
            project=self.project,
            ligation_state=1,
        )
        self.run = SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user)
        self.read = SangerRead.objects.create(run=self.run, name="forward")
        self.source_file = SangerReadFile.objects.create(
            read=self.read,
            format="ab1",
            original_name="forward.ab1",
            sha256="b" * 64,
            size=7,
        )
        self.source_file.file.save("forward.ab1", ContentFile(b"AB1DATA"), save=True)
        self.client.force_login(self.user)
        self.client.cookies["current_project_id"] = str(self.project.id)

    def test_validation_can_mark_construct_as_having_no_colony(self):
        response = self.client.post(
            reverse("plasmid_validation_edit", kwargs={"plasmid_id": self.plasmid.id}),
            {
                "ligation_state": "1",
                "working_colony": "27",
                "no_colony": "on",
                "colonypcr_state": "0",
                "digestion_state": "0",
                "sequencing_state": "0",
            },
        )

        self.assertEqual(response.status_code, 302)
        self.plasmid.refresh_from_db()
        self.assertTrue(self.plasmid.no_colony)
        self.assertIsNone(self.plasmid.working_colony)
        self.assertEqual(self.plasmid.colony_source_text(), "No colony")
        self.assertEqual(self.plasmid.working_colony_text_short(), "NC")

    def test_validation_link_accepts_no_colony_token(self):
        initial, error = plasmid_validation_initial_from_payload(
            "503_none_2026-08-20_pcr"
        )

        self.assertIsNone(error)
        self.assertIsNone(initial["working_colony"])
        self.assertTrue(initial["no_colony"])

        response = self.client.get(
            reverse(
                "plasmid_validation_from_link",
                kwargs={"validation_payload": "503_none_2026-08-20_pcr"},
            )
        )

        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.context["form"].initial["no_colony"])

    def test_numeric_validation_link_keeps_colony_behavior(self):
        initial, error = plasmid_validation_initial_from_payload(
            "503_7_2026-08-20_pcr"
        )

        self.assertIsNone(error)
        self.assertEqual(initial["working_colony"], 7)
        self.assertFalse(initial["no_colony"])

    def test_validation_form_lists_sequencing_file_and_review_links(self):
        response = self.client.get(reverse("plasmid_validation_edit", kwargs={"plasmid_id": self.plasmid.id}))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Sequencing Files")
        self.assertNotContains(response, "Sequencing clustal file")
        self.assertContains(response, "forward.ab1")
        self.assertContains(response, "Download")
        self.assertContains(response, "Review trace")
        self.assertContains(response, "Review alignment")

    def test_sequencing_file_download_is_scoped_to_its_plasmid_run_and_read(self):
        response = self.client.get(reverse("sanger_read_file_download", kwargs={
            "plasmid_id": self.plasmid.id,
            "run_id": self.run.id,
            "read_id": self.read.id,
            "file_id": self.source_file.id,
        }))

        self.assertEqual(response.status_code, 200)
        self.assertIn("forward.ab1", response["Content-Disposition"])
        self.assertEqual(b"".join(response.streaming_content), b"AB1DATA")


class SangerDecisionTests(TestCase):
    def setUp(self):
        self.editor = User.objects.create_user(username="seq-editor", password="pw")
        self.reader = User.objects.create_user(username="seq-reader", password="pw")
        self.project = Project.objects.create(name="DecisionProject", public=False)
        Membership.objects.create(member=self.editor, project=self.project, access_policies="w")
        Membership.objects.create(member=self.reader, project=self.project, access_policies="r")
        self.plasmid = Plasmid.objects.create(
            idx=501,
            name="Decision plasmid",
            intended_use="test",
            project=self.project,
            ligation_state=1,
            colonypcr_state=0,
            digestion_state=0,
            sequencing_state=1,
        )
        self.run = SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.editor)

    def test_verified_decision_updates_plasmid_and_evidence(self):
        self.client.force_login(self.editor)

        response = self.client.post(reverse("sanger_run_decision", kwargs={"plasmid_id": self.plasmid.id, "run_id": self.run.id}), {
            "manual_decision": "VERIFIED",
            "manual_decision_effective_date": "2026-07-30",
            "manual_decision_comment": "Clean overlapping reads.",
        })

        self.assertEqual(response.status_code, 302)
        self.run.refresh_from_db()
        self.plasmid.refresh_from_db()
        self.assertEqual(self.run.manual_decision, "VERIFIED")
        self.assertEqual(self.run.manual_decision_by, self.editor)
        self.assertEqual(str(self.run.manual_decision_effective_date), "2026-07-30")
        self.assertEqual(self.run.manual_decision_comment, "Clean overlapping reads.")
        self.assertEqual(self.plasmid.sequencing_state, 2)
        self.assertEqual(str(self.plasmid.sequencing_date), "2026-07-30")

    def test_inconclusive_decision_does_not_mark_plasmid_verified(self):
        self.client.force_login(self.editor)

        self.client.post(reverse("sanger_run_decision", kwargs={"plasmid_id": self.plasmid.id, "run_id": self.run.id}), {
            "manual_decision": "INCONCLUSIVE",
            "manual_decision_effective_date": "2026-07-31",
            "manual_decision_comment": "Coverage gap across junction.",
        })

        self.run.refresh_from_db()
        self.plasmid.refresh_from_db()
        self.assertEqual(self.run.manual_decision, "INCONCLUSIVE")
        self.assertEqual(self.plasmid.sequencing_state, 1)
        self.assertIsNone(self.plasmid.sequencing_date)

    def test_read_only_user_cannot_save_decision(self):
        self.client.force_login(self.reader)

        response = self.client.post(reverse("sanger_run_decision", kwargs={"plasmid_id": self.plasmid.id, "run_id": self.run.id}), {
            "manual_decision": "VERIFIED",
        })

        self.assertNotEqual(response.status_code, 302)
        self.run.refresh_from_db()
        self.plasmid.refresh_from_db()
        self.assertEqual(self.run.manual_decision, "")
        self.assertEqual(self.plasmid.sequencing_state, 1)


class SangerServicesListTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="sanger-services-user", password="pw")
        self.visible_project = Project.objects.create(name="Visible sequencing project", public=False)
        self.hidden_project = Project.objects.create(name="Hidden sequencing project", public=False)
        Membership.objects.create(member=self.user, project=self.visible_project, access_policies="r")
        self.visible_plasmid = Plasmid.objects.create(
            idx=505,
            name="Visible Sanger plasmid",
            intended_use="test",
            project=self.visible_project,
        )
        self.hidden_plasmid = Plasmid.objects.create(
            idx=506,
            name="Hidden Sanger plasmid",
            intended_use="test",
            project=self.hidden_project,
        )
        self.approved_run = SangerVerificationRun.objects.create(
            plasmid=self.visible_plasmid,
            created_by=self.user,
            automated_state="PASS",
            manual_decision="VERIFIED",
            manual_decision_by=self.user,
        )
        self.primer = Primer.objects.create(
            name="8001-Visible-Seq-F",
            sequence_3="ACGTACGTACGT",
            fwd_or_rev="f",
        )
        read = SangerRead.objects.create(
            run=self.approved_run,
            name="forward",
            detected_orientation="forward",
            is_usable=True,
        )
        self.read = read
        source_file = SangerReadFile.objects.create(
            read=read,
            primer=self.primer,
            format="ab1",
            original_name="visible-forward.ab1",
            sha256="c" * 64,
            size=123,
            metadata={"run_start_date": "2025-10-15", "run_start_time": "14:54:11"},
        )
        source_file.file.save("visible-forward.ab1", ContentFile(b"AB1DATA"), save=True)
        self.source_file = source_file
        reverse_read = SangerRead.objects.create(
            run=self.approved_run,
            name="reverse",
            detected_orientation="reverse",
            is_usable=True,
        )
        reverse_file = SangerReadFile.objects.create(
            read=reverse_read,
            format="ab1",
            original_name="visible-reverse.ab1",
            sha256="d" * 64,
            size=123,
        )
        reverse_file.file.save("visible-reverse.ab1", ContentFile(b"AB1DATA"), save=True)
        invalid_read = SangerRead.objects.create(
            run=self.approved_run,
            name="poor-quality-503",
            detected_orientation="forward",
            is_usable=False,
        )
        invalid_file = SangerReadFile.objects.create(
            read=invalid_read,
            format="ab1",
            original_name="poor-quality-503.ab1",
            sha256="e" * 64,
            size=123,
        )
        invalid_file.file.save("poor-quality-503.ab1", ContentFile(b"AB1DATA"), save=True)
        SangerVerificationRun.objects.create(
            plasmid=self.visible_plasmid,
            created_by=self.user,
            automated_state="FAIL",
            manual_decision="REJECTED",
        )
        SangerVerificationRun.objects.create(
            plasmid=self.hidden_plasmid,
            automated_state="PASS",
            manual_decision="VERIFIED",
        )
        self.client.force_login(self.user)

    def test_lists_accessible_runs_with_files_and_approval_status(self):
        response = self.client.get(reverse("services-sanger"))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Services")
        self.assertContains(response, "Sanger sequencing")
        self.assertNotContains(response, "Sanger sequencing runs")
        self.assertContains(response, 'id="table_search-input"')
        self.assertContains(response, "Search Sanger files ...")
        self.assertContains(response, "Visible Sanger plasmid")
        self.assertContains(response, "visible-forward.ab1")
        self.assertContains(response, "15 Oct 2025")
        self.assertContains(response, 'title="Date uploaded:')
        self.assertContains(response, "8001-Visible-Seq-F")
        self.assertContains(response, reverse("primer", kwargs={"primer_id": self.primer.id}))
        self.assertContains(response, "visible-reverse.ab1")
        self.assertContains(response, 'data-sanger-file-name="visible-forward.ab1"')
        self.assertContains(response, "FWD")
        self.assertContains(response, "REV")
        self.assertNotContains(response, "UNK")
        self.assertNotContains(response, 'class="btn btn-secondary sanger-read-direction"')
        self.assertNotContains(response, "poor-quality-503.ab1")
        self.assertContains(response, "Download original file: visible-forward.ab1")
        self.assertContains(response, 'download="visible-forward.ab1"')
        self.assertContains(response, reverse("sanger_read_file_download", kwargs={
            "plasmid_id": self.approved_run.plasmid_id,
            "run_id": self.approved_run.id,
            "read_id": self.read.id,
            "file_id": self.source_file.id,
        }))
        self.assertContains(response, "Approved")
        self.assertContains(response, "Authorized by sanger-services-user")
        self.assertContains(response, "All runs (2)")
        self.assertNotContains(response, 'aria-label="Select run for Visible Sanger plasmid"')
        self.assertNotContains(response, "Not approved")
        self.assertNotContains(response, "Hidden Sanger plasmid")
        self.assertContains(response, reverse("sanger_run_detail", kwargs={
            "plasmid_id": self.approved_run.plasmid_id,
            "run_id": self.approved_run.id,
        }))

    @patch("inventory.views.parse_ab1", side_effect=AssertionError("AB1 parsing must not run while listing saved runs"))
    def test_services_list_uses_persisted_metadata_without_reparsing_ab1(self, parse_ab1):
        response = self.client.get(reverse("services-sanger"))

        self.assertEqual(response.status_code, 200)
        parse_ab1.assert_not_called()

    def test_services_list_falls_back_to_newest_uploaded_run(self):
        fallback_plasmid = Plasmid.objects.create(
            idx=507,
            name="Fallback Sanger plasmid",
            intended_use="test",
            project=self.visible_project,
        )
        SangerVerificationRun.objects.create(
            plasmid=fallback_plasmid,
            created_by=self.user,
            label="older",
        )
        SangerVerificationRun.objects.create(
            plasmid=fallback_plasmid,
            created_by=self.user,
            label="newest",
        )

        response = self.client.get(reverse("services-sanger"))

        self.assertContains(response, "Fallback Sanger plasmid")
        self.assertContains(response, "newest")
        newest_url = reverse("sanger_run_detail", kwargs={
            "plasmid_id": fallback_plasmid.id,
            "run_id": SangerVerificationRun.objects.get(plasmid=fallback_plasmid, label="newest").id,
        })
        self.assertContains(response, newest_url)
        self.assertContains(response, "Selected: latest uploaded")
        self.assertContains(response, "All runs (2)")

    def test_reader_does_not_see_batch_upload_button(self):
        response = self.client.get(reverse("services-sanger"))

        self.assertNotContains(response, reverse("sanger_batch_upload"))
        self.assertNotContains(response, 'title="Re-run all Sanger analyses"')

    def test_uploader_sees_icon_only_batch_upload_control(self):
        membership = Membership.objects.get(member=self.user, project=self.visible_project)
        membership.access_policies = "w"
        membership.save(update_fields=["access_policies"])

        response = self.client.get(reverse("services-sanger"))

        self.assertContains(response, 'title="Upload AB1 batch"')
        self.assertContains(response, 'aria-label="Upload AB1 batch"')
        self.assertContains(response, 'title="Re-run all Sanger analyses"')
        self.assertContains(response, 'title="Refresh"')
        self.assertContains(response, 'class="btn btn-secondary sanger-header-action-button"')
        self.assertContains(response, 'class="btn btn-outline-primary sanger-batch-upload-button sanger-header-action-button"')
        self.assertContains(response, 'sanger-header-action-button')
        self.assertNotContains(response, "plasmids</span>")

    def test_writer_can_request_sanger_rerun_progress_manifest(self):
        membership = Membership.objects.get(member=self.user, project=self.visible_project)
        membership.access_policies = "w"
        membership.save(update_fields=["access_policies"])

        response = self.client.post(
            reverse("sanger_rerun_all"),
            {"progress": "1"},
            HTTP_ACCEPT="application/json",
        )

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload["total"], 2)
        self.assertEqual(set(payload["run_ids"]), {
            str(self.approved_run.id),
            str(SangerVerificationRun.objects.filter(
                plasmid=self.visible_plasmid,
            ).exclude(id=self.approved_run.id).get().id),
        })

    @patch("inventory.views.process_sanger_files")
    def test_writer_can_rerun_all_accessible_sanger_runs(self, process_sanger_files):
        membership = Membership.objects.get(member=self.user, project=self.visible_project)
        membership.access_policies = "w"
        membership.save(update_fields=["access_policies"])
        process_sanger_files.return_value = {
            "parameters": {},
            "uploaded_files": [],
            "reads": [{
                "name": "forward",
                "formats": ["ab1"],
                "selected_source": "ab1",
                "files": [UploadedSangerFile(
                    original_name="visible-forward.ab1",
                    data=b"AB1DATA",
                    size=7,
                    sha256="c" * 64,
                    format="ab1",
                    group_name="forward",
                )],
                "alignment": {},
                "display_alignment": {},
                "raw_sequence": "",
                "trimmed_sequence": "",
                "quality_metrics": {},
                "chromatogram": {},
                "warnings": [],
                "errors": [],
                "is_usable": False,
            }],
            "combined": {"variants": []},
            "classification": {"state": "NO_DATA", "reasons": []},
        }
        with patch("inventory.views.grab_seq", return_value=(True, Seq("ACGT" * 20))):
            response = self.client.post(reverse("sanger_rerun_all"))

        self.assertEqual(response.status_code, 302)
        self.assertIn("rerun_all=1", response.url)
        self.assertIn("processed=1", response.url)
        self.assertIn("skipped=1", response.url)
        self.assertIn("failed=0", response.url)
        self.assertTrue(SangerVerificationRun.objects.filter(id=self.approved_run.id).exists())
        self.assertEqual(SangerVerificationRun.objects.filter(plasmid=self.visible_plasmid).count(), 2)
        self.assertTrue(SangerReadFile.objects.filter(
            read__run_id=self.approved_run.id,
            original_name="visible-forward.ab1",
        ).exists())

    def test_services_list_shows_only_top_50_plasmids_by_id(self):
        for idx in range(600, 655):
            plasmid = Plasmid.objects.create(
                idx=idx,
                name="Sanger limit {}".format(idx),
                intended_use="test",
                project=self.visible_project,
            )
            SangerVerificationRun.objects.create(plasmid=plasmid, created_by=self.user)

        response = self.client.get(reverse("services-sanger"))

        self.assertEqual(response.content.count(b'class="sanger-services-row"'), 50)
        self.assertContains(response, "Sanger limit 654")
        self.assertNotContains(response, "Sanger limit 604")
        self.assertContains(response, 'aria-label="Next page"')
        self.assertContains(response, 'aria-label="Last page"')
        self.assertNotContains(response, 'aria-label="Previous page"')
        self.assertNotContains(response, 'aria-label="First page"')
        self.assertEqual(response.content.count(b'aria-label="Next page"'), 2)
        self.assertEqual(response.content.count(b'aria-label="Last page"'), 2)

        second_page = self.client.get(reverse("services-sanger") + "?page=2")

        self.assertEqual(second_page.content.count(b'class="sanger-services-row"'), 6)
        self.assertContains(second_page, "Sanger limit 604")
        self.assertNotContains(second_page, "Sanger limit 654")
        self.assertContains(second_page, 'aria-label="Previous page"')
        self.assertContains(second_page, 'aria-label="First page"')
        self.assertNotContains(second_page, 'aria-label="Next page"')
        self.assertNotContains(second_page, 'aria-label="Last page"')
        self.assertEqual(second_page.content.count(b'aria-label="Previous page"'), 2)
        self.assertEqual(second_page.content.count(b'aria-label="First page"'), 2)

    def test_user_without_project_membership_sees_no_runs(self):
        outsider = User.objects.create_user(username="sanger-services-outsider", password="pw")
        self.client.force_login(outsider)

        response = self.client.get(reverse("services-sanger"))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "No Sanger sequencing runs are available")
        self.assertNotContains(response, "Visible Sanger plasmid")


class SangerAb1MetadataTests(SimpleTestCase):
    def test_formats_ab1_run_start_datetime(self):
        self.assertEqual(
            ab1_run_datetime_label({"run_start_date": "2025-10-15", "run_start_time": "14:54:11"}),
            "15 Oct 2025 14:54",
        )

    def test_returns_empty_label_when_ab1_has_no_run_datetime(self):
        self.assertEqual(ab1_run_datetime_label({}), "")

    def test_extracts_run_datetime_from_filename(self):
        filename = "6098-Diego-503_Hr1_Chr4_Dw_c4-457-pL0-chech-F_2026-07-07-13-15-37_E02.ab1"
        self.assertEqual(filename_run_datetime(filename), datetime.datetime(2026, 7, 7, 13, 15, 37))
        self.assertEqual(ab1_run_datetime_label({}, filename), "07 Jul 2026 13:15")

    def test_prefers_filename_datetime_over_ab1_metadata(self):
        self.assertEqual(
            ab1_run_datetime_label(
                {"run_start_date": "2025-10-15", "run_start_time": "14:54:11"},
                "trace_2026-07-07-13-15-37.ab1",
            ),
            "07 Jul 2026 13:15",
        )

    def test_prefers_latest_forward_when_no_reverse_exists(self):
        old_forward = SimpleNamespace(name="4998-old-forward")
        latest_forward = SimpleNamespace(name="6098-latest-forward")
        selected = latest_confirmation_pair([
            {
                "read": old_forward,
                "orientation": "forward",
                "run_date": datetime.date(2026, 5, 4),
                "run_datetime": datetime.datetime(2026, 5, 4, 18, 35, 26),
                "registration_number": 4998,
            },
            {
                "read": latest_forward,
                "orientation": "forward",
                "run_date": datetime.date(2026, 7, 7),
                "run_datetime": datetime.datetime(2026, 7, 7, 13, 15, 37),
                "registration_number": 6098,
            },
        ])

        self.assertEqual(selected, [latest_forward])

    def test_prefers_latest_complete_pair_and_highest_close_registration_numbers(self):
        old_forward = SimpleNamespace(name="100-old-F")
        old_reverse = SimpleNamespace(name="101-old-R")
        latest_forward = SimpleNamespace(name="16837-latest-F")
        latest_forward_older = SimpleNamespace(name="16835-latest-F")
        latest_reverse = SimpleNamespace(name="16838-latest-R")
        latest_reverse_older = SimpleNamespace(name="16836-latest-R")
        selected = latest_confirmation_pair([
            {"read": old_forward, "orientation": "forward", "run_date": datetime.date(2024, 8, 29), "registration_number": 100},
            {"read": old_reverse, "orientation": "reverse", "run_date": datetime.date(2024, 8, 29), "registration_number": 101},
            {"read": latest_forward_older, "orientation": "forward", "run_date": datetime.date(2024, 9, 5), "registration_number": 16835},
            {"read": latest_reverse_older, "orientation": "reverse", "run_date": datetime.date(2024, 9, 5), "registration_number": 16836},
            {"read": latest_forward, "orientation": "forward", "run_date": datetime.date(2024, 9, 5), "registration_number": 16837},
            {"read": latest_reverse, "orientation": "reverse", "run_date": datetime.date(2024, 9, 5), "registration_number": 16838},
        ])

        self.assertEqual(selected, [latest_forward, latest_reverse])

    def test_omits_unknown_reads_when_latest_complete_pair_exists(self):
        old_forward = SimpleNamespace(name="100-old-F")
        old_reverse = SimpleNamespace(name="101-old-R")
        old_unknown = SimpleNamespace(name="102-old-unknown")
        latest_forward = SimpleNamespace(name="200-latest-F")
        latest_reverse = SimpleNamespace(name="201-latest-R")
        latest_unknown = SimpleNamespace(name="202-latest-unknown")

        selected = latest_confirmation_pair([
            {"read": old_forward, "orientation": "forward", "run_date": datetime.date(2024, 8, 29), "registration_number": 100},
            {"read": old_reverse, "orientation": "reverse", "run_date": datetime.date(2024, 8, 29), "registration_number": 101},
            {"read": old_unknown, "orientation": "unknown", "run_date": datetime.date(2024, 8, 29), "registration_number": 102},
            {"read": latest_forward, "orientation": "forward", "run_date": datetime.date(2024, 9, 5), "registration_number": 200},
            {"read": latest_reverse, "orientation": "reverse", "run_date": datetime.date(2024, 9, 5), "registration_number": 201},
            {"read": latest_unknown, "orientation": "unknown", "run_date": datetime.date(2024, 9, 5), "registration_number": 202},
        ])

        self.assertEqual(selected, [latest_forward, latest_reverse])

    def test_omits_unknown_reads_when_only_one_orientation_exists(self):
        latest_forward = SimpleNamespace(name="latest-forward")
        old_unknown = SimpleNamespace(name="old-unknown")
        latest_unknown = SimpleNamespace(name="latest-unknown")

        selected = latest_confirmation_pair([
            {
                "read": old_unknown,
                "orientation": "unknown",
                "run_date": datetime.date(2024, 8, 29),
                "registration_number": 100,
            },
            {
                "read": latest_forward,
                "orientation": "forward",
                "run_date": datetime.date(2024, 9, 5),
                "run_datetime": datetime.datetime(2024, 9, 5, 12, 0),
                "registration_number": 200,
            },
            {
                "read": latest_unknown,
                "orientation": "unknown",
                "run_date": datetime.date(2024, 9, 5),
                "registration_number": 201,
            },
        ])

        self.assertEqual(selected, [latest_forward])

    def test_review_uses_latest_complete_date_group_before_region_selection(self):
        old_forward = {"name": "old-forward", "is_usable": True, "alignment": {"orientation": "forward"}}
        latest_forward = {"name": "latest-forward", "is_usable": True, "alignment": {"orientation": "forward"}}
        latest_reverse = {"name": "latest-reverse", "is_usable": True, "alignment": {"orientation": "reverse"}}

        selected = select_sanger_review_candidates([
            {"read": old_forward, "orientation": "forward", "covered_positions": [1, 2], "run_date": datetime.date(2024, 8, 29)},
            {"read": latest_forward, "orientation": "forward", "covered_positions": [3, 4], "run_date": datetime.date(2024, 9, 5)},
            {"read": latest_reverse, "orientation": "reverse", "covered_positions": [5, 6], "run_date": datetime.date(2024, 9, 5)},
        ])

        self.assertEqual([item["read"] for item in selected], [latest_forward, latest_reverse])

    def test_groups_failed_reads_with_consecutive_registration_numbers(self):
        first = SimpleNamespace(id="first", files=SimpleNamespace(all=lambda: [SimpleNamespace(
            original_name="36720-read_2025-10-15-12-34-25.ab1", metadata={}
        )]))
        second = SimpleNamespace(id="second", files=SimpleNamespace(all=lambda: [SimpleNamespace(
            original_name="36721-read_2025-10-15-12-34-25.ab1", metadata={}
        )]))

        groups = sanger_failed_read_groups([first, second])

        self.assertEqual([[read.id for read in group["reads"]] for group in groups], [["second", "first"]])

    def test_collapses_overlapping_reads_per_orientation_to_latest(self):
        old_forward = SimpleNamespace(name="old-forward")
        latest_forward = SimpleNamespace(name="latest-forward")
        selected = select_nonredundant_read_candidates([
            {
                "read": old_forward,
                "orientation": "forward",
                "covered_positions": list(range(100, 200)),
                "run_datetime": datetime.datetime(2026, 5, 4, 18, 35, 26),
                "registration_number": 4998,
            },
            {
                "read": latest_forward,
                "orientation": "forward",
                "covered_positions": list(range(110, 210)),
                "run_datetime": datetime.datetime(2026, 7, 7, 13, 15, 37),
                "registration_number": 6098,
            },
        ])

        self.assertEqual([item["read"] for item in selected], [latest_forward])

    def test_keeps_same_orientation_reads_covering_distinct_regions(self):
        first_forward = SimpleNamespace(name="first-forward")
        second_forward = SimpleNamespace(name="second-forward")
        selected = select_nonredundant_read_candidates([
            {
                "read": first_forward,
                "orientation": "forward",
                "covered_positions": list(range(100, 200)),
                "run_datetime": datetime.datetime(2026, 7, 7, 13, 15, 37),
            },
            {
                "read": second_forward,
                "orientation": "forward",
                "covered_positions": list(range(500, 600)),
                "run_datetime": datetime.datetime(2026, 7, 7, 13, 20, 37),
            },
        ])

        self.assertEqual([item["read"] for item in selected], [first_forward, second_forward])

    def test_saved_run_review_omits_failed_reads_when_usable_evidence_exists(self):
        failed = {"name": "old-failed", "is_usable": False, "alignment": {}}
        latest = {"name": "latest-forward", "is_usable": True, "alignment": {"orientation": "forward"}}

        selected = select_sanger_review_candidates([
            {"read": failed, "orientation": "unknown", "covered_positions": []},
            {"read": latest, "orientation": "forward", "covered_positions": [100, 101]},
        ])

        self.assertEqual([item["read"] for item in selected], [latest])

    def test_prefers_manual_acceptance_then_automatic_pass(self):
        runs = [
            SimpleNamespace(id="accepted", created_at=datetime.datetime(2026, 1, 1), automated_state="PASS", manual_decision="VERIFIED"),
            SimpleNamespace(id="automatic", created_at=datetime.datetime(2026, 1, 2), automated_state="PASS", manual_decision=""),
            SimpleNamespace(id="review", created_at=datetime.datetime(2026, 1, 4), automated_state="REVIEW", manual_decision=""),
            SimpleNamespace(id="latest", created_at=datetime.datetime(2026, 1, 5), automated_state="FAIL", manual_decision=""),
        ]

        self.assertEqual(preferred_sanger_run(runs).id, "accepted")

        runs[0].manual_decision = ""
        self.assertEqual(preferred_sanger_run(runs).id, "automatic")

        runs[0].automated_state = "FAIL"
        runs[1].automated_state = "FAIL"
        self.assertEqual(preferred_sanger_run(runs).id, "review")

        runs[2].automated_state = "FAIL"
        self.assertEqual(preferred_sanger_run(runs).id, "latest")


class PrimerDimerAnalysisTests(SimpleTestCase):
    def test_normalizes_lowercase_and_spaces(self):
        primers, warnings = validate_primers([
            PrimerInput("p1", "aa aa cccc\n"),
        ])

        self.assertEqual(warnings, [])
        self.assertEqual(primers[0].sequence, "AAAACCCC")

    def test_rejects_invalid_sequence(self):
        with self.assertRaises(PrimerDimerInputError):
            validate_primers([PrimerInput("p1", "AAAN")])

    def test_rejects_duplicate_names(self):
        with self.assertRaises(PrimerDimerInputError):
            validate_primers([PrimerInput("p1", "AAAA"), PrimerInput("p1", "CCCC")])

    def test_reports_identical_sequences_with_different_names(self):
        primers, warnings = validate_primers([
            PrimerInput("p1", "AAAA"),
            PrimerInput("p2", "AAAA"),
        ])

        self.assertEqual(len(primers), 2)
        self.assertEqual(len(warnings), 1)

    def test_multiplex_combination_count(self):
        primers = [PrimerInput("p1", "AAAA"), PrimerInput("p2", "CCCC"), PrimerInput("p3", "TTTT")]

        combinations = list(primer_combinations(primers))

        self.assertEqual(len(combinations), 6)

    def test_reads_fasta(self):
        primers = read_primers(StringIO(">p1\nAAAA\n>p2\nCCCC\n"), "fasta")

        self.assertEqual([primer.name for primer in primers], ["p1", "p2"])

    def test_reads_csv(self):
        primers = read_primers(StringIO("name,sequence\np1,AAAA\np2,CCCC\n"), "csv")

        self.assertEqual([primer.sequence for primer in primers], ["AAAA", "CCCC"])

    def test_high_for_both_3prime_extendable_heterodimer(self):
        conditions = PrimerDimerConditions(annealing_temp_c=55)
        thresholds = PrimerDimerThresholds(high_delta_g_kcal_mol=-999, moderate_delta_g_kcal_mol=-999)

        result = analyze_pair(
            PrimerInput("a", "AAAACCCC"),
            PrimerInput("b", "GGGGTTTT"),
            "heterodimer",
            conditions,
            thresholds,
        )

        self.assertEqual(result["risk"], "HIGH")
        self.assertTrue(result["extendable_3prime"])
        self.assertEqual(result["terminal_3prime_run"], 8)
        self.assertLess(result["delta_g_kcal_mol"], 0)

    def test_four_base_3prime_run_with_weak_tm_is_moderate_for_pcr(self):
        result = analyze_pair(
            PrimerInput("a", "GCCGGAGACCCAGGCGCGGC"),
            PrimerInput("b", "TCGCCGAAGCACCGGAGAGT"),
            "heterodimer",
            PrimerDimerConditions(annealing_temp_c=60),
            PrimerDimerThresholds(high_delta_g_kcal_mol=-999, moderate_delta_g_kcal_mol=-999),
        )

        self.assertEqual(result["risk"], "MODERATE")
        self.assertEqual(result["terminal_3prime_run"], 4)
        self.assertLess(result["dimer_tm_c"], 55)

    def test_category_is_symmetric_for_heterodimer(self):
        conditions = PrimerDimerConditions(annealing_temp_c=55)
        thresholds = PrimerDimerThresholds(high_delta_g_kcal_mol=-999, moderate_delta_g_kcal_mol=-999)
        primer_a = PrimerInput("a", "AAAACCCC")
        primer_b = PrimerInput("b", "GGGGTTTT")

        result_ab = analyze_pair(primer_a, primer_b, "heterodimer", conditions, thresholds)
        result_ba = analyze_pair(primer_b, primer_a, "heterodimer", conditions, thresholds)

        self.assertEqual(result_ab["risk"], result_ba["risk"])

    def test_conditions_change_thermodynamic_result_reproducibly(self):
        primer_a = PrimerInput("a", "AAAACCCC")
        primer_b = PrimerInput("b", "GGGGTTTT")
        low_mg = analyze_pair(
            primer_a,
            primer_b,
            "heterodimer",
            PrimerDimerConditions(dv_conc_mM=0.5, annealing_temp_c=55),
            PrimerDimerThresholds(),
        )
        high_mg = analyze_pair(
            primer_a,
            primer_b,
            "heterodimer",
            PrimerDimerConditions(dv_conc_mM=3.0, annealing_temp_c=55),
            PrimerDimerThresholds(),
        )

        self.assertNotEqual(low_mg["dimer_tm_c"], high_mg["dimer_tm_c"])

    def test_internal_run_contributes_to_risk_even_if_preview_prefers_3prime(self):
        result = analyze_pair(
            PrimerInput("a", "AGGAGATGACTTAAGGCA"),
            PrimerInput("b", "AAAATCTCCTGAATACAA"),
            "heterodimer",
            PrimerDimerConditions(annealing_temp_c=55),
            PrimerDimerThresholds(high_delta_g_kcal_mol=-999, moderate_delta_g_kcal_mol=-999),
        )

        self.assertEqual(result["risk"], "MODERATE")
        self.assertEqual(result["longest_complementary_run"], 7)

    def test_analyze_primers_returns_expected_number_of_rows(self):
        analysis = analyze_primers([
            PrimerInput("p1", "AAAA"),
            PrimerInput("p2", "CCCC"),
            PrimerInput("p3", "TTTT"),
        ], conditions=PrimerDimerConditions(annealing_temp_c=55))

        self.assertEqual(len(analysis["results"]), 6)


class GenBankImportTests(TestCase):
    def setUp(self):
        self.project = Project.objects.create(name="YTK Import", public=False, assembly_standard="ytk")

    def write_genbank(self, directory, filename, labels, sequence="ATGC" * 40, date="18-DEC-2014"):
        record = SeqRecord(Seq(sequence), id=filename.replace(".gb", ""), name=filename.replace(".gb", ""), description=".")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "circular"
        record.annotations["date"] = date
        record.features = [
            SeqFeature(FeatureLocation(index * 4, (index * 4) + 4), type="misc_feature", qualifiers={"label": [label]})
            for index, label in enumerate(labels)
        ]
        path = Path(directory) / filename
        with path.open("w", encoding="utf-8") as handle:
            SeqIO.write(record, handle, "genbank")
        return path

    def write_ytk_part_genbank(
            self,
            directory,
            filename,
            upstream_overhang,
            downstream_overhang,
            payload="ATGC" * 40,
            feature_type="misc_feature",
            feature_label="Part",
            prefix="AA",
            suffix="TT",
            extra_features=None,
            date="18-DEC-2014"):
        sequence = (
            prefix
            + "GGTCTC"
            + "A"
            + upstream_overhang
            + payload
            + downstream_overhang
            + "A"
            + "GAGACC"
            + suffix
        )
        record = SeqRecord(Seq(sequence), id=filename.replace(".gb", ""), name=filename.replace(".gb", ""), description=".")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "circular"
        record.annotations["date"] = date

        payload_start = len(prefix) + 6 + 1 + 4
        payload_end = payload_start + len(payload)
        record.features = [
            SeqFeature(
                FeatureLocation(payload_start, payload_end),
                type=feature_type,
                qualifiers={"label": [feature_label]},
            )
        ]
        for feature in extra_features or []:
            record.features.append(feature)

        path = Path(directory) / filename
        with path.open("w", encoding="utf-8") as handle:
            SeqIO.write(record, handle, "genbank")
        return path

    def write_ytk_receiver_dropout_genbank(
            self,
            directory,
            filename,
            upstream_overhang,
            downstream_overhang,
            dropout_payload="ATGC" * 20,
            backbone_payload="GCTA" * 20,
            backbone_feature_label="CamR",
            backbone_feature_type="CDS",
            date="18-DEC-2014"):
        sequence = (
            downstream_overhang
            + "A"
            + "GAGACC"
            + dropout_payload
            + "GGTCTC"
            + "A"
            + upstream_overhang
            + backbone_payload
        )
        record = SeqRecord(Seq(sequence), id=filename.replace(".gb", ""), name=filename.replace(".gb", ""), description=".")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "circular"
        record.annotations["date"] = date

        backbone_start = len(downstream_overhang) + 1 + 6 + len(dropout_payload) + 6 + 1 + len(upstream_overhang)
        backbone_end = backbone_start + len(backbone_payload)
        record.features = [
            SeqFeature(
                FeatureLocation(backbone_start, backbone_end),
                type=backbone_feature_type,
                qualifiers={"label": [backbone_feature_label]},
            )
        ]

        path = Path(directory) / filename
        with path.open("w", encoding="utf-8") as handle:
            SeqIO.write(record, handle, "genbank")
        return path

    def test_import_genbank_creates_reference_plasmid_with_metadata(self):
        with tempfile.TemporaryDirectory() as tempdir:
            with override_settings(MEDIA_ROOT=tempdir):
                extra_features = [
                    SeqFeature(FeatureLocation(0, 7), type="protein_bind", qualifiers={"label": ["BsmBI"]}),
                    SeqFeature(FeatureLocation(20, 40), type="CDS", qualifiers={"label": ["CamR"]}),
                ]
                self.write_ytk_part_genbank(
                    tempdir,
                    "pYTK999.gb",
                    "GCTG",
                    "TACA",
                    feature_label="ConR1",
                    prefix="AACGTCTCA",
                    extra_features=extra_features,
                )

                result = import_plasmids_from_genbank_dir(tempdir, self.project)

        plasmid = Plasmid.objects.get(name="pYTK999", project=self.project)
        expected_size = len("AACGTCTCA" + "GGTCTC" + "A" + "GCTG" + ("ATGC" * 40) + "TACA" + "A" + "GAGACC" + "TT")

        self.assertEqual(result["created"], 1)
        self.assertEqual(plasmid.computed_size, expected_size)
        self.assertEqual(plasmid.level, 0)
        self.assertEqual(plasmid.type.name, "Insert")
        self.assertTrue(plasmid.reference_sequence)
        self.assertEqual(plasmid.created_on, datetime.date(2014, 12, 18))
        self.assertTrue(RestrictionEnzyme.objects.filter(name="BsaI", hf_version=False).exists())
        bsai_hf = RestrictionEnzyme.objects.get(name="BsaI", hf_version=True)
        self.assertTrue(RestrictionEnzyme.objects.filter(name="BsmBI", hf_version=False).exists())
        bsmbi = RestrictionEnzyme.objects.get(name="BsmBI", hf_version=False)
        self.assertEqual(bsai_hf.link_datasheet, "https://www.neb.com/en/products/r3733-bsai-hf-v2")
        self.assertEqual(bsai_hf.buffer_activity_map(), {
            "NEB 1.1": 100,
            "NEB 2.1": 100,
            "NEB 3.1": 100,
            "NEB CutSmart": 100,
        })
        self.assertEqual(bsmbi.link_datasheet, "https://www.neb.com/en/products/r0739-bsmbi-v2")
        self.assertEqual(bsmbi.buffer_activity_map(), {
            "NEB 1.1": 10,
            "NEB 2.1": 50,
            "NEB 3.1": 100,
            "NEB CutSmart": 25,
        })
        self.assertIn("<10%", bsmbi.description)
        self.assertEqual(list(plasmid.selectable_markers.values_list("three_letter_code", flat=True)), ["CLM"])
        self.assertIn("ConR1", plasmid.description)
        self.assertTrue(plasmid.sequence.name.endswith("pYTK999.gb"))
        self.assertEqual(plasmid.assembly_metadata["detected"]["part_type_key"], "ytk_5")
        self.assertEqual(plasmid.assembly_metadata["detected"]["source"], "digest")
        self.assertEqual(plasmid.assembly_metadata["confirmed"]["type_name"], "Insert")
        self.assertEqual(plasmid.assembly_metadata["confirmed"]["level"], 0)

    def test_plasmid_update_computed_size_without_restriction_enzyme_records(self):
        with tempfile.TemporaryDirectory() as tempdir:
            with override_settings(MEDIA_ROOT=tempdir):
                self.write_genbank(tempdir, "pYTK998.gb", ["BsaI", "BsaI(1)", "CamR", "sfGFP"])
                import_plasmids_from_genbank_dir(tempdir, self.project)

                plasmid = Plasmid.objects.get(name="pYTK998", project=self.project)
                plasmid.computed_size = None
                plasmid.insert_computed_size = None
                plasmid.save(update_fields=["computed_size", "insert_computed_size"])
                RestrictionEnzyme.objects.all().delete()

                result = plasmid_update_computed_size(plasmid)
                plasmid.refresh_from_db()

        self.assertTrue(result)
        self.assertEqual(plasmid.computed_size, 160)

    def test_import_genbank_updates_existing_plasmid_when_requested(self):
        with tempfile.TemporaryDirectory() as tempdir:
            with override_settings(MEDIA_ROOT=tempdir):
                self.write_genbank(tempdir, "pYTK001.gb", ["BsmBI", "BsmBI(1)", "CamR", "sfGFP"])
                Plasmid.objects.create(
                    name="pYTK001",
                    project=self.project,
                    intended_use="Legacy",
                    description="Old description",
                )

                result = import_plasmids_from_genbank_dir(
                    tempdir,
                    self.project,
                    update_existing=True,
                )

        plasmid = Plasmid.objects.get(name="pYTK001", project=self.project)

        self.assertEqual(result["updated"], 1)
        self.assertEqual(plasmid.type.name, "Receiver")
        self.assertEqual(plasmid.level, 0)
        self.assertEqual(plasmid.intended_use, "Imported from GenBank")
        self.assertEqual(list(plasmid.selectable_markers.values_list("three_letter_code", flat=True)), ["CLM"])
        self.assertEqual(plasmid.assembly_metadata["detected"]["source"], "legacy_labels")

    def test_classifier_is_invariant_to_rotation_and_reverse_complement(self):
        payload = "ATGC" * 20
        sequence = "AA" + "GGTCTC" + "A" + "GCTG" + payload + "TACA" + "A" + "GAGACC" + "TT"
        rotated_sequence = sequence[25:] + sequence[:25]
        reverse_complement_sequence = str(Seq(sequence).reverse_complement())

        records = []
        for index, sequence_variant in enumerate((sequence, rotated_sequence, reverse_complement_sequence), start=1):
            record = SeqRecord(Seq(sequence_variant), id=f"variant{index}", name=f"variant{index}", description=".")
            record.annotations["molecule_type"] = "DNA"
            record.annotations["topology"] = "circular"
            records.append(record)

        results = [classify_assembly_record(record, standard_id="ytk", allow_legacy_fallback=False) for record in records]

        self.assertTrue(all(result is not None for result in results))
        self.assertEqual({result.part_type_key for result in results}, {"ytk_5"})
        self.assertEqual({result.model_type_name for result in results}, {"Insert"})
        self.assertEqual({result.model_level for result in results}, {0})

    def test_import_prefers_digest_based_classification_even_with_extra_bsmbi_site(self):
        with tempfile.TemporaryDirectory() as tempdir:
            with override_settings(MEDIA_ROOT=tempdir):
                self.write_ytk_part_genbank(
                    tempdir,
                    "pYTK067_like.gb",
                    "GCTG",
                    "TACA",
                    feature_label="Con1",
                    prefix="AACGTCTCA",
                    extra_features=[SeqFeature(FeatureLocation(2, 8), type="protein_bind", qualifiers={"label": ["BsmBI"]})],
                )

                result = import_plasmids_from_genbank_dir(tempdir, self.project)

        plasmid = Plasmid.objects.get(name="pYTK067_like", project=self.project)

        self.assertEqual(result["created"], 1)
        self.assertEqual(plasmid.type.name, "Insert")
        self.assertEqual(plasmid.level, 0)
        self.assertEqual(plasmid.assembly_metadata["detected"]["part_type_key"], "ytk_5")

    def test_classifier_detects_ytk_234r_receiver_dropout(self):
        with tempfile.TemporaryDirectory() as tempdir:
            genbank_path = self.write_ytk_receiver_dropout_genbank(
                tempdir,
                "pYTK047_like.gb",
                "GCTG",
                "AACG",
            )
            record = SeqIO.read(str(genbank_path), "genbank")

        result = classify_assembly_record(record, standard_id="ytk", allow_legacy_fallback=False)

        self.assertIsNotNone(result)
        self.assertEqual(result.part_type_key, "ytk_234r")
        self.assertEqual(result.model_type_name, "Receiver")
        self.assertEqual(result.model_level, 1)

    def test_import_detects_ytk_234r_receiver_dropout_as_level_one_receiver(self):
        with tempfile.TemporaryDirectory() as tempdir:
            with override_settings(MEDIA_ROOT=tempdir):
                self.write_ytk_receiver_dropout_genbank(
                    tempdir,
                    "pYTK047_like.gb",
                    "GCTG",
                    "AACG",
                )

                result = import_plasmids_from_genbank_dir(tempdir, self.project)

        plasmid = Plasmid.objects.get(name="pYTK047_like", project=self.project)

        self.assertEqual(result["created"], 1)
        self.assertEqual(plasmid.type.name, "Receiver")
        self.assertEqual(plasmid.level, 1)
        self.assertEqual(plasmid.assembly_metadata["detected"]["part_type_key"], "ytk_234r")
        self.assertEqual(plasmid.assembly_metadata["confirmed"]["type_name"], "Receiver")
        self.assertEqual(plasmid.assembly_metadata["confirmed"]["level"], 1)

    def test_import_uploaded_genbanks_supports_multiple_files(self):
        with tempfile.TemporaryDirectory() as tempdir:
            first_path = self.write_genbank(tempdir, "pYTK111.gb", ["BsaI", "BsaI(1)", "CamR", "sfGFP"])
            second_path = self.write_genbank(tempdir, "pYTK112.gbk", ["BsmBI", "BsmBI(1)", "AmpR", "mRuby2"])
            uploads = [
                SimpleUploadedFile(first_path.name, first_path.read_bytes(), content_type="text/plain"),
                SimpleUploadedFile(second_path.name, second_path.read_bytes(), content_type="text/plain"),
            ]

            with override_settings(MEDIA_ROOT=tempdir):
                result = import_plasmids_from_uploaded_genbanks(uploads, self.project)

        self.assertEqual(result["created"], 2)
        self.assertEqual(Plasmid.objects.filter(project=self.project, name="pYTK111").count(), 1)
        self.assertEqual(Plasmid.objects.filter(project=self.project, name="pYTK112").count(), 1)

    def test_import_uploaded_genbanks_sorts_files_alphabetically_before_creation(self):
        with tempfile.TemporaryDirectory() as tempdir:
            first_path = self.write_genbank(tempdir, "pYTK002.gb", ["BsaI", "BsaI(1)", "CamR", "sfGFP"])
            second_path = self.write_genbank(tempdir, "pYTK001.gbk", ["BsmBI", "BsmBI(1)", "AmpR", "mRuby2"])
            uploads = [
                SimpleUploadedFile(first_path.name, first_path.read_bytes(), content_type="text/plain"),
                SimpleUploadedFile(second_path.name, second_path.read_bytes(), content_type="text/plain"),
            ]

            with override_settings(MEDIA_ROOT=tempdir):
                result = import_plasmids_from_uploaded_genbanks(uploads, self.project)

        plasmids_by_idx = list(Plasmid.objects.filter(project=self.project).order_by("idx").values_list("name", flat=True))

        self.assertEqual(result["created"], 2)
        self.assertEqual(plasmids_by_idx, ["pYTK001", "pYTK002"])

    def test_plasmid_import_view_imports_into_current_project(self):
        user = User.objects.create_user(username="genbank-user", password="pw")
        Membership.objects.create(member=user, project=self.project, access_policies="w")
        self.client.force_login(user)
        self.client.cookies["current_project_id"] = str(self.project.id)

        with tempfile.TemporaryDirectory() as tempdir:
            genbank_path = self.write_genbank(tempdir, "pYTK201.gb", ["BsaI", "BsaI(1)", "CamR", "Venus"])
            upload = SimpleUploadedFile(genbank_path.name, genbank_path.read_bytes(), content_type="text/plain")

            with override_settings(MEDIA_ROOT=tempdir):
                response = self.client.post(reverse("plasmid_import"), {
                    "target_mode": "existing",
                    "project": str(self.project.id),
                    "genbank_files": [upload],
                    "name_source": "filename",
                    "update_existing": "",
                    "public_visibility": "on",
                    "reference_sequence": "on",
                    "infer_ytk_metadata": "on",
                })

        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.url, reverse("plasmids") + "?form_result_plasmid_import_success=true")
        plasmid = Plasmid.objects.get(name="pYTK201", project=self.project)
        self.assertTrue(plasmid.public_visibility)
        self.assertTrue(plasmid.reference_sequence)
        self.assertEqual(response.cookies["current_project_id"].value, str(self.project.id))

    def test_plasmid_detail_view_shows_persisted_assembly_classification(self):
        user = User.objects.create_user(username="assembly-detail-user", password="pw")
        Membership.objects.create(member=user, project=self.project, access_policies="w")
        self.client.force_login(user)

        with tempfile.TemporaryDirectory() as tempdir:
            with override_settings(MEDIA_ROOT=tempdir):
                self.write_ytk_part_genbank(
                    tempdir,
                    "pYTK067_detail.gb",
                    "GCTG",
                    "TACA",
                    feature_label="Con1",
                    prefix="AACGTCTCA",
                )
                import_plasmids_from_genbank_dir(tempdir, self.project)

                plasmid = Plasmid.objects.get(name="pYTK067_detail", project=self.project)
                response = self.client.get(reverse("plasmid", args=(plasmid.id,)))

                self.assertEqual(response.status_code, 200)
                self.assertContains(response, "Assembly classification")
                self.assertContains(response, "YTK Part 5")
                self.assertContains(response, "GCTG")
                self.assertContains(response, "TACA")

    def test_api_plasmids_exposes_part_metadata_and_filters(self):
        user = User.objects.create_user(username="assembly-api-user", password="pw")
        Membership.objects.create(member=user, project=self.project, access_policies="w")
        self.client.force_login(user)
        self.client.cookies["current_project_id"] = str(self.project.id)

        with tempfile.TemporaryDirectory() as tempdir:
            with override_settings(MEDIA_ROOT=tempdir):
                self.write_ytk_part_genbank(
                    tempdir,
                    "pYTK067_api.gb",
                    "GCTG",
                    "TACA",
                    feature_label="Con1",
                    prefix="AACGTCTCA",
                )
                import_plasmids_from_genbank_dir(tempdir, self.project)

        imported = Plasmid.objects.get(name="pYTK067_api", project=self.project)
        imported.created_on = datetime.date(2026, 8, 20)
        imported.save(update_fields=["created_on"])
        newest = Plasmid.objects.create(
            name="newest-plasmid",
            intended_use="",
            project=self.project,
            created_on=datetime.date(2026, 8, 21),
        )

        response = self.client.get(reverse("api-plasmids"))

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload["page"], 1)
        self.assertEqual(payload["page_size"], 50)
        self.assertEqual(payload["plasmids"][0]["i"], str(newest.id))
        self.assertEqual(payload["plasmids"][1]["pk"], "ytk_5")
        self.assertEqual(payload["plasmids"][1]["pnm"], "YTK Part 5")
        self.assertIn(
            ["part", "Part", [["P5", "ap-ytk-5", "info"]]],
            payload["table_filters"],
        )

        response = self.client.get(reverse("api-plasmids"), {
            "recent_ids": f"{imported.id},{newest.id}",
        })
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["page_size"], 50)
        self.assertEqual(response.json()["main_page_size"], 48)
        self.assertEqual(
            [row["i"] for row in response.json()["recently_viewed"]],
            [str(imported.id), str(newest.id)],
        )

        response = self.client.get(reverse("api-plasmids"), {
            "q": "pYTK067_api",
            "search": "name",
        })
        self.assertEqual(response.status_code, 200)
        self.assertEqual(len(response.json()["plasmids"]), 1)

    def test_plasmid_import_view_can_create_destination_project(self):
        user = User.objects.create_user(username="new-project-user", password="pw")
        self.client.force_login(user)

        with tempfile.TemporaryDirectory() as tempdir:
            genbank_path = self.write_genbank(tempdir, "pYTK301.gb", ["BsaI", "BsaI(1)", "CamR", "Venus"])
            upload = SimpleUploadedFile(genbank_path.name, genbank_path.read_bytes(), content_type="text/plain")

            with override_settings(MEDIA_ROOT=tempdir):
                response = self.client.post(reverse("plasmid_import"), {
                    "target_mode": "new",
                    "project": "",
                    "new_project_name": "Imported From Form",
                    "new_project_public": "",
                    "new_project_assembly_standard": "ytk",
                    "genbank_files": [upload],
                    "name_source": "filename",
                    "update_existing": "",
                    "public_visibility": "",
                    "reference_sequence": "on",
                    "infer_ytk_metadata": "on",
                })

        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.url, reverse("plasmids") + "?form_result_plasmid_import_success=true")
        project = Project.objects.get(name="Imported From Form")
        self.assertEqual(project.assembly_standard, "ytk")
        self.assertTrue(Membership.objects.filter(member=user, project=project, access_policies="a").exists())
        self.assertTrue(Plasmid.objects.filter(name="pYTK301", project=project).exists())
        self.assertEqual(response.cookies["current_project_id"].value, str(project.id))


class ExperimentManagementTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="experiment-user", password="pw")
        self.project = Project.objects.create(name="Experiment Project", public=False)
        self.other_project = Project.objects.create(name="Other Project", public=False)
        Membership.objects.create(member=self.user, project=self.project, access_policies="w")
        Membership.objects.create(member=self.user, project=self.other_project, access_policies="r")
        self.plasmid = Plasmid.objects.create(
            name="pRoot", intended_use="Experiment", project=self.project
        )
        self.other_plasmid = Plasmid.objects.create(
            name="pOther", intended_use="Experiment", project=self.other_project
        )
        self.client.force_login(self.user)

    def test_create_assigns_plasmids_from_selected_project(self):
        response = self.client.post(reverse("experiment_create"), {
            "name": "Assembly run",
            "description": "Test assembly",
            "project": self.project.id,
            "plasmids": [str(self.plasmid.id)],
        })

        self.assertRedirects(response, reverse("experiments"))
        experiment = Experiment.objects.get(name="Assembly run")
        self.assertEqual(experiment.project, self.project)
        self.assertEqual(list(experiment.plasmids.all()), [self.plasmid])

    def test_create_form_renders_target_dependency_preview(self):
        response = self.client.get(reverse("experiment_create"))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "experiment-form-plasmids")
        self.assertContains(response, "Automatically included components")
        self.assertContains(response, "plasmid-dependencies")
        self.assertContains(response, "plasmid-levels")
        self.assertContains(response, "plasmid-urls")
        self.assertContains(response, "bi-box-arrow-up-right")
        self.assertContains(response, "Saltar a la página del plasmidio")
        self.assertContains(response, "Unassigned project")
        self.assertContains(response, "<details", html=False)

    def test_create_accepts_visible_plasmid_from_another_project(self):
        response = self.client.post(reverse("experiment_create"), {
            "name": "Cross-project run",
            "project": self.project.id,
            "plasmids": [str(self.other_plasmid.id)],
        })

        self.assertRedirects(response, reverse("experiments"))
        experiment = Experiment.objects.get(name="Cross-project run")
        self.assertEqual(experiment.project, self.project)
        self.assertEqual(list(experiment.plasmids.all()), [self.other_plasmid])

    def test_create_allows_unassigned_project(self):
        response = self.client.post(reverse("experiment_create"), {
            "name": "Unassigned run",
            "description": "No project yet",
            "plasmids": [str(self.plasmid.id)],
        })

        self.assertRedirects(response, reverse("experiments"))
        experiment = Experiment.objects.get(name="Unassigned run")
        self.assertIsNone(experiment.project)

        response = self.client.get(reverse("experiments"))
        self.assertContains(response, "Unassigned project")
        self.assertContains(response, "Unassigned run")

    def test_experiment_is_archived_when_all_plasmids_are_final(self):
        self.plasmid.ligation_state = 1
        self.plasmid.colonypcr_state = 0
        self.plasmid.digestion_state = 0
        self.plasmid.sequencing_state = 0
        self.plasmid.save()
        experiment = Experiment.objects.create(name="Finished run", project=self.project)
        experiment.plasmids.add(self.plasmid)

        response = self.client.get(reverse("experiments"))
        self.assertContains(response, "Archived")

        response = self.client.get(reverse("api-experiments-map"))
        payload = response.json()
        mapped = payload["projects"][0]["experiments"][0]
        self.assertTrue(mapped["archived"])
        self.assertEqual(mapped["stats"]["progress"], 100)

    def test_edit_updates_assignments_and_delete_does_not_delete_plasmids(self):
        experiment = Experiment.objects.create(
            name="Existing run", project=self.project
        )
        response = self.client.post(reverse("experiment_edit", args=(experiment.id,)), {
            "name": "Updated run",
            "description": "Updated",
            "project": self.project.id,
            "plasmids": [str(self.plasmid.id)],
        })
        self.assertRedirects(response, reverse("experiments"))
        experiment.refresh_from_db()
        self.assertEqual(experiment.name, "Updated run")
        self.assertEqual(list(experiment.plasmids.all()), [self.plasmid])

        response = self.client.post(reverse("experiment_delete", args=(experiment.id,)))
        self.assertRedirects(response, reverse("experiments"))
        self.assertFalse(Experiment.objects.filter(pk=experiment.id).exists())
        self.assertTrue(Plasmid.objects.filter(pk=self.plasmid.id).exists())


class RestrictionEnzymeCreateTests(TestCase):
    def test_recommended_enzyme_ignores_hf_duplicate(self):
        project = Project.objects.create(name="YTK", public=False, assembly_standard="ytk")
        plasmid = Plasmid.objects.create(
            name="L1 test",
            intended_use="test",
            level=1,
            project=project,
        )
        RestrictionEnzyme.objects.create(name="BsaI", hf_version=True)
        RestrictionEnzyme.objects.create(name="BsaI", hf_version=False)

        self.assertEqual(plasmid.recommended_enzyme_for_create(return_name=True), "BsaI")

    def test_build_restriction_enzymes_deduplicates_hf_variants(self):
        RestrictionEnzyme.objects.filter(name="BsaI").delete()
        regular = RestrictionEnzyme.objects.create(name="BsaI", hf_version=False)
        RestrictionEnzyme.objects.create(name="BsaI", hf_version=True)
        other = RestrictionEnzyme.objects.create(name="SapI", hf_version=True)

        enzymes = build_restriction_enzymes()

        enzyme_by_name = {enzyme.name: enzyme for enzyme in enzymes}
        self.assertEqual(list(enzyme_by_name), ["BsaI", "BsmBI", "SapI"])
        self.assertEqual(enzyme_by_name["BsaI"].pk, regular.pk)
        self.assertEqual(enzyme_by_name["SapI"].pk, other.pk)

    def test_restrictionenzyme_edit_updates_metadata_and_buffers(self):
        user = User.objects.create_user(username="enzyme-edit-user", password="pw")
        self.client.force_login(user)
        RestrictionEnzyme.objects.filter(name="BsaI", hf_version=False).delete()
        enzyme = RestrictionEnzyme.objects.create(
            name="BsaI",
            hf_version=False,
            description="Old description",
        )
        old_buffer = RestrictionBuffer.objects.create(name="Old buffer")
        new_buffer = RestrictionBuffer.objects.create(name="New buffer")
        RestrictionEnzymeBuffer.objects.create(
            restriction_enzyme=enzyme,
            buffer=old_buffer,
            activity_percent=60,
        )

        response = self.client.post(reverse("restrictionenzyme_edit", args=(enzyme.id,)), {
            "name": "BsaI",
            "hf_version": "",
            "link_datasheet": "https://example.com/bsai",
            "description": "Updated description",
            "buffers-TOTAL_FORMS": "1",
            "buffers-INITIAL_FORMS": "0",
            "buffers-MIN_NUM_FORMS": "0",
            "buffers-MAX_NUM_FORMS": "1000",
            "buffers-0-existing_buffer": str(new_buffer.id),
            "buffers-0-new_buffer_name": "",
            "buffers-0-activity_percent": "95",
            "buffers-0-DELETE": "",
        })

        enzyme.refresh_from_db()
        self.assertEqual(response.status_code, 302)
        self.assertIn("form_result_restrictionenzyme_edit_success=true", response.url)
        self.assertEqual(enzyme.description, "Updated description")
        self.assertEqual(
            list(enzyme.buffer_activity_map().items()),
            [("New buffer", 95)],
        )

    def test_restrictionenzyme_create_view_creates_record(self):
        user = User.objects.create_user(username="enzyme-user", password="pw")
        self.client.force_login(user)

        response = self.client.post(reverse("restrictionenzyme_create"), {
            "name": "SapI",
            "hf_version": "",
            "link_datasheet": "",
            "description": "Loaded from frontend",
            "buffers-TOTAL_FORMS": "1",
            "buffers-INITIAL_FORMS": "0",
            "buffers-MIN_NUM_FORMS": "0",
            "buffers-MAX_NUM_FORMS": "1000",
            "buffers-0-existing_buffer": "__new__",
            "buffers-0-new_buffer_name": "NEB CutSmart",
            "buffers-0-activity_percent": "100",
            "buffers-0-DELETE": "",
        })

        enzyme = RestrictionEnzyme.objects.get(name="SapI", hf_version=False)
        self.assertEqual(response.status_code, 302)
        self.assertEqual(
            response.url,
            reverse("restrictionenzyme", args=(enzyme.id,)) + "?form_result_restrictionenzyme_create_success=true",
        )
        self.assertEqual(enzyme.buffer_activity_map(), {"NEB CutSmart": 100})
        self.assertEqual(enzyme.description, "Loaded from frontend")
        self.assertTrue(RestrictionBuffer.objects.filter(name="NEB CutSmart").exists())

    def test_restrictionenzyme_create_view_reuses_existing_buffer(self):
        user = User.objects.create_user(username="enzyme-buffer-user", password="pw")
        self.client.force_login(user)
        buffer = RestrictionBuffer.objects.create(name="Buffer R")

        response = self.client.post(reverse("restrictionenzyme_create"), {
            "name": "SapI",
            "hf_version": "",
            "link_datasheet": "",
            "description": "",
            "buffers-TOTAL_FORMS": "1",
            "buffers-INITIAL_FORMS": "0",
            "buffers-MIN_NUM_FORMS": "0",
            "buffers-MAX_NUM_FORMS": "1000",
            "buffers-0-existing_buffer": str(buffer.id),
            "buffers-0-new_buffer_name": "",
            "buffers-0-activity_percent": "75",
            "buffers-0-DELETE": "",
        })

        enzyme = RestrictionEnzyme.objects.get(name="SapI", hf_version=False)
        self.assertEqual(response.status_code, 302)
        self.assertEqual(enzyme.buffer_activity_map(), {"Buffer R": 75})
        self.assertEqual(RestrictionBuffer.objects.filter(name="Buffer R").count(), 1)

    def test_restrictionenzyme_create_view_blocks_duplicate_name_and_hf(self):
        user = User.objects.create_user(username="enzyme-duplicate-user", password="pw")
        self.client.force_login(user)
        RestrictionEnzyme.objects.create(name="BsmBI", hf_version=True)

        response = self.client.post(reverse("restrictionenzyme_create"), {
            "name": "BsmBI",
            "hf_version": "on",
            "link_datasheet": "",
            "description": "",
            "buffers-TOTAL_FORMS": "1",
            "buffers-INITIAL_FORMS": "0",
            "buffers-MIN_NUM_FORMS": "0",
            "buffers-MAX_NUM_FORMS": "1000",
            "buffers-0-existing_buffer": "",
            "buffers-0-new_buffer_name": "",
            "buffers-0-activity_percent": "",
            "buffers-0-DELETE": "",
        })

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "BsmBI-HF is already loaded in Weaver.")

    def test_restrictionenzyme_delete_view_removes_record(self):
        user = User.objects.create_user(username="enzyme-delete-user", password="pw")
        self.client.force_login(user)
        enzyme = RestrictionEnzyme.objects.create(name="SapI", hf_version=False)

        response = self.client.post(reverse("restrictionenzyme_delete", args=(enzyme.id,)))

        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.url, reverse("restrictionenzymes") + "?form_result_object_deleted=true")
        self.assertFalse(RestrictionEnzyme.objects.filter(id=enzyme.id).exists())
