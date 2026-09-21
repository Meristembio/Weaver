from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from django.core.management.base import CommandError

from inventory.custom.genbank_import import GenBankImportError
from inventory.custom.genbank_import import import_plasmids_from_genbank_dir
from inventory.custom.standards import assembly_standards
from organization.models import Membership
from organization.models import Project


class Command(BaseCommand):
    help = "Import GenBank plasmids into Weaver's Plasmid inventory."

    def add_arguments(self, parser):
        parser.add_argument("genbank_dir", help="Directory containing .gb or .gbk files.")
        parser.add_argument(
            "--project",
            required=True,
            help="Project name or numeric project ID for imported plasmids.",
        )
        parser.add_argument(
            "--create-project",
            action="store_true",
            help="Create the project when it does not exist yet.",
        )
        parser.add_argument(
            "--project-public",
            action="store_true",
            help="Create the project as public when using --create-project.",
        )
        parser.add_argument(
            "--assembly-standard",
            choices=sorted(assembly_standards.keys()),
            default="ytk",
            help="Assembly standard to assign when creating a project.",
        )
        parser.add_argument(
            "--grant-user",
            action="append",
            default=[],
            help="Grant write access to a username on the target project. Repeatable.",
        )
        parser.add_argument(
            "--grant-superusers",
            action="store_true",
            help="Grant admin access to all superusers on the target project.",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Parse and report changes without writing plasmids to the database.",
        )
        parser.add_argument(
            "--update-existing",
            action="store_true",
            help="Update existing plasmids with the same name and project.",
        )
        parser.add_argument(
            "--public-visibility",
            action="store_true",
            help="Mark imported plasmids as publicly visible.",
        )
        parser.add_argument(
            "--constructed-plasmids",
            action="store_true",
            help="Import plasmids as constructed records instead of reference sequences.",
        )
        parser.add_argument(
            "--name-source",
            choices=("filename", "record"),
            default="filename",
            help="Use the filename stem or GenBank record name for Weaver plasmid names.",
        )
        parser.add_argument(
            "--no-infer-ytk-metadata",
            action="store_true",
            help="Disable conservative assembly resistance/type/level inference during GenBank import.",
        )

    def handle(self, *args, **options):
        project = self.get_project(
            options["project"],
            create_missing=options["create_project"],
            project_public=options["project_public"],
            assembly_standard=options["assembly_standard"],
            dry_run=options["dry_run"],
        )

        if not options["dry_run"]:
            self.grant_memberships(
                project,
                usernames=options["grant_user"],
                include_superusers=options["grant_superusers"],
            )

        try:
            result = import_plasmids_from_genbank_dir(
                options["genbank_dir"],
                project,
                dry_run=options["dry_run"],
                update_existing=options["update_existing"],
                public_visibility=options["public_visibility"],
                reference_sequence=not options["constructed_plasmids"],
                infer_ytk_metadata=not options["no_infer_ytk_metadata"],
                name_source=options["name_source"],
            )
        except GenBankImportError as error:
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

    def get_project(self, project_value, create_missing=False, project_public=False, assembly_standard="ytk", dry_run=False):
        if project_value.isdigit():
            project = Project.objects.filter(id=int(project_value)).first()
        else:
            project = Project.objects.filter(name__iexact=project_value).first()

        if project:
            return project

        if not create_missing:
            raise CommandError(f"Project not found: {project_value}")

        if dry_run:
            return Project(
                name=project_value,
                public=project_public,
                assembly_standard=assembly_standard,
            )

        return Project.objects.create(
            name=project_value,
            public=project_public,
            assembly_standard=assembly_standard,
        )

    def grant_memberships(self, project, usernames=(), include_superusers=False):
        seen = set()

        for username in usernames:
            user = User.objects.filter(username=username).first()
            if not user:
                raise CommandError(f"User not found: {username}")
            if user.pk not in seen:
                Membership.objects.get_or_create(
                    member=user,
                    project=project,
                    defaults={"access_policies": "w"},
                )
                seen.add(user.pk)

        if include_superusers:
            for user in User.objects.filter(is_superuser=True):
                if user.pk in seen:
                    continue
                membership, created = Membership.objects.get_or_create(
                    member=user,
                    project=project,
                    defaults={"access_policies": "a"},
                )
                if not created and membership.access_policies != "a":
                    membership.access_policies = "a"
                    membership.save(update_fields=["access_policies"])
