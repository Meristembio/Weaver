<div align="center">
  <img src="Django/static/img/logo-original-gh.png" width="100%" alt="Weaver">
</div>

# Weaver

Weaver is an open-source Django web application for DNA design, plasmid inventory, construct verification, and laboratory workflow tracking. It helps teams manage biological constructs, inspect plasmid sequences, plan molecular biology workflows, and keep validation evidence connected to inventory records.

The current Sticta fork combines a Django backend, server-rendered inventory views, bundled sequence visualization assets, and several React modules compiled into Django static files.

## Documentation

- [Current Weaver documentation](docs/current.md)
- [Feature documentation](docs/features/README.md)
- [Project history](docs/history/README.md)

## Main Features

- Project-based inventory for plasmids, glycerol stocks, primers, restriction enzymes, strains, boxes, locations, and protocol resources.
- Plasmid creation and design workflows, including manual records, assembly wizard flows, Level 0 design, and Golden Gate-style standards.
- Plasmid map and sequence visualization through the bundled Open Vector Editor integration.
- Sanger sequence verification with multiple-read uploads, AB1/PHD.1/SEQ processing, read alignment, chromatogram viewing, quality/low-confidence visualization, and manual verification decisions.
- PCR and primer tools, including plasmid PCR prediction, selected-region PCR suggestions, amplicon find, primer import, global PCR search, and primer dimer analysis.
- Restriction digest tools, including cut-site review, fragment prediction, buffer activity display, and OVE-integrated digest planning.
- Lab support services such as local BLAST search, stats, batch label printing, and an experiments map.

See [docs/features/README.md](docs/features/README.md) for detailed explanations and GIF walkthrough placeholders for each module.

## Architecture

```text
.
├── Django/
│   ├── manage.py
│   ├── Weaver/              # settings, URLs, ASGI/WSGI
│   ├── organization/        # projects and memberships
│   ├── inventory/           # plasmids, stocks, primers, Sanger, PCR, digest, services
│   ├── protocols/           # reactives, components, recipes, variants
│   └── static/              # CSS, JavaScript, images, compiled frontend assets
├── ReactModules/            # React source packages
├── docs/                    # feature and project-history documentation
├── requirements.txt
└── LICENSE.md
```

The configured database backend is SQLite. By default, Weaver stores local data, uploads, and media under `Django/local_data`, controlled by `WEAVER_DATA_DIR`.

## Requirements

- Python 3.12 is the verified local runtime.
- Django 6.0.7 is pinned in `requirements.txt`.
- `pip` and Python `venv` for backend dependencies.
- SQLite for the default database.
- Node.js and npm only when rebuilding or testing React modules.
- A modern browser for Open Vector Editor, chromatograms, and React-powered views.

## Installation

Clone the Sticta fork:

```bash
git clone https://github.com/Sticta-Biologicals/Weaver.git
cd Weaver
```

Create and activate a Python environment:

```bash
cd Django
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r ../requirements.txt
```

Configure local environment variables if needed:

```bash
export WEAVER_SECRET_KEY="change-this-for-your-installation"
export WEAVER_ALLOWED_HOSTS="127.0.0.1,localhost,weaver.local"
export WEAVER_DATA_DIR="$PWD/local_data"
```

Create the database schema and an admin user:

```bash
python manage.py migrate
python manage.py createsuperuser
```

Check and run Weaver:

```bash
python manage.py check
python manage.py runserver 127.0.0.1:8000
```

Open [http://127.0.0.1:8000/](http://127.0.0.1:8000/) and sign in with the superuser you created.

## Configuration

| Name | Required | Purpose | Safe example |
| --- | --- | --- | --- |
| `WEAVER_SECRET_KEY` | Required for production | Django secret key. If unset, settings use an insecure local-development fallback. | `change-this-for-production` |
| `WEAVER_ALLOWED_HOSTS` | Required outside default hosts | Comma-separated Django `ALLOWED_HOSTS`. Defaults to `weaver.sticta.com,weaver.local`. | `127.0.0.1,localhost,weaver.local` |
| `WEAVER_DATA_DIR` | Optional | Directory for SQLite database and media uploads. Defaults to `Django/local_data`. | `/srv/weaver/data` |

`DEBUG` is currently hard-coded to `True` in `Django/Weaver/settings.py`; review this before using Weaver as a public production service.

## React Modules

React modules are independent packages under `ReactModules/`. Install and build only the module you are working on.

Example:

```bash
cd ReactModules/experiments_reactflow
npm install
npm run build
```

The Vite-based experiments module builds into Django's served static assets at `Django/static/experiments-reactflow`.

## Quality Checks

Backend:

```bash
cd Django
source .venv/bin/activate
python manage.py check
python manage.py test
```

React modules:

```bash
cd ReactModules/experiments_reactflow
npm run lint
npm run build
```

There is no repository-level command that runs all backend and frontend checks.

## License

Weaver is licensed under the GNU Affero General Public License v3. See [LICENSE.md](LICENSE.md).
