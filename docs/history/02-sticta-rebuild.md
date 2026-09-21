# Sticta Rebuild

Historical notes for the Sticta rebuild stage.

## Scope

The rebuild restored local execution, aligned dependencies, consolidated repository sources, organized React modules, and updated runtime and configuration conventions.

## Documented Changes

- React module sources were consolidated under `ReactModules/`.
- The Django application runs from `Django/manage.py`.
- Local data and uploads use `Django/local_data` through `WEAVER_DATA_DIR`.
- Backend dependencies are recorded in `requirements.txt`.
- Python 3.12 and Django 6 are the verified local runtime combination in the repository context.

These are historical repository conventions, not deployment instructions. See the root README for the current installation procedure and the feature pages for current user workflows.
