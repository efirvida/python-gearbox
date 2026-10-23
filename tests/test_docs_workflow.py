"""Public documentation delivery configuration (requires the optional docs extra)."""

from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml", reason="documentation workflow checks need the docs extra")
ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = yaml.load((ROOT / ".github/workflows/docs.yml").read_text(), Loader=yaml.BaseLoader)
DOCS_URL = "https://efirvida.github.io/python-gearbox/"


def test_every_master_push_runs_documentation() -> None:
    push = WORKFLOW["on"]["push"]
    assert push["branches"] == ["master"]
    assert "paths" not in push
    assert "paths-ignore" not in push


def test_master_push_deploys_and_manual_publication_remains_explicit() -> None:
    deploy = WORKFLOW["jobs"]["deploy"]
    assert deploy["needs"] == "build"
    assert deploy["if"] == (
        "github.ref == 'refs/heads/master' && (github.event_name == 'push' || "
        "(github.event_name == 'workflow_dispatch' && inputs.publish))"
    )
    assert WORKFLOW["on"]["workflow_dispatch"]["inputs"]["publish"]["default"] == "false"


def test_pull_requests_keep_build_only_permissions_and_strict_checks() -> None:
    assert "pull_request" in WORKFLOW["on"]
    assert WORKFLOW["permissions"] == {"contents": "read"}
    build = WORKFLOW["jobs"]["build"]
    assert "permissions" not in build
    commands = [step.get("run", "") for step in build["steps"]]
    assert 'pip install -e ".[docs]"' in commands
    assert "sphinx-build -W --keep-going -b html docs docs/_build/html" in commands
    assert WORKFLOW["jobs"]["deploy"]["permissions"] == {
        "pages": "write", "id-token": "write"
    }
    assert WORKFLOW["jobs"]["deploy"]["environment"]["name"] == "github-pages"


def test_readme_links_to_the_online_documentation() -> None:
    readme = (ROOT / "README.md").read_text()
    assert f"[Documentation]({DOCS_URL})" in readme
