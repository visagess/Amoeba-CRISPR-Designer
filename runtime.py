"""Shared storage and optional tool locations for the distributable app."""
import json
import os
from pathlib import Path


def settings():
    path = Path(os.environ.get('AMOEBA_CRISPR_CONFIG', Path(__file__).with_name('config.json')))
    return json.loads(path.read_text(encoding='utf-8-sig')) if path.exists() else {}


def storage_path(*parts):
    root = os.environ.get('AMOEBA_CRISPR_HOME') or settings().get('storage', 'F:/AmoebaCrisp')
    return Path(root).joinpath(*parts)
