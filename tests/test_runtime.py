import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import runtime
import amoebacrisp


class RuntimeTests(unittest.TestCase):
    def test_shared_config_and_storage_override(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / 'settings.json'
            config.write_text(json.dumps({'storage': directory, 'linux_env_bin': '/custom/bin'}), encoding='utf-8')
            with patch.dict(os.environ, {'AMOEBA_CRISPR_CONFIG': str(config), 'AMOEBA_CRISPR_HOME': str(Path(directory) / 'override')}):
                self.assertEqual(runtime.storage_path('outputs'), Path(directory) / 'override' / 'outputs')
                self.assertEqual(amoebacrisp.default_config()['linux_env_bin'], '/custom/bin')
                self.assertEqual(amoebacrisp.default_config()['storage'], str(Path(directory) / 'override'))

    def test_missing_config_uses_defaults(self):
        with patch.dict(os.environ, {'AMOEBA_CRISPR_CONFIG': str(Path(tempfile.gettempdir()) / 'missing-amoeba-config.json'), 'AMOEBA_CRISPR_HOME': 'F:/example'}):
            self.assertEqual(runtime.settings(), {})
            self.assertEqual(runtime.storage_path('outputs'), Path('F:/example/outputs'))
