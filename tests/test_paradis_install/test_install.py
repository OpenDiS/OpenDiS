"""Exercise the ParaDiS installer without compiling the external simulator."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


MAKE = os.environ.get("MAKE", "make")
MAKEFILE = Path(__file__).resolve().parents[2] / "extensions/paradis/Makefile"
ARTIFACTS = ("lib/libparadis.so", "lib/Home.py", "python/paradis_util.py")


@unittest.skipUnless(shutil.which(MAKE), "GNU Make is required")
class InstallTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.extension = self.root / "extension"
        self.source = self.root / "source"
        for tree in (self.extension, self.source):
            (tree / "lib").mkdir(parents=True)
            (tree / "python").mkdir()
        shutil.copyfile(MAKEFILE, self.extension / "Makefile")
        for artifact in ARTIFACTS:
            (self.source / artifact).write_text(artifact, encoding="utf-8")

    def make(self, *args):
        env = os.environ.copy()
        env.pop("PARADIS_DIR", None)
        return subprocess.run(
            [MAKE, *args], cwd=self.extension, env=env,
            capture_output=True, text=True,
        )

    def assert_no_links(self):
        for artifact in ARTIFACTS:
            self.assertFalse(os.path.lexists(self.extension / artifact))

    def assert_installed(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        for artifact in ARTIFACTS:
            destination = self.extension / artifact
            self.assertEqual(destination.read_text(encoding="utf-8"), artifact)
            if os.name != "nt":
                self.assertTrue(destination.is_symlink())
                self.assertEqual(destination.resolve(), self.source / artifact)

    def test_install_requires_source_directory(self):
        result = self.make("install")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PARADIS_DIR is undefined", result.stderr)
        self.assert_no_links()

    def test_relative_source_directory(self):
        self.assert_installed(self.make("install", "PARADIS_DIR=../source"))

    def test_absolute_source_directory(self):
        self.assert_installed(self.make("install", f"PARADIS_DIR={self.source.as_posix()}"))

    def test_missing_wrapper_stops_parallel_install(self):
        (self.source / "python/paradis_util.py").unlink()
        result = self.make("-j4", "install", f"PARADIS_DIR={self.source.as_posix()}")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("python/paradis_util.py", result.stderr)
        self.assert_no_links()

    def test_direct_file_target_validates_source(self):
        result = self.make("lib/Home.py")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PARADIS_DIR is undefined", result.stderr)
        self.assert_no_links()

    def test_repeated_parallel_install(self):
        self.assert_installed(self.make("-j4", "install", "PARADIS_DIR=../source"))
        self.assert_installed(self.make("-j4", "install", "PARADIS_DIR=../source"))


if __name__ == "__main__":
    unittest.main()
