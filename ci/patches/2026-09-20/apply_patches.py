"""Apply this directory's patches to the core/exadis submodule.

Three fixes sent to the exadis developer on 2026-09-20 and not yet upstream.
The pinned submodule commit does not carry them, so a freshly cloned exadis
needs them before the tests can pass. Drop a patch from this directory once
its fix lands upstream.

    python3 ci/patches/2026-09-20/apply_patches.py [path/to/exadis]

Re-running is a no-op: a patch already in the tree is reported and skipped.
Exits nonzero if any patch fails to apply, which on a clean clone means the
patch no longer matches the pinned commit.
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXADIS = HERE.parents[2] / 'core' / 'exadis'   # ci/patches/<date>/ -> root


def git_apply(repo, patch, *flags):
    """git_apply: run git apply in repo, returning (ok, stderr)"""
    cmd = ['git', '-C', str(repo), 'apply', *flags, str(patch)]
    done = subprocess.run(cmd, capture_output=True, text=True)
    return done.returncode == 0, done.stderr.strip()


def apply_one(repo, patch):
    """apply_one: apply patch unless the tree already carries it"""
    if git_apply(repo, patch, '--reverse', '--check')[0]:
        print("already applied: %s" % patch.name)
        return True
    ok, err = git_apply(repo, patch)
    print("%s: %s" % ("applied" if ok else "FAILED", patch.name))
    if not ok:
        print(err)
    return ok


def main(argv):
    repo = Path(argv[1]) if len(argv) > 1 else EXADIS
    patches = sorted(HERE.glob('*.patch'))
    print("applying %d patch(es) to %s" % (len(patches), repo))
    if not patches:
        return 1
    return 0 if all([apply_one(repo, p) for p in patches]) else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv))
