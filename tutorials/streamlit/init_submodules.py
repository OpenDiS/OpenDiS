import os
import sys
import urllib.request

# 1. Anchor paths to the writeable directory
repo_root = "/mount/src/opendis"
pydis_root = os.path.join(repo_root, "core", "pydis", "python")
target_dir = os.path.join(pydis_root, "pydis", "graph")
target_file = os.path.join(target_dir, "graph.py")

# 2. Re-create file structure locally
os.makedirs(target_dir, exist_ok=True)
with open(os.path.join(pydis_root, "pydis", "__init__.py"), "a"): pass
with open(os.path.join(target_dir, "__init__.py"), "a"): pass

# 3. Pull down the missing source module code from GitHub
if not os.path.exists(target_file) or os.path.getsize(target_file) == 0:
    raw_url = "https://githubusercontent.com"
    try:
        urllib.request.urlretrieve(raw_url, target_file)
    except Exception as e:
        print(f"Submodule fetch fail: {e}")

# 4. Inject into the active running process path memory layout
if pydis_root not in sys.path:
    sys.path.insert(0, pydis_root)

# 5. GLOBAL FIX FOR SUBPROCESSES: Force the path into the OS environment track!
# This ensures background workers like frank_read_src_pydis.py inherit this variable.
existing_pythonpath = os.environ.get("PYTHONPATH", "")
if pydis_root not in existing_pythonpath:
    os.environ["PYTHONPATH"] = f"{pydis_root}{os.pathsep}{existing_pythonpath}".strip(os.pathsep)
