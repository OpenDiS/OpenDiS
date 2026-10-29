import os
import sys
import urllib.request
import site

# 1. Use the absolute base path provided by Streamlit Cloud
repo_root = "/mount/src/opendis"
pydis_root = os.path.join(repo_root, "core", "pydis", "python")

# 2. Try to locate where 'pydis' lives in the active virtual environment
pydis_env_path = None
for path in site.getsitepackages():
    test_path = os.path.join(path, "pydis")
    if os.path.exists(test_path) or "venv" in path:
        pydis_env_path = test_path
        break

# Fallback to local source folder if virtual env check isn't ready
if not pydis_env_path:
    pydis_env_path = os.path.join(pydis_root, "pydis")

# 3. Establish the target graph submodule folder paths
target_dir = os.path.join(pydis_env_path, "graph")
target_file = os.path.join(target_dir, "graph.py")

# 4. Safely create directory tree structural nodes
os.makedirs(target_dir, exist_ok=True)
with open(os.path.join(pydis_env_path, "__init__.py"), "a"): pass
with open(os.path.join(target_dir, "__init__.py"), "a"): pass

# 5. Fetch the missing graph module code over HTTPS from GitHub
if not os.path.exists(target_file) or os.path.getsize(target_file) == 0:
    raw_url = "https://githubusercontent.com"
    try:
        urllib.request.urlretrieve(raw_url, target_file)
    except Exception as e:
        print(f"Submodule fetch fail: {e}")

# 6. Ensure the local path is appended as a fallback for the main process
if pydis_root not in sys.path:
    sys.path.insert(0, pydis_root)
