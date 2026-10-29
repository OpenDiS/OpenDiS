import os
import sys
import urllib.request

# 1. Use the absolute base path provided by Streamlit Cloud
repo_root = "/mount/src/opendis"

# 2. Map out exactly where pydis and the submodule belong based on your traceback layout
pydis_root = os.path.join(repo_root, "core", "pydis", "python")
target_dir = os.path.join(pydis_root, "pydis", "graph")
target_file = os.path.join(target_dir, "graph.py")

# 3. Re-create the folder structure if needed
os.makedirs(target_dir, exist_ok=True)

# 4. Ensure initialization files exist inside the package structure
with open(os.path.join(pydis_root, "pydis", "__init__.py"), "a"): pass
with open(os.path.join(target_dir, "__init__.py"), "a"): pass

# 5. Download the script file over HTTPS directly into the package structure
if not os.path.exists(target_file) or os.path.getsize(target_file) == 0:
    raw_url = "https://githubusercontent.com"
    try:
        urllib.request.urlretrieve(raw_url, target_file)
    except Exception as e:
        print(f"Submodule fetch fail: {e}")

# 6. Inject the true Python package folder route into sys.path
if pydis_root not in sys.path:
    sys.path.insert(0, pydis_root)
