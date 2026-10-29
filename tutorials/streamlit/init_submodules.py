import os
import sys
import urllib.request

# 1. Map out where the submodule should live relative to this setup script
current_dir = os.path.dirname(__file__)
repo_root = os.path.abspath(os.path.join(current_dir, "../.."))
target_dir = os.path.join(repo_root, "pydis", "graph")
target_file = os.path.join(target_dir, "graph.py")

# 2. Re-create the folder tree structure
os.makedirs(target_dir, exist_ok=True)

# 3. Create blank __init__.py files so Python treats them as packages
with open(os.path.join(repo_root, "pydis", "__init__.py"), "a"): pass
with open(os.path.join(target_dir, "__init__.py"), "a"): pass

# 4. Automatically download the file directly from GitHub over HTTPS
if not os.path.exists(target_file) or os.path.getsize(target_file) == 0:
    raw_url = "https://githubusercontent.com"
    try:
        urllib.request.urlretrieve(raw_url, target_file)
    except Exception as e:
        print(f"Submodule fetch fail: {e}")

# 5. Inject the repository path context so your imports work seamlessly
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)
