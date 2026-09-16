# Eukaryoma PPI Webserver — Setup for Collaborators

You need two things: the **code** (this GitHub repo) and the **data**
(`eukaryoma_ppi_data_20260916.tar.gz`, ~4.8 GB — shared separately, too large
for git). The archive already contains the pre-built structures and score
indexes, so you do **not** need to run any data-conversion step or have the
custom `foldcomp` binary (it's macOS-only anyway).

Requires **Python 3.10+** and **git**.

## 1. Get the code

```bash
git clone https://github.com/evocellnet/Eukaryoma-PPI-webserver.git
cd Eukaryoma-PPI-webserver
```

## 2. Get the data

Ask for `eukaryoma_ppi_data_20260916.tar.gz` and extract it anywhere on disk —
it does not need to be near the repo. It unpacks into a single `data/`
folder.

**macOS / Linux:**
```bash
mkdir -p ~/eukaryoma-ppi && tar -xzf eukaryoma_ppi_data_20260916.tar.gz -C ~/eukaryoma-ppi
# -> ~/eukaryoma-ppi/data
```

**Windows (PowerShell — `tar` is built into Windows 10/11):**
```powershell
mkdir $HOME\eukaryoma-ppi
tar -xzf eukaryoma_ppi_data_20260916.tar.gz -C $HOME\eukaryoma-ppi
```

## 3. Point the app at the data

Set the `EUKARYOMA_DATA_DIR` environment variable to that extracted `data/`
folder (do this every time in a new terminal, or add it to your shell
profile / a permanent environment variable).

**macOS / Linux (bash/zsh):**
```bash
export EUKARYOMA_DATA_DIR=~/eukaryoma-ppi/data
```

**Windows (PowerShell):**
```powershell
$env:EUKARYOMA_DATA_DIR = "$HOME\eukaryoma-ppi\data"
```

## 4. Install and run

**macOS / Linux:**
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install -e .
streamlit run app.py
```

**Windows (PowerShell):**
```powershell
python -m venv venv
venv\Scripts\Activate.ps1
pip install -r requirements.txt
pip install -e .
streamlit run app.py
```

Then open **http://localhost:8501**.

## If a package manager is missing

- **macOS**: install [Homebrew](https://brew.sh), then `brew install python git`.
- **Linux (Debian/Ubuntu)**: `sudo apt install python3 python3-venv git`.
- **Windows**: install [Python](https://www.python.org/downloads/) (check
  "Add python.exe to PATH" during install) and
  [Git for Windows](https://git-scm.com/download/win).

## Data size

- Compressed archive: **~4.8 GB**.
- Uncompressed on disk after extracting: **~14.5 GB**, almost all of it
  `data/structures/` (the decompressed 3D structures).
