# Conda setup: one profile for Windows and WSL

This project uses a single **environment.yml** as the shared conda profile. You create one conda environment on Windows and one in WSL from the same file; they stay in sync by reusing this definition.

## 1. WSL: install Conda (one-time)

Conda on Windows does not run inside WSL. Install Linux Miniconda inside WSL:

```bash
# From project root in WSL
bash scripts/install_miniconda_wsl.sh
source ~/.bashrc
conda --version
```

## 2. Create the environment on each side

### Windows (PowerShell or CMD)

From the project root (`g:\repos\bio` or `\\wsl$\...\bio`):

```bash
conda env create -f environment.yml
conda activate bio
```

### WSL

From the project root in WSL:

```bash
conda env create -f environment.yml
conda activate bio
```

## 3. Updating the shared profile

After changing packages on one side, export and commit so the other side can recreate the same env:

```bash
conda env export --from-history > environment.yml
```

Then on the other OS:

```bash
conda env remove -n bio
conda env create -f environment.yml
conda activate bio
```

## 4. (Optional) Use Cursor with WSL

If you open this folder in Cursor via **WSL** (e.g. "Reopen Folder in WSL"), the integrated terminal and run/debug use WSL. Only the WSL conda env is used there; keep **environment.yml** as the single source of truth and create the `bio` env in WSL as in step 2.
