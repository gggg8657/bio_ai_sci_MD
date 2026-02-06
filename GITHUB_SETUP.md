# GitHub Repository Setup

## Repository Info
- **Name**: `bio_ai_sci_MD`
- **Owner**: `gggg8657`
- **URL**: https://github.com/gggg8657/bio_ai_sci_MD

---

## Option 1: Using GitHub CLI (gh)

```bash
# Install gh (if not installed)
# Ubuntu/WSL:
sudo apt install gh

# Authenticate
gh auth login

# Create repo and push (from project root)
cd /mnt/g/repos/bio
gh repo create bio_ai_sci_MD --public --source=. --push
```

---

## Option 2: Manual Setup

### Step 1: Create Repository on GitHub
1. Go to https://github.com/new
2. Repository name: `bio_ai_sci_MD`
3. Public/Private: choose as desired
4. Do NOT initialize with README (we already have files)
5. Click "Create repository"

### Step 2: Push from WSL

```bash
cd /mnt/g/repos/bio

# Add remote (after creating repo on GitHub)
git remote add origin https://github.com/gggg8657/bio_ai_sci_MD.git

# Push
git push -u origin feature/bio-tools-conda-env

# Or push all branches
git push -u origin --all
```

---

## Current Branch
The commit is on `feature/bio-tools-conda-env` branch.

If you want to push to `main`:
```bash
git checkout main
git merge feature/bio-tools-conda-env
git push -u origin main
```

---

## Quick Script
Run the helper script after creating the repo on GitHub:
```bash
bash scripts/github_push.sh
```
