# 📊 Survival Analysis for Data Scientists (2AMS11 - TU/e)

Welcome to the repository for **2AMS11 Survival Analysis** at Eindhoven University of Technology (TU/e).

[![Canvas Course Page](https://img.shields.io/badge/Canvas-Course%20Page-%23E05A47?style=for-the-badge&logo=instructure&logoColor=white)](https://canvas.tue.nl/courses/34519)

---

## 📌 Important Announcements & Deadlines

> **First Assignment Notice**
> * **Release Date:** Friday, September 11, 2026
> * **Deadline:** Friday, September 18, 2026

### 👥 Group Formation
To complete the assignments, you **must** join a group on Canvas:
* Please ensure you do this **before next Friday**.
* If you are not yet part of a group, please coordinate to form or join a group of **three students**.

---
*This repository will house code, notebooks, and deliverables for the course.*


## Assignments 2 and 3: simulation and results

Use Python 3.11 or newer and [uv](https://docs.astral.sh/uv/):

```sh
cd assignment_2
uv sync --locked --no-dev
uv run --locked --no-dev python data_generation.py
uv run --locked --no-dev python part3_results.py
```

The simulation runs 2,000 repetitions per dropout scenario and writes results
to `data/`. The second script creates `tables/coefficients.tex` and
`figures/recurrence_risk.pdf` (also saved as PNG) for Part III.
