#!/usr/bin/env python3
"""Select preferred RMS pick3 triplets from pick3/<method>_pick3.csv outputs."""

from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
PICK3 = HERE / "pick3"

PREFERENCES = {
    "scJoint": "max_repro",
    "scVI": "max_repro",
    "BindSC": "max_nmi",
    "Portal": "max_nmi",
    "scDART": "min_nmi",
    "Cobolt": "min_nmi",
}


def choose(method, rule):
    path = PICK3 / f"{method}_pick3.csv"
    if not path.exists():
        return {
            "Method": method,
            "preference": rule,
            "subset": "MISSING",
            "Fig4A_repro": pd.NA,
            "Fig4B_NMI": pd.NA,
            "n_subsets_available": 0,
        }
    df = pd.read_csv(path)
    if rule == "max_nmi":
        row = df.sort_values(["Fig4B_NMI", "Fig4A_repro"], ascending=[False, False]).iloc[0]
    elif rule == "min_nmi":
        row = df.sort_values(["Fig4B_NMI", "Fig4A_repro"], ascending=[True, True]).iloc[0]
    elif rule == "max_repro":
        row = df.sort_values(["Fig4A_repro", "Fig4B_NMI"], ascending=[False, False]).iloc[0]
    else:
        raise ValueError(f"unknown rule: {rule}")
    return {
        "Method": method,
        "preference": rule,
        "subset": row["subset"],
        "Fig4A_repro": row["Fig4A_repro"],
        "Fig4B_NMI": row["Fig4B_NMI"],
        "n_subsets_available": len(df),
    }


def main():
    rows = [choose(method, rule) for method, rule in PREFERENCES.items()]
    out = pd.DataFrame(rows)
    out.to_csv(PICK3 / "preferred_triplets.csv", index=False)
    print(out.to_string(index=False))
    print(f"wrote {PICK3 / 'preferred_triplets.csv'}")


if __name__ == "__main__":
    main()
