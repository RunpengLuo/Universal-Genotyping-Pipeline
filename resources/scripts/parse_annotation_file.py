"""Reformat a barcode/label Visium annotation file.

Last update: 2026-07-26

Inputs:
- argv: ann_file, out_file
- ann_file: barcode/label annotation table
Outputs:
- out_file: the same table with simplified labels
"""

import sys

import pandas as pd


_, ann_file, out_file = sys.argv


def simplify_label(v):
    """Simplify a Visium annotation label string.

    Labels starting with ``"T"`` without an underscore become ``"tumor"``.
    Labels starting with ``"T"`` with an underscore keep the part after the
    first underscore. All other labels are returned unchanged.

    Parameters
    ----------
    v : str
        Raw annotation label.

    Returns
    -------
    str
        Simplified label.
    """
    if v[0] == "T":
        if "_" not in v:
            return "tumor"
        else:
            return v[str(v).find("_") + 1 :]
    else:
        return v


# change accordingly
barcode_col = "Barcode"
raw_label = "Microregion_annotation"
target_label = "path_label"

anns = pd.read_table(ann_file, sep="\t", keep_default_na=True)
print(anns.head())

anns = anns.rename(columns={barcode_col: "BARCODE"})

# NA labels are non-tumor
anns[raw_label] = anns[raw_label].fillna("normal")
anns[target_label] = anns.apply(
    func=lambda r: simplify_label(r[raw_label]), axis=1
).astype("str")

print(anns.head())

anns.to_csv(out_file, sep="\t", header=True, index=False)
