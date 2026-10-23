#!/usr/bin/env python3
"""
Plot a 2D lookup table made by LookupTableGenerator (e.g. via TestMakeALookupTable).

The figure has two panels:
 * left: the lookup table's boxes, each drawn as a surface patch through its four corner points,
   i.e. the values the table interpolates between;
 * right: the evaluated points, with any that returned an error code marked in red.

The box tree isn't stored in the .dat file, so it is reconstructed from the points: starting from
[0,1]^2, a box was split along a dimension if both points half way along that dimension (on its edges)
were evaluated, otherwise it is a leaf. This works for both bisection and the older 2^DIM subdivision.

The channels, model and pacing frequency are read from the table's file name, which TestMakeALookupTable
writes as <model>_<N>d_<channel>_<channel>_<freq>Hz.dat, e.g. tentusscher_model_2006_epi_2d_hERG_INa_1Hz.dat

Usage:
    python3 plot_lookup_table_2d.py <table .dat> [-o out.png] [-q <QoI index>]

Requires numpy and matplotlib.
"""
import argparse
import os
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

# The QuantityOfInterest enumeration in src/lookup/QuantityOfInterest.hpp, with units.
QOI_NAMES = {0: ("APD50", "ms"), 1: ("APD90", "ms"), 2: ("Upstroke velocity", "mV/ms"),
             3: ("Peak voltage", "mV"), 4: ("qNet", "C/F")}


def parse_file_name(path):
    """Return (model, [channel names], frequency) from a lookup table's file name."""
    name = os.path.splitext(os.path.basename(path))[0]
    match = re.fullmatch(r"(?P<model>.+)_(?P<dim>\d+)d_(?P<channels>.+)_(?P<hz>[\d.]+)Hz", name)
    if not match:
        sys.exit(f"Can't read the model, channels and frequency from '{name}', "
                 "expected <model>_<N>d_<channels>_<freq>Hz")
    channels = match["channels"].split("_")
    if int(match["dim"]) != 2 or len(channels) != 2:
        sys.exit(f"'{name}' is a {match['dim']}D table of {channels}, this script only plots 2D tables")
    return match["model"], channels, match["hz"]


def read_table(path, qoi_index):
    """
    Read a lookup table .dat file.

    The header line is: number of parameters, number of QoIs, parameter names, QoI enumeration values.
    Each following line is: parameter values, error code, QoIs, and optionally the number of error
    estimates, then for each QoI its error estimate and (on most lines) the % of space meeting tolerance.

    Returns (points (N x 2), error codes, QoI values, QoI enum value, last reported % meeting tolerance).
    """
    with open(path) as f:
        header = f.readline().split()
        num_params, num_qois = int(header[0]), int(header[1])
        if num_params != 2:
            sys.exit(f"{path} has {num_params} parameters, this script only plots 2D tables")
        if not 0 <= qoi_index < num_qois:
            sys.exit(f"QoI index {qoi_index} out of range, this table has {num_qois} QoI(s)")
        qoi_enum = int(header[2 + num_params + qoi_index])

        points, codes, values, percent_met = [], [], [], None
        percent_column = 2 + 1 + num_qois + 1 + 2 * qoi_index + 1
        for line in f:
            parts = line.split()
            if not parts:
                continue
            points.append((float(parts[0]), float(parts[1])))
            codes.append(int(parts[2]))
            values.append(float(parts[3 + qoi_index]))
            if len(parts) > percent_column:
                percent_met = float(parts[percent_column])
    return np.array(points), np.array(codes), np.array(values), qoi_enum, percent_met


def leaf_boxes(value, xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0, out=None):
    """Reconstruct the leaf boxes, splitting wherever both points on a midline were evaluated."""
    if out is None:
        out = []
    xm, ym = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
    if key(xm, ymin) in value and key(xm, ymax) in value:
        leaf_boxes(value, xmin, xm, ymin, ymax, out)
        leaf_boxes(value, xm, xmax, ymin, ymax, out)
    elif key(xmin, ym) in value and key(xmax, ym) in value:
        leaf_boxes(value, xmin, xmax, ymin, ym, out)
        leaf_boxes(value, xmin, xmax, ym, ymax, out)
    else:
        out.append((xmin, xmax, ymin, ymax))
    return out


def key(x, y):
    """Points are dyadic fractions, so rounding gives a reliable dictionary key."""
    return (round(x, 12), round(y, 12))


def main():
    parser = argparse.ArgumentParser(description="Plot a 2D ApPredict lookup table (.dat file).")
    parser.add_argument("table", help="the lookup table .dat file")
    parser.add_argument("-o", "--output", help="the image to write (default: <table name>.png here)")
    parser.add_argument("-q", "--qoi", type=int, default=0,
                        help="which quantity of interest to plot, if the table has several (default 0)")
    args = parser.parse_args()

    model, (xlabel, ylabel), hz = parse_file_name(args.table)
    points, codes, values, qoi_enum, percent_met = read_table(args.table, args.qoi)
    qoi_name, qoi_units = QOI_NAMES.get(qoi_enum, (f"QoI {qoi_enum}", ""))
    zlabel = f"{qoi_name} ({qoi_units})" if qoi_units else qoi_name
    out_png = args.output or os.path.splitext(os.path.basename(args.table))[0] + ".png"

    value = {key(x, y): v for (x, y), v in zip(points, values)}
    boxes = leaf_boxes(value)
    quads = [[(x0, y0, value[key(x0, y0)]), (x1, y0, value[key(x1, y0)]),
              (x1, y1, value[key(x1, y1)]), (x0, y1, value[key(x0, y1)])] for x0, x1, y0, y1 in boxes]

    ok = codes == 0
    zmin, zmax = values.min(), values.max()  # Include error-coded values, as boxes pass through them.
    pad = 0.05 * (zmax - zmin) if zmax > zmin else 1.0

    fig = plt.figure(figsize=(16, 7.5))
    cmap = plt.get_cmap("jet")
    for col in (1, 2):
        ax = fig.add_subplot(1, 2, col, projection="3d")
        if col == 1:
            colours = [cmap((np.mean([c[2] for c in q]) - zmin) / ((zmax - zmin) or 1.0)) for q in quads]
            ax.add_collection3d(Poly3DCollection(quads, facecolors=colours, edgecolors="k", linewidths=0.3))
            ax.set_title(f"Lookup table boxes ({len(boxes)}), each through its corner points")
        else:
            ax.scatter(points[ok, 0], points[ok, 1], values[ok], c="b", s=8, marker="*",
                       depthshade=False, label="evaluated")
            if (~ok).any():
                ax.scatter(points[~ok, 0], points[~ok, 1], values[~ok], c="r", s=10, marker="x",
                           depthshade=False, label="error code")
            ax.legend(loc="upper left")
            pct = "" if percent_met is None else f", {percent_met:.1f}% of space meets tolerance"
            ax.set_title(f"{len(values)} evaluations ({int((~ok).sum())} error-coded){pct}")
        ax.view_init(elev=28, azim=45)  # Puts (1, 1), i.e. no block, at the front corner.
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_zlim(zmin - pad, zmax + pad)

    fig.suptitle(f"{model}, {hz} Hz: {qoi_name} lookup table for {xlabel} x {ylabel} block")
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.4 / fig.get_figheight()))
    fig.savefig(out_png, dpi=110)
    print(f"Wrote {out_png}: {len(values)} points ({int((~ok).sum())} error-coded) in {len(boxes)} boxes, "
          f"{qoi_name} {values.min():.4g} to {values.max():.4g}")


if __name__ == "__main__":
    main()
