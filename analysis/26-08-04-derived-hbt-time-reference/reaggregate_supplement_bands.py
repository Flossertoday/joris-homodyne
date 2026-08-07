"""Re-aggregate the supplement's subwindow traces into band traces with the
supplement-specific band windows, replacing the (statistics-starved) auto band
detection of the fig2a pipeline.

Windows were localized on the 72000-frame aligned subwindow profiles (analog and
photon agree): B1 350-440, B2 540-680, B3 770-860, B4 970-1100 — same four-band
geometry as the original batch (centers 395/610/815/1035 vs 415/625/845/1060),
and arm-resolved L/R ratios support the unchanged pairing L = B1+B3, R = B2+B4.
Writes a new npz with the original schema plus a `reaggregation` metadata key;
downstream render/check scripts consume it unchanged (--bands-l 0,2
--bands-r 1,3). Run from the repository/worktree root; the script's own
directory supplies extract_hbt_traces on sys.path.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from extract_hbt_traces import aggregate_band_traces, band_correlation

WINDOWS = np.array([[350, 440], [540, 680], [770, 860], [970, 1100]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", type=Path,
                        default=Path("drafts/26-08-06_fig2a_supplement/hbt_traces.npz"))
    parser.add_argument("--dst", type=Path, default=None,
                        help="default: <src stem>_bands.npz next to --src")
    args = parser.parse_args()
    dst = args.dst or args.src.with_name(f"{args.src.stem}_bands.npz")

    d = np.load(args.src)
    out = {key: d[key] for key in d.files}

    edges = d["subwindow_edges"]
    for prefix, mean_key, sem_key in (
        ("band_traces_analog", "traces_analog", "sem_analog"),
        ("band_traces_photon", "traces_photon", "sem_photon"),
    ):
        means, sems = aggregate_band_traces(d[mean_key], d[sem_key], WINDOWS, edges)
        out[prefix] = means
        out[prefix.replace("traces", "sem")] = sems
    out["band_windows"] = WINDOWS
    corr, _ = band_correlation(out["band_traces_analog"])
    out["band_corr"] = corr

    meta = json.loads(str(d["meta_json"]))
    meta["reaggregation"] = {
        "script": Path(__file__).name,
        "band_windows": WINDOWS.tolist(),
        "reason": "auto band detection is statistics-limited on the supplement "
                  "(lumpy off-band glow); windows localized on the 72000-frame "
                  "aligned subwindow profiles; pairing L={B1,B3} R={B2,B4} "
                  "supported by arm-resolved L/R ratios",
        "pairing": {"L": [0, 2], "R": [1, 3]},
    }
    out["meta_json"] = np.array(json.dumps(meta))

    np.savez(dst, **out)
    print(f"wrote {dst}")
    print("band windows:", WINDOWS.tolist())
    print("band corr (analog):", np.round(corr, 3).tolist())

    t = d["time_ps"]
    win = (t > 20) & (t < 800)
    sums = {src: out[f"band_traces_{src}"][:, win].sum(axis=1) for src in ("analog", "photon")}
    for src, s in sums.items():
        print(f"{src} window sums per frame: B1 {s[0]:.1f} B2 {s[1]:.1f} B3 {s[2]:.1f} B4 {s[3]:.1f}")
        print(f"  L=B1+B3 / R=B2+B4 = {(s[0] + s[2]) / (s[1] + s[3]):.3f}; "
              f"arm A B1/B2 = {s[0] / s[1]:.3f}, arm B B3/B4 = {s[2] / s[3]:.3f}")


if __name__ == "__main__":
    main()
