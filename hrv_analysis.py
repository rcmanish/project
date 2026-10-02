"""Reproducible PPG pulse-rate variability analysis for the BP35 study.

Run: python hrv_analysis.py
Requires: numpy, pandas, scipy, matplotlib.
All intervals are PPG peak-to-peak intervals (pulse-rate variability), not ECG RR.
"""

from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.signal import butter, find_peaks, sosfiltfilt, welch


ROOT = Path(__file__).resolve().parent
PHASES = {"Pre": ("pre", 300), "During": ("during", 600), "Post": ("post", 300)}
FS_PSD = 4.0
COLORS = {"Pre": "#246A73", "During": "#C96B3C", "Post": "#6B5B95"}


def extract_archive(phase: str, folder: str) -> Path:
    """Unpack one matching archive, checking paths before extraction."""
    matches = sorted(ROOT.glob(f"BP35_{phase}-*.zip"))
    if len(matches) != 1:
        raise ValueError(f"Expected one {phase} archive; found {len(matches)}")
    target = ROOT / folder
    target.mkdir(exist_ok=True)
    with zipfile.ZipFile(matches[0]) as archive:
        for member in archive.infolist():
            # Archives include macOS metadata and temporary Office lock files.
            if member.is_dir() or Path(member.filename).name.startswith((".", "~$")):
                continue
            if Path(member.filename).suffix.lower() not in {".csv", ".json"}:
                continue
            destination = (target / Path(member.filename).name).resolve()
            if destination.parent != target.resolve():
                raise ValueError(f"Unsafe archive member: {member.filename}")
            destination.write_bytes(archive.read(member))
    return target


def load_ppg(folder: Path) -> tuple[np.ndarray, np.ndarray, float, dict]:
    """Use device milliseconds for relative timing and JSON nominal PPG rate."""
    info_file, = folder.glob("*_info.json")
    metadata = json.loads(info_file.read_text(encoding="utf-8-sig"))
    ppg_info = next(row["info"] for row in metadata if row["info"].get("name") == "PPG")
    recording_info = next(row["info"] for row in metadata if row["info"].get("name") == "EmotiBitData")
    fs = float(ppg_info["nominal_srate"])
    pg_file, = folder.glob("*_PG.csv")
    df = pd.read_csv(pg_file, usecols=["EmotiBitTimestamp", "PG"])
    df = df.dropna().groupby("EmotiBitTimestamp", as_index=False)["PG"].mean()
    timestamp = df["EmotiBitTimestamp"].to_numpy(dtype=float) / 1000.0
    timestamp -= timestamp[0]
    return timestamp, df["PG"].to_numpy(dtype=float), fs, {
        "ppg": ppg_info, "recording": recording_info}


def band_power(freq: np.ndarray, psd: np.ndarray, low: float, high: float) -> float:
    """Trapezoidal integral, interpolating exact band limits."""
    internal = (freq > low) & (freq < high)
    x = np.r_[low, freq[internal], high]
    y = np.interp(x, freq, psd)
    return float(np.trapezoid(y, x))


def analyze_phase(name: str, folder: Path, requested_s: int) -> dict:
    timestamp, signal, fs, metadata = load_ppg(folder)
    span = float(timestamp[-1])
    # Start at the first recorded PPG sample. A short capture cannot be extended.
    end = min(float(requested_s), span)
    mask = timestamp <= end
    timestamp, signal = timestamp[mask], signal[mask]

    # Timestamps jitter slightly around 25 Hz. Interpolate to the nominal grid
    # before zero-phase filtering; retain real elapsed time as the grid axis.
    grid = np.arange(0, end, 1 / fs)
    uniform = np.interp(grid, timestamp, signal)
    sos = butter(4, [0.5, 3.5], btype="bandpass", fs=fs, output="sos")
    filtered = sosfiltfilt(sos, uniform)

    # A 300 ms refractory period sets the upper detectable rate at 200 bpm.
    # Positive height rejects small negative secondary waves seen in the Pre trace;
    # prominence scales to the recording's filtered signal amplitude.
    prominence = 0.35 * float(np.std(filtered))
    peaks, _ = find_peaks(filtered, distance=int(np.ceil(0.3 * fs)),
                          prominence=prominence, height=0)
    beat_t = grid[peaks]
    raw_ppi = np.diff(beat_t) * 1000.0
    interval_t = beat_t[1:]
    if len(raw_ppi) < 20:
        raise ValueError(f"{name}: too few detected intervals")

    # Mark implausible intervals and those >20% from a centered local median.
    # The local reference uses a broad 11-beat window, robust to isolated peaks.
    local = pd.Series(raw_ppi).rolling(11, center=True, min_periods=3).median().to_numpy()
    good = (raw_ppi >= 300) & (raw_ppi <= 1500) & (np.abs(raw_ppi - local) <= 0.2 * local)
    if good.sum() < 20:
        raise ValueError(f"{name}: too few clean intervals after artifact screening")
    # Cubic interpolation on elapsed beat time repairs only flagged intervals.
    # Keep all originally valid values unchanged for time-domain statistics.
    corrected = raw_ppi.copy()
    corrected[~good] = CubicSpline(interval_t[good], raw_ppi[good], extrapolate=True)(interval_t[~good])
    corrected = np.clip(corrected, 300, 1500)

    sdnn = float(np.std(corrected, ddof=1))
    rmssd = float(np.sqrt(np.mean(np.diff(corrected) ** 2)))

    # A 4 Hz uniformly sampled tachogram, then mean-centered Welch PSD.
    # Hanning is the conventional name for NumPy/SciPy's Hann window.
    regular_t = np.arange(interval_t[0], interval_t[-1], 1 / FS_PSD)
    tachogram = CubicSpline(interval_t, corrected)(regular_t)
    # Approximately 300 s windows resolve the 0.04 Hz LF boundary better than
    # short windows; for the 10-minute phase, Welch averages overlapping windows.
    segment_samples = min(int(300 * FS_PSD), len(tachogram))
    freq, psd = welch(tachogram - tachogram.mean(), fs=FS_PSD, window="hann",
                      nperseg=segment_samples, noverlap=segment_samples // 2,
                      detrend="constant", scaling="density")
    lf = band_power(freq, psd, 0.04, 0.15)
    hf = band_power(freq, psd, 0.15, 0.40)
    return {
        "name": name, "requested_s": requested_s, "recorded_span_s": span,
        "analyzed_span_s": end, "sampling_hz": fs,
        "device_id": metadata["recording"].get("device_id"),
        "firmware_version": metadata["recording"].get("firmware_version"),
        "recording_created_at": metadata["recording"].get("created_at"),
        "samples": len(signal), "detected_peaks": len(peaks),
        "intervals": len(raw_ppi), "corrected_intervals": int((~good).sum()),
        "artifact_percent": 100 * float((~good).mean()),
        "SDNN_ms": sdnn, "RMSSD_ms": rmssd, "LF_ms2": lf, "HF_ms2": hf,
        "interval_t": interval_t, "raw_ppi": raw_ppi,
        "corrected_ppi": corrected, "artifact": ~good,
        "freq": freq, "psd": psd, "metadata": metadata,
    }


def make_figure(results: list[dict]) -> None:
    """Save one vector and one 300 DPI raster four-panel figure."""
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.facecolor": "white"})
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), layout="constrained")
    ax = axes[0, 0]
    offset = 0.0
    for result in results:
        ax.plot((result["interval_t"] + offset) / 60, result["corrected_ppi"],
                lw=0.75, color=COLORS[result["name"]], label=result["name"])
        offset += result["requested_s"]
        if offset < 1200:
            ax.axvline(offset / 60, ls="--", lw=0.9, color="#404040")
    # These task lines are nominal only; no task event markers are present.
    ax.axvline(8, ls=":", lw=0.8, color="#888888")
    ax.axvline(12, ls=":", lw=0.8, color="#888888")
    ax.axvspan(8, 12, color="#888888", alpha=0.055)
    ax.text(10, 0.96, "Stroop (nominal)", ha="center", va="top",
            transform=ax.get_xaxis_transform(), fontsize=8, color="#666666")
    ax.set(xlabel="Protocol time (min)", ylabel="Corrected pulse interval (ms)",
           title="A  Beat-to-beat pulse intervals", xlim=(0, 20))
    ax.legend(frameon=False, ncol=3, loc="upper right")

    names = [r["name"] for r in results]
    x = np.arange(3)
    width = 0.35
    ax = axes[0, 1]
    ax.bar(x - width / 2, [r["SDNN_ms"] for r in results], width, label="SDNN", color="#246A73")
    ax.bar(x + width / 2, [r["RMSSD_ms"] for r in results], width, label="RMSSD", color="#C96B3C")
    ax.set(xticks=x, xticklabels=names, ylabel="Variability (ms)", title="B  Time-domain variability")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.axvspan(0.04, 0.15, color="#4A90A4", alpha=0.14, label="LF band")
    ax.axvspan(0.15, 0.40, color="#E8A46B", alpha=0.16, label="HF band")
    for result in results:
        m = result["freq"] <= 0.5
        ax.plot(result["freq"][m], result["psd"][m], lw=1.5,
                color=COLORS[result["name"]], label=result["name"])
    ax.set(xlabel="Frequency (Hz)", ylabel="PSD (ms²/Hz)",
           title="C  Pulse-interval power spectrum", xlim=(0, 0.5))
    ax.legend(frameon=False, ncol=3, fontsize=8)

    ax = axes[1, 1]
    ax.bar(x - width / 2, [r["LF_ms2"] for r in results], width, label="LF (0.04–0.15 Hz)", color="#4A90A4")
    ax.bar(x + width / 2, [r["HF_ms2"] for r in results], width, label="HF (0.15–0.40 Hz)", color="#E8A46B")
    ax.set(xticks=x, xticklabels=names, ylabel="Power (ms²)", title="D  Frequency-domain power")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(ROOT / "hrv_four_panel.png", dpi=300)
    fig.savefig(ROOT / "hrv_four_panel.pdf")
    plt.close(fig)


def markdown_table(frame: pd.DataFrame, signed: bool = False, index: bool = True) -> str:
    """Format a compact Markdown table without pandas' optional tabulate package."""
    display = frame.reset_index() if index else frame
    headers = list(display.columns)
    rows = []
    for row in display.itertuples(index=False, name=None):
        cells = []
        for value in row:
            if isinstance(value, (float, np.floating)):
                cells.append(f"{value:+.2f}" if signed else f"{value:.2f}")
            else:
                cells.append(str(value))
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join("---" for _ in headers) + " |", *rows])


def main() -> None:
    # Windows consoles may default to cp1252, which cannot print arrows or ms².
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    results = []
    for phase, (folder, requested_s) in PHASES.items():
        directory = extract_archive(phase, folder)
        results.append(analyze_phase(phase, directory, requested_s))

    metrics = ["SDNN_ms", "RMSSD_ms", "LF_ms2", "HF_ms2"]
    # One tidy CSV contains raw estimates and both requested percent changes.
    table = pd.DataFrame([{k: r[k] for k in ["name", "requested_s", "recorded_span_s",
                 "analyzed_span_s", "sampling_hz", "device_id", "firmware_version",
                 "recording_created_at", "samples", "detected_peaks", "intervals",
                 "corrected_intervals", "artifact_percent", *metrics]} for r in results])
    table = table.rename(columns={"name": "phase"})
    for index, r in enumerate(results):
        for metric in metrics:
            table.loc[index, f"{metric}_change_from_previous_pct"] = (
                np.nan if index == 0 else 100 * (r[metric] - results[index - 1][metric]) / results[index - 1][metric]
            )
    table.to_csv(ROOT / "hrv_analysis_results.csv", index=False, float_format="%.6f")

    beats = pd.concat([pd.DataFrame({"phase": r["name"], "time_s": r["interval_t"],
                        "nominal_subphase": np.where(r["interval_t"] < 180, "Listening 1",
                            np.where(r["interval_t"] < 420, "Stroop", "Listening 2"))
                            if r["name"] == "During" else r["name"],
                        "raw_ppi_ms": r["raw_ppi"], "corrected_ppi_ms": r["corrected_ppi"],
                        "artifact_corrected": r["artifact"]}) for r in results], ignore_index=True)
    beats.to_csv(ROOT / "hrv_beat_intervals.csv", index=False, float_format="%.6f")
    make_figure(results)

    raw = table.set_index("phase")[metrics].rename(columns={
        "SDNN_ms": "SDNN (ms)", "RMSSD_ms": "RMSSD (ms)",
        "LF_ms2": "LF (ms²)", "HF_ms2": "HF (ms²)"})
    changes = pd.DataFrame({
        "Pre → During (%)": 100 * (raw.loc["During"] - raw.loc["Pre"]) / raw.loc["Pre"],
        "During → Post (%)": 100 * (raw.loc["Post"] - raw.loc["During"]) / raw.loc["During"],
    })
    raw_md = markdown_table(raw)
    changes_md = markdown_table(changes, signed=True)
    during_beats = beats[beats["phase"] == "During"]
    subphase_table = pd.DataFrame([
        {"Nominal activity": label, "Start (s)": start, "End (s)": stop,
         "Intervals": int(((during_beats["time_s"] >= start) &
                           (during_beats["time_s"] < stop)).sum())}
        for label, start, stop in [("Listening 1", 0, 180),
                                   ("Task while listening", 180, 420),
                                   ("Listening 2", 420, 600)]
    ])
    subphase_md = markdown_table(subphase_table, index=False)
    quality_md = markdown_table(table[["phase", "recorded_span_s", "analyzed_span_s",
                                       "detected_peaks", "corrected_intervals", "artifact_percent"]], index=False)
    print("\nRAW METRICS\n", raw_md)
    print("\nPERCENT CHANGES\n", changes_md)
    print("\nQUALITY\n", quality_md)
    lines = ["# BP35 PPG pulse-rate variability results", "", "## Raw metrics", "",
             raw_md, "", "## Percentage changes", "", changes_md, "",
             "## During-intervention segmentation", "", subphase_md, "",
             "## Data quality and interpretation", "", quality_md,
             "", "* All three archives report device MD-V5-0000471, firmware 1.14.0, and a nominal PPG sampling rate of 25 Hz.",
             "", "* Pre and Post recordings are shorter than 300 s; no samples were fabricated.",
             "* During uses the first 600 s of a 706.7 s recording. The nominal task boundaries at 180 and 420 s are not verified by event markers.",
             "* Metrics are PPG pulse-rate variability estimates, not ECG-derived HRV. Interpret comparisons cautiously, especially because During is twice as long as Pre and Post.",
             "* The 25 Hz PPG sampling grid quantizes peak times in approximately 40 ms steps, limiting interval precision, especially for RMSSD.",
             "* Artifact correction uses a 300–1500 ms physiologic range and a 20% centered 11-interval local-median threshold, followed by cubic interpolation.",
             "* PSD uses 4 Hz cubic-spline interpolation and approximately 300 s Hann-window Welch segments with 50% overlap; LF and HF use exact-edge trapezoidal integration.",
             "", "## Method references", "",
             "* [Kubios HRV Scientific User's Guide](https://www.kubios.com/downloads/HRV-Scientific-Users-Guide.pdf): standard 4 Hz interpolation, LF/HF bands, and Welch spectral settings.",
             "* [Task Force HRV standards (1996)](https://www.escardio.org/static-file/Escardio/Guidelines/Scientific-Statements/guidelines-Heart-Rate-Variability-FT-1996.pdf): conventional definitions of short-term HRV bands and metrics.",
             ""]
    (ROOT / "hrv_analysis_report.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
