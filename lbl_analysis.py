#!/usr/bin/env python3
r"""
Wendelstein LBL Analysis Toolkit
================================

A lightweight interactive tool for inspecting one LBL radial-velocity RDB
file or every RDB file in a directory.

Expected LBL columns
--------------------
rjd, vrad, svrad

Available analyses
------------------
1. Relative-RV time series
2. Lomb-Scargle periodogram
3. Phase-folded RV curve

Usage
-----
    python lbl_analysis_v1_1.py

The user may enter either:

    ...\lblrdb\lbl_TOI5786_TOI5786.rdb

or a directory:

    ...\lblrdb\

When a directory is supplied, every ``*.rdb`` file directly inside that
directory is processed.

Dependencies
------------
numpy
pandas
matplotlib
astropy
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle


REQUIRED_COLUMNS = ("rjd", "vrad", "svrad")
DEFAULT_MIN_PERIOD = 1.0
DEFAULT_MAX_PERIOD = 100.0
SAMPLES_PER_PEAK = 10

# Used only for the strongest periodogram peak so it remains clearly visible
# against both the periodogram curve and the light figure background.
PEAK_COLOR = "crimson"


##############################################################
# Input helpers
##############################################################


def ask_float(
    prompt: str,
    default: float | None = None,
    allow_skip: bool = False,
) -> float | None:
    """Request a positive floating-point number from the user."""
    while True:
        default_text = f" [{default:g}]" if default is not None else ""
        skip_text = " (or s to skip)" if allow_skip else ""
        value = input(f"{prompt}{default_text}{skip_text}: ").strip()

        if allow_skip and value.lower() in {"s", "skip"}:
            return None

        if not value and default is not None:
            return float(default)

        try:
            number = float(value)
        except ValueError:
            print("Please enter a valid number.")
            continue

        if not np.isfinite(number) or number <= 0:
            print("Please enter a positive finite number.")
            continue

        return number


def ask_yes_no(prompt: str, default: bool = False) -> bool:
    """Request a yes/no answer from the user."""
    suffix = " [Y/n]" if default else " [y/N]"

    while True:
        answer = input(f"{prompt}{suffix}: ").strip().lower()

        if not answer:
            return default
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False

        print("Please answer y or n.")


##############################################################
# File discovery and RDB loading
##############################################################


def discover_rdb_files(input_path: str | Path) -> list[Path]:
    """
    Return one RDB file or all RDB files directly inside a directory.

    Directory searching is intentionally non-recursive so that unrelated
    files in nested folders are not processed accidentally.
    """
    path = Path(input_path).expanduser()

    if path.is_file():
        if path.suffix.lower() != ".rdb":
            raise ValueError(f"The selected file is not an RDB file: {path}")
        return [path]

    if path.is_dir():
        files = sorted(
            candidate
            for candidate in path.iterdir()
            if candidate.is_file() and candidate.suffix.lower() == ".rdb"
        )

        if not files:
            raise FileNotFoundError(
                f"No .rdb files were found directly inside: {path}"
            )

        return files

    raise FileNotFoundError(f"File or directory not found: {path}")


def load_rdb(filename: str | Path) -> pd.DataFrame:
    """
    Load an LBL RDB file and return clean RV measurements.

    Non-numeric rows, including the standard dashed RDB separator row, are
    removed automatically.
    """
    path = Path(filename).expanduser()

    try:
        table = pd.read_csv(
            path,
            sep=r"\s+",
            comment="#",
            engine="python",
        )
    except Exception as exc:
        raise ValueError(f"Could not read the RDB file: {exc}") from exc

    table.columns = [str(column).strip().lower() for column in table.columns]

    missing = [column for column in REQUIRED_COLUMNS if column not in table.columns]
    if missing:
        available = ", ".join(table.columns)
        raise ValueError(
            "The RDB file does not contain the required columns "
            f"{missing}. Available columns: {available}"
        )

    data = table.loc[:, REQUIRED_COLUMNS].copy()

    for column in REQUIRED_COLUMNS:
        data[column] = pd.to_numeric(data[column], errors="coerce")

    data = data.replace([np.inf, -np.inf], np.nan)
    data = data.dropna(subset=list(REQUIRED_COLUMNS))
    data = data[data["svrad"] > 0]
    data = data.sort_values("rjd").reset_index(drop=True)

    if data.empty:
        raise ValueError(
            "No valid measurements remained after removing non-numeric rows, "
            "missing values, and non-positive uncertainties."
        )

    return data


def infer_target_name(filename: str | Path) -> str:
    """
    Infer a target name from a typical LBL filename.

    Example
    -------
    lbl_TOI5786_TOI5786.rdb -> TOI5786
    """
    stem = Path(filename).stem
    cleaned = re.sub(r"^lbl2?_", "", stem, flags=re.IGNORECASE)
    parts = [part for part in cleaned.split("_") if part]

    if not parts:
        return stem

    return parts[0]


def relative_rv(data: pd.DataFrame) -> np.ndarray:
    """Return RV measurements after subtracting their arithmetic mean."""
    rv = data["vrad"].to_numpy(dtype=float)
    return rv - np.mean(rv)


def print_data_summary(data: pd.DataFrame, target: str) -> None:
    """Print a short summary for one loaded target."""
    time = data["rjd"].to_numpy(dtype=float)
    rv = data["vrad"].to_numpy(dtype=float)
    err = data["svrad"].to_numpy(dtype=float)
    rel_rv = rv - np.mean(rv)

    print(f"\nTarget              : {target}")
    print(f"Observations        : {len(data)}")
    print(f"Time span           : {np.ptp(time):.3f} days")
    print(f"Mean RV             : {np.mean(rv):.3f} m/s")
    print(f"Relative-RV RMS     : {np.sqrt(np.mean(rel_rv**2)):.3f} m/s")
    print(f"Median uncertainty  : {np.median(err):.3f} m/s")


##############################################################
# Plotting helpers
##############################################################


def configure_plotting() -> None:
    """Set a clean, publication-oriented Matplotlib configuration."""
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.labelsize": 13,
            "axes.titlesize": 14,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 9,
            "axes.linewidth": 1.0,
            "lines.linewidth": 1.3,
            "savefig.dpi": 300,
        }
    )


def setup_figure() -> tuple[plt.Figure, plt.Axes]:
    """Create a consistently formatted figure."""
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.grid(alpha=0.22)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")
    return fig, ax


def safe_filename(text: str) -> str:
    """Convert text into a safe output filename component."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")


def save_figure(
    fig: plt.Figure,
    source_file: Path,
    target: str,
    plot_name: str,
) -> tuple[Path, Path]:
    """Save PNG and PDF copies next to the source RDB file."""
    base_name = f"{safe_filename(target)}_{plot_name}"
    png_path = source_file.parent / f"{base_name}.png"
    pdf_path = source_file.parent / f"{base_name}.pdf"

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")

    return png_path, pdf_path


def finish_figure(
    fig: plt.Figure,
    source_file: Path,
    target: str,
    plot_name: str,
    save: bool,
    batch_mode: bool,
) -> None:
    """Save or display a completed figure."""
    if save:
        png_path, pdf_path = save_figure(
            fig,
            source_file,
            target,
            plot_name,
        )
        print(f"Saved: {png_path}")
        print(f"Saved: {pdf_path}")
        plt.close(fig)
        return

    if batch_mode:
        plt.show(block=False)
        plt.pause(0.1)
    else:
        plt.show()


##############################################################
# RV time series
##############################################################


def plot_time_series(
    data: pd.DataFrame,
    target: str,
) -> plt.Figure:
    """Create a relative-RV time-series figure."""
    time = data["rjd"].to_numpy(dtype=float)
    rv = relative_rv(data)
    err = data["svrad"].to_numpy(dtype=float)

    fig, ax = setup_figure()

    ax.errorbar(
        time,
        rv,
        yerr=err,
        fmt="o",
        markersize=6,
        capsize=3,
        elinewidth=1.1,
    )

    ax.axhline(0.0, linewidth=1.0, alpha=0.5)
    ax.set_xlabel("RJD")
    ax.set_ylabel("Relative RV (m/s)")
    ax.set_title(target)

    fig.tight_layout()
    print("Mean RV removed from the displayed measurements.")

    return fig


##############################################################
# Lomb-Scargle periodogram
##############################################################


def plot_periodogram(
    data: pd.DataFrame,
    target: str,
    min_period: float,
    max_period: float,
) -> plt.Figure:
    """Create a weighted Lomb-Scargle periodogram figure."""
    if min_period >= max_period:
        raise ValueError("Minimum period must be smaller than maximum period.")

    time = data["rjd"].to_numpy(dtype=float)
    rv = relative_rv(data)
    err = data["svrad"].to_numpy(dtype=float)

    if len(time) < 3:
        raise ValueError("At least three observations are required.")

    lomb_scargle = LombScargle(
        time,
        rv,
        dy=err,
        fit_mean=True,
        center_data=True,
    )

    frequency, power = lomb_scargle.autopower(
        minimum_frequency=1.0 / max_period,
        maximum_frequency=1.0 / min_period,
        samples_per_peak=SAMPLES_PER_PEAK,
    )

    if frequency.size == 0:
        raise ValueError("The selected period range produced no frequency grid.")

    period = 1.0 / frequency
    best_index = int(np.nanargmax(power))
    best_period = float(period[best_index])
    best_power = float(power[best_index])

    try:
        best_fap = float(
            lomb_scargle.false_alarm_probability(
                best_power,
                minimum_frequency=1.0 / max_period,
                maximum_frequency=1.0 / min_period,
            )
        )
    except Exception:
        best_fap = float("nan")

    probabilities = np.array([0.10, 0.01, 0.001])

    try:
        fap_levels = lomb_scargle.false_alarm_level(
            probabilities,
            minimum_frequency=1.0 / max_period,
            maximum_frequency=1.0 / min_period,
        )
    except Exception:
        fap_levels = np.full(probabilities.shape, np.nan)

    order = np.argsort(period)
    period_sorted = period[order]
    power_sorted = power[order]

    fig, ax = setup_figure()

    ax.plot(period_sorted, power_sorted, linewidth=1.3)

    # The strongest peak uses a contrasting color and both a line and point.
    ax.axvline(
        best_period,
        color=PEAK_COLOR,
        linestyle="--",
        linewidth=1.6,
        label=f"Best period: {best_period:.4f} d",
        zorder=4,
    )
    ax.scatter(
        [best_period],
        [best_power],
        color=PEAK_COLOR,
        marker="v",
        s=65,
        zorder=5,
    )

    for probability, level in zip(probabilities, fap_levels):
        if np.isfinite(level):
            label = f"{100 * probability:g}% FAP"
            ax.axhline(level, linestyle=":", linewidth=1.0, label=label)

    ax.set_xlim(min_period, max_period)
    ax.set_xlabel("Period (days)")
    ax.set_ylabel("Lomb–Scargle Power")
    ax.set_title(target)
    ax.legend()

    fig.tight_layout()

    print("Strongest periodogram peak")
    print(f"Period : {best_period:.6f} days")
    print(f"Power  : {best_power:.6f}")
    if np.isfinite(best_fap):
        print(f"FAP    : {best_fap:.6g}")
    else:
        print("FAP    : unavailable")

    return fig


##############################################################
# Phase-folded RV curve
##############################################################


def plot_phase_fold(
    data: pd.DataFrame,
    target: str,
    period: float,
) -> plt.Figure:
    """Create a relative-RV phase-folded figure over two cycles."""
    time = data["rjd"].to_numpy(dtype=float)
    rv = relative_rv(data)
    err = data["svrad"].to_numpy(dtype=float)

    reference_time = np.min(time)
    phase = ((time - reference_time) / period) % 1.0

    order = np.argsort(phase)
    phase = phase[order]
    rv = rv[order]
    err = err[order]

    phase_two_cycles = np.concatenate((phase, phase + 1.0))
    rv_two_cycles = np.concatenate((rv, rv))
    err_two_cycles = np.concatenate((err, err))

    fig, ax = setup_figure()

    ax.errorbar(
        phase_two_cycles,
        rv_two_cycles,
        yerr=err_two_cycles,
        fmt="o",
        markersize=6,
        capsize=3,
        elinewidth=1.1,
    )

    ax.axhline(0.0, linewidth=1.0, alpha=0.5)
    ax.set_xlim(0.0, 2.0)
    ax.set_xlabel("Orbital Phase")
    ax.set_ylabel("Relative RV (m/s)")
    ax.set_title(f"{target} — P = {period:g} d")

    fig.tight_layout()

    print(f"Phase folded with P = {period:g} days.")
    print(f"Reference epoch: RJD {reference_time:.6f}")

    return fig


##############################################################
# Main interface
##############################################################


def print_header() -> None:
    """Print the application header."""
    print("=" * 57)
    print("          Wendelstein LBL Analysis Toolkit")
    print("=" * 57)


def choose_analysis() -> str:
    """Ask the user which analysis should be run."""
    print("\nSelect analysis")
    print("1. RV Time Series")
    print("2. Lomb–Scargle Periodogram")
    print("3. Phase Fold")

    while True:
        choice = input("\nChoice [1-3]: ").strip()
        if choice in {"1", "2", "3"}:
            return choice
        print("Please enter 1, 2, or 3.")


def process_one_file(
    rdb_file: Path,
    choice: str,
    save: bool,
    batch_mode: bool,
    min_period: float | None = None,
    max_period: float | None = None,
) -> bool:
    """
    Process one RDB file.

    Returns True when a figure was generated and False when the file was
    skipped or failed.
    """
    target = infer_target_name(rdb_file)

    print("\n" + "-" * 57)
    print(f"File: {rdb_file.name}")

    try:
        data = load_rdb(rdb_file)
        print_data_summary(data, target)

        if choice == "1":
            fig = plot_time_series(data, target)
            plot_name = "timeseries"

        elif choice == "2":
            if min_period is None or max_period is None:
                raise ValueError("Periodogram limits were not provided.")

            fig = plot_periodogram(
                data,
                target,
                min_period=min_period,
                max_period=max_period,
            )
            plot_name = "periodogram"

        else:
            period = ask_float(
                f"Orbital period for {target} in days",
                allow_skip=batch_mode,
            )

            if period is None:
                print(f"Skipped phase fold for {target}.")
                return False

            fig = plot_phase_fold(data, target, period)
            plot_name = "phasefold"

        finish_figure(
            fig,
            source_file=rdb_file,
            target=target,
            plot_name=plot_name,
            save=save,
            batch_mode=batch_mode,
        )
        return True

    except (ValueError, OSError) as exc:
        print(f"Could not process {rdb_file.name}: {exc}", file=sys.stderr)
        return False


def main() -> int:
    """Run the interactive analysis workflow."""
    configure_plotting()
    print_header()

    raw_path = input(
        "\nEnter the path to an LBL RDB file or a directory containing RDB files:\n> "
    ).strip()
    raw_path = raw_path.strip('"').strip("'")

    if not raw_path:
        print("No file or directory was provided.")
        return 1

    try:
        rdb_files = discover_rdb_files(raw_path)
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        return 1

    batch_mode = len(rdb_files) > 1

    if batch_mode:
        print(f"\nFound {len(rdb_files)} RDB files:")
        for rdb_file in rdb_files:
            print(f"  - {rdb_file.name}")
    else:
        print(f"\nSelected: {rdb_files[0].name}")

    choice = choose_analysis()

    min_period: float | None = None
    max_period: float | None = None

    if choice == "2":
        min_period = ask_float(
            "Minimum period in days",
            default=DEFAULT_MIN_PERIOD,
        )
        max_period = ask_float(
            "Maximum period in days",
            default=DEFAULT_MAX_PERIOD,
        )

        if min_period is None or max_period is None:
            return 1

        if min_period >= max_period:
            print(
                "\nError: Minimum period must be smaller than maximum period.",
                file=sys.stderr,
            )
            return 1

    if batch_mode:
        save = ask_yes_no(
            "\nSave PNG and PDF figures for all successfully processed files?",
            default=True,
        )
    else:
        save = ask_yes_no(
            "\nSave the figure as PNG and PDF?",
            default=False,
        )

    generated = 0

    try:
        for rdb_file in rdb_files:
            success = process_one_file(
                rdb_file=rdb_file,
                choice=choice,
                save=save,
                batch_mode=batch_mode,
                min_period=min_period,
                max_period=max_period,
            )
            generated += int(success)

        if batch_mode and not save and generated:
            print(
                "\nAll figures are open. Close the figure windows to finish."
            )
            plt.show()

    except KeyboardInterrupt:
        print("\nOperation cancelled.")
        plt.close("all")
        return 130

    print("\n" + "=" * 57)
    print(f"Finished. Figures generated: {generated}/{len(rdb_files)}")
    print("=" * 57)

    return 0 if generated else 1


if __name__ == "__main__":
    raise SystemExit(main())
