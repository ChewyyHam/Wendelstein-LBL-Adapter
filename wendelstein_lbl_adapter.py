#!/usr/bin/env python3
"""
Wendelstein -> LBL adapter

Convert Wendelstein Observatory *_ods.fits spectra into FITS files that are
compatible with the LBL generic-instrument workflow.

Notes:
    - Both physical-order limits are inclusive.
      Example: order_lower=84 and order_upper=114 keeps orders 84...114.
    - snr_ext is an array index after the selected orders have been sorted in
      descending order.
    - Target coordinates are read from PERA2000 / PEDE2000 when available.
      If valid coordinates cannot be read from the FITS header, the user is
      asked to enter them manually. No default target coordinates are used.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.constants import c
from astropy.coordinates import SkyCoord
from astropy.coordinates.solar_system import get_body_barycentric
from astropy.io import fits
from astropy.table import MaskedColumn, Table, join
from astropy.time import Time
import astropy.units as u


DEFAULT_ORDER_LOWER = 84
DEFAULT_ORDER_UPPER = 114
DEFAULT_SNR_EXT = 10
DEFAULT_AIRMASS = 1.325
DEFAULT_BERV = 0.0


def der_snr(flux: np.ndarray) -> float:
    """Compute signal-to-noise ratio using the DER_SNR estimator."""
    flux = np.asarray(flux)
    flux = flux[flux != 0.0]
    flux = flux[np.isfinite(flux)]
    n = len(flux)

    if n <= 4:
        return 0.0

    signal = np.median(flux)
    noise = 0.6052697 * np.median(
        np.abs(2.0 * flux[2 : n - 2] - flux[0 : n - 4] - flux[4:n])
    )

    if noise == 0 or not np.isfinite(noise):
        return 0.0

    return float(signal / noise)


def ask_path(prompt: str, default: Optional[str] = None, must_exist: bool = False) -> Path:
    """Ask interactively for a folder path."""
    while True:
        suffix = f" [{default}]" if default else ""
        value = input(f"{prompt}{suffix}: ").strip()
        if not value and default:
            value = default

        path = Path(value).expanduser()
        if must_exist and not path.exists():
            print(f"  Folder does not exist: {path}")
            continue
        return path


def ask_int(prompt: str, default: int) -> int:
    """Ask interactively for an integer value."""
    while True:
        value = input(f"{prompt} [{default}]: ").strip()
        if not value:
            return default
        try:
            return int(value)
        except ValueError:
            print("  Please enter an integer.")


def ask_target_coordinates() -> SkyCoord:
    """Ask the user for target coordinates and validate them with Astropy."""
    print(
        "\nTarget coordinates could not be read from the FITS header.\n"
        "Please enter them manually. No default coordinates will be used.\n\n"
        "Accepted examples:\n"
        "  RA  = 01:34:32.18     (hour:minute:second)\n"
        "  DEC = +68:56:55.07    (degree:arcminute:arcsecond)\n\n"
        "Decimal degrees are also accepted:\n"
        "  RA  = 23.634083\n"
        "  DEC = 68.948631\n"
    )

    while True:
        ra_text = input("Target RA: ").strip()
        dec_text = input("Target DEC: ").strip()

        if not ra_text or not dec_text:
            print("  Both RA and DEC are required. Please try again.\n")
            continue

        try:
            ra_unit = u.hourangle if ":" in ra_text else u.deg
            target = SkyCoord(
                ra=ra_text,
                dec=dec_text,
                unit=(ra_unit, u.deg),
                frame="icrs",
            )
        except (ValueError, TypeError) as exc:
            print(
                "\n  Invalid coordinate format.\n"
                f"  Astropy message: {exc}\n"
                "  Please use one of the example formats shown above.\n"
            )
            continue

        print(
            "\nCoordinates accepted:\n"
            f"  RA  = {target.ra.deg:.8f} deg\n"
            f"  DEC = {target.dec.deg:.8f} deg\n"
        )
        return target


def get_target_from_header(header: fits.Header) -> SkyCoord:
    """
    Return target coordinates from the FITS header.

    PERA2000 and PEDE2000 are used when both are present and valid. If either
    keyword is missing, empty, or invalid, the user is asked to enter the
    coordinates manually. The adapter never substitutes default coordinates.
    """
    ra = header.get("PERA2000")
    dec = header.get("PEDE2000")

    if ra not in (None, "") and dec not in (None, ""):
        try:
            target = SkyCoord(
                ra=float(ra) * u.deg,
                dec=float(dec) * u.deg,
                frame="icrs",
            )
            return target
        except (ValueError, TypeError) as exc:
            print(
                "  PERA2000 / PEDE2000 were found but could not be interpreted.\n"
                f"  Header values: RA={ra!r}, DEC={dec!r}\n"
                f"  Astropy message: {exc}"
            )
    else:
        missing = []
        if ra in (None, ""):
            missing.append("PERA2000")
        if dec in (None, ""):
            missing.append("PEDE2000")
        print(f"  Missing coordinate keyword(s): {', '.join(missing)}")

    return ask_target_coordinates()


def compute_bjd(date_obs: str, exptime_seconds: float, header: fits.Header) -> tuple[float, float, float]:
    """Compute MJDMID, RJD and BJD from DATE-OBS and exposure time."""
    mjd_mid = Time(date_obs, format="isot", scale="utc").jd + exptime_seconds / 2 / 86400
    rjd = mjd_mid - 2400000.5

    target = get_target_from_header(header)
    time_tt = Time(mjd_mid, format="jd", scale="tt")
    barycentric_offset = get_body_barycentric("earth", time_tt).dot(target.cartesian) / c
    barycentric_offset_days = barycentric_offset.to(u.s) / (86400 * u.s)
    bjd = mjd_mid + barycentric_offset_days.value

    return float(mjd_mid), float(rjd), float(bjd)


def align_orders(data, order_lower: int, order_upper: int) -> Table:
    """Keep the requested physical orders; both limits are inclusive."""
    filtered_data = data[data["fiber"] != "B"]
    filtered_data["order"] = filtered_data["order"].astype(int)

    full_order = np.arange(order_lower, order_upper + 1, 1)
    template = Table({"order": full_order})

    aligned = join(template, filtered_data, keys="order", join_type="left")
    aligned.sort("order")
    aligned.reverse()

    for col in aligned.colnames:
        if isinstance(aligned[col], MaskedColumn) and np.issubdtype(aligned[col].dtype, np.floating):
            aligned[col] = aligned[col].filled(np.nan)

    return aligned


def convert_one_file(
    input_file: Path,
    output_folder: Path,
    blaze_folder: Path,
    order_lower: int,
    order_upper: int,
    snr_ext: int,
    overwrite: bool = True,
) -> None:
    """Convert one Wendelstein *_ods.fits file to LBL-compatible products."""
    file_name = input_file.name
    output_file_flux = output_folder / f"new_{file_name}"
    output_file_blaze = blaze_folder / f"blaze_{file_name}"

    with fits.open(input_file) as hdul:
        header = hdul[0].header
        data = hdul[1].data

        date_obs = header["DATE-OBS"]
        exptime = float(header["EXPOSURE"])
        mjd_mid, rjd, bjd = compute_bjd(date_obs, exptime, header)

        objname = header.get("OBJECT", "UNKNOWN")
        airmass = float(header.get("AIRMASS", DEFAULT_AIRMASS))

        aligned = align_orders(data, order_lower, order_upper)

        flux = np.asarray(aligned["flux_sum"])
        wavelength = np.asarray(aligned["wavelength"]) / 10.0
        blaze = np.asarray(aligned["blaze"])
        blaze = np.nan_to_num(blaze, nan=0.0)

        if snr_ext < 0 or snr_ext >= len(flux):
            raise IndexError(
                f"snr_ext={snr_ext} is outside the selected order array "
                f"with length {len(flux)}. Choose a value from 0 to {len(flux) - 1}."
            )

        snr = der_snr(flux[snr_ext])
        blaze_header = f"blaze_{file_name}"

        primary_hdu = fits.PrimaryHDU(header=fits.Header())
        primary_hdu.header["MJDMID"] = mjd_mid
        primary_hdu.header["MJSTART"] = mjd_mid - exptime / 2 / 86400
        primary_hdu.header["RJD"] = rjd
        primary_hdu.header["BJD"] = bjd
        primary_hdu.header["BERV"] = DEFAULT_BERV
        primary_hdu.header["EXPTIME"] = exptime
        primary_hdu.header["DATE"] = date_obs
        primary_hdu.header["OBJNAME"] = objname
        primary_hdu.header["AIRMASS"] = airmass
        primary_hdu.header["SNR"] = snr
        primary_hdu.header["EXT_SNR"] = snr_ext
        primary_hdu.header["BLAZE"] = blaze_header
        primary_hdu.header["ORD_LOW"] = order_lower
        primary_hdu.header["ORD_UP"] = order_upper

        flux_hdu = fits.ImageHDU(data=flux, name="FLUX")
        wavelength_hdu = fits.ImageHDU(data=wavelength, name="WAVELENGTH")
        blaze_hdu = fits.ImageHDU(data=blaze, name="BLAZE")

        fits.HDUList([primary_hdu, flux_hdu, wavelength_hdu, blaze_hdu]).writeto(
            output_file_flux, overwrite=overwrite
        )
        fits.HDUList([primary_hdu, blaze_hdu, wavelength_hdu]).writeto(
            output_file_blaze, overwrite=overwrite
        )

    print(f"  Created science file: {output_file_flux}")
    print(f"  Created blaze file:   {output_file_blaze}")


def convert_wendelstein_to_lbl(
    input_folder: Path,
    output_folder: Path,
    blaze_folder: Path,
    order_lower: int = DEFAULT_ORDER_LOWER,
    order_upper: int = DEFAULT_ORDER_UPPER,
    snr_ext: int = DEFAULT_SNR_EXT,
    overwrite: bool = True,
) -> None:
    """Convert all Wendelstein *_ods.fits files in input_folder."""
    input_folder = Path(input_folder).expanduser()
    output_folder = Path(output_folder).expanduser()
    blaze_folder = Path(blaze_folder).expanduser()

    if not input_folder.exists():
        raise FileNotFoundError(f"Input folder does not exist: {input_folder}")

    if order_upper < order_lower:
        raise ValueError("order_upper must be greater than or equal to order_lower.")

    output_folder.mkdir(parents=True, exist_ok=True)
    blaze_folder.mkdir(parents=True, exist_ok=True)

    fits_files = sorted(input_folder.rglob("*_ods.fits"))
    fits_files = [path for path in fits_files if "new" not in path.name]

    if not fits_files:
        print(f"No *_ods.fits files found in {input_folder}")
        return

    print("\nConversion settings")
    print("-------------------")
    print(f"Input folder:   {input_folder}")
    print(f"Science output: {output_folder}")
    print(f"Blaze output:   {blaze_folder}")
    print(f"Orders:         {order_lower} to {order_upper} (inclusive)")
    print(f"snr_ext index:  {snr_ext}")
    print(f"Files found:    {len(fits_files)}\n")

    success = 0
    failed = 0
    for input_file in fits_files:
        print(f"Processing {input_file.name}")
        try:
            convert_one_file(
                input_file=input_file,
                output_folder=output_folder,
                blaze_folder=blaze_folder,
                order_lower=order_lower,
                order_upper=order_upper,
                snr_ext=snr_ext,
                overwrite=overwrite,
            )
            success += 1
        except Exception as exc:
            failed += 1
            print(f"  Failed: {exc}")

    print("\nDone")
    print("----")
    print(f"Successfully converted: {success}")
    print(f"Failed:                 {failed}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Wendelstein *_ods.fits spectra into LBL-compatible FITS files."
    )
    parser.add_argument("--input-folder", type=Path, help="Folder containing Wendelstein *_ods.fits files.")
    parser.add_argument("--output-folder", type=Path, help="LBL science output folder.")
    parser.add_argument("--blaze-folder", type=Path, help="LBL calib/blaze output folder.")
    parser.add_argument("--order-lower", type=int, default=None, help="Lower physical order limit, included.")
    parser.add_argument("--order-upper", type=int, default=None, help="Upper physical order limit, included.")
    parser.add_argument("--snr-ext", type=int, default=None, help="Index of selected order used by DER_SNR.")
    parser.add_argument("--no-overwrite", action="store_true", help="Do not overwrite existing output files.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("\nWendelstein -> LBL adapter")
    print("==========================")
    print("Press Enter to accept the default shown in brackets.\n")

    input_folder = args.input_folder or ask_path(
        "Input folder containing Wendelstein *_ods.fits files", must_exist=True
    )
    output_folder = args.output_folder or ask_path("Output folder for LBL science files")
    blaze_folder = args.blaze_folder or ask_path("Output folder for LBL blaze/calib files")

    order_lower = args.order_lower if args.order_lower is not None else ask_int(
        "Lower order limit, included", DEFAULT_ORDER_LOWER
    )
    order_upper = args.order_upper if args.order_upper is not None else ask_int(
        "Upper order limit, included", DEFAULT_ORDER_UPPER
    )
    snr_ext = args.snr_ext if args.snr_ext is not None else ask_int(
        "snr_ext: order index used by DER_SNR", DEFAULT_SNR_EXT
    )

    convert_wendelstein_to_lbl(
        input_folder=input_folder,
        output_folder=output_folder,
        blaze_folder=blaze_folder,
        order_lower=order_lower,
        order_upper=order_upper,
        snr_ext=snr_ext,
        overwrite=not args.no_overwrite,
    )


if __name__ == "__main__":
    main()
