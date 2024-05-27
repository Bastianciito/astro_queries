from ..queries import get_detections_and_force_photometry_by_oid
from psql_connector import PgConnector
import matplotlib.pyplot as plt
import logging
import numpy as np


class AstroFrame:

    def __init__(self, oid: str = None, **kwargs) -> None:

        self.oid = oid
        self.detections = None
        self.forced_photometry = None
        self.detections_filtered = None
        self.forced_photometry_filtered = None
        self.psql_conn = PgConnector()

        # From Masci+2023 sections 6.4 and 6.5
        self.SNT = 3.0
        self.SNU = 5.0

    def _connect_to_database(self, **kwargs):
        self.psql_conn.create_conn(**kwargs)

    def magdiff2flux_uJy(self, df=None, col_mag=None, col_isdiffpos="isdiffpos"):
        return 10.0 ** (-0.4 * (df[col_mag] - 23.9)) * df[col_isdiffpos]

    def magtot2flux_uJy(self, df=None, col_mag=None):
        return 10.0 ** (-0.4 * (df[col_mag] - 23.9))

    def fluxerr(self, df=None, col_magerr=None, col_flux=None):
        return df[col_magerr].abs() * df[col_flux].abs()

    def flux_uJy2magupperlim(self, df=None, col_fluxerr=None):
        return -2.5 * np.log10(self.SNU * df[col_fluxerr].abs()) + 23.9

    def flux_dets(self, df_dets):
        df = df_dets.copy()

        df["fluxdiff_uJy"] = self.magdiff2flux_uJy(df=df, col_mag="magpsf")
        df["fluxerrdiff_uJy"] = self.fluxerr(
            df=df, col_magerr="sigmapsf", col_flux="fluxdiff_uJy"
        )

        if "magpsf_corr" in df.columns:
            df["fluxtot_uJy"] = self.magtot2flux_uJy(df=df, col_mag="magpsf_corr")
            mask = ~df["sigmapsf_corr_ext"].isna() & ~df["sigmapsf_corr_ext"].isna()
            df.loc[mask, "fluxerrtot_uJy"] = self.fluxerr(
                df=df[mask], col_magerr="sigmapsf_corr_ext", col_flux="fluxtot_uJy"
            )
            df.loc[~mask, "fluxerrtot_uJy"] = np.nan

        return df

    def flux_forced(self, df_forced):
        df = df_forced.copy()

        df["fluxdiff_uJy_forced"] = self.magdiff2flux_uJy(df=df, col_mag="mag")
        df["fluxerrdiff_uJy_forced"] = self.fluxerr(
            df=df, col_magerr="e_mag", col_flux="fluxdiff_uJy_forced"
        )

        if "mag_corr" in df.columns:
            df["fluxtot_uJy_forced"] = self.magtot2flux_uJy(df=df, col_mag="mag_corr")
            mask = ~df["e_mag_corr_ext"].isna() & ~df["fluxtot_uJy_forced"].isna()
            df.loc[mask, "fluxerrtot_uJy_forced"] = self.fluxerr(
                df=df[mask], col_magerr="e_mag_corr_ext", col_flux="fluxtot_uJy_forced"
            )
            df.loc[~mask, "fluxerrtot_uJy_forced"] = np.nan

        return df

    def _get_flux_uJd(self):
        pass

    def _fetch_lightcurve(self):

        if self.psql_conn.conn is None:
            logging.warning(
                "PSQL connector was not initializateed, try to create db connection with '_connect_to_database' method."
            )
        if self.oid is None:
            logging.warning("oid was not initializateed.")

        if self.oid is not None and self.psql_conn is not None:
            logging.info("Fetching raw lightcurves")
            astro_query_dd, astro_query_ff = get_detections_and_force_photometry_by_oid(
                oid=self.oid
            )
            df_dd = self.psql_conn.sqlio_query(astro_query_dd)
            df_ff = self.psql_conn.sqlio_query(astro_query_ff)

            self.detections = df_dd
            self.forced_photometry = df_ff

    def deduplicate(self):
        self.detections_filtered = None
        self.forced_photometry_filtered = None

    def plot_curves_ztf_gr(self, **kwargs):

        plt.style.use("bmh")
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))

        if self.detections is not None:

            g_band = self.detections.loc[self.detections.fid == 1]
            r_band = self.detections.loc[self.detections.fid == 2]

            if len(g_band) > 0:
                ax[0].errorbar(
                    g_band.mjd,
                    g_band.magpsf,
                    yerr=g_band.sigmapsf,
                    color="g",
                    fmt="o",
                    alpha=0.5,
                    label="Detections g-band",
                )
            if len(r_band) > 0:
                ax[0].errorbar(
                    r_band.mjd,
                    r_band.magpsf,
                    yerr=r_band.sigmapsf,
                    color="r",
                    fmt="o",
                    alpha=0.5,
                    label="Detections r-band",
                )

        if self.forced_photometry is not None:

            g_band = self.forced_photometry.loc[self.forced_photometry.fid == 1]
            r_band = self.forced_photometry.loc[self.forced_photometry.fid == 2]

            if len(g_band) > 0:
                ax[0].plot(
                    g_band.mjd,
                    g_band.mag,
                    "ro",
                    marker="s",
                    alpha=0.5,
                    label="Forced Photometry g-band",
                )
            if len(r_band) > 0:
                ax[0].plot(
                    r_band.mjd,
                    r_band.mag,
                    "ro",
                    marker="s",
                    alpha=0.5,
                    label="Forced Photometry r-band",
                )

        if "time_lapse" in kwargs.keys():
            try:
                ax[0].set_xlim(kwargs["time_lapse"])
            except:
                pass
        ax[0].grid(False)
        ax[1].grid(False)
        ax[0].set_title("Lightcurve without deduplication")
        ax[1].set_title("Lightcurve with deduplication")
        ax[0].set_xlabel("[MJD]")
        ax[1].set_xlabel("[MJD]")
        ax[0].legend()
        ax[1].legend()
        return fig, ax
