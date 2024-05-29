from ..queries import get_detections_and_force_photometry_by_oid
from psql_connector import PgConnector
import matplotlib.pyplot as plt
import logging
import numpy as np
import pandas as pd


class AstroFrame:

    def __init__(self, oid: str = None, credentials: dict = None, **kwargs) -> None:

        self.oid = oid
        self.detections = None
        self.forced_photometry = None
        self.detections_filtered = None
        self.forced_photometry_filtered = None
        self.psql_conn = PgConnector(**{"credentials": credentials})
        # From Masci+2023 sections 6.4 and 6.5
        self.SNT = 3.0
        self.SNU = 5.0

    def _connect_to_database(self, **kwargs):
        self.psql_conn.create_conn(**kwargs)
        self.psql_conn.create_sqlalchemy_engine(**kwargs)

    # Transformations extracted from https://github.com/alercebroker/usecases/blob/master/notebooks/CosmicStreams_Dec2023_tutorial/ALeRCE_Basics_CosmicStreams.ipynb
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
        logging.info(
            "Creating Flux diff columns to detections and forced photometry dataframes"
        )
        self.detections = self.flux_dets(df_dets=self.detections)
        self.forced_photometry = self.flux_forced(df_forced=self.forced_photometry)
        logging.info("Columns created")

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
            # get fluxes
            self._get_flux_uJd()
            self.deduplicate()

    def deduplicate_by_mjd(self):

        g_band = self.detections.loc[self.detections.fid == 1]
        r_band = self.detections.loc[self.detections.fid == 2]

        g_band_forced = self.forced_photometry.loc[self.forced_photometry.fid == 1]
        r_band_forced = self.forced_photometry.loc[self.forced_photometry.fid == 2]
        logging.info("Before deduplication by mjd ....")
        logging.info(f"num detections g:{len(g_band)} r:{len(r_band)}")
        logging.info(
            f"num forced photometry g:{len(g_band_forced)} r:{len(r_band_forced)}"
        )
        g_band = g_band.drop_duplicates(subset=["mjd"])
        r_band = r_band.drop_duplicates(subset=["mjd"])

        g_band_forced = g_band_forced.drop_duplicates(subset=["mjd"])
        r_band_forced = r_band_forced.drop_duplicates(subset=["mjd"])

        self.detections_filtered = pd.concat([g_band, r_band], ignore_index=True)
        self.forced_photometry_filtered = pd.concat(
            [g_band_forced, r_band_forced], ignore_index=True
        )
        logging.info("After deduplication by mjd ....")
        logging.info(f"num detections g:{len(g_band)} r:{len(r_band)}")
        logging.info(
            f"num forced photometry g:{len(g_band_forced)} r:{len(r_band_forced)}"
        )

    def deduplicate_by_oid_pid(self):
        # code provided by ALE

        g_band = self.detections_filtered.loc[self.detections_filtered.fid == 1]
        r_band = self.detections_filtered.loc[self.detections_filtered.fid == 2]

        g_band_forced = self.forced_photometry_filtered.loc[
            self.forced_photometry_filtered.fid == 1
        ]
        r_band_forced = self.forced_photometry_filtered.loc[
            self.forced_photometry_filtered.fid == 2
        ]
        logging.info("Before deduplication by oid-pid ....")
        logging.info(f"num detections g:{len(g_band)} r:{len(r_band)}")
        logging.info(
            f"num forced photometry g:{len(g_band_forced)} r:{len(r_band_forced)}"
        )

        objs_epochs = self.detections_filtered.sort_values(by="mjd").copy()
        objs_epochs.rename(columns={"magpsf": "mag"}, inplace=True)
        objs_epochs["detected"] = True
        aux = self.forced_photometry_filtered.sort_values(by="mjd").copy()
        aux["detected"] = False
        aux = aux[~aux["pid"].isin(objs_epochs["pid"])].copy()
        # addapted to work with astroframe
        self.forced_photometry_filtered = aux
        self.detections_filtered = objs_epochs

        g_band = self.detections_filtered.loc[self.detections_filtered.fid == 1]
        r_band = self.detections_filtered.loc[self.detections_filtered.fid == 2]

        g_band_forced = self.forced_photometry_filtered.loc[
            self.forced_photometry_filtered.fid == 1
        ]
        r_band_forced = self.forced_photometry_filtered.loc[
            self.forced_photometry_filtered.fid == 2
        ]
        logging.info("Before deduplication by oid-pid ....")
        logging.info(f"num detections g:{len(g_band)} r:{len(r_band)}")
        logging.info(
            f"num forced photometry g:{len(g_band_forced)} r:{len(r_band_forced)}"
        )

    def deduplicate(self):

        self.deduplicate_by_mjd()
        self.deduplicate_by_oid_pid()

    def define_cadence(self, delta):
        def find_slot(value, arr):
            return np.argmax(value < arr)

        g_band = self.detections_filtered.loc[self.detections_filtered.fid == 1].copy()
        r_band = self.detections_filtered.loc[self.detections_filtered.fid == 2].copy()

        g_band_forced = self.forced_photometry_filtered.loc[
            self.forced_photometry_filtered.fid == 1
        ].copy()
        r_band_forced = self.forced_photometry_filtered.loc[
            self.forced_photometry_filtered.fid == 2
        ].copy()

        grid_g = np.arange(g_band.mjd.min(), g_band.mjd.max(), delta)
        grid_r = np.arange(r_band.mjd.min(), r_band.mjd.max(), delta)
        grid_forced_g = np.arange(
            g_band_forced.mjd.min(), g_band_forced.mjd.max(), delta
        )
        grid_forced_r = np.arange(
            r_band_forced.mjd.min(), r_band_forced.mjd.max(), delta
        )

        g_band["time_slot"] = g_band.apply(
            lambda x: find_slot(x["mjd"], grid_g), axis=1
        )
        r_band["time_slot"] = g_band.apply(
            lambda x: find_slot(x["mjd"], grid_r), axis=1
        )
        g_band_forced["time_slot"] = g_band.apply(
            lambda x: find_slot(x["mjd"], grid_forced_g), axis=1
        )
        r_band_forced["time_slot"] = g_band.apply(
            lambda x: find_slot(x["mjd"], grid_forced_r), axis=1
        )

        g_band.drop_duplicates(subset=["time_slot"], inplace=True)
        r_band.drop_duplicates(subset=["time_slot"], inplace=True)
        g_band_forced.drop_duplicates(subset=["time_slot"], inplace=True)
        r_band_forced.drop_duplicates(subset=["time_slot"], inplace=True)

        return g_band, r_band, g_band_forced, r_band_forced

    def plot_curves_ztf_gr_cadence(self, **kwargs):

        plt.style.use("bmh")
        fig, ax = plt.subplots(1, 2, figsize=(16, 8))
        meta = {}

        g_band = self.detections_filtered[self.detections_filtered.fid == 1]
        r_band = self.detections_filtered[self.detections_filtered.fid == 2]

        g_band_forced = self.forced_photometry_filtered[
            self.forced_photometry_filtered.fid == 1
        ]
        r_band_forced = self.forced_photometry_filtered[
            self.forced_photometry_filtered.fid == 2
        ]

        g_band_filter, r_band_filter, g_band_forced_filter, r_band_forced_filter = (
            self.define_cadence(**kwargs)
        )

        meta["n_det_g"] = len(g_band)
        meta["n_det_r"] = len(r_band)
        meta["n_det_g_filter"] = len(g_band_filter)
        meta["n_det_r_filter"] = len(r_band_filter)

        meta["n_forced_g"] = len(g_band_forced)
        meta["n_forced_r"] = len(r_band_forced)
        meta["n_forced_g_filter"] = len(g_band_forced_filter)
        meta["n_forced_r_filter"] = len(r_band_forced_filter)

        for detections_band, color in zip([g_band, r_band], ["g", "r"]):

            if len(detections_band) > 0:
                ax[0].errorbar(
                    detections_band.mjd,
                    detections_band.fluxdiff_uJy,
                    yerr=detections_band.fluxerrdiff_uJy,
                    color=color,
                    fmt="o",
                    alpha=0.5,
                    label=f"Detections {color}-band",
                )

        for detections_band, color in zip([g_band_filter, r_band_filter], ["g", "r"]):
            if len(detections_band) > 0:
                ax[1].errorbar(
                    detections_band.mjd,
                    detections_band.fluxdiff_uJy,
                    yerr=detections_band.fluxerrdiff_uJy,
                    color=color,
                    fmt="o",
                    alpha=0.5,
                    label=f"Detections {color}-band",
                )

        for detections_band, color in zip([g_band_forced, r_band_forced], ["g", "r"]):
            if len(detections_band) > 0:
                ax[0].plot(
                    g_band.mjd,
                    g_band.fluxdiff_uJy_forced,
                    f"{color}s",
                    # marker="s",
                    alpha=0.5,
                    label=f"Forced Photometry {color}-band",
                )

        for detections_band, color in zip(
            [g_band_forced_filter, r_band_forced_filter], ["g", "r"]
        ):
            if len(detections_band) > 0:
                ax[0].plot(
                    g_band.mjd,
                    g_band.fluxdiff_uJy_forced,
                    f"{color}s",
                    # marker="s",
                    alpha=0.5,
                    label=f"Forced Photometry {color}-band",
                )

        if "time_lapse" in kwargs.keys():
            try:
                ax[0].set_xlim(kwargs["time_lapse"])
                ax[1].set_xlim(kwargs["time_lapse"])
            except:
                pass

        ax[0].grid(False)
        ax[1].grid(False)
        ax[0].set_title(
            f"Lightcurve without deduplication\n[n_det_g={meta['n_det_g']}]-[n_forced_g={meta['n_forced_g']}]\n[n_det_r={meta['n_det_r']}]-[n_forced_r={meta['n_det_r']}]"
        )
        ax[1].set_title(
            f"Lightcurve with deduplication\n[n_det_g={meta['n_det_g_filter']}]-[n_forced_g={meta['n_forced_g_filter']}]\n[n_det_r={meta['n_det_r_filter']}]-[n_forced_r={meta['n_det_r_filter']}]"
        )
        ax[0].set_xlabel("[MJD]")
        ax[1].set_xlabel("[MJD]")
        ax[0].set_ylabel("Flux difference [uJy]")
        ax[1].set_ylabel("Flux difference [uJy]")
        ax[0].legend()
        ax[1].legend()
        return fig, ax

    def plot_curves_ztf_gr(self, **kwargs):

        plt.style.use("bmh")
        fig, ax = plt.subplots(1, 2, figsize=(16, 8))

        meta = {}
        if self.detections is not None:
            g_band = self.detections.loc[self.detections.fid == 1]
            r_band = self.detections.loc[self.detections.fid == 2]

            g_band_filter = self.detections_filtered[self.detections_filtered.fid == 1]
            r_band_filter = self.detections_filtered[self.detections_filtered.fid == 2]

            meta["n_det_g"] = len(g_band)
            meta["n_det_r"] = len(r_band)
            meta["n_det_g_filter"] = len(g_band_filter)
            meta["n_det_r_filter"] = len(r_band_filter)

            if len(g_band) > 0:
                ax[0].errorbar(
                    g_band.mjd,
                    g_band.fluxdiff_uJy,
                    yerr=g_band.fluxerrdiff_uJy,
                    color="g",
                    fmt="o",
                    alpha=0.5,
                    label="Detections g-band",
                )
            if len(r_band) > 0:
                ax[0].errorbar(
                    r_band.mjd,
                    r_band.fluxdiff_uJy,
                    yerr=r_band.fluxerrdiff_uJy,
                    color="r",
                    fmt="o",
                    alpha=0.5,
                    label="Detections r-band",
                )

            if len(g_band_filter) > 0:
                ax[1].errorbar(
                    g_band_filter.mjd,
                    g_band_filter.fluxdiff_uJy,
                    yerr=g_band_filter.fluxerrdiff_uJy,
                    color="g",
                    fmt="o",
                    alpha=0.5,
                    label="Detections g-band (filter)",
                )

            if len(r_band_filter) > 0:
                ax[1].errorbar(
                    r_band_filter.mjd,
                    r_band_filter.fluxdiff_uJy,
                    yerr=r_band_filter.fluxerrdiff_uJy,
                    color="r",
                    fmt="o",
                    alpha=0.5,
                    label="Detections r-band (filter)",
                )

        if self.forced_photometry is not None:

            g_band = self.forced_photometry.loc[self.forced_photometry.fid == 1]
            r_band = self.forced_photometry.loc[self.forced_photometry.fid == 2]

            g_band_filter = self.forced_photometry_filtered[
                self.forced_photometry_filtered.fid == 1
            ]
            r_band_filter = self.forced_photometry_filtered[
                self.forced_photometry_filtered.fid == 2
            ]

            meta["n_forced_g"] = len(g_band)
            meta["n_forced_r"] = len(r_band)
            meta["n_forced_g_filter"] = len(g_band_filter)
            meta["n_forced_r_filter"] = len(r_band_filter)

            if len(g_band) > 0:
                ax[0].plot(
                    g_band.mjd,
                    g_band.fluxdiff_uJy_forced,
                    "gs",
                    # marker="s",
                    alpha=0.5,
                    label="Forced Photometry g-band",
                )
            if len(r_band) > 0:
                ax[0].plot(
                    r_band.mjd,
                    r_band.fluxdiff_uJy_forced,
                    "rs",
                    # marker="s",
                    alpha=0.5,
                    label="Forced Photometry r-band",
                )

            if len(g_band) > 0:
                ax[1].plot(
                    g_band_filter.mjd,
                    g_band_filter.fluxdiff_uJy_forced,
                    "gs",
                    # marker="s",
                    alpha=0.5,
                    label="Forced Photometry g-band",
                )
            if len(r_band) > 0:
                ax[1].plot(
                    r_band_filter.mjd,
                    r_band_filter.fluxdiff_uJy_forced,
                    "rs",
                    # marker="s",
                    alpha=0.5,
                    label="Forced Photometry r-band",
                )

        if "time_lapse" in kwargs.keys():
            try:
                ax[0].set_xlim(kwargs["time_lapse"])
                ax[1].set_xlim(kwargs["time_lapse"])
            except:
                pass

        ax[0].grid(False)
        ax[1].grid(False)
        ax[0].set_title(
            f"Lightcurve without deduplication\n[n_det_g={meta['n_det_g']}]-[n_forced_g={meta['n_forced_g']}]\n[n_det_r={meta['n_det_r']}]-[n_forced_r={meta['n_det_r']}]"
        )
        ax[1].set_title(
            f"Lightcurve with deduplication\n[n_det_g={meta['n_det_g_filter']}]-[n_forced_g={meta['n_forced_g_filter']}]\n[n_det_r={meta['n_det_r_filter']}]-[n_forced_r={meta['n_det_r_filter']}]"
        )
        ax[0].set_xlabel("[MJD]")
        ax[1].set_xlabel("[MJD]")
        ax[0].set_ylabel("Flux difference [uJy]")
        ax[1].set_ylabel("Flux difference [uJy]")
        ax[0].legend()
        ax[1].legend()
        return fig, ax
