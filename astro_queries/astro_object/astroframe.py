from ..queries import get_detections_and_force_photometry_by_oid
from psql_connector import PgConnector
import matplotlib.pyplot as plt
import logging


class AstroFrame:

    def __init__(self, oid: str = None, **kwargs) -> None:

        self.oid = oid
        self.detections = None
        self.forced_photometry = None
        self.detections_filtered = None
        self.forced_photometry_filtered = None
        self.psql_conn = PgConnector()

    def _connect_to_database(self, **kwargs):
        self.psql_conn.create_conn(**kwargs)

    def _fetch_lightcurve(self):

        if self.psql_conn.conn is not None:
            logging.info("Fetching raw lightcurves")
            astro_query_dd, astro_query_ff = get_detections_and_force_photometry_by_oid(
                oid=self.oid
            )
            df_dd = self.psql_conn.sqlio_query(astro_query_dd)
            df_ff = self.psql_conn.sqlio_query(astro_query_ff)

            self.detections = df_dd
            self.forced_photometry = df_ff
        else:
            logging.warning(
                "PSQL connector was not initializateed, try to create db connection with '_connect_to_database' method."
            )

    def deduplicate(self):
        self.detections_filtered = None
        self.forced_photometry_filtered = None

    def plot_curves(self, **kwargs):
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        plt.close()
